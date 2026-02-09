#!/usr/bin/env python3
"""
Harvest Failed Pharmaceutical Trials from ClinicalTrials.gov API v2.

Criteria:
  - lastSummaryUpdateDate (LastUpdatePostDate) within last 24 months
  - overallStatus = TERMINATED or WITHDRAWN, hasResults = True
  - Target major sponsors: Pfizer, AbbVie, Merck

Output: Olympus_Research/{NCTId}/
  - raw_api_results.json   (full study from API)
  - recovery_pitch_draft.md (Architecture Review template + recoupment pillars)

Rate limit: 10 requests per second (configurable).
Data integrity: Pulls PrimaryOutcomeMeasure and StatisticalAnalysis from API.
"""

import argparse
import json
import time
from datetime import datetime, UTC, timedelta
from pathlib import Path
from typing import Any

import requests

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
API_BASE = "https://clinicaltrials.gov/api/v2"
DEFAULT_OUTPUT_DIR = Path("Olympus_Research")
REQUESTS_PER_SECOND = 10
MIN_INTERVAL = 1.0 / REQUESTS_PER_SECOND  # 0.1 seconds between requests

MAJOR_SPONSORS = [
    "Pfizer",
    "AbbVie",
    "Merck",
]


class RateLimiter:
    """Enforce max N requests per second."""

    def __init__(self, requests_per_second: float = REQUESTS_PER_SECOND):
        self.interval = 1.0 / requests_per_second
        self.last_call = 0.0

    def wait(self) -> None:
        elapsed = time.monotonic() - self.last_call
        if elapsed < self.interval:
            time.sleep(self.interval - elapsed)
        self.last_call = time.monotonic()


def date_24_months_ago() -> str:
    """Return YYYY-MM-DD for 24 months ago (UTC)."""
    d = datetime.now(UTC) - timedelta(days=730)
    return d.strftime("%Y-%m-%d")


def search_studies(
    sponsor: str,
    from_date: str,
    limiter: RateLimiter,
    page_size: int = 100,
) -> list[dict[str, Any]]:
    """
    Search API v2 for studies: TERMINATED/WITHDRAWN, hasResults, sponsor, updated since from_date.
    Returns list of study summary objects (may include hasResults; full fetch needed for outcomes).
    """
    all_studies: list[dict[str, Any]] = []
    params: dict[str, Any] = {
        "query.spons": sponsor,
        "filter.overallStatus": "TERMINATED,WITHDRAWN",
        "query.term": f"AREA[LastUpdatePostDate]RANGE[{from_date},MAX]",
        "pageSize": page_size,
        "format": "json",
    }
    url = f"{API_BASE}/studies"

    while True:
        limiter.wait()
        resp = requests.get(url, params=params, timeout=30)
        resp.raise_for_status()
        data = resp.json()
        studies = data.get("studies") or []
        for s in studies:
            # Prefer to filter by hasResults here if API returns it; otherwise filter after full fetch
            if "hasResults" in s and s.get("hasResults") is not True:
                continue
            all_studies.append(s)
        next_token = data.get("nextPageToken")
        if not next_token:
            break
        params["pageToken"] = next_token

    return all_studies


def get_full_study(nct_id: str, limiter: RateLimiter) -> dict[str, Any] | None:
    """Fetch full study by NCT ID. Returns None on failure."""
    limiter.wait()
    url = f"{API_BASE}/studies/{nct_id}"
    try:
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
        return resp.json()
    except Exception:
        return None


def extract_primary_outcomes(study: dict[str, Any]) -> list[dict[str, Any]]:
    """Extract PrimaryOutcome measures from protocolSection (data integrity)."""
    try:
        proto = study.get("protocolSection") or {}
        outcomes = proto.get("outcomesModule") or {}
        return list(outcomes.get("primaryOutcomes") or [])
    except Exception:
        return []


def extract_results_and_analysis(study: dict[str, Any]) -> dict[str, Any]:
    """Extract resultsSection: outcomes, p-values, confidence intervals (data integrity)."""
    try:
        results = study.get("resultsSection") or {}
        # Preserve key outcome and statistical analysis blocks
        out = {}
        if "outcomeMeasuresModule" in results:
            out["outcomeMeasuresModule"] = results["outcomeMeasuresModule"]
        if "adverseEventsModule" in results:
            out["adverseEventsModule"] = results["adverseEventsModule"]
        return out
    except Exception:
        return {}


def get_nct_id(study: dict[str, Any]) -> str | None:
    """Get NCT ID from study (search result or full study)."""
    if "protocolSection" in study:
        ident = (study["protocolSection"] or {}).get("identificationModule") or {}
        return ident.get("nctId")
    return study.get("nctId")


def get_brief_title(study: dict[str, Any]) -> str:
    """Get brief title from study."""
    try:
        ident = (study.get("protocolSection") or {}).get("identificationModule") or {}
        return ident.get("briefTitle") or "Unknown"
    except Exception:
        return "Unknown"


def get_lead_sponsor(study: dict[str, Any]) -> str:
    """Get lead sponsor name."""
    try:
        spons = (study.get("protocolSection") or {}).get("sponsorCollaboratorsModule") or {}
        lead = spons.get("leadSponsor") or {}
        return lead.get("name") or "Unknown"
    except Exception:
        return "Unknown"


def get_conditions(study: dict[str, Any]) -> list[str]:
    """Get conditions list."""
    try:
        cond = (study.get("protocolSection") or {}).get("conditionsModule") or {}
        return list(cond.get("conditions") or [])
    except Exception:
        return []


def write_recovery_pitch_draft(
    folder: Path,
    study: dict[str, Any],
    primary_outcomes: list[dict[str, Any]],
    results_section: dict[str, Any],
) -> None:
    """Write recovery_pitch_draft.md with Architecture Review format and recoupment pillars."""
    nct_id = get_nct_id(study) or "Unknown"
    title = get_brief_title(study)
    sponsor = get_lead_sponsor(study)
    conditions = get_conditions(study)
    conditions_str = ", ".join(conditions) if conditions else "Not specified"

    # Format primary outcomes for display
    outcomes_text = ""
    for i, o in enumerate(primary_outcomes, 1):
        measure = o.get("measure") or o.get("title") or "—"
        outcomes_text += f"{i}. {measure}\n"

    # Results/statistical summary for "where the architecture broke"
    results_text = json.dumps(results_section, indent=2) if results_section else "(No results section in API response.)"

    md = f"""# Recovery Pitch Draft — {nct_id}

**Study:** {title}  
**Sponsor:** {sponsor}  
**Condition(s):** {conditions_str}  
**NCT ID:** {nct_id}

---

## Architecture Review

### Observed Failure Point

- **Trial status:** Terminated/Withdrawn with posted results.
- **Primary outcome measures (from API):**
{outcomes_text}
- **Results / statistical analysis (excerpt):**  
  *(Use p-values and confidence intervals below to pinpoint where the architecture broke.)*

```
{results_text}
```

---

### Olympus Law Deviation

- Identify which **Law of 51** (or internal Olympus rule) was violated (e.g., binding geometry, exposure–response alignment, population selection).
- *To be completed from 35D analysis and manual review.*

---

### Recoupment Strategy

Frame outreach using these three pillars:

**Pillar 1 — Precision Pivot**  
"Our analysis shows the drug didn't fail because of the mechanism, but because of a 3D geometric misalignment in the [Protein] binding site. We have the coordinates to fix it."

**Pillar 2 — Indication Expansion**  
"The architecture that failed in [{conditions_str}] is perfectly suited for [Indication B]. We can demonstrate why with a 90%+ confidence interval."

**Pillar 3 — R&D De-risking**  
"We have identified the specific 'Law of 51' that was violated. Applying this to your current pipeline will prevent the next $500M write-down."

---

## Adversarial Lab Partner Check

- **Data integrity:** Primary outcome measures and results/statistical analysis sections above are pulled directly from the ClinicalTrials.gov API v2 for this study.
- **Genius guardrail:** Do not send these results to the general "Contact Us" email. Use LinkedIn or professional databases to find the **Head of External Innovation** or **VP of Translational Medicine** for this therapeutic area. They are the ones whose bonuses depend on fixing these failures.

---

*Placeholder: attach your manual **olympus_analysis.pdf** (35D output) to this folder when ready.*
"""
    (folder / "recovery_pitch_draft.md").write_text(md, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Harvest failed trials (TERMINATED/WITHDRAWN, hasResults) from ClinicalTrials.gov API v2."
    )
    parser.add_argument(
        "-o", "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Output root directory (default: {DEFAULT_OUTPUT_DIR})",
    )
    parser.add_argument(
        "--sponsors",
        nargs="*",
        default=MAJOR_SPONSORS,
        help="Sponsor names to search (default: Pfizer AbbVie Merck)",
    )
    parser.add_argument(
        "--months",
        type=int,
        default=24,
        help="Only trials updated in the last N months (default: 24)",
    )
    parser.add_argument(
        "--rate",
        type=float,
        default=REQUESTS_PER_SECOND,
        help=f"Max requests per second (default: {REQUESTS_PER_SECOND})",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Only search and print NCT IDs; do not create folders or fetch full studies",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Max number of trials to harvest (default: all); use e.g. 3 for a quick test",
    )
    args = parser.parse_args()

    from_date = (datetime.now(UTC) - timedelta(days=args.months * 30)).strftime("%Y-%m-%d")
    limiter = RateLimiter(requests_per_second=args.rate)
    output_root = args.output_dir
    output_root.mkdir(parents=True, exist_ok=True)

    seen_nct: set[str] = set()
    hit_list: list[dict[str, Any]] = []

    print("Searching ClinicalTrials.gov API v2...")
    print(f"  Sponsors: {', '.join(args.sponsors)}")
    print(f"  Updated since: {from_date}")
    print(f"  Status: TERMINATED, WITHDRAWN; hasResults = True")
    print()

    for sponsor in args.sponsors:
        try:
            studies = search_studies(sponsor, from_date, limiter)
            for s in studies:
                nct_id = get_nct_id(s)
                if not nct_id or nct_id in seen_nct:
                    continue
                # Re-check hasResults from summary if present
                # hasResults filter: if present in summary, require True; else allow and filter after full fetch
                if s.get("hasResults") is False:
                    continue
                seen_nct.add(nct_id)
                hit_list.append({"nctId": nct_id, "sponsor": sponsor, "study_summary": s})
        except Exception as e:
            print(f"  Warning: search for sponsor '{sponsor}' failed: {e}")

    print(f"Found {len(hit_list)} unique trials with results.")
    if not hit_list:
        print("Nothing to save. Exiting.")
        return

    if args.dry_run:
        for h in hit_list:
            print(f"  {h['nctId']}  ({h['sponsor']})")
        return

    to_harvest = hit_list[: args.limit] if args.limit else hit_list
    if args.limit:
        print(f"Harvesting first {len(to_harvest)} trials (--limit {args.limit}).")

    for i, hit in enumerate(to_harvest, 1):
        nct_id = hit["nctId"]
        folder = output_root / nct_id
        folder.mkdir(parents=True, exist_ok=True)

        full = get_full_study(nct_id, limiter)
        if not full:
            print(f"  [{i}/{len(to_harvest)}] {nct_id} — full study fetch failed, skipping.")
            continue
        if full.get("hasResults") is not True:
            print(f"  [{i}/{len(to_harvest)}] {nct_id} — hasResults=False, skipping.")
            continue

        # Save raw API response
        (folder / "raw_api_results.json").write_text(
            json.dumps(full, indent=2), encoding="utf-8"
        )

        primary_outcomes = extract_primary_outcomes(full)
        results_section = extract_results_and_analysis(full)
        write_recovery_pitch_draft(folder, full, primary_outcomes, results_section)

        print(f"  [{i}/{len(to_harvest)}] {nct_id} — {get_brief_title(full)[:50]}...")

    print()
    print(f"Done. Hit list under: {output_root.absolute()}")
    print("  Each folder contains: raw_api_results.json, recovery_pitch_draft.md")
    print("  Add olympus_analysis.pdf manually after 35D review.")


if __name__ == "__main__":
    main()
