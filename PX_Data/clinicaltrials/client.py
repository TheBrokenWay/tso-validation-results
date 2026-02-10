"""
ClinicalTrials.gov API v2 Client

Fetches interventional trial data for PRV-eligible diseases.
Uses the ClinicalTrials.gov v2 REST API (public, no key required).

Respects net_policy and uses retry_transient_http for resilience.
Results are written to PX_Data/clinicaltrials/ as CSV and JSONL.

API docs: https://clinicaltrials.gov/data-api/about-api

Uses only: requests, json, csv, pathlib, datetime, hashlib (+ PX integrations).
"""
from __future__ import annotations

import csv
import hashlib
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

try:
    import requests
except ImportError:
    requests = None  # type: ignore[assignment]

from PX_System.foundation.integrations.net_policy import (
    check_external_access_allowed,
    online_sync_enabled,
)
from PX_System.foundation.integrations.retry import retry_transient_http

# ── Constants ──────────────────────────────────────────────────────────

BASE_URL = "https://clinicaltrials.gov/api/v2"
DATA_DIR = Path(__file__).resolve().parent
SERVICE_NAME = "clinicaltrials"

# PRV-eligible disease search terms (aligned with OLE.py and PX_Domain)
PRV_DISEASE_QUERIES = [
    "nipah", "ebola", "marburg", "lassa fever", "dengue",
    "chikungunya", "zika", "rabies", "tuberculosis", "malaria",
    "leishmaniasis", "chagas disease", "onchocerciasis", "schistosomiasis",
    "lymphatic filariasis", "dracunculiasis", "trachoma", "buruli ulcer",
    "yaws", "leprosy", "cholera", "yellow fever",
]

# Fields to request from the v2 API
FIELDS = [
    "NCTId", "BriefTitle", "Condition", "InterventionName",
    "InterventionType", "Phase", "OverallStatus", "StartDate",
    "PrimaryCompletionDate", "EnrollmentCount", "StudyType",
    "LeadSponsorName",
]


# ── Core Functions ─────────────────────────────────────────────────────

@retry_transient_http(max_attempts=3, delay=2.0, backoff=2.0)
def _fetch_studies(query: str, page_size: int = 100, page_token: Optional[str] = None) -> Dict[str, Any]:
    """Fetch a page of studies from ClinicalTrials.gov v2 API."""
    if requests is None:
        raise ImportError("requests library required for ClinicalTrials.gov API")

    params: Dict[str, Any] = {
        "query.cond": query,
        "filter.overallStatus": "RECRUITING,ACTIVE_NOT_RECRUITING,COMPLETED",
        "fields": ",".join(FIELDS),
        "pageSize": page_size,
    }
    if page_token:
        params["pageToken"] = page_token

    resp = requests.get(f"{BASE_URL}/studies", params=params, timeout=30)
    resp.raise_for_status()
    return resp.json()


def fetch_disease_trials(disease: str, max_pages: int = 5) -> List[Dict[str, Any]]:
    """
    Fetch all trials for a disease (up to max_pages * 100 results).

    Args:
        disease: Disease search term (e.g. "nipah", "malaria")
        max_pages: Maximum number of pages to fetch

    Returns:
        List of study dicts with FIELDS
    """
    if not online_sync_enabled():
        print(f"  [clinicaltrials] Network sync disabled, skipping {disease}", file=sys.stderr)
        return []
    if not check_external_access_allowed(SERVICE_NAME):
        print(f"  [clinicaltrials] Access denied by net_policy for {SERVICE_NAME}", file=sys.stderr)
        return []

    all_studies: List[Dict[str, Any]] = []
    page_token = None

    for page in range(max_pages):
        data = _fetch_studies(disease, page_token=page_token)
        studies = data.get("studies", [])
        if not studies:
            break
        for study in studies:
            proto = study.get("protocolSection", {})
            ident = proto.get("identificationModule", {})
            status_mod = proto.get("statusModule", {})
            design = proto.get("designModule", {})
            conditions = proto.get("conditionsModule", {})
            interventions = proto.get("armsInterventionsModule", {})
            sponsor = proto.get("sponsorCollaboratorsModule", {})

            row = {
                "nct_id": ident.get("nctId", ""),
                "brief_title": ident.get("briefTitle", ""),
                "conditions": "|".join(conditions.get("conditions", [])),
                "interventions": "|".join(
                    i.get("name", "") for i in interventions.get("interventions", [])
                ),
                "intervention_types": "|".join(
                    i.get("type", "") for i in interventions.get("interventions", [])
                ),
                "phase": "|".join((design.get("phases") or [])),
                "overall_status": status_mod.get("overallStatus", ""),
                "start_date": (status_mod.get("startDateStruct") or {}).get("date", ""),
                "primary_completion_date": (status_mod.get("primaryCompletionDateStruct") or {}).get("date", ""),
                "enrollment": (design.get("enrollmentInfo") or {}).get("count", ""),
                "study_type": design.get("studyType", ""),
                "lead_sponsor": (sponsor.get("leadSponsor") or {}).get("name", ""),
                "disease_query": disease,
            }
            all_studies.append(row)

        page_token = data.get("nextPageToken")
        if not page_token:
            break
        time.sleep(0.5)  # Rate limiting

    return all_studies


def refresh_all_diseases(max_pages_per_disease: int = 3) -> Path:
    """
    Fetch trials for all PRV-eligible diseases and write to CSV + metadata.

    Returns:
        Path to the output CSV file
    """
    print(f"[clinicaltrials] Refreshing trials for {len(PRV_DISEASE_QUERIES)} diseases...")
    all_trials: List[Dict[str, Any]] = []

    for disease in PRV_DISEASE_QUERIES:
        try:
            trials = fetch_disease_trials(disease, max_pages=max_pages_per_disease)
            print(f"  {disease}: {len(trials)} trials")
            all_trials.extend(trials)
        except Exception as e:
            print(f"  {disease}: FAILED ({e})", file=sys.stderr)

    # Deduplicate by NCT ID
    seen: set = set()
    unique: List[Dict[str, Any]] = []
    for t in all_trials:
        nct = t.get("nct_id", "")
        if nct and nct not in seen:
            seen.add(nct)
            unique.append(t)

    # Write CSV
    csv_path = DATA_DIR / "clinicaltrials_interventions.csv"
    if unique:
        fieldnames = list(unique[0].keys())
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(unique)

    # Write metadata
    meta = {
        "last_refresh_utc": datetime.now(timezone.utc).isoformat(),
        "total_trials": len(unique),
        "diseases_queried": len(PRV_DISEASE_QUERIES),
        "source": "ClinicalTrials.gov v2 API",
        "checksum": hashlib.sha256(csv_path.read_bytes()).hexdigest()[:16] if csv_path.exists() else "",
    }
    meta_path = DATA_DIR / "metadata.json"
    meta_path.write_text(json.dumps(meta, indent=2), encoding="utf-8")

    print(f"[clinicaltrials] Done: {len(unique)} unique trials -> {csv_path}")
    return csv_path


def search_trials(compound_name: str, disease: Optional[str] = None) -> List[Dict[str, Any]]:
    """
    Search the local CSV for trials involving a compound.

    Args:
        compound_name: Drug/compound name to search for
        disease: Optional disease filter

    Returns:
        List of matching trial rows
    """
    csv_path = DATA_DIR / "clinicaltrials_interventions.csv"
    if not csv_path.exists():
        return []

    matches: List[Dict[str, Any]] = []
    query_lower = compound_name.lower()

    with open(csv_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            interventions = (row.get("interventions") or "").lower()
            title = (row.get("brief_title") or "").lower()
            if query_lower in interventions or query_lower in title:
                if disease:
                    conditions = (row.get("conditions") or "").lower()
                    disease_q = (row.get("disease_query") or "").lower()
                    if disease.lower() not in conditions and disease.lower() not in disease_q:
                        continue
                matches.append(dict(row))

    return matches


# ── CLI Entry Point ────────────────────────────────────────────────────

def main() -> int:
    import argparse
    parser = argparse.ArgumentParser(description="ClinicalTrials.gov API client")
    sub = parser.add_subparsers(dest="command")

    refresh_p = sub.add_parser("refresh", help="Refresh trial data for all PRV diseases")
    refresh_p.add_argument("--max-pages", type=int, default=3, help="Max pages per disease")

    search_p = sub.add_parser("search", help="Search local trials by compound name")
    search_p.add_argument("compound", help="Compound name to search for")
    search_p.add_argument("--disease", help="Filter by disease")

    args = parser.parse_args()
    if args.command == "refresh":
        refresh_all_diseases(max_pages_per_disease=args.max_pages)
        return 0
    elif args.command == "search":
        results = search_trials(args.compound, disease=args.disease)
        print(f"Found {len(results)} matching trials")
        for r in results[:20]:
            print(f"  {r.get('nct_id')}: {r.get('brief_title', '')[:80]}")
        return 0
    else:
        parser.print_help()
        return 0


if __name__ == "__main__":
    sys.exit(main())
