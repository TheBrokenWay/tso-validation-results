#!/usr/bin/env python3
"""
Fetch intervention records from ClinicalTrials.gov API v2 (all statuses by default).
Writes PX_Data/clinicaltrials/clinicaltrials_interventions.csv. Default target 100k rows for ~24h intake.
Rate-friendly pagination with 0.25s delay between pages.
"""
from __future__ import annotations

import argparse
import csv
import time
import sys
from pathlib import Path

try:
    import requests
except ImportError:
    print("pip install requests", file=sys.stderr)
    sys.exit(1)

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = REPO_ROOT / "PX_Data" / "clinicaltrials" / "clinicaltrials_interventions.csv"
BASE = "https://clinicaltrials.gov/api/v2/studies"
PAGE_SIZE = 100
PAGE_DELAY = 0.25  # seconds between API pages


def main() -> int:
    ap = argparse.ArgumentParser(description="Fetch ClinicalTrials.gov interventions to CSV")
    ap.add_argument("--output", type=Path, default=DEFAULT_OUT, help="Output CSV path")
    ap.add_argument("--limit", type=int, default=100000, help="Max intervention rows (default 100k for ~24h intake)")
    ap.add_argument("--status", type=str, default="", help="Filter by overallStatus (comma-separated); empty = all statuses")
    args = ap.parse_args()
    rows = []
    next_token = None
    statuses = [s.strip() for s in (args.status or "").split(",") if s.strip()]
    while len(rows) < args.limit:
        params = {"pageSize": PAGE_SIZE}
        if next_token:
            params["pageToken"] = next_token
        if statuses:
            params["filter.overallStatus"] = statuses[0]
        try:
            r = requests.get(BASE, params=params, timeout=30)
            r.raise_for_status()
            data = r.json()
        except Exception as e:
            print(f"API error: {e}", file=sys.stderr)
            break
        studies = data.get("studies") or []
        next_token = data.get("nextPageToken")
        for s in studies:
            if len(rows) >= args.limit:
                break
            protocol = s.get("protocolSection") or {}
            ident = protocol.get("identificationModule") or {}
            nct_id = ident.get("nctId") or ""
            status_mod = protocol.get("statusModule") or {}
            overall = status_mod.get("overallStatus") or ""
            cond_mod = protocol.get("conditionsModule") or {}
            conditions = cond_mod.get("conditions") or []
            condition = "; ".join(conditions)[:200] if conditions else ""
            sponsor_mod = protocol.get("sponsorCollaboratorsModule") or {}
            lead_sponsor = sponsor_mod.get("leadSponsor") or {}
            sponsor = (lead_sponsor.get("name") or "")[:100]
            design = protocol.get("designModule") or {}
            phases = design.get("phases") or []
            phase = "; ".join(phases).lower() if phases else "na"
            arms = protocol.get("armsInterventionsModule") or {}
            interventions = arms.get("interventions") or []
            for intr in interventions:
                if len(rows) >= args.limit:
                    break
                name = (intr.get("name") or "").strip()
                if not name:
                    continue
                itype = (intr.get("type") or "").upper()
                if itype not in ("DRUG", "BIOLOGICAL", "DEVICE", "RADIATION", "BEHAVIORAL", "PROCEDURE", "OTHER"):
                    itype = "OTHER"
                rows.append({
                    "nct_id": nct_id,
                    "intervention_name": name,
                    "intervention_type": itype,
                    "condition": condition,
                    "phase": phase,
                    "overall_status": overall,
                    "sponsor": sponsor,
                    "smiles": "",
                })
        if not next_token or not studies:
            break
        time.sleep(PAGE_DELAY)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["nct_id", "intervention_name", "intervention_type", "condition", "phase", "overall_status", "sponsor", "smiles"])
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {len(rows)} intervention rows to {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
