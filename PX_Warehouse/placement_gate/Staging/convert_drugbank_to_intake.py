#!/usr/bin/env python3
"""
Convert a DrugBank CSV export to the intake harvester format.
Usage: python convert_drugbank_to_intake.py [path/to/drugbank_export.csv]
Output: PX_Data/drugbank/drugbank.csv (or --output path).
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = REPO_ROOT / "PX_Data" / "drugbank" / "drugbank.csv"


def normalize_status(groups: str) -> str:
    """Map DrugBank 'groups' (e.g. 'approved; investigational') to approved/investigational."""
    g = (groups or "").lower()
    if "approved" in g:
        return "approved"
    if "investigational" in g or "experimental" in g:
        return "investigational"
    return "other"


def main() -> int:
    ap = argparse.ArgumentParser(description="Convert DrugBank CSV to intake format")
    ap.add_argument("input", type=Path, nargs="?", help="DrugBank export CSV (default: stdin)")
    ap.add_argument("--output", type=Path, default=DEFAULT_OUT, help="Output CSV path")
    args = ap.parse_args()
    if args.input and not args.input.exists():
        print(f"File not found: {args.input}", file=sys.stderr)
        return 1
    fin = args.input.open("r", encoding="utf-8", errors="replace") if args.input else sys.stdin
    reader = csv.DictReader(fin)
    rows = list(reader)
    if args.input:
        fin.close()
    if not rows:
        print("No rows in input.", file=sys.stderr)
        return 1
    # Map common DrugBank column names to our schema
    out_rows = []
    for r in rows:
        dbid = r.get("drugbank_id") or r.get("drugbank-id") or r.get("id") or ""
        name = r.get("name") or r.get("title") or ""
        smiles = r.get("smiles") or r.get("canonical_smiles") or r.get("canonical-smiles") or ""
        groups = r.get("groups") or r.get("state") or ""
        status = normalize_status(groups)
        atc = r.get("atc_codes") or r.get("atc-codes") or ""
        targets = r.get("targets") or r.get("polypeptides") or ""
        if not dbid or not smiles:
            continue
        out_rows.append({
            "drugbank_id": dbid,
            "name": name,
            "status": status,
            "smiles": smiles,
            "atc_codes": atc,
            "targets": targets,
        })
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["drugbank_id", "name", "status", "smiles", "atc_codes", "targets"])
        w.writeheader()
        w.writerows(out_rows)
    print(f"Wrote {len(out_rows)} rows to {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
