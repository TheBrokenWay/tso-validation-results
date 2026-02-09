#!/usr/bin/env python3
"""
Fetch real molecule data from ChEMBL REST API and write
PX_Data/chembl/chembl_molecules.csv (or path from --output).

Requires: requests
Usage: python fetch_chembl_csv.py [--output path] [--limit N] [--phase 0|1|2|3|4]  # 0=all phases (for 10k)
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

try:
    import requests
except ImportError:
    print("pip install requests", file=sys.stderr)
    sys.exit(1)

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = REPO_ROOT / "PX_Data" / "chembl" / "chembl_molecules.csv"
BASE_URL = "https://www.ebi.ac.uk/chembl/api/data/molecule.json"


def fetch_molecules(max_phase: int | None = 4, limit: int = 2000) -> list[dict]:
    """Fetch molecules from ChEMBL API. max_phase 4 = approved, 1-3 = investigational, None = all phases (for 10k)."""
    out = []
    offset = 0
    page_size = min(500, limit)
    while len(out) < limit:
        params = {"limit": page_size, "offset": offset}
        if max_phase is not None and max_phase >= 1:
            params["max_phase"] = max_phase
        try:
            r = requests.get(BASE_URL, params=params, timeout=60)
            r.raise_for_status()
            data = r.json()
        except Exception as e:
            print(f"ChEMBL API error: {e}", file=sys.stderr)
            break
        molecules = data.get("molecules") or data.get("molecule") or []
        if not molecules:
            break
        for m in molecules:
            if len(out) >= limit:
                break
            chembl_id = m.get("molecule_chembl_id") or m.get("chembl_id")
            pref_name = m.get("pref_name") or ""
            structures = m.get("molecule_structures") or {}
            if isinstance(structures, list) and structures:
                structures = structures[0] if isinstance(structures[0], dict) else {}
            canonical_smiles = (structures.get("canonical_smiles") or "").strip()
            if not chembl_id or not canonical_smiles:
                continue
            phase = m.get("max_phase", 0)
            status = "approved" if phase == 4 else "investigational"
            out.append({
                "chembl_id": chembl_id,
                "pref_name": pref_name,
                "canonical_smiles": canonical_smiles,
                "status": status,
                "targets": "",
                "mechanism_of_action": "",
            })
        offset += len(molecules)
        if len(molecules) < page_size:
            break
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description="Fetch ChEMBL molecules to CSV")
    ap.add_argument("--output", type=Path, default=DEFAULT_OUT, help="Output CSV path")
    ap.add_argument("--limit", type=int, default=10000, help="Max molecules to fetch")
    ap.add_argument("--phase", type=int, default=0, help="Max phase: 1-4 or 0 for all (to get 10k)")
    args = ap.parse_args()
    phase_arg = None if args.phase == 0 else args.phase
    if phase_arg is not None and phase_arg not in (1, 2, 3, 4):
        phase_arg = 4
    rows = fetch_molecules(max_phase=phase_arg, limit=args.limit)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["chembl_id", "pref_name", "canonical_smiles", "status", "targets", "mechanism_of_action"])
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {len(rows)} molecules to {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
