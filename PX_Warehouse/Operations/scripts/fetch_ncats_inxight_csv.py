#!/usr/bin/env python3
"""
Fetch NCATS Inxight Drugs (FDA + NIH curated) and write PX_Data/ncats_inxight/ncats_inxight.csv.
Downloads FRDB zip from drugs.ncats.io, extracts frdb-drugs.tsv. ~4.5k in current FRDB; default limit 50k for future bulk.
Optional: merge multiple FRDB versions via --extra-urls for larger intake.
"""
from __future__ import annotations

import argparse
import csv
import io
import sys
import zipfile
from pathlib import Path

try:
    import requests
except ImportError:
    print("pip install requests", file=sys.stderr)
    sys.exit(1)

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = REPO_ROOT / "PX_Data" / "ncats_inxight" / "ncats_inxight.csv"
ZIP_URL = "https://drugs.ncats.io/downloads-public/frdb-v2024-12-30.zip"
# Optional extra zips to merge (e.g. older FRDB for more rows)
EXTRA_ZIP_URLS = [
    "https://drugs.ncats.io/downloads-public/frdb-v2023-07-05.zip",
]


def _parse_frdb_tsv(content: bytes, limit: int, seen_ids: set[str], rows: list[dict]) -> None:
    text = content.decode("utf-8", errors="replace")
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    for row in reader:
        if len(rows) >= limit:
            return
        ncats_id = row.get("compound_unii") or row.get("unii") or row.get("compound_id") or row.get("id") or ""
        name_val = row.get("compound_name") or row.get("name") or row.get("preferred_name") or row.get("title") or ""
        smiles = row.get("smiles") or row.get("canonical_smiles") or row.get("structure") or ""
        status = row.get("status") or row.get("development_status") or "approved"
        atc = row.get("atc") or row.get("atc_codes") or ""
        mechanism = row.get("mechanism") or row.get("mechanism_of_action") or ""
        if not ncats_id and not name_val:
            continue
        key = (ncats_id or "").strip() or (name_val or "").strip()
        if key and key in seen_ids:
            continue
        if key:
            seen_ids.add(key)
        rows.append({
            "ncats_id": ncats_id,
            "name": name_val,
            "smiles": smiles,
            "status": status,
            "atc": atc,
            "mechanism": mechanism,
        })


def main() -> int:
    ap = argparse.ArgumentParser(description="Fetch NCATS Inxight Drugs to CSV")
    ap.add_argument("--output", type=Path, default=DEFAULT_OUT, help="Output CSV path")
    ap.add_argument("--limit", type=int, default=50000, help="Max rows (default 50k; FRDB has ~4.5k)")
    ap.add_argument("--extra", action="store_true", help="Also fetch --extra-urls zips and merge (dedupe by id/name)")
    args = ap.parse_args()
    seen: set[str] = set()
    rows: list[dict] = []
    urls = [ZIP_URL]
    if args.extra:
        urls = [ZIP_URL] + EXTRA_ZIP_URLS
    for url in urls:
        if len(rows) >= args.limit:
            break
        print(f"Downloading {url}...")
        r = requests.get(url, timeout=120)
        r.raise_for_status()
        with zipfile.ZipFile(io.BytesIO(r.content), "r") as z:
            names = [n for n in z.namelist() if "frdb-drugs" in n and n.endswith(".tsv")]
            if not names:
                names = [n for n in z.namelist() if n.endswith(".tsv")]
            for name in names[:1]:
                with z.open(name) as f:
                    _parse_frdb_tsv(f.read(), args.limit, seen, rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["ncats_id", "name", "smiles", "status", "atc", "mechanism"])
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {len(rows)} rows to {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
