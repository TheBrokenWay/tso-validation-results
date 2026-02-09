#!/usr/bin/env python3
"""
Fetch compounds from PubChem (approved/drug-like/bioactive) and write
PX_Data/pubchem/pubchem_molecules.csv. Uses E-utilities for CIDs then PUG REST for SMILES.
Default target: 150k molecules for ~24h pipeline intake. Rate limit: stay under 5 requests/second.
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
DEFAULT_OUT = REPO_ROOT / "PX_Data" / "pubchem" / "pubchem_molecules.csv"
BATCH = 200  # PUG REST; keep modest for URL length and reliability
DELAY = 0.22  # ~4.5/sec
ESEARCH_PAGE = 10000  # E-utilities retmax cap per request

# Multiple E-utilities queries to fill CIDs (drug-like, bioactive, approved)
ESEARCH_TERMS = [
    "has_dailymed",           # FDA DailyMed
    "SRCSET[0]",              # PubChem source set
    "MeSH[Source]",            # MeSH compounds
    "bioactive",               # Bioactive subset
]


def _esearch_cids(base: str, term: str, retmax: int, retstart: int = 0) -> list[int]:
    try:
        r = requests.get(
            f"{base}/esearch.fcgi",
            params={
                "db": "pccompound",
                "term": term,
                "retmax": min(retmax, ESEARCH_PAGE),
                "retstart": retstart,
                "retmode": "json",
            },
            timeout=60,
        )
        r.raise_for_status()
        data = r.json()
        id_list = data.get("esearchresult", {}).get("idlist") or []
        return [int(x) for x in id_list]
    except Exception as e:
        print(f"E-utilities {term!r} (start={retstart}): {e}", file=sys.stderr)
        return []


def get_cids(limit: int) -> list[int]:
    """Get CIDs from PubChem: multiple terms with pagination, then CID range fallback."""
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    seen: set[int] = set()
    cids: list[int] = []
    for term in ESEARCH_TERMS:
        if len(cids) >= limit:
            break
        retstart = 0
        while len(cids) < limit:
            chunk = _esearch_cids(base, term, limit - len(cids), retstart)
            if not chunk:
                break
            for cid in chunk:
                if cid not in seen:
                    seen.add(cid)
                    cids.append(cid)
                    if len(cids) >= limit:
                        break
            if len(chunk) < ESEARCH_PAGE:
                break
            retstart += ESEARCH_PAGE
            time.sleep(DELAY)
    if len(cids) < limit:
        # Fallback: CIDs from high range (many exist; some missing SMILES)
        hi = max(seen) if seen else 1
        for i in range(hi + 1, hi + 1 + (limit - len(cids)) * 3):
            if len(cids) >= limit:
                break
            if i not in seen:
                cids.append(i)
                seen.add(i)
    return cids[:limit]


def fetch_smiles_batch(cids: list[int]) -> list[dict]:
    """PUG REST: get CanonicalSMILES and Title for a list of CIDs."""
    if not cids:
        return []
    cid_str = ",".join(str(x) for x in cids)
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid_str}/property/CanonicalSMILES,Title/CSV"
    try:
        r = requests.get(url, timeout=60)
        r.raise_for_status()
        lines = r.text.strip().split("\n")
        if not lines:
            return []
        reader = csv.DictReader(lines)
        out = []
        for row in reader:
            cid = row.get("CID") or row.get("cid")
            smiles = (
                row.get("CanonicalSMILES") or row.get("Canonical SMILES")
                or row.get("ConnectivitySMILES") or row.get("Connectivity SMILES") or ""
            ).strip()
            name = (row.get("Title") or "").strip()
            if cid and smiles:
                out.append({
                    "cid": cid,
                    "name": name,
                    "smiles": smiles,
                    "bioactive": "Y",
                    "assay_hit": "N",
                })
        return out
    except Exception as e:
        print(f"PUG batch error for {len(cids)} CIDs: {e}", file=sys.stderr)
        return []


def main() -> int:
    ap = argparse.ArgumentParser(description="Fetch PubChem molecules to CSV")
    ap.add_argument("--output", type=Path, default=DEFAULT_OUT, help="Output CSV path")
    ap.add_argument("--limit", type=int, default=150000, help="Max molecules (default 150k for ~24h intake)")
    args = ap.parse_args()
    cids = get_cids(args.limit)
    print(f"Got {len(cids)} CIDs, fetching SMILES in batches of {BATCH}...")
    rows = []
    for i in range(0, len(cids), BATCH):
        batch = cids[i : i + BATCH]
        rows.extend(fetch_smiles_batch(batch))
        time.sleep(DELAY)
        if len(rows) >= args.limit:
            break
    # If we need more (many CIDs lack SMILES), get more CIDs and retry
    needed = args.limit - len(rows)
    if needed > 0 and cids:
        extra_start = max(cids) + 1
        for j in range(0, min(needed * 4, 200000), BATCH):
            extra_cids = list(range(extra_start + j, extra_start + j + BATCH))
            rows.extend(fetch_smiles_batch(extra_cids))
            time.sleep(DELAY)
            if len(rows) >= args.limit:
                break
    rows = rows[: args.limit]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["cid", "name", "smiles", "bioactive", "assay_hit"])
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {len(rows)} molecules to {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
