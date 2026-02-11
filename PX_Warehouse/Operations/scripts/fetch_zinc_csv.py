#!/usr/bin/env python3
"""
Fetch commercially purchasable molecules from ZINC (ZINC20 2D).
Downloads 2D tranche .smi files from files.docking.org. Default target 200k molecules for ~24h intake.
"""
from __future__ import annotations

import argparse
import csv
import sys
import time
from pathlib import Path

try:
    import requests
except ImportError:
    print("pip install requests", file=sys.stderr)
    sys.exit(1)

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = REPO_ROOT / "PX_Data" / "zinc" / "zinc_molecules.csv"
BASE_2D = "https://files.docking.org/2D"
# ZINC20 2D: 121 first-level tranches (AA..KK); within each, 4-letter .smi files (e.g. AAAA, AAAB, ...)
FIRST_LETTERS = "ABCDEFGHIJK"
SECOND_LETTERS = "ABCDEFGHIJK"
REACTIVITY = "ABCDEF"   # anodyne..annotated
PURCHASABILITY = "AB"   # in stock (A/B) for repurposing-friendly set; expand to "ABCDEF" for more
# Build many tranche URLs: 2D/{T2}/{T2}{R}{P}.smi
def _build_zinc_2d_urls(max_urls: int = 400) -> list[str]:
    urls: list[str] = []
    for c1 in FIRST_LETTERS:
        for c2 in SECOND_LETTERS:
            t2 = c1 + c2
            for r in REACTIVITY:
                for p in PURCHASABILITY:
                    if len(urls) >= max_urls:
                        return urls
                    four = t2 + r + p
                    urls.append(f"{BASE_2D}/{t2}/{four}.smi")
    return urls

# Up to 600 URLs: in-stock (A/B) across reactivity and many tranches for 200k+ molecules
TRANCH_URLS = _build_zinc_2d_urls(600)
HEADERS = {"User-Agent": "PredatorX-Intake/1.0 (research; rate-limited)"}
ZINC_DELAY = 0.15  # seconds between tranche requests


def _parse_smi_text(text: str, limit: int, rows: list) -> None:
    for line in text.splitlines():
        if len(rows) >= limit:
            return
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t") if "\t" in line else line.split()
        if len(parts) >= 1 and parts[0]:
            smiles = parts[0]
            if len(smiles) < 3:
                continue
            name_or_id = parts[1] if len(parts) > 1 else ""
            zinc_id = parts[2] if len(parts) > 2 else name_or_id or f"ZINC{len(rows)}"
            rows.append({
                "zinc_id": zinc_id,
                "name": name_or_id if name_or_id != zinc_id else "",
                "smiles": smiles,
                "catalog": "ZINC20",
            })


def main() -> int:
    ap = argparse.ArgumentParser(description="Fetch ZINC molecules to CSV")
    ap.add_argument("--output", type=Path, default=DEFAULT_OUT, help="Output CSV path")
    ap.add_argument("--limit", type=int, default=200000, help="Max molecules (default 200k for ~24h intake)")
    ap.add_argument("--url", type=str, default="", help="Override: single ZINC SMI URL (ignores multi-tranche)")
    args = ap.parse_args()
    rows = []
    if args.url:
        urls = [args.url]
    else:
        urls = TRANCH_URLS
    print(f"Fetching from ZINC (target {args.limit} molecules)...")
    for url in urls:
        if len(rows) >= args.limit:
            break
        time.sleep(ZINC_DELAY)
        try:
            r = requests.get(url, timeout=120, stream=True, headers=HEADERS)
            r.raise_for_status()
            content = r.content
            if url.endswith(".gz"):
                import gzip
                content = gzip.decompress(content)
            text = content.decode("utf-8", errors="replace")
            _parse_smi_text(text, args.limit, rows)
        except Exception as e:
            print(f"Skip {url}: {e}", file=sys.stderr)
    if not rows:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("w", newline="", encoding="utf-8") as f:
            f.write("zinc_id,name,smiles,catalog\n")
        print("Wrote header-only CSV. Set --url to a valid ZINC SMI URL.", file=sys.stderr)
        return 0
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["zinc_id", "name", "smiles", "catalog"])
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {len(rows)} molecules to {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
