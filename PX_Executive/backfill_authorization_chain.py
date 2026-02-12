#!/usr/bin/env python3
"""
backfill_authorization_chain.py — Re-run legacy dossiers through the 12-engine pipeline.

Legacy dossiers (pre-v3.0.0) only have 4 engines (OPE, ADMET, OBE, OLE) and lack
the authorization_chain required by Finalization_Spec v3.0.0.

This script:
  1. Scans Novel_Dossiers and Prv_Dossiers for dossiers missing authorization_chain
  2. Extracts SMILES and metadata from each
  3. Re-runs through the full 12-engine pipeline (px_prv.run_full_pipeline)
  4. Saves updated dossiers (with authorization_chain) back to warehouse

After running this, use `px_finalize.py --reprocess` to finalize them.

Usage:
    python PX_Executive/backfill_authorization_chain.py
    python PX_Executive/backfill_authorization_chain.py --dry-run
    python PX_Executive/backfill_authorization_chain.py --limit 5
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def main() -> int:
    parser = argparse.ArgumentParser(description="Backfill authorization_chain for legacy dossiers")
    parser.add_argument("--dry-run", action="store_true", help="Report only, no writes")
    parser.add_argument("--limit", type=int, default=0, help="Max dossiers to process (0 = all)")
    parser.add_argument("--finalized-only", action="store_true",
                        help="Only backfill dossiers that have corresponding finalized versions")
    parser.add_argument("--tiers", nargs="+", default=None,
                        help="Only backfill specific tiers (e.g. Diamond Gold)")
    args = parser.parse_args()

    from PX_Warehouse.warehouse_layout import TIERS

    wh = REPO_ROOT / "PX_Warehouse"

    # If --finalized-only, collect stems of existing finalized dossiers
    finalized_stems: set[str] | None = None
    if args.finalized_only:
        finalized_stems = set()
        finalized_dir = wh / "Finalized_Dossiers"
        if finalized_dir.exists():
            for tier in TIERS:
                td = finalized_dir / tier
                if td.exists():
                    for f in td.glob("*.json"):
                        if not f.name.startswith("PATH_TEST"):
                            finalized_stems.add(f.stem)

    # Tier filter
    tiers_filter = set(args.tiers) if args.tiers else None

    # Collect dossiers missing authorization_chain
    to_backfill: list[tuple[Path, str, bool]] = []
    for is_novel, folder in [(True, "Novel_Dossiers"), (False, "Prv_Dossiers")]:
        for tier in TIERS:
            if tiers_filter and tier not in tiers_filter:
                continue
            d_dir = wh / folder / tier
            if not d_dir.exists():
                continue
            for f in d_dir.glob("*.json"):
                if f.name.startswith("PATH_TEST"):
                    continue
                if finalized_stems is not None and f.stem not in finalized_stems:
                    continue
                try:
                    dossier = json.loads(f.read_text(encoding="utf-8"))
                    if "authorization_chain" not in dossier:
                        to_backfill.append((f, f.stem, is_novel))
                except Exception:
                    continue

    if args.limit > 0:
        to_backfill = to_backfill[:args.limit]

    print("=" * 64)
    print("   PREDATOR X — AUTHORIZATION CHAIN BACKFILL")
    print("=" * 64)
    print(f"Dossiers needing backfill: {len(to_backfill)} | Dry-run: {args.dry_run}")
    print("-" * 64)

    if not to_backfill:
        print("All dossiers already have authorization_chain. Nothing to do.")
        return 0

    from PX_Executive.px_prv import run_full_pipeline, save_dossier

    ok = 0
    fail = 0

    for path, item_id, is_novel in to_backfill:
        try:
            old_dossier = json.loads(path.read_text(encoding="utf-8"))
        except Exception as e:
            print(f"  Skip {path.name}: read error {e}")
            fail += 1
            continue

        # Extract candidate info from existing dossier
        candidate = old_dossier.get("candidate", {})
        smiles = candidate.get("smiles", "")
        name = candidate.get("name", item_id)

        if not smiles:
            print(f"  Skip {item_id}: no SMILES in dossier")
            fail += 1
            continue

        # Determine indication from existing engine results
        indication = ""
        engines = old_dossier.get("engines", {})
        ole = engines.get("ole", {})
        if isinstance(ole, dict):
            indication = ole.get("indication", "")

        if args.dry_run:
            print(f"  [dry-run] Would re-run {item_id} (smiles={smiles[:40]}...)")
            ok += 1
            continue

        # Build item dict matching what px_prv expects
        item = {
            "id": item_id,
            "smiles": smiles,
            "name": name,
            "type": "N" if is_novel else "R",
            "indication": indication,
        }

        print(f"  Re-running {item_id} through 12-engine pipeline...", end="", flush=True)
        try:
            new_dossier, err = run_full_pipeline(item)
        except Exception as e:
            print(f" ERROR: {e}")
            fail += 1
            continue

        if err:
            print(f" FAILED: {err}")
            fail += 1
            continue

        # Verify authorization_chain is present
        auth = new_dossier.get("authorization_chain", {})
        auth_count = auth.get("authorization_count", "0/0")
        grade = (new_dossier.get("discovery_grading") or {}).get("grade", "?")

        # Save updated dossier (may land at different tier than original)
        out_path = save_dossier(new_dossier, item)

        # Remove orphaned original if tier changed (old copy has no auth_chain)
        if out_path and Path(out_path).resolve() != path.resolve():
            try:
                path.unlink()
            except OSError:
                pass

        print(f" OK -> {out_path} (auth={auth_count}, grade={grade})")
        ok += 1

    print("\n" + "=" * 64)
    print("BACKFILL SUMMARY")
    print("=" * 64)
    print(f"Updated:     {ok}")
    print(f"Failed/Skip: {fail}")
    print("=" * 64)

    if ok > 0 and not args.dry_run:
        print("\nNext step: python PX_Executive/px_finalize.py --reprocess -v")

    return 0 if fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
