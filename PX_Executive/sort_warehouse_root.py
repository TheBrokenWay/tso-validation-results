#!/usr/bin/env python3
"""
Sort all loose dossier files from PX_Warehouse root into sanctioned tier directories.

Routes:
  PRV_NOV_*.json  → Novel_Dossiers/<Diamond|Gold|Silver|Bronze>/
  PRV_REP_*.json  → Prv_Dossiers/<Diamond|Gold|Silver|Bronze>/
  TRIAL_SIMULATION_DOSSIER*.json → Calibration_Molecules/

Uses get_tier() from warehouse_layout for deterministic placement based on toxicity.
After placement, removes the original from warehouse root.

Run from repo root:
  python PX_Executive/sort_warehouse_root.py --dry-run     # Report only
  python PX_Executive/sort_warehouse_root.py               # Execute sort
  python PX_Executive/sort_warehouse_root.py --limit 100   # Process max 100
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from PX_System.foundation.ZeusLaws import run_zeus_gate
from PX_Warehouse.warehouse_layout import (
    ensure_structure,
    get_tier,
    get_prv_dossier_dir,
    get_calibration_molecules_dir,
)


def sort_root_files(dry_run: bool = False, limit: int = 0) -> dict:
    """Sort all loose dossier files from warehouse root into tier directories."""
    ensure_structure(REPO_ROOT)
    wh = REPO_ROOT / "PX_Warehouse"

    stats = {
        "total": 0,
        "placed": 0,
        "skipped": 0,
        "rejected": 0,
        "errors": 0,
        "tiers": {"Diamond": 0, "Gold": 0, "Silver": 0, "Bronze": 0},
        "calibration": 0,
        "by_type": {"PRV_NOV": 0, "PRV_REP": 0, "TRIAL_SIM": 0},
    }

    # Collect all root dossier files
    to_sort: list[tuple[Path, str]] = []  # (path, category)
    for f in sorted(wh.glob("PRV_NOV_*.json")):
        to_sort.append((f, "PRV_NOV"))
    for f in sorted(wh.glob("PRV_REP_*.json")):
        to_sort.append((f, "PRV_REP"))
    for f in sorted(wh.glob("TRIAL_SIMULATION_DOSSIER*.json")):
        to_sort.append((f, "TRIAL_SIM"))

    stats["total"] = len(to_sort)
    if limit > 0:
        to_sort = to_sort[:limit]

    print(f"Sorting {len(to_sort)} root files (dry_run={dry_run})")
    if stats["total"] > len(to_sort):
        print(f"  (limited from {stats['total']} total)")

    for idx, (src, category) in enumerate(to_sort, start=1):
        if idx % 500 == 0:
            print(f"  Progress: {idx}/{len(to_sort)} (placed={stats['placed']})")

        try:
            if category == "TRIAL_SIM":
                dest_dir = get_calibration_molecules_dir(REPO_ROOT)
                dest_dir.mkdir(parents=True, exist_ok=True)
                dest = dest_dir / src.name
                if not dry_run:
                    shutil.move(str(src), str(dest))
                stats["calibration"] += 1
                stats["placed"] += 1
                stats["by_type"]["TRIAL_SIM"] += 1
                continue

            dossier = json.loads(src.read_text(encoding="utf-8"))

            # Zeus gate: constitutional governance check before tier placement
            zeus_verdict = run_zeus_gate(dossier)
            if not zeus_verdict.get("authorized", False):
                rationale = zeus_verdict.get("rationale", "governance failure")
                print(f"  ZEUS REJECTED: {src.name}: {rationale}")
                try:
                    from PX_System.foundation.Sovereign_Log_Chain import append as slc_append
                    slc_append("ZEUS_GATE_REJECTION", {
                        "file": src.name,
                        "rationale": rationale,
                        "source": "sort_warehouse_root",
                    }, context="warehouse_sort")
                except Exception:
                    pass  # Logging failure must not block rejection
                stats["rejected"] += 1
                continue

            tier = get_tier(dossier)
            is_novel = category == "PRV_NOV"
            dest_dir = get_prv_dossier_dir(is_novel, tier, REPO_ROOT)
            dest_dir.mkdir(parents=True, exist_ok=True)
            dest = dest_dir / src.name

            if dest.exists():
                # Already placed — just remove the root copy
                if not dry_run:
                    src.unlink()
                stats["skipped"] += 1
                continue

            if not dry_run:
                shutil.move(str(src), str(dest))

            stats["placed"] += 1
            stats["tiers"][tier] += 1
            stats["by_type"][category] += 1

        except Exception as e:
            print(f"  ERROR: {src.name}: {e}")
            stats["errors"] += 1

    return stats


def main() -> int:
    parser = argparse.ArgumentParser(description="Sort warehouse root dossiers into tier directories.")
    parser.add_argument("--dry-run", action="store_true", help="Report what would be sorted; no writes.")
    parser.add_argument("--limit", type=int, default=0, help="Max files to process (0 = all).")
    args = parser.parse_args()

    stats = sort_root_files(dry_run=args.dry_run, limit=args.limit)

    print(f"\n{'='*60}")
    print("WAREHOUSE ROOT SORT COMPLETE")
    print(f"{'='*60}")
    print(f"Total scanned:   {stats['total']}")
    print(f"Placed:          {stats['placed']}")
    print(f"Zeus rejected:   {stats['rejected']}")
    print(f"Already present: {stats['skipped']}")
    print(f"Errors:          {stats['errors']}")
    print(f"\nBy tier:")
    for tier, count in stats["tiers"].items():
        print(f"  {tier}: {count}")
    print(f"  Calibration_Molecules: {stats['calibration']}")
    print(f"\nBy type:")
    for t, count in stats["by_type"].items():
        print(f"  {t}: {count}")
    print(f"{'='*60}")

    return 0 if stats["errors"] == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
