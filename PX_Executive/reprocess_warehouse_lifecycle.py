#!/usr/bin/env python3
"""
One-time batch lifecycle run: push all unfinalized Prv/Novel dossiers through
grading + trial binding + Zeus gate, then two-tier placement.

Two-tier bar (re-ordered grading):
  1. Zeus authorized (tox/harm < 0.021) → Finalized_Dossiers/<tier> (Diamond/Gold/Silver/Bronze).
  2. Else discovery bar (tox/harm < 0.05) → Finalized_Dossiers/Discovery_Accepted.
  3. Else rejected (not written).

Pipeline (inside Finalization_Pipeline.run_finalization):
  get_tier(dossier) → GradingEngine (GradingSchema_Discovery.json) →
  SoC/novelty/SA metrics → WorldLine_Database.record_discovery →
  ZeusLaws.run_zeus_gate → finalize_and_place (Zeus → tier; else discovery bar → Discovery_Accepted).

Constraints:
  - Does not modify or delete raw dossier data in Prv_Dossiers/Novel_Dossiers.
  - Uses only existing engines and PX_Warehouse APIs (warehouse_layout, WorldLine_Database).
  - After batch, you can run run_all_tests.py and run_e2e_layers.py to confirm.

Run from repo root:
  python PX_Executive/reprocess_warehouse_lifecycle.py
  python PX_Executive/reprocess_warehouse_lifecycle.py --dry-run
  python PX_Executive/reprocess_warehouse_lifecycle.py --limit 100
  python PX_Executive/reprocess_warehouse_lifecycle.py --validate-after
  python PX_Executive/reprocess_warehouse_lifecycle.py --validate-before
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def run_validation() -> tuple[bool, str]:
    """Run run_all_tests.py and run_e2e_layers.py. Return (both_ok, message)."""
    try:
        r1 = subprocess.run(
            [sys.executable, "PX_Validation/tests/run_all_tests.py"],
            cwd=str(REPO_ROOT),
            capture_output=True,
            text=True,
            timeout=300,
        )
        if r1.returncode != 0:
            return False, f"run_all_tests.py failed (exit {r1.returncode})"
        r2 = subprocess.run(
            [sys.executable, "run_e2e_layers.py"],
            cwd=str(REPO_ROOT),
            capture_output=True,
            text=True,
            timeout=120,
        )
        if r2.returncode != 0:
            return False, f"run_e2e_layers.py failed (exit {r2.returncode})"
        return True, "run_all_tests.py and run_e2e_layers.py both passed"
    except subprocess.TimeoutExpired:
        return False, "Validation timed out"
    except Exception as e:
        return False, str(e)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Batch reprocess unfinalized dossiers through grading + finalization → Finalized_Dossiers/<tier>."
    )
    parser.add_argument("--dry-run", action="store_true", help="Do not write; report what would be finalized.")
    parser.add_argument("--limit", type=int, default=0, help="Max number to process (0 = all).")
    parser.add_argument("--validate-before", action="store_true", help="Run run_all_tests + run_e2e_layers before batch; abort if either fails.")
    parser.add_argument("--validate-after", action="store_true", help="Run run_all_tests + run_e2e_layers after batch; exit 1 if either fails.")
    parser.add_argument("--verbose", "-v", action="store_true", help="Print Zeus failure reason for each rejected dossier.")
    args = parser.parse_args()

    if args.validate_before:
        print("Validation before batch...")
        ok, msg = run_validation()
        if not ok:
            print(f"  ABORT: {msg}")
            return 1
        print(f"  {msg}\n")

    from PX_Warehouse.warehouse_layout import TIERS, get_finalized_dossier_dir, get_discovery_accepted_dir
    from PX_Warehouse.Finalization_Pipeline import run_finalization, finalize_and_place, passes_discovery_bar

    wh = REPO_ROOT / "PX_Warehouse"
    discovery_dir = get_discovery_accepted_dir(REPO_ROOT)
    to_process: list[tuple[Path, str, bool]] = []
    for is_novel, folder in [(True, "Novel_Dossiers"), (False, "Prv_Dossiers")]:
        for tier in TIERS:
            d_dir = wh / folder / tier
            if not d_dir.exists():
                continue
            for f in d_dir.glob("*.json"):
                if f.name.startswith("PATH_TEST"):
                    continue
                item_id = f.stem
                to_process.append((f, item_id, is_novel))

    finalized_dir = wh / "Finalized_Dossiers"
    seen_finalized: set[str] = set()
    if finalized_dir.exists():
        for t in TIERS:
            td = finalized_dir / t
            if td.exists():
                for f in td.glob("*.json"):
                    if not f.name.startswith("PATH_TEST"):
                        seen_finalized.add(f.stem)
        # Two-tier: also skip if already in Discovery_Accepted
        discovery_dir = finalized_dir / "Discovery_Accepted"
        if discovery_dir.exists():
            for f in discovery_dir.glob("*.json"):
                if not f.name.startswith("PATH_TEST"):
                    seen_finalized.add(f.stem)
    to_process = [(p, iid, nov) for p, iid, nov in to_process if iid not in seen_finalized]

    if args.limit > 0:
        to_process = to_process[: args.limit]

    total = len(to_process)
    print(f"Reprocess warehouse lifecycle: {total} unfinalized dossiers (dry_run={args.dry_run}, validate_after={args.validate_after})")
    if total == 0:
        print("Nothing to process. Exiting.")
        return 0

    # Two-tier bar (re-ordered grading): Zeus 0.021 → Finalized_Dossiers/<tier>; discovery bar 0.05 → Discovery_Accepted; else reject
    ok_finalized = 0
    ok_discovery = 0
    fail = 0
    for idx, (path, item_id, is_novel) in enumerate(to_process, start=1):
        if total >= 50 and idx % 100 == 0:
            print(f"  Progress: {idx}/{total} (finalized={ok_finalized}, discovery={ok_discovery}, rejected={fail})")
        try:
            dossier = json.loads(path.read_text(encoding="utf-8"))
        except Exception as e:
            print(f"  Skip {path.name}: read error {e}")
            from PX_System.finalization_log import log_finalization_failure
            log_finalization_failure(
                source_file="reprocess_warehouse_lifecycle.py",
                candidate_id=item_id,
                error=str(e),
                context="dossier JSON read during warehouse lifecycle reprocessing",
            )
            fail += 1
            continue
        if args.dry_run:
            finalized, tier, err = run_finalization(dossier, item_id, is_novel, REPO_ROOT)
            if err:
                print(f"  {item_id}: {err}")
                fail += 1
                continue
            zv = (finalized.get("finalization") or {}).get("zeus_verdict") or {}
            if zv.get("authorized"):
                print(f"  [dry-run] Would finalize {item_id} → Finalized_Dossiers/{tier}/")
                ok_finalized += 1
            elif passes_discovery_bar(finalized):
                print(f"  [dry-run] Would place {item_id} → Finalized_Dossiers/Discovery_Accepted/")
                ok_discovery += 1
            else:
                if args.verbose:
                    failed = [k for k, v in (zv.get("laws_results") or {}).items() if isinstance(v, dict) and not v.get("passed")]
                    print(f"  {item_id}: would reject (Zeus + discovery bar) — failed: {failed}")
                fail += 1
            continue
        out_path = finalize_and_place(dossier, item_id, is_novel, REPO_ROOT)
        if out_path is None:
            if args.verbose:
                print(f"  {item_id}: rejected (did not pass Zeus or discovery bar)")
            fail += 1
            continue
        if discovery_dir is not None and str(discovery_dir) in str(out_path):
            ok_discovery += 1
        else:
            ok_finalized += 1

    print(f"\nDone: {ok_finalized} finalized (Zeus 0.021 → Finalized_Dossiers/<tier>), {ok_discovery} discovery-accepted (→ Discovery_Accepted), {fail} rejected.")

    if args.validate_after and (ok_finalized + ok_discovery) > 0:
        print("Validation after batch...")
        ok_val, msg = run_validation()
        if not ok_val:
            print(f"  FAIL: {msg}")
            return 1
        print(f"  {msg}")

    return 0 if fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
