#!/usr/bin/env python3
"""
DEPRECATED: Use px_finalize.py instead.

  python PX_Executive/px_finalize.py
  python PX_Executive/px_finalize.py --dry-run
  python PX_Executive/px_finalize.py --limit 100
  python PX_Executive/px_finalize.py --reprocess

This file is kept for backward compatibility and still works.
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
    parser = argparse.ArgumentParser(description="Finalize dossiers from Prv/Novel into Finalized_Dossiers.")
    parser.add_argument("--dry-run", action="store_true", help="Do not write; only report what would be finalized.")
    parser.add_argument("--limit", type=int, default=0, help="Max number to process (0 = all).")
    parser.add_argument("--reprocess", action="store_true", help="Process ALL dossiers including already-finalized (add new metrics, version, compliance report; Zeus at end).")
    parser.add_argument("--verbose", "-v", action="store_true", help="Print Zeus failure reason for each rejected dossier (toxicity/harm vs threshold 0.0210).")
    args = parser.parse_args()

    from PX_Warehouse.warehouse_layout import (
        TIERS,
        get_prv_dossier_dir,
        get_finalized_dossier_dir,
    )
    from PX_Warehouse.Finalization_Pipeline import run_finalization, write_finalized_dossier

    wh = REPO_ROOT / "PX_Warehouse"
    to_process: list[tuple[Path, str, bool]] = []  # (path, item_id, is_novel)
    # Process Novel first so we hit passing (low-tox) dossiers sooner when backfilling
    for is_novel, folder in [(True, "Novel_Dossiers"), (False, "Prv_Dossiers")]:
        for tier in TIERS:
            d_dir = wh / folder / tier
            if not d_dir.exists():
                continue
            for f in d_dir.glob("*.json"):
                if f.name.startswith("PATH_TEST"):
                    continue
                item_id = f.stem  # e.g. PRV_NOV_xxx or PRV_REP_xxx
                to_process.append((f, item_id, is_novel))

    # Skip already finalized unless --reprocess (reprocess adds new metrics, version, compliance report)
    if not args.reprocess:
        finalized_dir = wh / "Finalized_Dossiers"
        seen_finalized = set()
        if finalized_dir.exists():
            for t in TIERS:
                td = finalized_dir / t
                if td.exists():
                    for f in td.glob("*.json"):
                        if not f.name.startswith("PATH_TEST"):
                            seen_finalized.add(f.stem)
        to_process = [(p, iid, nov) for p, iid, nov in to_process if iid not in seen_finalized]

    if args.limit > 0:
        to_process = to_process[: args.limit]

    print(f"Finalize run: {len(to_process)} dossiers to process (dry_run={args.dry_run}, reprocess={args.reprocess}, verbose={args.verbose})")
    ok = 0
    fail = 0
    for path, item_id, is_novel in to_process:
        try:
            dossier = json.loads(path.read_text(encoding="utf-8"))
        except Exception as e:
            print(f"  Skip {path.name}: read error {e}")
            from PX_System.finalization_log import log_finalization_failure
            log_finalization_failure(
                source_file="run_finalize_dossiers.py",
                candidate_id=item_id,
                error=str(e),
                context="dossier JSON read before finalization",
            )
            fail += 1
            continue
        finalized, tier, err = run_finalization(dossier, item_id, is_novel, REPO_ROOT)
        if err:
            print(f"  {item_id}: {err}")
            from PX_System.finalization_log import log_finalization_failure
            log_finalization_failure(
                source_file="run_finalize_dossiers.py",
                candidate_id=item_id,
                error=err,
                context="run_finalization pipeline returned error",
            )
            fail += 1
            continue
        # Only the WRITE is gated by Zeus; full pipeline already ran for every dossier
        zeus_authorized = (finalized.get("finalization") or {}).get("zeus_verdict") or {}
        if not zeus_authorized.get("authorized"):
            if not args.dry_run:
                print(f"  {item_id}: full pipeline run, Zeus not authorized (not written)", end="")
                if args.verbose:
                    rationale = zeus_authorized.get("rationale", "")
                    laws = zeus_authorized.get("laws_results") or {}
                    failed = [k for k, v in laws.items() if isinstance(v, dict) and not v.get("passed")]
                    if failed:
                        print(f" — {rationale}; failed: {failed}")
                        for k in failed:
                            r = laws.get(k, {})
                            if isinstance(r, dict) and r.get("reason"):
                                print(f"      {k}: {r['reason']}")
                    else:
                        print(f" — {rationale}")
                else:
                    print()
            fail += 1
            continue
        if args.dry_run:
            print(f"  [dry-run] Would finalize {item_id} → Finalized_Dossiers/{tier}/")
            ok += 1
            continue
        out_path = write_finalized_dossier(finalized, item_id, tier, REPO_ROOT)
        print(f"  {item_id} → {out_path}")
        ok += 1

    print(f"Done: {ok} finalized (written), {fail} failed/skipped or Zeus-rejected (full pipeline run, not written).")
    return 0 if fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
