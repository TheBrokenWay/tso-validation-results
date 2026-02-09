#!/usr/bin/env python3
"""
Purge error logs and clear failed FTO flags from the last batch so the 4,605
processed items (and full backlog) are re-evaluated with clean metrics.

- Archives RUN_LOG_*.md and FTO-blocked artifacts so Command Center yield
  reflects the next run (expect 15–20% after FTS5 quoting fix).
- Queue (prv_24h_queue.json) is unchanged; orchestrator re-run will process
  all candidates from the top, re-evaluating the previously failed batch.

Usage (from repo root):
  python PX_Executive/purge_fto_error_logs.py
  python PX_Executive/purge_fto_error_logs.py --dry-run
"""
from __future__ import annotations

import argparse
from pathlib import Path
from datetime import datetime, timezone

REPO_ROOT = Path(__file__).resolve().parents[1]
PX_LOGS = REPO_ROOT / "PX_LOGS"


def main() -> int:
    ap = argparse.ArgumentParser(description="Purge FTO error logs and clear failed-FTO flags.")
    ap.add_argument("--dry-run", action="store_true", help="Only list what would be archived.")
    args = ap.parse_args()

    stamp = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
    archive_dir = PX_LOGS / f"archive_fto_purge_{stamp}"

    to_archive: list[Path] = []
    # RUN_LOG files (drive yield / DONE counts in Command Center)
    for f in sorted(PX_LOGS.glob("RUN_LOG_*.md")):
        to_archive.append(f)
    # FTO-blocked artifacts
    for name in ("fto_audit.log", "fto_recheck_candidates.json"):
        p = PX_LOGS / name
        if p.exists():
            to_archive.append(p)
    for f in PX_LOGS.glob("FTO_BLOCKED_*.md"):
        to_archive.append(f)
    # Optional: monitoring summaries from same period (so yield is from fresh run only)
    for f in sorted(PX_LOGS.glob("REALTIME_MONITORING_SUMMARY_*.md")):
        to_archive.append(f)
    for f in sorted(PX_LOGS.glob("WAREHOUSE_DRIFT_REPORT_*.md")):
        to_archive.append(f)
    for f in sorted(PX_LOGS.glob("PRV_DISCOVERY_SUMMARY_*.md")):
        to_archive.append(f)

    if not to_archive:
        print("Nothing to purge (no RUN_LOG or FTO artifacts found).")
        return 0

    if args.dry_run:
        print(f"[DRY-RUN] Would archive {len(to_archive)} items to {archive_dir.name}/")
        for p in to_archive[:30]:
            print(f"  {p.name}")
        if len(to_archive) > 30:
            print(f"  ... and {len(to_archive) - 30} more")
        return 0

    archive_dir.mkdir(parents=True, exist_ok=True)
    for p in to_archive:
        try:
            dest = archive_dir / p.name
            if dest.exists():
                dest = archive_dir / f"{p.stem}_{id(p)}{p.suffix}"
            p.rename(dest)
            print(f"Archived: {p.name} -> {archive_dir.name}/")
        except Exception as e:
            print(f"Skip {p.name}: {e}")

    print(f"\nPurge complete. {len(to_archive)} items moved to {archive_dir}.")
    print("Re-run PRV_24H_Orchestrator to process backlog with FTS5 quoting; expect yield 15–20%.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
