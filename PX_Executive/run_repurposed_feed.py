#!/usr/bin/env python3
"""
DEPRECATED: Use px_feed.py --mode repurpose instead.

  python PX_Executive/px_feed.py --mode repurpose
  python PX_Executive/px_feed.py --mode repurpose --replace

This file is kept for backward compatibility.
"""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
HARVESTER = REPO_ROOT / "PX_Warehouse" / "Operations" / "scripts" / "run_intake_harvester.py"


def main() -> int:
    replace = "--replace" in sys.argv
    argv = [sys.executable, str(HARVESTER)]
    if not replace:
        argv.append("--merge")
    env = os.environ.copy()
    env["PYTHONUNBUFFERED"] = "1"
    if str(REPO_ROOT) not in os.environ.get("PYTHONPATH", ""):
        env["PYTHONPATH"] = str(REPO_ROOT) + (os.pathsep + os.environ.get("PYTHONPATH", "") if os.environ.get("PYTHONPATH") else "")

    print("Repurposed feed: intake harvester → PX_Warehouse/Feeder/prv_24h_queue.json")
    print("Mode:", "replace (overwrite)" if replace else "merge (keep existing novel)")
    r = subprocess.run(argv, cwd=str(REPO_ROOT), env=env)
    if r.returncode == 0:
        print("Done. Run repurposed orchestrator: python PX_Executive/run_prv_repurposed.py")
    return r.returncode


if __name__ == "__main__":
    sys.exit(main())
