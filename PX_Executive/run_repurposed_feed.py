#!/usr/bin/env python3
"""
Repurposed feed — populate the PRV queue with repurposed (R) candidates from seed libraries.

Sources (from intake_policy.json): ChEMBL, DrugBank, ClinicalTrials, FDA_OrangeBook,
PubChem, NCATS_Repurposing, ZINC. Data is read from PX_Data/<source>/ and written to
PX_Warehouse/Feeder/prv_24h_queue.json (canonical queue for the 24H orchestrator).

Usage:
  python PX_Executive/run_repurposed_feed.py         # merge into existing queue (keeps novel N)
  python PX_Executive/run_repurposed_feed.py --replace   # overwrite queue with harvest only

The repurposed orchestrator (run_prv_repurposed.py) reads from the same Feeder queue
and processes only type "R". Run this feed before or alongside the orchestrator so
the repurposed pipeline has candidates.
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
