#!/usr/bin/env python3
"""
E2E test: Genesis Feed → Novel PRV Orchestrator.

1. Generate novel molecules via Genesis Engine (Vector Core gate).
2. Write them to a dedicated test queue (e2e_test_queue.json).
3. Run the novel PRV orchestrator on that queue (PRV_MAX_ITEMS = count).

Run from repo root: python tests/run_genesis_e2e_test.py
"""
from __future__ import annotations

import os
import sys
import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

E2E_QUEUE_NAME = "e2e_test_queue.json"
GENESIS_COUNT = 3
RNG_SEED = 999


def main() -> int:
    print("=" * 60)
    print("E2E: Genesis → Novel PRV Orchestrator")
    print("=" * 60)

    # 1. Generate novel candidates via Genesis
    from PX_Engine.Genesis_Engine import generate_novel_candidates
    candidates = generate_novel_candidates(n=GENESIS_COUNT, use_vector_core=True, rng_seed=RNG_SEED)
    if not candidates:
        print("FAIL: Genesis produced no candidates.")
        return 1
    print(f"\n1. Genesis generated {len(candidates)} novel candidates (Vector Core authorized).")

    # 2. Write test queue to Feeder (canonical queue location)
    from PX_Warehouse.warehouse_layout import get_feeder_dir
    feeder_dir = get_feeder_dir(REPO_ROOT)
    feeder_dir.mkdir(parents=True, exist_ok=True)
    queue_path = feeder_dir / E2E_QUEUE_NAME
    queue_path.write_text(json.dumps({"candidates": candidates}, indent=2), encoding="utf-8")
    print(f"2. Wrote {E2E_QUEUE_NAME} with {len(candidates)} candidates.")

    # 3. Run novel orchestrator on this queue
    os.environ["PRV_QUEUE_FILE"] = E2E_QUEUE_NAME
    os.environ["PRV_MODE"] = "novel"
    os.environ["PRV_MAX_ITEMS"] = str(GENESIS_COUNT)
    os.environ["PRV_24H_DURATION_SEC"] = "600"
    if "PRV_API_PACING_SEC" not in os.environ:
        os.environ["PRV_API_PACING_SEC"] = "2"

    print(f"3. Running novel PRV orchestrator (PRV_MAX_ITEMS={GENESIS_COUNT}, pacing 2s)...\n")
    from PX_Executive.PRV_24H_Orchestrator import main as orchestrator_main
    code = orchestrator_main()
    if code != 0:
        print(f"\nE2E: Orchestrator exited with code {code}.")
        return code
    print("\n" + "=" * 60)
    print("E2E PASS: Genesis → Queue → Novel PRV pipeline completed.")
    print("=" * 60)
    return 0


if __name__ == "__main__":
    sys.exit(main())
