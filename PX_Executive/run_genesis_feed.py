#!/usr/bin/env python3
"""
Genesis Feed (Continuous) — The Perpetual Creation Engine.

Purpose:
  Runs continuously to invent novel molecules (SMILES) via the Genesis Engine.
  Appends them to the PRV queue so the Novel Pipeline never starves.

Configuration (Env Vars):
  GENESIS_COUNT=20        Molecules to invent per batch (default 50 in high-chaos)
  GENESIS_INTERVAL=60     Seconds to wait between batches (omit or 0 = run once and exit)
  GENESIS_QUEUE           Target queue file (default: prv_24h_queue.json)
  GENESIS_SKIP_VECTOR_CORE  Set to "1" to bypass Vector Core physics checks (not recommended)
  GENESIS_SEED            Optional fixed random seed (if set, use only for first batch in continuous mode)
  GENESIS_HIGH_CHAOS=1    High-entropy mode: expanded chemical space, wilder molecules, more per batch
  GENESIS_ENTROPY=high    Same as GENESIS_HIGH_CHAOS=1

Run (perpetual):
  python PX_Executive/run_genesis_feed.py

Run (single batch then exit):
  GENESIS_INTERVAL=0 python PX_Executive/run_genesis_feed.py
"""
from __future__ import annotations

import os
import sys
import json
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from PX_Engine.Genesis_Engine import generate_novel_candidates
from PX_Warehouse.warehouse_layout import get_feeder_dir, get_queue_path


def main() -> int:
    # --- CONFIGURATION ---
    high_chaos = (
        os.environ.get("GENESIS_HIGH_CHAOS", "").strip().lower() in ("1", "true", "yes")
        or os.environ.get("GENESIS_ENTROPY", "").strip().lower() == "high"
    )
    default_n = "50" if high_chaos else "20"
    n = int(os.environ.get("GENESIS_COUNT", default_n))
    interval = int(os.environ.get("GENESIS_INTERVAL", "60"))
    seed_env = os.environ.get("GENESIS_SEED")
    rng_seed_first = int(seed_env) if seed_env not in (None, "") else None

    queue_name = os.environ.get("GENESIS_QUEUE", "prv_24h_queue.json")
    read_path = get_queue_path(queue_name, REPO_ROOT)
    feeder_dir = get_feeder_dir(REPO_ROOT)

    print("================================================================")
    print("   GENESIS DRIVE :: " + ("HIGH ENTROPY MODE" if high_chaos else "CONTINUOUS FEED MODE"))
    print("================================================================")
    print(f"TARGET QUEUE: {queue_name}")
    print(f"BATCH SIZE:   {n} candidates")
    if high_chaos:
        print("   (High chaos: expanded chemistry, wilder molecules)")
    print(f"INTERVAL:     {interval} seconds" + (" (single run)" if interval <= 0 else ""))
    print("----------------------------------------------------------------\n")

    use_vector_core = os.environ.get("GENESIS_SKIP_VECTOR_CORE", "").strip().lower() not in ("1", "true", "yes")
    # Production: Vector Core remains mandatory physics gate even in high-entropy (OPE + Law U34 for every novel SMILES)
    batch_count = 0

    while True:
        batch_count += 1
        timestamp = time.strftime("%H:%M:%S")
        print(f"[{timestamp}] Batch #{batch_count}: Engaging " + ("High-Entropy Genesis..." if high_chaos else "Vector Core..."))

        rng_seed = rng_seed_first if (batch_count == 1 and rng_seed_first is not None) else None
        candidates = generate_novel_candidates(
            n=n,
            use_vector_core=use_vector_core,
            rng_seed=rng_seed,
            high_entropy=high_chaos,
        )

        if not candidates:
            print(f"[{timestamp}] ⚠️  No candidates generated. Retrying in {max(1, interval)}s...")
            if interval <= 0:
                return 1
            time.sleep(max(1, interval))
            continue

        # Load queue (read-modify-write); dedupe by id so queue doesn't grow with duplicates
        existing: list = []
        if read_path.exists():
            try:
                data = json.loads(read_path.read_text(encoding="utf-8"))
                existing = data.get("candidates", []) if isinstance(data, dict) else []
                if not isinstance(existing, list):
                    existing = []
            except Exception:
                pass

        seen_ids = {str(c.get("id") or "").strip() for c in existing if c.get("id")}
        added = 0
        for c in candidates:
            kid = str(c.get("id") or "").strip()
            if kid and kid not in seen_ids:
                seen_ids.add(kid)
                existing.append(c)
                added += 1
        payload = {"candidates": existing}

        try:
            feeder_dir.mkdir(parents=True, exist_ok=True)
            out_path = feeder_dir / queue_name
            tmp_path = out_path.with_suffix(out_path.suffix + ".tmp")
            tmp_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
            tmp_path.replace(out_path)
            print(f"[{timestamp}] ✅ INJECTED: {added} new (queue has {len(existing)} unique). Queue size: {len(payload['candidates'])}")
        except Exception as e:
            print(f"[{timestamp}] ❌ Write failed: {e}")
            if interval <= 0:
                return 1

        if interval <= 0:
            return 0

        print(f"[{timestamp}] 💤 Cooling down for {interval}s...")
        time.sleep(interval)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("\n[STOP] Genesis Drive shutdown safely.")
        sys.exit(0)
