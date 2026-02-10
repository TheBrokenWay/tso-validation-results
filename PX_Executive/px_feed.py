#!/usr/bin/env python3
"""
px_feed.py — Unified Feed Orchestrator

Single entry point for both novel (Genesis) and repurposed (intake harvest) feeds.
Writes candidates to PX_Warehouse/Feeder/prv_24h_queue.json.

Usage:
    python PX_Executive/px_feed.py --mode novel          # Genesis_Engine → new SMILES
    python PX_Executive/px_feed.py --mode repurpose      # Existing drugs → queue
    python PX_Executive/px_feed.py --mode novel --count 50
    python PX_Executive/px_feed.py --mode novel --high-chaos
    python PX_Executive/px_feed.py --mode novel --interval 0   # Single batch, exit

Environment variables (novel mode):
    GENESIS_COUNT       Molecules per batch (default 20, 50 in high-chaos)
    GENESIS_INTERVAL    Seconds between batches (0 = single run)
    GENESIS_SEED        Fixed random seed for first batch
    GENESIS_HIGH_CHAOS  High-entropy mode (expanded chemistry)
    GENESIS_SKIP_VECTOR_CORE  Bypass Vector Core physics checks (not recommended)
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def run_novel_feed(count: int, interval: int, high_chaos: bool, seed: int | None) -> int:
    """Run Genesis Engine to invent novel SMILES and write to queue."""
    from PX_Engine.Genesis_Engine import generate_novel_candidates
    from PX_Warehouse.warehouse_layout import get_feeder_dir, get_queue_path

    queue_name = os.environ.get("GENESIS_QUEUE", "prv_24h_queue.json")
    read_path = get_queue_path(queue_name, REPO_ROOT)
    feeder_dir = get_feeder_dir(REPO_ROOT)

    use_vector_core = os.environ.get("GENESIS_SKIP_VECTOR_CORE", "").strip().lower() not in ("1", "true", "yes")

    print("=" * 64)
    print("   GENESIS FEED :: " + ("HIGH ENTROPY" if high_chaos else "STANDARD"))
    print("=" * 64)
    print(f"BATCH SIZE: {count} | INTERVAL: {interval}s | VECTOR_CORE: {use_vector_core}")
    print("-" * 64)

    batch_count = 0
    while True:
        batch_count += 1
        ts = time.strftime("%H:%M:%S")
        print(f"[{ts}] Batch #{batch_count}: generating {count} candidates...")

        rng_seed = seed if (batch_count == 1 and seed is not None) else None
        candidates = generate_novel_candidates(
            n=count,
            use_vector_core=use_vector_core,
            rng_seed=rng_seed,
            high_entropy=high_chaos,
        )

        if not candidates:
            print(f"[{ts}] No candidates generated.")
            if interval <= 0:
                return 1
            time.sleep(max(1, interval))
            continue

        # Load existing queue (read-modify-write, dedupe by id)
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
            print(f"[{ts}] Injected {added} new (queue: {len(existing)} unique)")
        except Exception as e:
            print(f"[{ts}] Write failed: {e}", file=sys.stderr)
            if interval <= 0:
                return 1

        if interval <= 0:
            return 0

        time.sleep(interval)


def run_repurposed_feed(replace: bool) -> int:
    """Run intake harvester to populate queue with repurposed candidates."""
    harvester = REPO_ROOT / "PX_Warehouse" / "Operations" / "scripts" / "run_intake_harvester.py"
    if not harvester.exists():
        print(f"Harvester not found: {harvester}", file=sys.stderr)
        return 1

    argv = [sys.executable, str(harvester)]
    if not replace:
        argv.append("--merge")
    env = os.environ.copy()
    env["PYTHONUNBUFFERED"] = "1"
    if str(REPO_ROOT) not in os.environ.get("PYTHONPATH", ""):
        env["PYTHONPATH"] = str(REPO_ROOT) + (os.pathsep + os.environ.get("PYTHONPATH", "") if os.environ.get("PYTHONPATH") else "")

    print("=" * 64)
    print("   REPURPOSED FEED :: INTAKE HARVESTER")
    print("=" * 64)
    print(f"Mode: {'replace' if replace else 'merge (keep existing novel)'}")
    print("-" * 64)

    r = subprocess.run(argv, cwd=str(REPO_ROOT), env=env)
    if r.returncode == 0:
        print("Done. Queue ready for: python PX_Executive/px_prv.py --type repurpose")
    return r.returncode


def main() -> int:
    parser = argparse.ArgumentParser(description="Unified feed orchestrator (novel + repurposed)")
    parser.add_argument("--mode", choices=["novel", "repurpose"], required=True, help="Feed mode")
    parser.add_argument("--count", type=int, default=None, help="Candidates per batch (novel mode)")
    parser.add_argument("--interval", type=int, default=None, help="Seconds between batches (0 = single run)")
    parser.add_argument("--high-chaos", action="store_true", help="High-entropy Genesis mode")
    parser.add_argument("--seed", type=int, default=None, help="Fixed random seed for first batch")
    parser.add_argument("--replace", action="store_true", help="Replace queue (repurpose mode)")
    args = parser.parse_args()

    if args.mode == "novel":
        high_chaos = args.high_chaos or os.environ.get("GENESIS_HIGH_CHAOS", "").strip().lower() in ("1", "true", "yes")
        default_n = 50 if high_chaos else 20
        count = args.count or int(os.environ.get("GENESIS_COUNT", str(default_n)))
        interval = args.interval if args.interval is not None else int(os.environ.get("GENESIS_INTERVAL", "60"))
        seed = args.seed
        if seed is None:
            seed_env = os.environ.get("GENESIS_SEED")
            if seed_env not in (None, ""):
                seed = int(seed_env)
        return run_novel_feed(count, interval, high_chaos, seed)
    else:
        return run_repurposed_feed(args.replace)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("\n[STOP] Feed shutdown.")
        sys.exit(0)
