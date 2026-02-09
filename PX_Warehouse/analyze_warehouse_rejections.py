#!/usr/bin/env python3
"""
Analyze warehouse: all raw/processed dossiers that did NOT pass to Finalized_Dossiers.
Runs finalization (without writing) on each unfinalized dossier and aggregates:
  - Why they failed (Zeus: L11, L1, U34, U27; or early error)
  - Distribution of toxicity_index and harm_energy for rejections
  - Comparison to current thresholds (0.0210) and grading schema

Purpose: Recalibrate standards — ensure we're not too strict vs industry and not impeding process.

Usage (from repo root):
  python PX_Warehouse/analyze_warehouse_rejections.py
  python PX_Warehouse/analyze_warehouse_rejections.py --limit 50
  python PX_Warehouse/analyze_warehouse_rejections.py --limit 0   # all
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from collections import defaultdict

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

# Thresholds (governance lock)
ZEUS_TOXICITY_THRESHOLD = 0.0210
ZEUS_HARM_THRESHOLD = 0.0210
WAREHOUSE_GOLD_TOX = 0.021   # get_tier: Gold if tox < 0.021
WAREHOUSE_SILVER_TOX = 0.10  # Silver if tox < 0.10


def main() -> int:
    import argparse
    parser = argparse.ArgumentParser(description="Analyze why processed dossiers did not finalize")
    parser.add_argument("--limit", type=int, default=200, help="Max dossiers to run (0 = all)")
    parser.add_argument("--sample", type=int, default=0, help="If set, randomly sample this many instead of first N")
    args = parser.parse_args()

    from PX_Warehouse.warehouse_layout import TIERS
    from PX_Warehouse.Finalization_Pipeline import run_finalization

    wh = REPO_ROOT / "PX_Warehouse"
    to_process: list[tuple[Path, str, bool]] = []
    for is_novel, folder in [(True, "Novel_Dossiers"), (False, "Prv_Dossiers")]:
        for tier in TIERS:
            d_dir = wh / folder / tier
            if not d_dir.exists():
                continue
            for f in d_dir.glob("*.json"):
                if f.name.startswith("PATH_TEST"):
                    continue
                to_process.append((f, f.stem, is_novel))

    finalized_dir = wh / "Finalized_Dossiers"
    seen_finalized: set[str] = set()
    if finalized_dir.exists():
        for t in TIERS:
            td = finalized_dir / t
            if td.exists():
                for x in td.glob("*.json"):
                    if not x.name.startswith("PATH_TEST"):
                        seen_finalized.add(x.stem)

    to_process = [(p, iid, nov) for p, iid, nov in to_process if iid not in seen_finalized]
    if args.sample > 0 and len(to_process) > args.sample:
        import random
        to_process = random.sample(to_process, args.sample)
    elif args.limit > 0:
        to_process = to_process[: args.limit]

    total = len(to_process)
    print(f"Unfinalized dossiers (Prv+Novel, not in Finalized_Dossiers): {total}")
    if total == 0:
        print("Nothing to analyze. Exiting.")
        return 0

    # Run finalization (no write); collect outcomes
    early_errors = []
    zeus_rejected = []   # list of {item_id, toxicity_index, harm_energy, failed_laws, rationale, whole_profile_passed}
    passed = 0

    for idx, (path, item_id, is_novel) in enumerate(to_process, start=1):
        if total >= 20 and idx % 50 == 0:
            print(f"  Progress: {idx}/{total} (passed={passed}, zeus_rejected={len(zeus_rejected)}, errors={len(early_errors)})")
        try:
            dossier = json.loads(path.read_text(encoding="utf-8"))
        except Exception as e:
            early_errors.append({"item_id": item_id, "error": str(e)})
            continue
        finalized, tier, err = run_finalization(dossier, item_id, is_novel, REPO_ROOT)
        if err:
            early_errors.append({"item_id": item_id, "error": err})
            continue
        zv = (finalized.get("finalization") or {}).get("zeus_verdict") or {}
        if zv.get("authorized"):
            passed += 1
            continue
        laws = zv.get("laws_results") or {}
        failed = [k for k, v in laws.items() if isinstance(v, dict) and not v.get("passed")]
        zeus_rejected.append({
            "item_id": item_id,
            "tier_computed": tier,
            "toxicity_index": zv.get("toxicity_index"),
            "harm_energy": zv.get("harm_energy"),
            "failed_laws": failed,
            "rationale": zv.get("rationale", ""),
            "whole_profile_passed": zv.get("whole_profile_passed", False),
            "all_single_passed": zv.get("all_single_passed", False),
        })

    # --- Aggregation ---
    by_failed_law = defaultdict(int)
    tox_values = []
    harm_values = []
    for r in zeus_rejected:
        for law in r["failed_laws"]:
            by_failed_law[law] += 1
        if r.get("toxicity_index") is not None:
            try:
                tox_values.append(float(r["toxicity_index"]))
            except (TypeError, ValueError):
                pass
        if r.get("harm_energy") is not None:
            try:
                harm_values.append(float(r["harm_energy"]))
            except (TypeError, ValueError):
                pass

    # --- Report ---
    print()
    print("=" * 70)
    print("WAREHOUSE REJECTION ANALYSIS")
    print("=" * 70)
    print(f"  Total processed (unfinalized): {total}")
    print(f"  Would pass Zeus (authorized):  {passed}")
    print(f"  Zeus rejected:                 {len(zeus_rejected)}")
    print(f"  Early errors (load/finalize):  {len(early_errors)}")
    print()

    if early_errors:
        print("Early errors (first 10):")
        for e in early_errors[:10]:
            print(f"  {e['item_id']}: {e['error'][:80]}")
        print()

    print("Zeus rejection — failed laws (count):")
    for law in sorted(by_failed_law.keys()):
        print(f"  {law}: {by_failed_law[law]}")
    print()

    if tox_values:
        tox_values.sort()
        n = len(tox_values)
        print("Toxicity index (Zeus-rejected) distribution:")
        print(f"  min={min(tox_values):.4f}, max={max(tox_values):.4f}, median={tox_values[n//2]:.4f}")
        print(f"  Current Zeus threshold (L11): {ZEUS_TOXICITY_THRESHOLD}")
        above = sum(1 for t in tox_values if t >= ZEUS_TOXICITY_THRESHOLD)
        print(f"  Rejected with tox >= threshold: {above}/{n}")
        # Bins
        bins = [(0.021, 0.05), (0.05, 0.10), (0.10, 0.20), (0.20, 0.50), (0.50, 2.0)]
        for lo, hi in bins:
            c = sum(1 for t in tox_values if lo <= t < hi)
            if c:
                print(f"    [{lo:.2f}, {hi:.2f}): {c}")
    print()

    if harm_values:
        harm_values.sort()
        n = len(harm_values)
        print("Harm energy (Zeus-rejected) distribution:")
        print(f"  min={min(harm_values):.4f}, max={max(harm_values):.4f}, median={harm_values[n//2]:.4f}")
        print(f"  Current Zeus threshold (L1): {ZEUS_HARM_THRESHOLD}")
    print()

    # --- Industry comparison and recommendations ---
    print("=" * 70)
    print("THRESHOLDS & RECALIBRATION")
    print("=" * 70)
    print("Current (governance lock):")
    print("  TOXICITY_HARD_LIMIT (L11/L1): 0.0210  — reject if toxicity_index or harm_energy >= 0.0210")
    print("  Whole-profile override: TOXICITY_DIAMOND, or tox<0.01, or (tox>0.02 and safety_margin>50)")
    print("  Warehouse tier: Gold tox<0.021, Silver tox<0.10, Bronze else")
    print()
    print("Discovery grading schema (GradingSchema_Discovery.json):")
    print("  GOLD_TIER:  toxicity_max 0.2 (20%)   — much looser than Zeus 2.1%")
    print("  SILVER_TIER: toxicity_max 0.35 (35%)")
    print("  BRONZE_TIER: toxicity_max 0.5 (50%)")
    print()
    print("Industry context:")
    print("  - Preclinical toxicity filters often use 10–30% predicted risk as 'acceptable' for progression.")
    print("  - 2.1% (0.021) is very conservative and aligns with late-stage / high-safety bar.")
    print("  - Consequence: many discovery-stage candidates fail Zeus despite being Silver/Bronze by grading.")
    print()
    print("Recommendations:")
    if by_failed_law.get("L11_DETERMINISTIC_ENGINE", 0) + by_failed_law.get("L1_HARM_LAW", 0) > len(zeus_rejected) * 0.8:
        print("  1. Most rejections are L11/L1 (toxicity or harm >= 0.021). Consider:")
        print("     - Keeping 0.021 for WRITE to Finalized_Dossiers (compliance bar) but recording")
        print("       'Zeus-rejected' dossiers in a separate audit path (e.g. Operations/Zeus_Rejected)")
        print("       so they remain visible and can be re-reviewed if thresholds are ever recalibrated.")
        print("     - Or introduce a two-tier bar: 'Finalized' (0.021) vs 'Discovery_Accepted' (e.g. 0.05)")
        print("       with only Finalized going to Finalized_Dossiers.")
    print("  2. Align Zeus single-metric bar with grading: e.g. allow Zeus pass when dossier tier is")
    print("     Silver or better (tox < 0.10) and whole_profile or L11/L1 pass at a relaxed discovery")
    print("     threshold (e.g. 0.05) for discovery-only pipeline; keep 0.021 for regulatory/filing path.")
    print("  3. Do not loosen 0.021 for any regulatory-facing output without explicit governance change.")
    print("  4. Add an audit report (this script) to CI or weekly run to monitor rejection rate and")
    print("     distribution so recalibration is data-driven.")
    print()

    return 0


if __name__ == "__main__":
    sys.exit(main())
