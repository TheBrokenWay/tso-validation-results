"""
GOVERNANCE POISON PILL — Temporal paradox (outbreak year before first outbreak).

Expected: PX_Validation rejects with TEMPORAL_DRIFT (outbreak_year < first_outbreak_year).
Layer pinning: last_failure_layer == "PX_Validation".
"""
import os
import sys

ROOT_DIR = os.getcwd()
if ROOT_DIR not in sys.path:
    sys.path.append(ROOT_DIR)

try:
    from Nipah_Analysis.adapters.nipah_miner_adapter import NipahMinerAdapter
    from governance.layer_monotonicity import assert_layer_monotonicity
except ImportError as e:
    print("CRITICAL: Could not import adapter or layer_monotonicity. Run from E:\\foundation:", e)
    sys.exit(1)


def run_temporal_paradox_test():
    print("╔══════════════════════════════════════════════════════════════════════════════╗")
    print("║ GOVERNANCE STRESS TEST: TEMPORAL PARADOX (OUTBREAK YEAR BEFORE FIRST)       ║")
    print("╚══════════════════════════════════════════════════════════════════════════════╝")

    adapter = NipahMinerAdapter()
    toxic_file = os.path.join("Nipah_Analysis", "data", "raw", "temporal_paradox_nipah.json")

    if not os.path.exists(toxic_file):
        print(f"❌ TEST ERROR: Poison pill file not found at: {toxic_file}")
        sys.exit(1)

    print(f"TARGET: {toxic_file}")
    print("TOXIN:  outbreak_year 1998 for NiV_Bangladesh (first outbreak 2001)")
    print("-" * 80)

    result = adapter.mine_candidate(toxic_file)

    print("-" * 80)

    if result is not None:
        print("❌ GOVERNANCE_POISON_PILL canary broken: adapter accepted temporal paradox.")
        sys.exit(1)

    assert_layer_monotonicity(adapter, "PX_Validation")

    print("✅ TEST PASSED: DEFENSE SUCCESSFUL")
    print("   Failure layer: PX_Validation (expected). Temporal drift rejected.")
    print("   Layer monotonicity preserved.")

    try:
        from governance.coverage import generate_governance_coverage
        generate_governance_coverage()
    except Exception:
        pass


if __name__ == "__main__":
    run_temporal_paradox_test()
