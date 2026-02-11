"""
GOVERNANCE POISON PILL — Non-regression canary test.

IDE classification (recommended):
  type: GOVERNANCE_POISON_PILL
  violation: CFR_ENVELOPE
  expected_behavior: HARD_REJECT
  adapter_contract: returns None
  downstream_calls: NONE

IDE instructions (drop-in):
  test_type: GOVERNANCE_POISON_PILL
  layer_expected_to_fail: PX_Validation
  failure_class: ONTOLOGY_VIOLATION
  artifact_expected: false
  state_mutation_expected: false
  retry_allowed: false
  non_optional: true
  non_skippable: true
  pre_execution_gate: true
  pipeline_policy: If this ever passes incorrectly → block the entire pipeline.

This canary must scream forever: the Sovereign System MUST reject the poison pill.
If this test ever passes with result != None, governance has failed.
"""
import os
import sys

# Classification: canonical tag for this non-regression test
GOVERNANCE_POISON_PILL_CLASSIFICATION = {
    "type": "GOVERNANCE_POISON_PILL",
    "violation": "CFR_ENVELOPE",
    "expected_behavior": "HARD_REJECT",
    "adapter_contract": "returns None",
    "downstream_calls": "NONE",
}

# IDE instructions (drop-in): tell the IDE how to treat this test
IDE_INSTRUCTIONS = {
    "test_type": "GOVERNANCE_POISON_PILL",
    "layer_expected_to_fail": "PX_Validation",
    "failure_class": "ONTOLOGY_VIOLATION",
    "artifact_expected": False,
    "state_mutation_expected": False,
    "retry_allowed": False,
    "non_optional": True,
    "non_skippable": True,
    "pre_execution_gate": True,
    "pipeline_policy": "If this ever passes incorrectly → block the entire pipeline.",
}

# Governance signature: IDE can hash this and alert on semantic drift
GOVERNANCE_SIGNATURE = (
    IDE_INSTRUCTIONS["test_type"],
    IDE_INSTRUCTIONS["failure_class"],
    IDE_INSTRUCTIONS["pipeline_policy"],
)

# Ensure we can import from the root context
ROOT_DIR = os.getcwd()
if ROOT_DIR not in sys.path:
    sys.path.append(ROOT_DIR)

try:
    from Nipah_Analysis.adapters.nipah_miner_adapter import NipahMinerAdapter
except ImportError:
    print("CRITICAL: Could not import NipahMinerAdapter. Run this from E:\\foundation")
    sys.exit(1)

def run_governance_stress_test():
    print("╔══════════════════════════════════════════════════════════════════════════════╗")
    print("║ GOVERNANCE STRESS TEST: 'ZOMBIE VIRUS' (THERMODYNAMIC PARADOX)              ║")
    print("╚══════════════════════════════════════════════════════════════════════════════╝")
    
    # 1. Setup
    adapter = NipahMinerAdapter()
    toxic_file = os.path.join("Nipah_Analysis", "data", "raw", "toxic_nipah_candidate.json")
    
    if not os.path.exists(toxic_file):
        print(f"❌ TEST ERROR: Poison pill file not found at: {toxic_file}")
        print("   pre_execution_gate failed → blocking pipeline.")
        sys.exit(1)

    print(f"TARGET: {toxic_file}")
    print(f"TOXIN:  High Complexity (Variance) + Zero Energy (0% CFR)")
    print(f"CLASS:  {GOVERNANCE_POISON_PILL_CLASSIFICATION['type']} | {GOVERNANCE_POISON_PILL_CLASSIFICATION['violation']} → {GOVERNANCE_POISON_PILL_CLASSIFICATION['expected_behavior']}")
    print("-" * 80)

    # 2. Execution — adapter_contract: returns None; downstream_calls: NONE
    result = adapter.mine_candidate(toxic_file)
    
    print("-" * 80)

    # 3. Verification — canary: must HARD_REJECT (adapter_contract: returns None)
    # If this ever passes incorrectly (result is not None) → block the entire pipeline.
    if result is not None:
        print("❌ GOVERNANCE_POISON_PILL canary broken: adapter accepted poison pill.")
        print("   pipeline_policy: block the entire pipeline.")
        sys.exit(1)
    assert result is None  # redundant guard; exit(1) above blocks pipeline

    # 4. Layer pinning — failure must occur at PX_Validation (ontology), not Engine/Audit
    # Generalized monotonicity check (same contract for all adapters)
    from governance.layer_monotonicity import assert_layer_monotonicity
    assert_layer_monotonicity(adapter, "PX_Validation")
    print("✅ TEST PASSED: DEFENSE SUCCESSFUL")
    print("   The Sovereign System correctly IDENTIFIED and REJECTED the physics violation.")
    print(f"   Failure layer: {adapter.last_failure_layer} (expected).")
    print("   The Byzantium Council (or Vector Core) enforced the Constitution.")

    # Stamp for execution-order enforcement: gate passed this run
    _write_gate_stamp()


def _write_gate_stamp():
    """Write stamp file so pipeline runners can verify poison pill gate passed."""
    stamp_dir = os.path.join(ROOT_DIR, "governance")
    stamp_file = os.path.join(stamp_dir, ".poison_pill_gate_passed")
    os.makedirs(stamp_dir, exist_ok=True)
    signature_digest = str(hash(GOVERNANCE_SIGNATURE))
    policy_line = IDE_INSTRUCTIONS["pipeline_policy"]
    payload = (signature_digest + "\n" + policy_line).encode("utf-8")
    lines = [signature_digest, policy_line]
    try:
        from PX_Constitution.gate_stamp import sign_gate_stamp, get_gate_stamp_secret
        if get_gate_stamp_secret():
            hmac_hex = sign_gate_stamp(payload)
            lines.append("hmac_sha256=" + hmac_hex)
    except Exception:
        pass
    with open(stamp_file, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    # Auto-generate governance coverage metric for audits
    try:
        from governance.coverage import generate_governance_coverage
        generate_governance_coverage()
    except Exception:
        pass  # non-fatal; coverage is informational


if __name__ == "__main__":
    run_governance_stress_test()
    # Exit 0 only if gate passed; any failure above already sys.exit(1)
