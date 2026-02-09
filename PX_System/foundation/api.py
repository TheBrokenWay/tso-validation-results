# Olympus API Layer - Deterministic Implementation
import os
import json

_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

def run_boot_sequence():
    state_path = os.path.join(_REPO_ROOT, "PX_STATE", "current_state.json")
    if os.path.exists(state_path):
        return "BOOT SEQUENCE: SYSTEM_READY_V2"
    return "BOOT SEQUENCE: INITIALIZING_CORE"

def verify_integrity():
    required = [
        os.path.join(_REPO_ROOT, "PX_System", "foundation", "ZeusLaws.py"),
        os.path.join(_REPO_ROOT, "PX_Executive", "Gold_Rush_Miner.py"),
    ]
    missing = [f for f in required if not os.path.exists(f)]
    if not missing:
        return "INTEGRITY: VERIFIED_BY_CORE_V2"
    return f"INTEGRITY: FAILED - MISSING {missing}"

def adjudicate_package(evidence):
    # Deterministic adjudication based on ZeusLaws
    return {
        "verdict": "VERIFIED_GOLD" if evidence.get("harm_energy", 1.0) < 0.02 else "REJECTED",
        "input_hash": evidence.get("constitutional_seal", "UNKNOWN")
    }

def enforce_law_l1(internal_snapshot, external_snapshot):
    """
    Law L1: Internal 35D physics snapshot is primary.
    External data may provide context but cannot override internal physics.
    """
    if external_snapshot is None:
        return internal_snapshot, False
    if internal_snapshot != external_snapshot:
        return internal_snapshot, True
    return internal_snapshot, False

def append_audit_record(record):
    audit_path = os.path.join(_REPO_ROOT, "PX_Audit", "audit_trail.jsonl")
    os.makedirs(os.path.dirname(audit_path), exist_ok=True)
    with open(audit_path, "a", encoding="utf-8") as f:
        f.write(json.dumps(record) + "\n")

def get_system_status():
    audit_path = os.path.join(_REPO_ROOT, "PX_Audit", "audit_trail.jsonl")
    if os.path.exists(audit_path):
        with open(audit_path, "r") as f:
            lines = f.readlines()
            return {"audit_chain_length": len(lines), "last_record": lines[-1] if lines else None}
    return {"audit_chain": "EMPTY"}
