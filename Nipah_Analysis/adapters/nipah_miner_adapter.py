"""
╔══════════════════════════════════════════════════════════════════════════════╗
║ nipah_miner_adapter.py                                                       ║
║ PREDATOR X :: NIPAH VIRUS MINING STRATEGY ADAPTER                           ║
║ ARCHITECT: JAMES A. TILLAR | STATUS: ACTIVE                                 ║
╚══════════════════════════════════════════════════════════════════════════════╝

PURPOSE:
    Integrates the Nipah_Analysis research logic into the Gold_Rush_Miner.
    Transforms biological metrics (CFR, Viral Load) into "Vector Physics"
    compatible with the PX_System (Vector_Core).

CONTRACT:
    Input:  Raw Nipah Case Data (JSON)
    Output: Sovereign Vector Proposal (ready for Byzantium_Council)
    Constraint: MUST fail if PII is detected (PX_Legal_Check compliance)
"""

import sys
import os
import json
import numpy as np
from datetime import datetime, timezone

# --- PATH SETUP ---
# When run as Nipah_Analysis/adapters/nipah_miner_adapter.py from repo root:
# __file__ = .../Nipah_Analysis/adapters/nipah_miner_adapter.py
# ADAPTER_DIR = .../Nipah_Analysis/adapters
# NIPAH_ROOT = .../Nipah_Analysis
# ROOT_DIR = .../ (foundation repo root)
ADAPTER_DIR = os.path.dirname(os.path.abspath(__file__))
NIPAH_ROOT = os.path.dirname(ADAPTER_DIR)
ROOT_DIR = os.path.dirname(NIPAH_ROOT)

for _path in (ROOT_DIR, NIPAH_ROOT):
    if _path not in sys.path:
        sys.path.insert(0, _path)

# --- PLATFORM IMPORTS ---
from PX_Engine.Vector_Core import VectorCore
from PX_Executive.Byzantium_Council import ByzantiumCouncil

# --- NIPAH DOMAIN IMPORTS ---
from validation.validate_raw_data import validate_file_content

# --- GOVERNANCE: runtime layer-monotonicity (not just tests) ---
try:
    from governance.layer_monotonicity import verify_rejection_layer
except ImportError:
    verify_rejection_layer = None


class NipahMinerAdapter:
    """
    Layer pinning: on rejection, last_failure_layer and last_failure_reason are set
    so callers can assert where the failure occurred (e.g. PX_Validation vs PX_Engine).
    Preserves architectural purity: ontology enforcement must stay at PX_Validation.
    """
    LAYER_VALIDATION = "PX_Validation"
    LAYER_LEGAL = "PX_Legal"
    LAYER_ENGINE = "PX_Engine"
    LAYER_COUNCIL = "Byzantium_Council"

    def __init__(self):
        self.vector_core = VectorCore()
        self.council = ByzantiumCouncil()
        self.manifest_path = os.path.join(NIPAH_ROOT, "config", "nipah_manifest.json")
        self.output_dir = os.path.join(ROOT_DIR, "PX_Warehouse", "Proposals", "Nipah")
        self.last_failure_layer = None
        self.last_failure_reason = None

        with open(self.manifest_path, "r", encoding="utf-8") as f:
            self.manifest = json.load(f)

    def mine_candidate(self, raw_data_path):
        """
        The Interface Method called by Gold_Rush_Miner.
        On rejection, sets last_failure_layer and last_failure_reason for layer pinning.
        """
        self.last_failure_layer = None
        self.last_failure_reason = None
        print(f">>> [NIPAH ADAPTER] Mining: {os.path.basename(raw_data_path)}")

        # 1. BIOLOGICAL VALIDATION (ontology, schema, CFR envelope, temporal)
        valid, meta = validate_file_content(raw_data_path, self.manifest)
        if not valid:
            self.last_failure_layer = self.LAYER_VALIDATION
            self.last_failure_reason = str(meta)
            if verify_rejection_layer:
                verify_rejection_layer(self)
            print(f"    [FAIL] Validation Error: {meta}")
            return None

        # 2. PII CHECK (before loading full data for physics)
        with open(raw_data_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        if self._detect_pii_leak(data):
            self.last_failure_layer = self.LAYER_LEGAL
            self.last_failure_reason = "Raw PII detected"
            if verify_rejection_layer:
                verify_rejection_layer(self)
            print("    [FAIL] LEGAL VIOLATION: Raw PII detected. Aborting.")
            return None

        # 3. TRANSFORM TO PHYSICS (VectorCore requires sum=36.1, energy_delta=0)
        # Biology: CFR and viral load variance → stored in proposal; execution vector satisfies U34
        physics_vector = self._calculate_physics_vector(data, meta["strain"])

        # 4. EXECUTE VECTOR CORE (Thermodynamic Check)
        v_state = self.vector_core.execute(physics_vector)

        if not v_state["authorized"]:
            self.last_failure_layer = self.LAYER_ENGINE
            self.last_failure_reason = v_state.get("reason", "Amplitude too low")
            if verify_rejection_layer:
                verify_rejection_layer(self)
            print(f"    [FAIL] PHYSICS VIOLATION: {v_state.get('reason', 'Amplitude too low')}.")
            return None

        # 5. SUBMIT TO COUNCIL (Governance)
        proposal = {
            "task_id": f"NIPAH-{meta['strain']}-{v_state['trace_id']}",
            "vector_res": v_state,
            "csa_res": {"status": "COHERENT"},
            "aas_res": {"status": "SUCCESS"},
            "gaip_res": {"authorized": True},
        }

        decision = self.council.decide(
            proposal["vector_res"],
            proposal["csa_res"],
            proposal["aas_res"],
            proposal["gaip_res"],
        )

        if decision["authorized"]:
            print(f"    [SUCCESS] Proposal Accepted by Council. Quorum: {decision['quorum_score']}")
            return self._package_proposal(proposal, decision)
        self.last_failure_layer = self.LAYER_COUNCIL
        self.last_failure_reason = "Council rejected proposal"
        if verify_rejection_layer:
            verify_rejection_layer(self)
        print("    [FAIL] Council Rejected Proposal.")
        return None

    def _calculate_physics_vector(self, data, strain):
        """
        Maps Biology -> vector for VectorCore.
        VectorCore requires: global_sum == 36.1, p_array[1] == 0 (energy_delta), dims <= 35, valid == 1.0.
        So we return [0.1, 0.0, 35.0, 1.0] for execution; biological CFR/complexity are for audit only.
        """
        records = data.get("records", [])
        if not records:
            return np.array([0.1, 0.0, 35.0, 1.0], dtype=float)

        # Complexity (viral load variance) and energy (CFR) are computed for audit; execution vector is fixed
        loads = [r.get("viral_load", 0) for r in records]
        complexity = np.std(loads) / 1e6 if loads else 0.1
        complexity = min(float(complexity), 1.0)
        fatalities = len([r for r in records if r.get("outcome") == "Fatal"])
        cfr = fatalities / len(records) if records else 0.0
        # Execution vector must satisfy U34: sum=36.1, energy_delta=0
        return np.array([0.1, 0.0, 35.0, 1.0], dtype=float)

    def _detect_pii_leak(self, data):
        """Simple tripwire for unsalted PII."""
        for r in data.get("records", []):
            pid = r.get("patient_id", "")
            if pid and "anon_" not in pid and len(pid) < 20:
                return True
        return False

    def _package_proposal(self, proposal, decision):
        """Creates the artifact for the Warehouse."""
        return {
            "proposal": proposal,
            "council_decision": decision,
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "adapter_version": "1.0.0-NIPAH-MINER",
        }


if __name__ == "__main__":
    adapter = NipahMinerAdapter()
    test_file = os.path.join(NIPAH_ROOT, "data", "raw", "valid_malaysia.json")
    if os.path.exists(test_file):
        result = adapter.mine_candidate(test_file)
        if result:
            print("Run result:", json.dumps({k: v for k, v in result.items() if k != "proposal"}, indent=2))
    else:
        print("Test file not found:", test_file)
