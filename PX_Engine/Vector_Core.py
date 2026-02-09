import numpy as np
import os
from datetime import datetime, timezone
from PX_Constitution.Virtual_Machine import get_vm_fingerprint

_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
METADATA = {
    "role": "VectorCore",
    "version": "1.2.0-GAIP",
    "vm_fingerprint": get_vm_fingerprint()
}

class VectorCore:
    def __init__(self, threshold=0.95, dims_limit=35.0, global_sum_target=36.1):
        self.threshold = threshold
        self.dims_limit = dims_limit
        self.global_sum_target = global_sum_target
        self.W = np.array([0.2, 0.196, 0.2, 0.2, 0.198])
        os.makedirs(os.path.join(_REPO_ROOT, "PX_Audit"), exist_ok=True)

    def execute(self, p_array, physical_descriptors=None):
        """
        Calculates dynamic resonance from 35D manifold density and physical descriptors.
        """
        energy_delta = p_array[1]
        global_sum = float(np.sum(p_array))
        if abs(global_sum - self.global_sum_target) > 1e-6:
            return {
                "trace_id": f"PRD-V-{int(np.sum(p_array * 1000) % 1e8)}",
                "authorized": False,
                "amplitude": 0.0,
                "metadata": METADATA,
                "ts": datetime.now(timezone.utc).isoformat(),
                "status": "LAW_U34_VIOLATION",
                "reason": "Law U34 Violation: GLOBAL_SUM mismatch",
                "global_sum": global_sum,
                "energy_delta": float(energy_delta),
            }
        if energy_delta != 0:
            return {
                "trace_id": f"PRD-V-{int(np.sum(p_array * 1000) % 1e8)}",
                "authorized": False,
                "amplitude": 0.0,
                "metadata": METADATA,
                "ts": datetime.now(timezone.utc).isoformat(),
                "status": "LAW_U34_VIOLATION",
                "reason": "Law U34 Violation: Energy Delta non-zero",
                "global_sum": global_sum,
                "energy_delta": float(energy_delta),
            }
        # 1. 35D Manifold Density Check
        # p_array[2] is the dimension count (must be <= 35)
        dims = p_array[2]
        if dims > self.dims_limit:
            return {"authorized": False, "reason": "Manifold Dimension Overflow"}

        # 2. Dynamic Resonance Calculation
        # Resonance is a function of complexity (p_array[0]) and physical validity
        # If physical_descriptors (from RDKit) are provided, use them to adjust W
        if physical_descriptors:
            # Adjust weights based on molecular weight and logp
            mw = physical_descriptors.get("molecular_weight", 350.0)
            logp = physical_descriptors.get("logp", 2.5)
            # Higher MW or extreme LogP adds friction to resonance
            friction = (mw / 1000.0) + (abs(logp - 2.5) / 10.0)
            resonance_base = 1.0 - (friction * 0.05)
        else:
            resonance_base = 0.98

        # Complexity Gate
        bridge = max(0.6, 1.0 - (p_array[0] * 0.25))
        eval_v = np.array([1.0, resonance_base, 1.0, bridge, 0.99])
        v_unified = np.dot(eval_v, self.W)
        
        # 3. Invariant Mask (Law U34 & L11)
        invariants = np.array([
            p_array[1] == 0,             # Energy Delta (Zero-Sum Law U34)
            p_array[2] <= self.dims_limit, # Dimensions
            p_array[3] == 1.0            # Validation Status
        ])
        
        # Hard-fail if v_unified < 0.95 (Constitutional Requirement)
        authorized = (v_unified >= 0.95) and np.all(invariants)
        
        return {
            "trace_id": f"PRD-V-{int(np.sum(p_array * 1000) % 1e8)}",
            "authorized": bool(authorized),
            "amplitude": float(f"{v_unified:.8f}"),
            "metadata": METADATA,
            "ts": datetime.now(timezone.utc).isoformat(),
            "status": "RESONANCE_ACHIEVED" if authorized else "DETERMINISTIC_FAIL",
            "global_sum": global_sum,
            "energy_delta": float(energy_delta),
        }

