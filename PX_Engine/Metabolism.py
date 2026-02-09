import json
import os
from PX_Engine.operations import OCE
from PX_Security.AAS_CircuitBreaker import AASCircuitBreaker

_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

class Metabolism:
    def __init__(self, state_path=None):
        self.state_path = state_path or os.path.join(_REPO_ROOT, "PX_Warehouse", "system_state.json")
        os.makedirs(os.path.dirname(self.state_path), exist_ok=True)
        self.cycle = self._load_cycle()
        self.circuit_breaker = AASCircuitBreaker()
        self.drift_lock_threshold = 0.05

    def _load_cycle(self):
        if os.path.exists(self.state_path):
            with open(self.state_path, "r", encoding="utf-8") as f:
                data = json.load(f)
            return data.get("metabolic_cycle") or data.get("Metabolic_Cycle", 0) or 0
        return 0

    def pulse(self, task_id, p_vec=None, csa_s=None):
        """
        Pulses the metabolic cycle. Link to drift-detection audit.
        System age only increments if 35D manifold resonance is achieved.
        """
        # Verify 35D Manifold Coherence (via OCE)
        payload = {
            "p_vector": p_vec or [0.1, 0.0, 35.0, 1.0],
            "csa_scores": csa_s or [1.0, 1.0, 1.0, 1.0, 1.0]
        }
        oce_result = OCE.execute(payload)

        drift_score = oce_result.get("drift_score", 1.0)

        # AAS L2/L14 Circuit Breaker: lock on drift
        if drift_score > self.drift_lock_threshold:
            self.circuit_breaker.lock(
                "MANIFOLD_DRIFT_DETECTED",
                {"drift_score": drift_score, "task_id": task_id, "law": "L2/L14"},
            )

        # Self-healing loop: refuse execution while locked until normalized
        if self.circuit_breaker.is_locked():
            normalization = self._normalize_manifold()
            if not normalization["normalized"]:
                status = "LOCKED_AAS_CIRCUIT"
                self._persist_state(task_id, status, oce_result)
                return self.cycle, status
            self.circuit_breaker.mark_normalized(normalization)

        if oce_result["authorized"]:
            self.cycle += 1
            status = "RESONANCE_ACHIEVED"
            self.circuit_breaker.close(
                "RESONANCE_ACHIEVED",
                {"drift_score": drift_score, "task_id": task_id},
            )
        else:
            status = "DRIFT_DETECTED_OPEN"
            self.circuit_breaker.open(
                "COHERENCE_BELOW_THRESHOLD",
                {"drift_score": drift_score, "task_id": task_id},
            )
            # Do not increment cycle if resonance is not achieved

        self._persist_state(task_id, status, oce_result)
            
        return self.cycle, status

    def get_current_age(self):
        return self.cycle

    def _normalize_manifold(self):
        normalization_payload = {
            "p_vector": [0.1, 0.0, 35.0, 1.0],
            "csa_scores": [1.0, 1.0, 1.0, 1.0, 1.0],
            "security_score": 1.0,
        }
        result = OCE.execute(normalization_payload)
        drift_score = result.get("drift_score", 1.0)
        normalized = result.get("authorized", False) and drift_score <= self.drift_lock_threshold
        return {
            "normalized": normalized,
            "drift_score": drift_score,
            "coherence": result.get("coherence"),
            "engine": result.get("engine"),
        }

    def _persist_state(self, task_id, status, oce_result):
        with open(self.state_path, "w") as f:
            json.dump({
                "metabolic_cycle": self.cycle,
                "last_action": f"PULSE_{task_id}",
                "manifold_status": status,
                "coherence": oce_result["coherence"],
                "drift_score": oce_result.get("drift_score"),
                "aas_circuit": self.circuit_breaker.snapshot(),
            }, f)
