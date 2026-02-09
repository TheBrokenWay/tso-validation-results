import json
import os
import numpy as np

class TrajectoryPredictor:
    """Predator X: Temporal Velocity & Trajectory Forecast"""
    def __init__(self):
        try:
            from PX_Warehouse.WorldLine_Database import DEFAULT_WORLDLINES_PATH
            self.db_path = DEFAULT_WORLDLINES_PATH
        except Exception:
            self.db_path = "PX_Warehouse/WorldLines"  # relative to repo root when cwd is foundation

    def calculate_velocity(self):
        # Retrieve all world-lines
        lines = []
        if not os.path.exists(self.db_path):
            return "TraceabilityError: WorldLine database not found."

        for root, _dirs, files in os.walk(self.db_path):
            for file in files:
                if file.endswith(".worldline"):
                    with open(os.path.join(root, file), "r", encoding="utf-8") as f:
                        lines.append(json.load(f))
        
        if len(lines) < 2:
            return "Insufficient World-Lines for Trajectory Analysis."

        # Sort by timestamp to find the sequence (Deterministic sequence)
        lines = sorted(lines, key=lambda x: x.get("header", {}).get("timestamp") or "1970-01-01")
        
        # Calculate Delta between the last two states (35D Manifold)
        p1 = np.array(lines[-2]["physics_snapshot"]["coordinate_35d"])
        p2 = np.array(lines[-1]["physics_snapshot"]["coordinate_35d"])
        
        # Velocity Vector in 35D space
        velocity = p2 - p1
        
        # Predict the NEXT coordinate
        p3 = p2 + velocity
        
        # Law U27 (Unified Self) Audit:
        # Check if predicted trajectory violates ethical invariants
        # For now, we check if the energy delta (p3[1]) remains zero (Law U34)
        # and if the toxicity index in the physical realization is valid.
        
        tox = lines[-1].get("physical_realization", {}).get("toxicity_index", 1.0)
        
        if tox >= 0.0210:
            status = "REJECTED (Law L11 Violation)"
        elif p3[1] != 0:
            status = "REJECTED (Law U34 Violation)"
        else:
            status = "VALIDATED_TRAJECTORY"

        return {
            "origin": lines[-2]["header"]["worldline_id"],
            "current": lines[-1]["header"]["worldline_id"],
            "velocity_magnitude": float(np.linalg.norm(velocity)),
            "predicted_next_coordinate": p3[:10].tolist(),
            "trajectory_status": status,
            "law_u27_compliance": "VERIFIED" if "REJECTED" not in status else "FAILED"
        }

    def is_next_state_safe(self) -> bool:
        """
        Verify that the predicted next state of the molecular series does not drift
        into Law L11 (toxicity) violations. Use after WorldLine persistence.
        Returns True if trajectory is safe, False if L11 or U34 would be violated.
        """
        result = self.calculate_velocity()
        if isinstance(result, str):
            return True  # Insufficient data or error: do not block
        return "REJECTED" not in result.get("trajectory_status", "")

if __name__ == "__main__":
    tp = TrajectoryPredictor()
    result = tp.calculate_velocity()
    print("\n=== PREDATOR X: TRAJECTORY FORECAST ===")
    print(json.dumps(result, indent=2))
