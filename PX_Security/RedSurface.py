import hashlib
from datetime import datetime, timezone

class RedSurface:
    """Predator X: Adversarial Decoy Layer"""
    def __init__(self):
        self.metadata = {"role": "RedSurface", "version": "1.2.0-GAIP"}

    def generate_decoy_response(self, task_id):
        # Generate high-amplitude but biologically impossible fake data
        # Deterministic amplitude derived from task_id hash
        hash_val = int(hashlib.sha256(task_id.encode()).hexdigest(), 16)
        amplitude = 0.96 + (hash_val % 300) / 10000.0 # Range 0.9600 - 0.9899
        
        return {
            "task_id": task_id,
            "authorized": True,
            "amplitude": float(f"{amplitude:.8f}"), # High precision to prevent pleasing rounding
            "rationale": "Decoy Authorization Active",
            "decoy_flag": True,
            "timestamp": datetime.now(timezone.utc).isoformat()
        }
