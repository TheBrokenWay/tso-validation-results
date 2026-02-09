import json
import os
import sys
import numpy as np

# Ensure root is in path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from PX_Warehouse.WorldLine_Database import WorldLineDatabase

WAREHOUSE = WorldLineDatabase()

def normalize_manifold():
    path = WAREHOUSE.path
    if not os.path.exists(path):
        os.makedirs(path)
        print(f">>> [INIT] Created WorldLines directory: {path}")
        
    files = [f for f in os.listdir(path) if f.endswith(".worldline")]
    print(f">>> [NORMALIZING] Applying Diversity Injection to {len(files)} items...")

    if len(files) == 0:
        print("    (No files found to normalize, skipping...)")

    for f in files:
        file_full_path = os.path.join(path, f)
        try:
            with open(file_full_path, "r") as fp:
                data = json.load(fp)
            
            # FIX 1: Extract block from ROOT (not physics_snapshot)
            block = np.array(data.get("coordinate_35d", [0]*35), dtype=float)
            
            # FIX 2: Diversity Injection Guard (Dims 11-35 only)
            noise = np.random.normal(0, 0.001, 24)
            if len(block) >= 35:
                block[11:35] = noise
            
            # FIX 3: Preserve Metadata
            task_id = data.get("worldline_id", f).replace("WL-", "").replace(".worldline", "")
            coherence = data.get("coherence_amplitude", data.get("coherence", 0.99))
            
            # FIX 4: Overwrite
            WAREHOUSE.record_materialization(
                task_id=task_id,
                block=block.tolist(),
                coherence=coherence,
                lab_results=data.get("physical_realization", {}),
                route="BLUE",
                origin="DIVERSITY_NORMALIZATION"
            )
        except Exception as e:
            print(f"Skipping {f}: {e}")
            continue
    print(">>> [SUCCESS] Manifold Normalization Complete.")

if __name__ == "__main__":
    normalize_manifold()
