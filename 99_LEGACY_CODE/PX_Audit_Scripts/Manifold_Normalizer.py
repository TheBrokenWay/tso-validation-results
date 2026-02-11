import json
import os
import numpy as np
from PX_Warehouse.WorldLine_Database import WorldLineDatabase
from PX_Warehouse.Worldline_Indexer import WorldlineIndexer

WAREHOUSE = WorldLineDatabase()
INDEXER = WorldlineIndexer()

def normalize_manifold():
    path = WAREHOUSE.path
    files = [f for f in os.listdir(path) if f.endswith(".worldline")]
    print(f">>> [NORMALIZING] Applying Diversity Injection to {len(files)} items...")

    for f in files:
        file_full_path = os.path.join(path, f)
        with open(file_full_path, "r") as fp:
            try:
                data = json.load(fp)
            except (json.JSONDecodeError, TypeError):
                continue
            
        # FIX 1: Extract block from ROOT (not physics_snapshot)
        block = np.array(data.get("coordinate_35d", [0]*35), dtype=float)
        
        # FIX 2: Diversity Injection Guard (Dims 11-35 only)
        # Prevents infinite density in the KD-Tree
        noise = np.random.normal(0, 0.001, 24)
        block[11:35] = noise
        
        # FIX 3: Preserve Lineage and Metadata from real schema
        task_id = data.get("worldline_id", f).replace("WL-", "")
        coherence = data.get("coherence_amplitude", data.get("coherence", 0.0))
        route = data.get("route", "BLUE")
        cycle = data.get("metabolic_cycle", 0)
        
        # FIX 4: Overwrite with corrected structure
        WAREHOUSE.record_materialization(
            task_id=task_id,
            block=block.tolist(),
            coherence=coherence,
            lab_results=data.get("physical_realization", {}),
            route=route,
            cycle=cycle,
            origin="DIVERSITY_NORMALIZATION",
            parent_id=data.get("parent_id")
        )

    print(">>> [NORMALIZATION COMPLETE] Manifold Shaken.")
    print(">>> [INDEXER] Rebuilding KD-Tree with new diversity...")
    INDEXER.rebuild_index()

if __name__ == "__main__":
    normalize_manifold()
