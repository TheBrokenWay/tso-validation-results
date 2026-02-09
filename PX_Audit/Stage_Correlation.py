import json
import os
import numpy as np

WL_PATH = "E:/foundation/PX_Warehouse/WorldLines"

def load_worldlines():
    records = []
    files = [f for f in os.listdir(WL_PATH) if f.endswith(".worldline")]
    for f in files:
        with open(os.path.join(WL_PATH, f), "r") as fp:
            try:
                data = json.load(fp)
                # Dims 0-9 contain our legacy Stage 1-10 data (support both root and physics_snapshot)
                block = np.array(
                    data.get("coordinate_35d") or data.get("physics_snapshot", {}).get("coordinate_35d", [0]*35),
                    dtype=float
                )
                phys = data.get("physical_realization", {})
                affinity = phys.get("binding_affinity_kj", None)
                if affinity is not None:
                    records.append((block[:10], affinity))
            except Exception:
                continue
    return records

def compute_stage_correlations(records):
    if not records:
        return []
    stages = np.array([r[0] for r in records])   # (N, 10)
    affinities = np.array([r[1] for r in records])

    correlations = []
    for i in range(10):
        stage_vals = stages[:, i]
        # Calculate Pearson Correlation
        if np.std(stage_vals) == 0:
            corr = 0.0
        else:
            corr = np.corrcoef(stage_vals, affinities)[0, 1]
        correlations.append(corr)
    return correlations

def print_heatmap(correlations):
    print("\n=== LEGACY STAGE -> BINDING AFFINITY CORRELATION HEATMAP ===\n")
    print("Stage Index | Correlation | Influence")
    print("--------------------------------------")

    for i, c in enumerate(correlations):
        # A correlation > 0.6 is a strong positive driver; < -0.6 is a strong negative driver
        influence = "HIGH" if abs(c) > 0.6 else "MEDIUM" if abs(c) > 0.3 else "LOW"
        print(f"Stage {i:2d}    | {c:11.4f} | {influence}")

if __name__ == "__main__":
    print(">>> [CORRELATION] Loading 15,382 recovered world-lines...")
    records = load_worldlines()
    print(f">>> [CORRELATION] Loaded {len(records)} valid samples.")
    correlations = compute_stage_correlations(records)
    print_heatmap(correlations)
