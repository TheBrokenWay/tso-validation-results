import json
import os
import numpy as np
from collections import Counter

def find_common_denominator():
    path = "E:/foundation/PX_Warehouse/WorldLines"
    all_coords = []
    
    files = [f for f in os.listdir(path) if f.endswith(".worldline")]
    print(f">>> [PROBE] Extracting DNA from {len(files)} scrubbed items...")

    for file in files:
        with open(os.path.join(path, file), "r") as f:
            data = json.load(f)
            c = data.get("physics_snapshot", {}).get("coordinate_35d", [])
            if c: all_coords.append(c)

    if not all_coords:
        print("[!] No data found to analyze.")
        return

    data_matrix = np.array(all_coords)
    
    # 1. Coordinate Center (Where is the pile located?)
    mean_vec = np.mean(data_matrix, axis=0)
    
    # 2. Variance (Is the pile tight or scattered?)
    variance_vec = np.var(data_matrix, axis=0)
    
    # 3. Dimensional Dominance (Which dimensions are "locked"?)
    # We look for dimensions with near-zero variance
    locked_dims = [i for i, v in enumerate(variance_vec) if v < 0.001]

    print("\n=== SYSTEM PROBE: THE COMMON DENOMINATOR ===")
    print(f"Global Center (Dim 0-4): {mean_vec[:5]}")
    print(f"Average Variance: {np.mean(variance_vec):.6f}")
    print(f"Locked Dimensions (Static DNA): {locked_dims}")
    
    # Check for "Zero-Value" poisoning
    zero_count = np.count_nonzero(mean_vec == 0)
    print(f"Dead Dimensions (Zeroed): {zero_count} / 35")

if __name__ == "__main__":
    find_common_denominator()
