import json
import os

def explore_manifold():
    path = "PX_Warehouse/WorldLines"
    candidates = []
    for file in os.listdir(path):
        if file.endswith(".worldline"):
            with open(os.path.join(path, file), "r") as f:
                data = json.load(f)
                phys = data["physical_realization"]
                candidates.append({
                    "id": data["worldline_id"],
                    "coherence": data["coherence_amplitude"],
                    "affinity": phys["binding_affinity_kj"],
                    "tox": phys.get("toxicity_index", 0.02)
                })

    ranked = sorted(candidates, key=lambda x: x["affinity"], reverse=True)
    
    print("\n=== PREDATOR X: GLOBAL MANIFOLD RANKING ===")
    print(f"{'WORLD-LINE ID':<25} | {'COHERENCE':<10} | {'AFFINITY':<10} | {'TOXICITY':<10}")
    print("-" * 65)
    for c in ranked:
        print(f"{c['id']:<25} | {c['coherence']:<10.4f} | {c['affinity']:<10.2f} | {c['tox']:<10.4f}")

if __name__ == "__main__":
    explore_manifold()
