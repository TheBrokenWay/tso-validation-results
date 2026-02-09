import json
import os

def run_dossier_audit():
    path = "E:/foundation/PX_Warehouse/WorldLines"
    results = []
    
    files = [f for f in os.listdir(path) if f.endswith(".worldline")]
    print(f">>> [AUDIT] Scanning {len(files)} re-mapped world-lines for peak affinity...")

    for file in files:
        try:
            with open(os.path.join(path, file), "r") as f:
                data = json.load(f)
                
                # Extracting affinity from physical_realization
                # Samples show binding_affinity_kj as the key
                phys = data.get("physical_realization", {})
                affinity = phys.get("binding_affinity_kj", 0)
                
                if affinity > 0:
                    results.append({
                        "id": data["header"]["worldline_id"],
                        "affinity": affinity,
                        "coherence": data["physics_snapshot"]["coherence"],
                        "origin": data["physics_snapshot"]["origin_context"]
                    })
        except Exception:
            continue

    # Sort by affinity descending
    results.sort(key=lambda x: x["affinity"], reverse=True)

    print("\n=== TOP 10 RE-MAPPED CANDIDATES BY BINDING AFFINITY ===")
    print(f"{'WorldLine ID':<30} | {'Affinity (kJ)':<15} | {'Coherence'}")
    print("-" * 60)
    
    for item in results[:10]:
        print(f"{item['id']:<30} | {item['affinity']:<15.2f} | {item['coherence']}")

if __name__ == "__main__":
    run_dossier_audit()
