import os
import sys
import json
import time

ROOT_DIR = "E:/foundation"
if ROOT_DIR not in sys.path:
    sys.path.append(ROOT_DIR)

def load_warehouse_data():
    try:
        from PX_Warehouse.WorldLine_Database import DEFAULT_WORLDLINES_PATH
        path = DEFAULT_WORLDLINES_PATH
    except Exception:
        path = os.path.join(ROOT_DIR, "PX_Warehouse", "WorldLines")
    files = [f for f in os.listdir(path) if f.endswith(".worldline")]
    data = []
    for f in files:
        try:
            with open(os.path.join(path, f), "r") as fp:
                d = json.load(fp)
                # Handle legacy vs new schema
                phys = d.get("physics_snapshot", {})
                real = d.get("physical_realization", {})
                
                # Check for direct keys if schema is flat
                coh = phys.get("coherence", d.get("coherence_amplitude", 0))
                tox = real.get("toxicity_index", d.get("toxicity_index", 1.0))
                
                data.append({"id": f, "coh": coh, "tox": tox})
        except (OSError, json.JSONDecodeError, KeyError, TypeError):
            continue
    return data

def render_dashboard():
    print("\n=== OLYMPUS RESEARCH DASHBOARD ===")
    print(f"Time: {time.strftime('%H:%M:%S')} | FDA-GAIP-2026 Status: ACTIVE\n")
    
    data = load_warehouse_data()
    golden_threshold = 0.0200
    
    print(f"{'CANDIDATE ID':<30} | {'COHERENCE':<10} | {'TOXICITY':<10} | {'STATUS'}")
    print("-" * 70)
    
    # Sort by Toxicity (Best to Worst)
    sorted_data = sorted(data, key=lambda x: x["tox"])
    
    for item in sorted_data[:15]: # Top 15
        marker = " "
        status = "STANDARD"
        
        if item["tox"] < golden_threshold:
            marker = "*"
            status = "GOLDEN"
        elif item["coh"] > 1.0:
            status = "OVERDRIVE"
            
        print(f"{marker} {item['id']:<28} | {item['coh']:.4f}     | {item['tox']:.4f}     | {status}")
        
    print("-" * 70)
    print(f"Total Artifacts: {len(data)}")

if __name__ == "__main__":
    render_dashboard()
