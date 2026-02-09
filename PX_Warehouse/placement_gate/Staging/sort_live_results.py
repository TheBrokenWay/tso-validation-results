import os
import json
import shutil
from pathlib import Path

# Configuration
LIVE_RUNS_DIR = Path(r"E:\foundation\PX_Warehouse\TrialSimulations\LiveRuns")
WAREHOUSE_ROOT = Path(r"E:\foundation\PX_Warehouse\CommercialAssets")
GOLD_DIR = WAREHOUSE_ROOT / "Gold"
SILVER_DIR = WAREHOUSE_ROOT / "Silver"
LEARNING_DIR = WAREHOUSE_ROOT / "Learning_Material"

# Ensure directories exist
for d in [GOLD_DIR, SILVER_DIR, LEARNING_DIR]:
    d.mkdir(parents=True, exist_ok=True)

def grade_dossier(data):
    """Grade dossier based on PTA and Responder Rate."""
    efficacy = data.get("virtual_efficacy", {})
    responder_rate = efficacy.get("responder_rate", {}).get("response_ratio", 0.0)
    
    # Try different PTA paths
    pta = efficacy.get("pk_pta", {}).get("auc_mg_h_per_L", {}).get("pta", 0.0)
    
    # Grading logic from consolidate_warehouse.py
    if pta >= 30.0 and responder_rate >= 0.3:
        return "GOLD"
    elif pta >= 20.0 and responder_rate >= 0.2:
        return "SILVER"
    else:
        return "LEARNING"

def sort_live_runs():
    print(f"Scanning {LIVE_RUNS_DIR} for results...")
    
    count = 0
    for run_dir in LIVE_RUNS_DIR.iterdir():
        if not run_dir.is_dir():
            continue
            
        for dossier_file in run_dir.glob("TRIAL_SIMULATION_DOSSIER-*.json"):
            try:
                with open(dossier_file, 'r', encoding='utf-8') as f:
                    data = json.load(f)
                
                grade = grade_dossier(data)
                dest_dir = GOLD_DIR if grade == "GOLD" else SILVER_DIR if grade == "SILVER" else LEARNING_DIR
                
                # We copy instead of move to keep the LiveRuns record intact
                dest_path = dest_dir / dossier_file.name
                if dest_path.exists():
                    dest_path = dest_dir / f"{dossier_file.stem}_{run_dir.name}{dossier_file.suffix}"
                
                shutil.copy(str(dossier_file), str(dest_path))
                count += 1
                
            except Exception as e:
                print(f"Error sorting {dossier_file}: {e}")

    print(f"Sorted {count} new dossiers into pillars.")

if __name__ == "__main__":
    sort_live_runs()
