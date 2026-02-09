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
    """Grade dossier based on PTA, Responder Rate, and Toxicity."""
    # Try to find efficacy metrics
    # In Evidence Package v3 (Live Orchestrator output), the structure is different
    # virtual_efficacy.pk_pta.auc_mg_h_per_L.pta
    # virtual_efficacy.responder_rate.response_ratio
    
    ve = data.get("virtual_efficacy", {})
    pta = ve.get("pk_pta", {}).get("auc_mg_h_per_L", {}).get("pta", 0.0)
    responder_rate = ve.get("responder_rate", {}).get("response_ratio", 0.0)
    
    # Toxicity in admet.toxicity.toxicity_index
    admet = data.get("admet", {})
    tox = admet.get("toxicity", {}).get("toxicity_index", 1.0)

    # Grading logic
    if tox >= 0.0210:
        return "LEARNING"
    elif pta >= 30.0 and responder_rate >= 0.3 and tox < 0.0200:
        return "GOLD"
    elif pta >= 20.0 and responder_rate >= 0.2 and tox < 0.0210:
        return "SILVER"
    else:
        return "LEARNING"

def sort_live_runs():
    print(f"Sorting Live Runs from {LIVE_RUNS_DIR}...")
    
    count = 0
    # Live runs are in subfolders: run_YYYYMMDD_HHMMSS/TRIAL_SIMULATION_DOSSIER-runid.json
    for dossier_file in LIVE_RUNS_DIR.glob("**/TRIAL_SIMULATION_DOSSIER-*.json"):
        try:
            with open(dossier_file, 'r', encoding='utf-8') as f:
                data = json.load(f)
            
            grade = grade_dossier(data)
            
            dest_dir = GOLD_DIR if grade == "GOLD" else SILVER_DIR if grade == "SILVER" else LEARNING_DIR
            dest_path = dest_dir / dossier_file.name
            
            # Handle filename collisions
            if dest_path.exists():
                dest_path = dest_dir / f"{dossier_file.stem}_{dossier_file.parent.name}{dossier_file.suffix}"
            
            shutil.copy(str(dossier_file), str(dest_path))
            count += 1
            
            if count % 100 == 0:
                print(f"Sorted {count} dossiers...")
                
        except Exception as e:
            print(f"Error sorting {dossier_file}: {e}")

    print(f"Completed sorting {count} dossiers.")

if __name__ == "__main__":
    sort_live_runs()
