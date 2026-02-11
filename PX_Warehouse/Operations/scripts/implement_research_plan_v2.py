import os
import json
import shutil
from pathlib import Path
import subprocess

# Configuration — canonical warehouse (Prv_Dossiers, Novel_Dossiers, Learning_Material)
_REPO = Path(__file__).resolve().parents[1]
try:
    from PX_Warehouse.warehouse_layout import get_learning_material_dir, get_prv_dossier_dir
    WAREHOUSE_ROOT = _REPO / "PX_Warehouse"
    SOURCE_DIRS = [get_learning_material_dir(_REPO)]
    GOLD_DIR = get_prv_dossier_dir(False, "Gold", _REPO)
    SILVER_DIR = get_prv_dossier_dir(False, "Silver", _REPO)
    LEARNING_DIR = get_learning_material_dir(_REPO)
except Exception:
    WAREHOUSE_ROOT = _REPO / "PX_Warehouse"
    SOURCE_DIRS = [WAREHOUSE_ROOT / "Learning_Material"]
    GOLD_DIR = WAREHOUSE_ROOT / "Prv_Dossiers" / "Gold"
    SILVER_DIR = WAREHOUSE_ROOT / "Prv_Dossiers" / "Silver"
    LEARNING_DIR = WAREHOUSE_ROOT / "Learning_Material"

ORCHESTRATOR_PATH = _REPO / "PX_Executive" / "orchestrators" / "PX_Live_Orchestrator_v2.py"

# Ensure directories exist
for d in [GOLD_DIR, SILVER_DIR, LEARNING_DIR]:
    d.mkdir(parents=True, exist_ok=True)

def get_smiles(data):
    """Extract SMILES from various dossier formats."""
    # Format 1: metadata.smiles
    if "metadata" in data and "smiles" in data["metadata"]:
        return data["metadata"]["smiles"]
    # Format 2: prv_candidate.smiles
    if "prv_candidate" in data and "smiles" in data["prv_candidate"]:
        return data["prv_candidate"]["smiles"]
    # Format 3: molecular_identity.smiles
    if "molecular_identity" in data and "smiles" in data["molecular_identity"]:
        return data["molecular_identity"]["smiles"]
    # Format 4: activities[0].canonical_smiles
    if "data" in data and "activities" in data["data"] and len(data["data"]["activities"]) > 0:
        return data["data"]["activities"][0].get("canonical_smiles")
    return None

def get_name(data):
    """Extract Name from various dossier formats."""
    if "metadata" in data and "name" in data["metadata"]:
        return data["metadata"]["name"]
    if "prv_candidate" in data and "common_name" in data["prv_candidate"]:
        return data["prv_candidate"]["common_name"]
    if "dossier_header" in data and "candidate_id" in data["dossier_header"]:
        return data["dossier_header"]["candidate_id"]
    if "data" in data and "activities" in data["data"] and len(data["data"]["activities"]) > 0:
        return data["data"]["activities"][0].get("molecule_chembl_id")
    return "Unknown"

def grade_dossier(data):
    """Grade dossier based on PTA, Responder Rate, and Toxicity."""
    efficacy = data.get("virtual_efficacy", {})
    responder_rate = efficacy.get("responder_rate", {}).get("response_ratio")
    
    # Try different PTA paths
    pta = efficacy.get("pk_pta", {}).get("auc_mg_h_per_L", {}).get("pta")
    if pta is None:
        pta = data.get("pkpd_analysis", {}).get("pta")
    
    # Toxicity check from consolidate_warehouse.py
    admet = data.get("admet", {})
    tox = admet.get("toxicity", {}).get("toxicity_index", 1.0)

    if responder_rate is None: responder_rate = 0.0
    if pta is None: pta = 0.0

    # GRADING LOGIC from consolidate_warehouse.py
    if tox >= 0.0210:
        return "LEARNING" # REJECTED_TRASH
    elif pta >= 30.0 and responder_rate >= 0.3 and tox < 0.0200:
        return "GOLD"
    elif pta >= 20.0 and responder_rate >= 0.2 and tox < 0.0210:
        return "SILVER"
    else:
        return "LEARNING"

def main():
    candidates = {} # smiles -> name
    files_to_move = [] # (src_path, grade)

    print("--- Phase 1: Data Discovery and Extraction ---")
    for source_dir in SOURCE_DIRS:
        if not source_dir.exists():
            continue
        print(f"Scanning {source_dir}...")
        for json_file in source_dir.rglob("*.json"):
            try:
                if len(files_to_move) % 100 == 0:
                    print(f"Processed {len(files_to_move)} files...")
                with open(json_file, 'r', encoding='utf-8') as f:
                    data = json.load(f)
                
                smiles = get_smiles(data)
                name = get_name(data)
                
                if smiles:
                    candidates[smiles] = name
                
                grade = grade_dossier(data)
                files_to_move.append((json_file, grade))
                
            except Exception as e:
                print(f"Error processing {json_file}: {e}")
                files_to_move.append((json_file, "LEARNING"))

    print(f"Extracted {len(candidates)} unique SMILES candidates.")

    print("\n--- Phase 2: Grading and Sorting ---")
    for src_path, grade in files_to_move:
        dest_dir = GOLD_DIR if grade == "GOLD" else SILVER_DIR if grade == "SILVER" else LEARNING_DIR
        
        # Avoid moving to the same directory
        if src_path.parent == dest_dir:
            continue
            
        dest_path = dest_dir / src_path.name
        
        # Handle filename collisions
        if dest_path.exists():
            dest_path = dest_dir / f"{src_path.stem}_{src_path.parent.name}{src_path.suffix}"
            
        try:
            shutil.move(str(src_path), str(dest_path))
            # print(f"Moved {src_path.name} -> {grade}")
        except Exception as e:
            print(f"Failed to move {src_path}: {e}")
    print("Grading and sorting complete.")

    print("\n--- Phase 3: Pipeline Execution (Parallel) ---")
    # Limit to unique SMILES and run orchestrator in parallel batches
    batch_size = 10
    smiles_list = list(candidates.items())
    for i in range(0, len(smiles_list), batch_size):
        batch = smiles_list[i:i+batch_size]
        processes = []
        for smiles, name in batch:
            # print(f"Executing pipeline for: {name} ({smiles})")
            cmd = [
                "python",
                ORCHESTRATOR_PATH,
                "--smiles", smiles,
                "--name", name,
                "--quiet"
            ]
            processes.append(subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL))
        
        print(f"Waiting for batch {i//batch_size + 1}/{(len(smiles_list)-1)//batch_size + 1}...")
        for p in processes:
            p.wait()

    print("\n--- Phase 4: Validation and Reporting ---")
    print("All phases completed.")

if __name__ == "__main__":
    main()
