import os
import json
import shutil
from pathlib import Path

# Configuration
LEARNING_MATERIAL_DIR = Path(r"E:\foundation\PX_Warehouse\CommercialAssets\Learning_Material")
WAREHOUSE_ROOT = Path(r"E:\foundation\PX_Warehouse\CommercialAssets")
GOLD_DIR = WAREHOUSE_ROOT / "Gold"
SILVER_DIR = WAREHOUSE_ROOT / "Silver"
LEARNING_DIR = WAREHOUSE_ROOT / "Learning_Material"

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
    return None

def get_name(data):
    """Extract Name from various dossier formats."""
    if "metadata" in data and "name" in data["metadata"]:
        return data["metadata"]["name"]
    if "prv_candidate" in data and "common_name" in data["prv_candidate"]:
        return data["prv_candidate"]["common_name"]
    if "dossier_header" in data and "candidate_id" in data["dossier_header"]:
        return data["dossier_header"]["candidate_id"]
    return "Unknown"

def grade_dossier(data):
    """Grade dossier based on PTA, Responder Rate, and Toxicity."""
    # Try to find efficacy metrics
    efficacy = data.get("virtual_efficacy", {})
    responder_rate = efficacy.get("responder_rate", {}).get("response_ratio", 0.0)
    
    # Try different PTA paths
    pta = efficacy.get("pk_pta", {}).get("auc_mg_h_per_L", {}).get("pta")
    if pta is None:
        pta = data.get("pkpd_analysis", {}).get("pta", 0.0)

    # Robustly find toxicity (handle different schema versions)
    admet = data.get("admet", {})
    tox = admet.get("toxicity", {}).get("toxicity_index")
    if tox is None:
        # Check alternative path
        tox = data.get("stages", {}).get("admet", {}).get("toxicity", {}).get("toxicity_index", 0.0)

    # Grading logic from consolidate_warehouse.py
    # TRASH: Toxicity index >= 0.0210 MUST be discarded
    if tox >= 0.0210:
        return "LEARNING"
    # GOLD: High efficacy, commercially viable PTA (>30%), and SAFE (tox < 0.0200)
    elif pta >= 30.0 and responder_rate >= 0.3 and tox < 0.0200:
        return "GOLD"
    # SILVER: Moderate efficacy, acceptable PTA (>20%), and ACCEPTABLE (tox < 0.0210)
    elif pta >= 20.0 and responder_rate >= 0.2 and tox < 0.0210:
        return "SILVER"
    else:
        return "LEARNING"

def process_files():
    candidates = {} # smiles -> name
    files_to_move = [] # (src_path, grade)

    print(f"Scanning {LEARNING_MATERIAL_DIR}...")
    
    for json_file in LEARNING_MATERIAL_DIR.rglob("*.json"):
        try:
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

    # Sorting
    print(f"Sorting {len(files_to_move)} files...")
    for src_path, grade in files_to_move:
        dest_dir = GOLD_DIR if grade == "GOLD" else SILVER_DIR if grade == "SILVER" else LEARNING_DIR
        dest_path = dest_dir / src_path.name
        
        # Handle filename collisions
        if dest_path.exists():
            dest_path = dest_dir / f"{src_path.stem}_{src_path.parent.name}{src_path.suffix}"
            
        try:
            shutil.move(str(src_path), str(dest_path))
        except Exception as e:
            print(f"Failed to move {src_path}: {e}")

    # Save candidates for next phase
    with open("extracted_candidates.json", "w") as f:
        json.dump(candidates, f, indent=2)
    
    print(f"Extracted {len(candidates)} unique SMILES candidates.")
    print("Sorting complete.")

if __name__ == "__main__":
    process_files()
