import os
import json
import re

def extract_from_dir(directory, candidates):
    if not os.path.exists(directory):
        print(f"Directory not found: {directory}")
        return
    
    print(f"Scanning {directory}...")
    files = [f for f in os.listdir(directory) if f.endswith(".json")]
    for filename in files:
        path = os.path.join(directory, filename)
        try:
            with open(path, "r", encoding="utf-8") as f:
                data = json.load(f)
            
            smiles = None
            name = None
            
            # Try different common locations for SMILES and Name
            if "metadata" in data:
                smiles = data["metadata"].get("smiles")
                name = data["metadata"].get("name")
            elif "smiles" in data:
                smiles = data.get("smiles")
                name = data.get("name") or data.get("chembl_id")
            elif "molecule_structures" in data: # ChEMBL format
                smiles = data["molecule_structures"].get("canonical_smiles")
                name = data.get("pref_name") or data.get("molecule_chembl_id")
            
            if smiles and smiles not in candidates:
                candidates[smiles] = name or "Unknown_Candidate"
                
        except Exception as e:
            # Skip files that aren't valid JSON or don't have the data
            continue

def main():
    from pathlib import Path
    _root = Path(__file__).resolve().parent.parent.parent.parent
    base_path = _root / "PX_Warehouse" / "CommercialAssets"
    folders = ["Gold", "Silver", "Learning_Material"]
    
    all_candidates = {}
    
    for folder in folders:
        extract_from_dir(str(base_path / folder), all_candidates)
    
    output_path = _root / "PX_Warehouse" / "Feeder" / "reprocess_candidates.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(all_candidates, f, indent=2)
    
    print(f"\nExtracted {len(all_candidates)} unique candidates for reprocessing.")
    print(f"Saved to {output_path}")

if __name__ == "__main__":
    main()
