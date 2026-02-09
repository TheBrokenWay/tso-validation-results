import os
import json
import re
from pathlib import Path

def extract_from_dir(directory, candidates):
    if not os.path.exists(directory):
        print(f"Directory not found: {directory}")
        return
    
    print(f"Scanning {directory}...")
    # Handle both .json and .md (for summaries)
    files = [f for f in os.listdir(directory) if f.endswith(".json") or f.endswith(".md")]
    for filename in files:
        path = os.path.join(directory, filename)
        try:
            if filename.endswith(".json"):
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
            
            elif filename.endswith(".md"):
                # Extract from summary files
                with open(path, "r", encoding="utf-8") as f:
                    content = f.read()
                
                # Look for SMILES: [SMILES] or similar patterns
                smiles_match = re.search(r"SMILES:\s*([^\s\n]+)", content)
                name_match = re.search(r"#\s*(.+)\s*Summary", content) or re.search(r"Name:\s*(.+)", content)
                
                if smiles_match:
                    smiles = smiles_match.group(1).strip()
                    name = name_match.group(1).strip() if name_match else filename.replace("_Summary.md", "")
                    if smiles not in candidates:
                        candidates[smiles] = name
                
        except Exception as e:
            continue

def main():
    _root = Path(__file__).resolve().parent.parent.parent.parent
    base_path = _root / "PX_Warehouse" / "CommercialAssets"
    folders = ["Gold", "Silver", "Learning_Material", "Executive_Summary"]
    
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
