"""
PX_Executive/harvest_leads.py
Harvests 'Silver Tier' candidates from Discovery_Accepted for Lead Optimization.
"""
import json
import glob
import os
from pathlib import Path

# Paths
REPO_ROOT = Path(__file__).resolve().parents[1]
DISCOVERY_DIR = REPO_ROOT / "PX_Warehouse" / "Finalized_Dossiers" / "Discovery_Accepted"
FEEDER_FILE = REPO_ROOT / "PX_Warehouse" / "Feeder" / "reprocess_candidates.json"

def harvest():
    print(f"--- HARVESTING LEADS FROM: {DISCOVERY_DIR} ---")
    
    candidates_to_refine = []
    
    # 1. Scan the Discovery Folder
    files = list(DISCOVERY_DIR.glob("*.json"))
    print(f"Scanning {len(files)} discovery dossiers...")

    for f in files:
        try:
            data = json.loads(f.read_text(encoding='utf-8'))
            
            # Extract Key Metrics
            candidate = data.get("candidate", {})
            smiles = candidate.get("smiles")
            name = candidate.get("name", "Unknown")
            
            # Get Toxicity (The thing we need to fix)
            zeus = data.get("finalization", {}).get("zeus_verdict", {})
            tox = zeus.get("toxicity_index", 1.0)
            
            # FILTER: Only harvest the "Close Calls" (Toxicity 0.021 - 0.045)
            # If it's > 0.045, it might be too toxic to save.
            if 0.021 <= tox <= 0.05:
                print(f"  [HARVEST] {name} (Tox: {tox:.4f}) -> Queued for Optimization")
                
                candidates_to_refine.append({
                    "source_id": f.stem,
                    "name": f"{name}_v2", # Mark as version 2
                    "smiles": smiles,
                    "target_action": "OPTIMIZE_TOXICITY", # Instruction to the engine
                    "previous_tox": tox
                })
                
        except Exception as e:
            print(f"  [ERROR] Could not read {f.name}: {e}")

    # 2. Write to Feeder
    if candidates_to_refine:
        # Load existing queue if it exists
        current_queue = []
        if FEEDER_FILE.exists():
            try:
                current_queue = json.loads(FEEDER_FILE.read_text())
            except:
                current_queue = []
        
        # Merge (avoid duplicates based on SMILES)
        existing_smiles = {c.get("smiles") for c in current_queue}
        added_count = 0
        for c in candidates_to_refine:
            if c["smiles"] not in existing_smiles:
                current_queue.append(c)
                added_count += 1
        
        # Save
        FEEDER_FILE.parent.mkdir(parents=True, exist_ok=True)
        FEEDER_FILE.write_text(json.dumps(current_queue, indent=2))
        print(f"--- HARVEST COMPLETE ---")
        print(f"Added {added_count} candidates to the Reprocessing Queue.")
        print(f"Total Queue Size: {len(current_queue)}")
        print(f"Ready for: 'python PX_Executive/run_repurposing.py'")
    else:
        print("No candidates met the harvest criteria (Tox 0.021-0.05).")

if __name__ == "__main__":
    harvest()