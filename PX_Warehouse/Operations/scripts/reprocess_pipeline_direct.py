import os
import json
import time
from datetime import datetime, UTC
import sys
from pathlib import Path

# Repo root (this file is under PX_Warehouse/)
_root = Path(__file__).resolve().parents[1]
if str(_root) not in sys.path:
    sys.path.insert(0, str(_root))
from PX_Executive.orchestrators.PX_Live_Orchestrator_v2 import PredatorXOrchestratorV2
from PX_Warehouse.consolidate_warehouse import grade_and_sort_asset, ensure_architecture
from PX_Warehouse.warehouse_layout import get_queue_path

def run_reprocess_pipeline():
    # Ensure warehouse architecture exists
    ensure_architecture()

    _root = Path(__file__).resolve().parents[1]
    candidate_file = get_queue_path("reprocess_candidates.json", _root)
    log_file = _root / "PX_LOGS" / "archive" / "reprocess_pipeline_direct_log.txt"
    
    candidate_file = str(candidate_file)
    log_file = str(log_file)
    Path(log_file).parent.mkdir(parents=True, exist_ok=True)
    if not os.path.exists(candidate_file):
        print(f"Candidate file not found: {candidate_file}")
        return

    with open(candidate_file, "r") as f:
        candidates = json.load(f)

    print(f"Starting DIRECT reprocessing of {len(candidates)} candidates...")
    
    processed_count = 0
    
    # Initialize orchestrator once
    orchestrator = PredatorXOrchestratorV2(verbose=False)
    
    # Clear log for fresh start
    with open(log_file, "w") as log:
        log.write(f"--- DIRECT REPROCESS START: {datetime.now(UTC).isoformat()} ---\n")
    
    for smiles, name in candidates.items():
        processed_count += 1
        
        try:
            # 1. Run the pipeline directly
            metadata = {
                "name": name,
                "id": name # Use name as ID for now
            }
            
            # This is synchronous
            results = orchestrator.run_pipeline(smiles, metadata)
            
            dossier_path = results.get("dossier_path")
            
            with open(log_file, "a") as log:
                if dossier_path and os.path.exists(dossier_path):
                    # 2. Sort the result immediately
                    sort_result = grade_and_sort_asset(dossier_path)
                    log.write(f"SUCCESS: {name} | Sorted: {sort_result}\n")
                else:
                    log.write(f"ERROR: {name} | Pipeline finished but dossier missing\n")
            
        except Exception as e:
            import traceback
            with open(log_file, "a") as log:
                log.write(f"CRITICAL: {name} | Exception: {str(e)}\n")
                log.write(traceback.format_exc())
        
        # Small delay
        time.sleep(0.01)

if __name__ == "__main__":
    run_reprocess_pipeline()
