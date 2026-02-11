import os
import json
import subprocess
import time
from datetime import datetime, UTC
import sys
from pathlib import Path
import re
import shutil

# Add parent directory to path to import sorting logic
_root = Path(__file__).resolve().parent.parent.parent.parent
sys.path.insert(0, str(_root))
from PX_Warehouse.consolidate_warehouse import grade_and_sort_asset
from PX_Warehouse.warehouse_layout import get_queue_path

def run_reprocess_pipeline():
    candidate_file = get_queue_path("reprocess_candidates.json", _root)
    log_file = _root / "PX_LOGS" / "archive" / "reprocess_pipeline_final_log.txt"
    python_exe = sys.executable
    script_path = _root / "PX_Executive" / "orchestrators" / "PX_Live_Orchestrator_v2.py"
    
    candidate_file = str(candidate_file)
    log_file = str(log_file)
    if not os.path.exists(candidate_file):
        print(f"Candidate file not found: {candidate_file}")
        return
    Path(log_file).parent.mkdir(parents=True, exist_ok=True)

    with open(candidate_file, "r") as f:
        candidates = json.load(f)

    print(f"Starting reprocessing of {len(candidates)} candidates...")
    
    processed_count = 0
    
    # Clear log for fresh start
    with open(log_file, "w") as log:
        log.write(f"--- REPROCESS START: {datetime.now(UTC).isoformat()} ---\n")
        log.write(f"Python: {python_exe}\n")
        log.write(f"Script: {script_path}\n")
    
    for smiles, name in candidates.items():
        processed_count += 1
        
        try:
            # Use a string for the command and shell=True to let Windows handle it
            # Escape quotes in smiles and name
            safe_smiles = smiles.replace('"', '""')
            safe_name = name.replace('"', '""')
            cmd_str = f'"{python_exe}" "{script_path}" --smiles "{safe_smiles}" --name "{safe_name}"'
            
            result = subprocess.run(cmd_str, capture_output=True, text=True, shell=True)
            
            with open(log_file, "a") as log:
                if result.returncode == 0:
                    match = re.search(r"Dossier: (.*\.json)", result.stdout)
                    if match:
                        dossier_path = match.group(1).strip()
                        if os.path.exists(dossier_path):
                            sort_result = grade_and_sort_asset(dossier_path)
                            log.write(f"SUCCESS: {name} | Sorted: {sort_result}\n")
                        else:
                            log.write(f"ERROR: {name} | Dossier path found in output but file missing: {dossier_path}\n")
                    else:
                        log.write(f"ERROR: {name} | Pipeline succeeded but dossier path not found in output\n")
                else:
                    log.write(f"ERROR: {name} | Pipeline failed with exit code {result.returncode}\n")
                    log.write(f"STDERR: {result.stderr[-200:]}\n")
            
        except Exception as e:
            with open(log_file, "a") as log:
                log.write(f"CRITICAL: {name} | Exception: {str(e)}\n")
        
        time.sleep(0.1)

if __name__ == "__main__":
    run_reprocess_pipeline()
