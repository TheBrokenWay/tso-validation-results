import json
import subprocess
import sys
from pathlib import Path
from datetime import datetime, UTC

# Configuration
CANDIDATES_FILE = "extracted_candidates.json"
ORCHESTRATOR_PATH = r"E:\foundation\PX_Executive\orchestrators\PX_Live_Orchestrator_v2.py"
LOG_FILE = "pipeline_execution_log.txt"

def run_pipeline():
    if not Path(CANDIDATES_FILE).exists():
        print(f"Error: {CANDIDATES_FILE} not found.")
        return

    with open(CANDIDATES_FILE, 'r') as f:
        candidates = json.load(f)

    print(f"Starting pipeline for {len(candidates)} candidates...")
    
    # We will process in batches or just a subset if it's too many, 
    # but the plan says "all unique SMILES found".
    # 957 is a lot, let's do them all but log progress.

    with open(LOG_FILE, "a") as log:
        log.write(f"--- Pipeline Execution Started: {datetime.now(UTC).isoformat()} ---\n")

    for i, (smiles, name) in enumerate(candidates.items(), 1):
        print(f"[{i}/{len(candidates)}] Processing {name}...")
        
        try:
            # Execute orchestrator
            result = subprocess.run(
                [sys.executable, ORCHESTRATOR_PATH, "--smiles", smiles, "--name", name, "--quiet"],
                capture_output=True,
                text=True,
                timeout=300 # 5 min timeout per candidate
            )
            
            status = "SUCCESS" if result.returncode == 0 else "FAILED"
            log_msg = f"[{i}] {status}: {name} ({smiles}) - {result.stdout.strip() if status == 'SUCCESS' else result.stderr.strip()}\n"
            
            with open(LOG_FILE, "a") as log:
                log.write(log_msg)
                
            if result.returncode != 0:
                print(f"  FAILED: {result.stderr.strip()[:100]}...")
            else:
                print(f"  SUCCESS")

        except subprocess.TimeoutExpired:
            print(f"  TIMEOUT: {name}")
            with open(LOG_FILE, "a") as log:
                log.write(f"[{i}] TIMEOUT: {name} ({smiles})\n")
        except Exception as e:
            print(f"  ERROR: {e}")
            with open(LOG_FILE, "a") as log:
                log.write(f"[{i}] ERROR: {name} ({smiles}) - {str(e)}\n")

    print(f"Pipeline execution complete. Results logged to {LOG_FILE}")

if __name__ == "__main__":
    run_pipeline()
