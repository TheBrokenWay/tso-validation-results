import os
import json
import subprocess
import sys
import time
from pathlib import Path
from datetime import datetime, UTC

# Configuration
CANDIDATES_FILE = "extracted_candidates.json"
ORCHESTRATOR_PATH = r"E:\foundation\PX_Executive\orchestrators\PX_Live_Orchestrator_v2.py"
SORT_SCRIPT_PATH = r"E:\foundation\PX_Warehouse\Operations\scripts\sort_live_results.py"
LOG_FILE = "one_hour_pipeline_log.txt"
DURATION_SECONDS = 3600 # 1 hour

def run_one_hour_pipeline():
    start_time = time.time()
    end_time = start_time + DURATION_SECONDS
    
    if not Path(CANDIDATES_FILE).exists():
        print(f"Error: {CANDIDATES_FILE} not found.")
        return

    with open(CANDIDATES_FILE, 'r') as f:
        candidates = json.load(f)
    
    candidate_list = list(candidates.items())
    num_candidates = len(candidate_list)
    
    print(f"Starting one-hour pipeline execution at {datetime.now(UTC).isoformat()}")
    print(f"Total candidates available: {num_candidates}")

    with open(LOG_FILE, "a") as log:
        log.write(f"--- One-Hour Pipeline Started: {datetime.now(UTC).isoformat()} ---\n")

    processed_count = 0
    idx = 0
    
    while time.time() < end_time:
        if idx >= num_candidates:
            print("All candidates processed before time limit. Restarting list if time remains...")
            idx = 0 # Loop back if we finish all 957 in under an hour
            
        smiles, name = candidate_list[idx]
        elapsed = time.time() - start_time
        remaining = DURATION_SECONDS - elapsed
        
        print(f"[{processed_count+1}] ({int(elapsed)}s elapsed, {int(remaining)}s left) Processing {name}...")
        
        try:
            # Step 1: Run Orchestrator
            result = subprocess.run(
                [sys.executable, ORCHESTRATOR_PATH, "--smiles", smiles, "--name", name, "--quiet"],
                capture_output=True,
                text=True,
                timeout=300
            )
            
            status = "SUCCESS" if result.returncode == 0 else "FAILED"
            log_msg = f"[{processed_count+1}] {status}: {name} ({smiles}) - {result.stdout.strip() if status == 'SUCCESS' else result.stderr.strip()}\n"
            
            with open(LOG_FILE, "a") as log:
                log.write(log_msg)
            
            if result.returncode == 0:
                print("  SUCCESS. Sorting...")
                # Step 2: Sort the result immediately
                sort_result = subprocess.run(
                    [sys.executable, SORT_SCRIPT_PATH],
                    capture_output=True,
                    text=True
                )
                if sort_result.returncode == 0:
                    print("  SORTED.")
                else:
                    print(f"  SORT FAILED: {sort_result.stderr.strip()}")
            else:
                print(f"  PIPELINE FAILED: {result.stderr.strip()[:100]}...")

        except Exception as e:
            print(f"  ERROR: {e}")
            with open(LOG_FILE, "a") as log:
                log.write(f"[{processed_count+1}] ERROR: {name} ({smiles}) - {str(e)}\n")
        
        processed_count += 1
        idx += 1
        
        # Small sleep to prevent CPU hogging if it runs too fast
        time.sleep(0.5)

    print(f"One-hour pipeline complete. Processed {processed_count} candidates.")
    with open(LOG_FILE, "a") as log:
        log.write(f"--- One-Hour Pipeline Finished: {datetime.now(UTC).isoformat()}. Total processed: {processed_count} ---\n")

if __name__ == "__main__":
    run_one_hour_pipeline()
