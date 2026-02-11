import json
import subprocess
import os
from pathlib import Path
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

CANDIDATES_FILE = "extracted_candidates.json"
ORCHESTRATOR_PATH = r"E:\foundation\PX_Executive\orchestrators\PX_Live_Orchestrator_v2.py"
LOG_FILE = "batch_pipeline_execution.log"
PROGRESS_FILE = "batch_progress.json"
MAX_WORKERS = 4  # Adjust based on CPU cores

def load_candidates():
    with open(CANDIDATES_FILE, "r") as f:
        return json.load(f)

def load_progress():
    if os.path.exists(PROGRESS_FILE):
        with open(PROGRESS_FILE, "r") as f:
            return json.load(f)
    return {"completed": [], "failed": []}

def save_progress(progress):
    with open(PROGRESS_FILE, "w") as f:
        json.dump(progress, f, indent=2)

def run_single_candidate(smiles, name):
    try:
        cmd = [
            "python",
            ORCHESTRATOR_PATH,
            "--smiles", smiles,
            "--name", name,
            "--quiet"
        ]
        
        start_time = time.time()
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        duration = time.time() - start_time
        
        return smiles, name, True, duration, None
                
    except subprocess.CalledProcessError as e:
        return smiles, name, False, 0, f"STDOUT: {e.stdout}\nSTDERR: {e.stderr}"
    
    except Exception as e:
        return smiles, name, False, 0, str(e)

def run_batch():
    candidates = load_candidates()
    progress = load_progress()
    
    total = len(candidates)
    
    # Filter out already processed
    to_process = [(smiles, name) for smiles, name in candidates.items() 
                  if smiles not in progress["completed"] and smiles not in progress["failed"]]
    
    processed_count = len(progress["completed"]) + len(progress["failed"])
    
    print(f"Starting batch execution with {MAX_WORKERS} workers: {processed_count}/{total} already processed.")
    
    if not to_process:
        print("All candidates already processed.")
        return

    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        future_to_candidate = {executor.submit(run_single_candidate, smiles, name): (smiles, name) 
                               for smiles, name in to_process}
        
        for future in as_completed(future_to_candidate):
            smiles, name, success, duration, error = future.result()
            processed_count += 1
            
            if success:
                print(f"[{processed_count}/{total}] Success: {name} ({duration:.2f}s)")
                progress["completed"].append(smiles)
                with open(LOG_FILE, "a") as f:
                    f.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} - SUCCESS - {name} - {smiles} - {duration:.2f}s\n")
            else:
                print(f"[{processed_count}/{total}] FAILED: {name} - {error[:100]}...")
                progress["failed"].append(smiles)
                with open(LOG_FILE, "a") as f:
                    f.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} - FAILED - {name} - {smiles}\n{error}\n")
            
            # Save progress every 10 candidates to avoid losing much work
            if processed_count % 10 == 0:
                save_progress(progress)
                
            if processed_count % 100 == 0:
                print(f"--- Processed {processed_count}/{total} ---")

    save_progress(progress)
    print(f"Batch execution complete. Total processed: {processed_count}/{total}")

if __name__ == "__main__":
    run_batch()
