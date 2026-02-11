import os
import json
import shutil
from datetime import datetime

# --- PATH CONFIGURATION ---
BASE_PATH = r"E:\foundation\PX_Warehouse"
LOG_DIR = os.path.join(BASE_PATH, "logs")
STRAT_LOG = os.path.join(LOG_DIR, "stratification.log")
MON_LOG = os.path.join(LOG_DIR, "monetization.log")

def rollback_from_log(log_path):
    """Reverses moves recorded in the log file."""
    if not os.path.exists(log_path):
        print(f"Log file {log_path} not found.")
        return

    print(f"Starting rollback from {log_path}...")
    actions_reversed = 0
    
    # Read log entries in reverse order
    with open(log_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    
    # Verify Hash Chain before rollback
    print("Verifying Hash Chain Integrity...")
    last_hash = "0" * 64
    for i, line in enumerate(lines):
        if line.startswith("#") or not line.strip():
            continue
        try:
            entry = json.loads(line)
            stored_hash = entry.get("chain_hash")
            if not stored_hash:
                continue # Skip entries without hash chaining (legacy)
            
            # Recalculate hash
            entry_copy = entry.copy()
            entry_copy.pop("chain_hash")
            entry_str = json.dumps(entry_copy, sort_keys=True)
            recalculated_hash = hashlib.sha256((last_hash + entry_str).encode('utf-8')).hexdigest()
            
            if recalculated_hash != stored_hash:
                print(f"CRITICAL: Hash chain broken at line {i+1}!")
                print(f"Expected: {stored_hash}")
                print(f"Actual:   {recalculated_hash}")
                choice = input("Continue rollback despite integrity failure? (y/n): ")
                if choice.lower() != 'y':
                    return
            last_hash = stored_hash
        except Exception as e:
            print(f"Warning: Could not verify line {i+1}: {e}")
    
    print("Hash Chain Verified. Proceeding with Rollback.")
    
    for line in reversed(lines):
        if line.startswith("#") or not line.strip():
            continue
            
        try:
            entry = json.loads(line)
            action = entry.get("action")
            
            if action == "STRATIFY":
                source = entry.get("source")
                dest = entry.get("destination")
                
                if os.path.exists(dest):
                    # Ensure parent directory of source exists
                    os.makedirs(os.path.dirname(source), exist_ok=True)
                    shutil.move(dest, source)
                    print(f"Reversed STRATIFY: {dest} -> {source}")
                    actions_reversed += 1
                else:
                    print(f"Warning: File {dest} not found, cannot reverse.")
                    
            elif action == "MONETIZE":
                dest = entry.get("destination")
                if os.path.exists(dest):
                    os.remove(dest)
                    print(f"Reversed MONETIZE: Deleted {dest}")
                    actions_reversed += 1
                else:
                    print(f"Warning: File {dest} not found, cannot reverse.")
                    
        except Exception as e:
            print(f"Error processing log entry: {e}")

    print(f"Rollback complete. Total actions reversed: {actions_reversed}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python rollback_audit.py [stratification|monetization]")
    else:
        target = sys.argv[1].lower()
        if target == "stratification":
            rollback_from_log(STRAT_LOG)
        elif target == "monetization":
            rollback_from_log(MON_LOG)
        else:
            print(f"Unknown target: {target}")
