import os
import shutil
from datetime import datetime

# --- PATH CONFIGURATION ---
BASE_PATH = r"E:\foundation\PX_Warehouse"
COMMERCIAL_DIR = os.path.join(BASE_PATH, "CommercialAssets")
ARCHIVE_ROOT = os.path.join(BASE_PATH, "Archive")

def archive_state():
    """Creates a timestamped snapshot of the CommercialAssets directory."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    archive_dir = os.path.join(ARCHIVE_ROOT, f"Snapshot_{timestamp}")
    
    if not os.path.exists(COMMERCIAL_DIR):
        print(f"Commercial directory {COMMERCIAL_DIR} not found. Skipping archive.")
        return

    print(f"Creating state archive: {archive_dir}...")
    try:
        # Copy the entire CommercialAssets tree
        shutil.copytree(COMMERCIAL_DIR, archive_dir)
        
        # Also archive the current logs
        log_src = os.path.join(BASE_PATH, "logs")
        if os.path.exists(log_src):
            shutil.copytree(log_src, os.path.join(archive_dir, "logs_snapshot"))
            
        print("Archive complete.")
    except Exception as e:
        print(f"Error during archiving: {e}")

if __name__ == "__main__":
    archive_state()
