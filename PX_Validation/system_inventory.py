import os
import json
from datetime import datetime, timezone

def perform_system_inventory():
    """Project operates on repo drive only (E: when repo is E:\\foundation). No C: or D: paths."""
    print("=== PREDATOR X: SYSTEM INVENTORY & MANIFEST SCAN ===")
    ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    manifest = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "structure": {},
        "file_counts": {},
        "critical_organs": ["PX_Constitution", "PX_Security", "PX_Engine", "PX_Warehouse", "PX_Audit", "PX_Laboratory"]
    }

    for organ in os.listdir(ROOT):
        path = os.path.join(ROOT, organ)
        if os.path.isdir(path):
            files = os.listdir(path)
            manifest["structure"][organ] = files[:10] # Sample first 10 for review
            manifest["file_counts"][organ] = len(files)

    # Specific check for the World-Line Database (canonical path)
    try:
        from PX_Warehouse.WorldLine_Database import DEFAULT_WORLDLINES_PATH
        wl_path = DEFAULT_WORLDLINES_PATH
    except Exception:
        wl_path = os.path.join(ROOT, "PX_Warehouse", "WorldLines")
    if os.path.exists(wl_path):
        manifest["worldline_total"] = len(os.listdir(wl_path))
    else:
        manifest["worldline_total"] = 0

    print(f"\n>>> [SCAN] Total Departments/Organs: {len(manifest['structure'])}")
    print(f">>> [SCAN] World-Lines Materialized: {manifest['worldline_total']}")
    
    manifest_path = os.path.join(ROOT, "PX_Validation", "system_manifest.json")
    with open(manifest_path, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=4)
    print(f"\n[INVENTORY COMPLETE] Manifest saved to {manifest_path}")

if __name__ == "__main__":
    perform_system_inventory()
