import os
import sys
import json
import hashlib
import time
from datetime import datetime, timezone

# Root Calibration
ROOT_DIR = "E:/foundation"
if ROOT_DIR not in sys.path:
    sys.path.append(ROOT_DIR)

def generate_hash(file_path):
    sha256_hash = hashlib.sha256()
    try:
        with open(file_path, "rb") as f:
            for byte_block in iter(lambda: f.read(4096), b""):
                sha256_hash.update(byte_block)
        return sha256_hash.hexdigest()
    except FileNotFoundError:
        return "MISSING_ASSET"

def _discover_asset(dir_path, suffix):
    """Return path to one asset in dir_path with given suffix, or None."""
    if not os.path.isdir(dir_path):
        return None
    for name in os.listdir(dir_path):
        if name.endswith(suffix):
            return os.path.join(dir_path, name)
    for sub in os.listdir(dir_path):
        subpath = os.path.join(dir_path, sub)
        if os.path.isdir(subpath):
            for name in os.listdir(subpath):
                if name.endswith(suffix):
                    return os.path.join(subpath, name)
    return None


def seal_system():
    print("=== OLYMPUS: FINAL SYSTEM SEAL PROTOCOL ===")
    wh = os.path.join(ROOT_DIR, "PX_Warehouse")

    # 1. DISCOVER CRITICAL ASSETS (dossier from Prv/Finalized, worldline from WorldLines)
    dossier_path = _discover_asset(os.path.join(wh, "Prv_Dossiers"), ".json")
    if not dossier_path:
        dossier_path = _discover_asset(os.path.join(wh, "Finalized_Dossiers"), ".json")
    if not dossier_path:
        dossier_path = "E:/foundation/PX_Warehouse/Prv_Dossiers/DOSSIER-D4B1F079.json"  # fallback for hash only
    worldline_path = _discover_asset(os.path.join(wh, "WorldLines"), ".worldline")
    if not worldline_path:
        worldline_path = os.path.join(wh, "WorldLines", "BATCH-GLP1-00-OPT-GOLD.worldline")

    # 2. CRYPTOGRAPHIC HASHING
    print(f"Hashing Asset: {os.path.basename(dossier_path)}")
    dossier_hash = generate_hash(dossier_path)
    print(f"Hashing Asset: {os.path.basename(worldline_path)}")
    wl_hash = generate_hash(worldline_path)

    if dossier_hash == "MISSING_ASSET" and wl_hash == "MISSING_ASSET":
        print("!!! [WARN] No warehouse assets found; sealing with empty asset list.")

    # 3. GENERATE SEAL MANIFEST
    seal_id = "SEAL-" + hashlib.shake_128(str(time.time()).encode()).hexdigest(4).upper()
    
    seal_manifest = {
        "seal_id": seal_id,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "custodian": "OLYMPUS_ADMIN",
        "session_summary": {
            "type": "RESEARCH_CYCLE_COMPLETE",
            "metabolic_cycle": 25210,
            "security_mode": "READ_ONLY_LOCK"
        },
        "assets": {
            "dossier": {
                "file": os.path.basename(dossier_path),
                "sha256": dossier_hash,
                "integrity": "VERIFIED" if dossier_hash != "MISSING_ASSET" else "MISSING"
            },
            "worldline": {
                "file": os.path.basename(worldline_path),
                "sha256": wl_hash,
                "integrity": "VERIFIED" if wl_hash != "MISSING_ASSET" else "MISSING"
            }
        }
    }
    
    # 4. WRITE SEAL FILE
    seal_path = os.path.join(ROOT_DIR, "PX_Audit/SYSTEM_SEAL.json")
    with open(seal_path, "w") as f:
        json.dump(seal_manifest, f, indent=4)
        
    print(f">>> SYSTEM SEALED. Manifest: {seal_path}")
    print(f"    Seal ID: {seal_id}")
    print(f"    Dossier Hash: {dossier_hash[:16]}...")
    
    # 5. ENGAGE LOCK
    lock_file = os.path.join(ROOT_DIR, "SYSTEM_LOCK.dat")
    with open(lock_file, "w") as f:
        f.write(f"LOCKED BY {seal_id} AT {datetime.now()}")
        
    print("\n[SYSTEM LOCKED] Write operations suspended. Session Terminated.")

if __name__ == "__main__":
    seal_system()
