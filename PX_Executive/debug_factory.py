import sys
import os
import json
from pathlib import Path

print("🩺 STARTING SYSTEM DIAGNOSTIC...\n")

# 1. CHECK PATHS
ROOT = Path.cwd()
print(f"📂 Root Directory: {ROOT}")
if not (ROOT / "PX_Executive").exists():
    print("❌ CRITICAL: PX_Executive folder missing!")
    sys.exit(1)

# Add Root to Path so imports work
sys.path.insert(0, str(ROOT))

# 2. TEST ENGINE IMPORTS
print("🔧 Testing Engine Imports...")
try:
    from PX_Engine.operations.ADMET import run_admet
    from PX_Engine.operations.OPE import run_ope
    print("   ✅ ADMET Engine: LOADED")
except Exception as e:
    print(f"   ❌ ADMET Engine FAILURE: {e}")
    sys.exit(1)

try:
    from PX_System.governance.Placement_Gatekeeper import PlacementGatekeeper
    print("   ✅ Gatekeeper: LOADED")
except Exception as e:
    print(f"   ❌ Gatekeeper FAILURE: {e}")
    sys.exit(1)

# 3. TEST FEEDED
print("📥 Checking Feeder...")
queue_path = ROOT / "PX_Warehouse" / "Feeder" / "prv_24h_queue.json"
if queue_path.exists():
    print(f"   ✅ Queue Found: {queue_path.name}")
else:
    print("   ⚠️  Queue Missing in Feeder. (Creating Dummy...)")
    (ROOT / "PX_Warehouse" / "Feeder").mkdir(parents=True, exist_ok=True)
    dummy = {"candidates": [{"id": "DEBUG_001", "type": "R", "smiles": "CCO", "name": "Ethanol_Debug"}]}
    queue_path.write_text(json.dumps(dummy), encoding="utf-8")
    print("   ✅ Dummy Queue Created.")

# 4. DRY RUN (Generate 1 Asset)
print("⚙️  Attempting Single Production Run...")
try:
    # A. Run Physics (Mock)
    smiles = "CC(=O)Oc1ccccc1C(=O)O" # Aspirin
    print(f"   -> Processing: {smiles}")
    
    # B. Run ADMET (single API: run_admet with OPE)
    ope = run_ope(smiles)
    admet = run_admet(smiles, ope)
    tier = admet["toxicity"]["risk_level"]
    print(f"   -> ADMET Result: {tier}")
    
    # C. Mock Dossier
    safety_margin = admet["safety_margins"].get("safety_margin", 0)
    dossier = {
        "candidate": {"name": "DEBUG_ASPIRIN", "smiles": smiles},
        "engines": {"admet": admet, "ope": ope},
        "safety_margin": safety_margin
    }
    
    # D. Gatekeeper Placement
    root_folder, tier_name = PlacementGatekeeper.resolve_placement(dossier, "R")
    print(f"   -> Gatekeeper Routing: {root_folder} / {tier_name}")
    
    final_path = PlacementGatekeeper.get_secure_path(ROOT, "PRV_REP_DEBUG_001.json", dossier, "R")
    final_path.write_text(json.dumps(dossier, indent=2))
    print(f"   ✅ SUCCESS! File written to: {final_path}")

except Exception as e:
    print(f"\n❌ PRODUCTION FAILURE: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

print("\n🟢 SYSTEM IS HEALTHY. The issue is likely in the 'run_prv_repurposed.py' wrapper.")