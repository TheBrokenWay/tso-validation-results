"""
╔══════════════════════════════════════════════════════════════════════════════╗
║ Gold_Rush_Miner.py                                                           ║
║ PREDATOR X :: WAREHOUSE MINING & COMMERCIAL DOSSIER GENERATION              ║
║ ARCHITECT: JAMES A. TILLAR | STATUS: PRODUCTION READY (OFFSET-AWARE)        ║
╚══════════════════════════════════════════════════════════════════════════════╝

PURPOSE:
    Mine 30,000+ WorldLines for gold-tier drug candidates.
    
UPDATED LOGIC (SMART OFFSETS):
    1. CLASSIC GOLD: Toxicity Index < 0.021 (The perfect ones)
    2. OFFSET DIAMOND: Toxicity > 0.021 BUT Safety Margin > 50x (The real ones)
       (Saves drugs like Caffeine, Aspirin, Aticaprant)

EXECUTION:
    python PX_Executive/Gold_Rush_Miner.py
"""

import os
import sys
import json
import time
from pathlib import Path

# ROOT CALIBRATION
ROOT_DIR = Path(__file__).resolve().parents[1]
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

# IMPORT ORGANS
from PX_Warehouse.WorldLine_Database import WorldLineDatabase
from PX_Laboratory.Simulation_Engine import SimulationEngine
from PX_Executive.Sovereign_Commercial_Pipeline import generate_dossier
from PX_Executive.PX_Legal_Check import PX_Legal_Check

def main():
    print("\n" + "="*80)
    print("║ GOLD RUSH MINER :: OFFSET-AWARE PROTOCOL ║".center(80))
    print("="*80 + "\n")

    # 1. CONNECT TO DATABASE
    db = WorldLineDatabase()
    all_candidates = db.get_all_candidates() # Expects list of dicts or paths
    
    print(f"Loaded {len(all_candidates)} WorldLines from Warehouse.")
    print("Engaging Smart Filter (Tox < 0.021 OR SafetyMargin > 50x)...")

    processed = 0
    dossiers_created = 0
    
    # 2. ENGAGE LEGAL GATE
    legal_checker = PX_Legal_Check(mode="REGULATORY")

    for candidate in all_candidates:
        # Extract Physics
        tox = candidate.get("physics", {}).get("toxicity_index", 1.0)
        sm = candidate.get("physics", {}).get("safety_margin", 0.0)
        
        # --- THE SMART FILTER ---
        is_classic_gold = tox < 0.021
        is_offset_diamond = (tox >= 0.021) and (sm > 50.0)
        
        if not (is_classic_gold or is_offset_diamond):
            continue # TRASH

        # If we get here, it's valuable.
        processed += 1
        name = candidate.get("name", "Unknown")
        smiles = candidate.get("smiles")
        
        print(f"💎 FOUND GEM: {name} (Tox: {tox:.4f} | SM: {sm:.1f}x)")

        # 3. LEGAL CHECK
        fto = legal_checker.check_compound(smiles, name)
        if not fto.freedom_to_operate:
            print(f"   ⛔ FTO BLOCKED: Patent active.")
            continue

        # 4. GENERATE COMMERCIAL DOSSIER
        # Writes to canonical warehouse (Prv_Dossiers/Novel_Dossiers by tier via pipeline)
        try:
            dossier_path = generate_dossier(candidate, candidate.get("worldline_path"))
            if dossier_path:
                print(f"   ✅ DOSSIER MINTED: {os.path.basename(dossier_path)}")
                dossiers_created += 1
        except Exception as e:
            print(f"   ⚠️ GENERATION FAILED: {e}")

    print("\n" + "="*80)
    print(f"MINING COMPLETE. GENERATED {dossiers_created} NEW ASSETS.")
    print("="*80)

if __name__ == "__main__":
    main()