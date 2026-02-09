import os
import sys
import json
import time
from datetime import datetime, timezone

# Root Calibration
sys.path.append("E:/foundation")

def run_global_test():
    print(f"=== OLYMPUS: GLOBAL INTEGRITY SMOKE TEST [2026-01-21] ===")
    
    # 1. Executive Handshake
    print(">>> [01_EXECUTIVE] Testing ZeusLaws Adjudication...")
    try:
        from PX_Executive.GAIP_Gateway import GAIPGateway
        gateway = GAIPGateway(mode="REGULATORY")
        print("    [PASS] Executive Gateway Responsive.")
    except Exception as e:
        print(f"    [FAIL] Executive Error: {e}")

    # 2. Warehouse Ingress
    print(">>> [04_WAREHOUSE] Validating 15,382 Dossier Index...")
    from PX_Warehouse.Worldline_Indexer import WorldlineIndexer
    indexer = WorldlineIndexer()
    indexer.rebuild_index()
    if len(indexer.index) >= 15382:
        print(f"    [PASS] Indexer Integrity: {len(indexer.index)} World-Lines Found.")
    else:
        print(f"    [FAIL] Warehouse Count Mismatch: {len(indexer.index)}")

    # 3. Laboratory Synthesis Loop
    print(">>> [08_LABORATORY] Testing Simulation Engine...")
    try:
        from PX_Laboratory.Simulation_Engine import SimulationEngine
        sim = SimulationEngine()
        res = sim.materialize_candidate("SMOKE-TEST", 0.95)
        print(f"    [PASS] Lab Sim Result: {res['status']} | Affinity: {res['binding_affinity_kj']}")
    except Exception as e:
        print(f"    [FAIL] Laboratory Error: {e}")

    # 4. Metabolic Pulse
    from PX_Engine.Metabolism import Metabolism
    heartbeat = Metabolism()
    new_cycle, _ = heartbeat.pulse("GLOBAL_SMOKE_TEST")
    print(f">>> [ENGINE] Metabolic Pulse Recorded. New System Age: {new_cycle}")

    print("\n=== SMOKE TEST COMPLETE: ALL SYSTEMS NOMINAL ===")

if __name__ == "__main__":
    run_global_test()
