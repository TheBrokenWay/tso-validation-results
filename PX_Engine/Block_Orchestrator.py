import json
from PX_Security.PredatorImmune_Block import PredatorImmuneBlock
from PX_Warehouse.WorldLine_Database import WorldLineDatabase
from PX_Engine.Metabolism import Metabolism

IMMUNE = PredatorImmuneBlock()
WAREHOUSE = WorldLineDatabase()
HEARTBEAT = Metabolism()

def execute_living_pulse(task_id, source_ctx, p_vec, csa_s):
    # 1. Pulse the Heartbeat (Metabolic Cycle)
    cycle_age, status = HEARTBEAT.pulse(task_id, p_vec, csa_s)
    
    if status != "RESONANCE_ACHIEVED":
        print(f">>> [HEARTBEAT] Cycle: {cycle_age} | {task_id} FAILED (Manifold Drift).")
        return

    # 2. Process through Manifold Verification
    result = IMMUNE.handle_block_request(task_id, source_ctx, p_vec, csa_s)
    
    if result['authorized'] and result['route'] == "BLUE":
        # 3. Record with Metabolic Anchor
        WAREHOUSE.record_materialization(
            task_id=task_id, 
            block=result['block'], 
            coherence=result['coherence'], 
            lab_results={"status": "DETERMINISTIC_SYNTHESIS", "cycle": cycle_age},
            route="BLUE"
        )
        print(f">>> [HEARTBEAT] Cycle: {cycle_age} | {task_id} MATERIALIZED.")
    else:
        print(f">>> [HEARTBEAT] Cycle: {cycle_age} | {task_id} FAILED (Unauthorized).")

# Alias for PX_Laboratory.Synthetic_Expansion and other callers
execute_health_aware_pulse = execute_living_pulse

if __name__ == "__main__":
    ctx = {"source_id": "PX-KERNEL-01", "fingerprint": "PX-ROOT-Sovereign"}
    execute_living_pulse("METABOLIC-DRUG-01", ctx, [0.1, 0.0, 35.0, 1.0], [1.0, 1.0, 1.0, 1.0, 1.0])
