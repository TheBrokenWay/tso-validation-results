import time
import numpy as np
from PX_Engine.Engine_Orchestrator import predator_discovery_cycle
from PX_Audit.Mural_Network import get_mural

def run_stress_test(cycles=100):
    print(f">>> [PREDATOR X] Starting {cycles}-cycle Stress Test...")
    latencies = []
    
    for i in range(cycles):
        t0 = time.perf_counter()
        # Simulated GLP-1 Discovery Path
        predator_discovery_cycle(f"STRESS_{i}", 0.1, 0.0, 35.0, 1.0)
        latencies.append(time.perf_counter() - t0)
        
    avg_lat = sum(latencies) / cycles
    print(f"\n=== STRESS TEST COMPLETE ===")
    print(f"Avg Cycle Latency: {avg_lat:.6f}s")
    print(f"Total Mural Edges: {len(get_mural()['edges'])}")

if __name__ == "__main__":
    run_stress_test(100)
