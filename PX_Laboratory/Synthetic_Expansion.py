import numpy as np
from PX_Engine.Block_Orchestrator import execute_health_aware_pulse

def run_expansion(steps=5):
    print(f"\n>>> [EXPANSION] Sculpting Manifold: Generating {steps} Synthetic World-Lines...")
    
    # Interpolating between 0.10 (Alpha) and 0.11 (Beta)
    complexities = np.linspace(0.10, 0.11, steps)
    
    ctx = {"source_id": "PX-INTERNAL-LAB", "fingerprint": "PX-ROOT-Sovereign"}
    
    for i, comp in enumerate(complexities):
        task_id = f"SYNTH-EXP-{i:02d}"
        p_vec = [float(comp), 0.0, 35.0, 1.0]
        csa_s = [1.0, 1.0, 1.0, 1.0, 1.0]
        execute_health_aware_pulse(task_id, ctx, p_vec, csa_s)

if __name__ == "__main__":
    run_expansion(steps=10)
