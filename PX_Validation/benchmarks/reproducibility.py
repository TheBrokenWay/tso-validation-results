import sys
import os
import json
import hashlib
from pathlib import Path
from datetime import datetime, UTC

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from PX_Executive.orchestrators.PX_Live_Orchestrator_v2 import PredatorXOrchestratorV2

def test_reproducibility(smiles_list, iterations=5):
    """
    Test if Predator X yields bit-for-bit identical results for the same SMILES.
    """
    orchestrator = PredatorXOrchestratorV2(verbose=False)
    results = {}
    
    print(f"Starting Reproducibility Test ({iterations} iterations per molecule)...")
    
    for smiles in smiles_list:
        print(f"Testing SMILES: {smiles}")
        iteration_results = []
        for i in range(iterations):
            metadata = {"name": f"ReproTest_{i}", "run_id": f"repro_{i}"}
            # We need to mock run_id and timestamp to be constant for bit-for-bit comparison
            # Actually, let's just compare the 'stages' part of the results
            res = orchestrator.run_pipeline(smiles, metadata)
            iteration_results.append(res["stages"])
        
        # Compare all iterations to the first one
        first_res = json.dumps(iteration_results[0], sort_keys=True)
        is_reproducible = True
        for i in range(1, iterations):
            current_res = json.dumps(iteration_results[i], sort_keys=True)
            if current_res != first_res:
                print(f"  ❌ Iteration {i} failed reproducibility!")
                is_reproducible = False
                break
        
        if is_reproducible:
            print(f"  ✅ SMILES {smiles} is 100% reproducible.")
        
        results[smiles] = {
            "is_reproducible": is_reproducible,
            "iterations": iterations
        }
    
    return results

if __name__ == "__main__":
    test_smiles = [
        "CC(=O)OC1=CC=CC=C1C(=O)O",  # Aspirin
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",  # Caffeine
        "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"  # Ibuprofen
    ]
    
    repro_results = test_reproducibility(test_smiles)
    
    # Output structured JSON
    output_path = Path("PX_Warehouse/logs/reproducibility_metrics.json")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(repro_results, f, indent=2)
    
    print(f"\nReproducibility metrics saved to {output_path}")
