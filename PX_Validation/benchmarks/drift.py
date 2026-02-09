import sys
import os
import json
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from PX_Executive.orchestrators.PX_Live_Orchestrator_v2 import PredatorXOrchestratorV2

def track_grade_drift(smiles_list):
    """
    Track how grades change between engine versions or settings.
    """
    orchestrator = PredatorXOrchestratorV2(verbose=False)
    
    print("Starting Grade Drift Tracking...")
    
    # Mock historical grades (v2.0)
    historical_grades = {
        "CC(=O)OC1=CC=CC=C1C(=O)O": "SILVER",
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C": "BRONZE",
        "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O": "SILVER"
    }
    
    drift_data = []
    
    for smiles in smiles_list:
        res = orchestrator.run_pipeline(smiles, {"name": "DriftTrack"})
        stages = res["stages"]
        
        # Calculate v3 grade
        sm = stages.get("dose_optimization", {}).get("safety_margin", 0)
        rr = stages.get("virtual_efficacy", {}).get("responder_rate", {}).get("response_ratio", 0)
        
        v3_grade = "REJECT"
        if sm >= 50 and rr >= 0.10: v3_grade = "GOLD"
        elif sm >= 20 and rr >= 0.05: v3_grade = "SILVER"
        elif sm >= 10 and rr >= 0.02: v3_grade = "BRONZE"
        
        historical = historical_grades.get(smiles, "UNKNOWN")
        
        drift_data.append({
            "smiles": smiles,
            "historical_grade": historical,
            "current_grade": v3_grade,
            "drifted": historical != v3_grade
        })
        
        if historical != v3_grade:
            print(f"  ⚠️  DRIFT: {smiles} moved from {historical} to {v3_grade}")
            
    drift_count = sum(1 for d in drift_data if d["drifted"])
    print(f"  ✅ Grade Drift: {drift_count}/{len(smiles_list)} assets changed grade.")
    
    results = {
        "drift_data": drift_data,
        "total_drift_count": drift_count,
        "drift_percentage": drift_count / len(smiles_list)
    }
    
    return results

if __name__ == "__main__":
    test_smiles = [
        "CC(=O)OC1=CC=CC=C1C(=O)O",
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"
    ]
    
    drift_results = track_grade_drift(test_smiles)
    
    output_path = Path("PX_Warehouse/logs/drift_metrics.json")
    with open(output_path, "w") as f:
        json.dump(drift_results, f, indent=2)
    
    print(f"Drift metrics saved to {output_path}")
