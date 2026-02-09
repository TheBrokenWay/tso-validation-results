import sys
import os
import json
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from PX_Executive.orchestrators.PX_Live_Orchestrator_v2 import PredatorXOrchestratorV2

def test_accuracy():
    """
    Evaluate False Positive Rate (FPR) against known constraints.
    """
    orchestrator = PredatorXOrchestratorV2(verbose=False)
    
    # Curated list of compounds that SHOULD fail Predator X v3 gates
    # 1. High toxicity
    # 2. Low bioavailability
    # 3. Poor binding affinity
    negative_controls = [
        {"smiles": "C1=CC=C(C=C1)N=NC2=CC=C(C=C2)N", "reason": "Potential mutagen"},
        {"smiles": "CCCCCCCCCCCCCCCCCCCC", "reason": "Extremely low bioavailability"},
        {"smiles": "C1CCCCC1", "reason": "No functional groups for binding"}
    ]
    
    print("Starting Accuracy / FPR Evaluation...")
    
    false_positives = 0
    total = len(negative_controls)
    
    for control in negative_controls:
        smiles = control["smiles"]
        res = orchestrator.run_pipeline(smiles, {"name": "NegativeControl"})
        stages = res["stages"]
        
        # Check if it passes DIAMOND or GOLD gates
        sm = stages.get("dose_optimization", {}).get("safety_margin", 0)
        rr = stages.get("virtual_efficacy", {}).get("responder_rate", {}).get("response_ratio", 0)
        
        # v3 GOLD gate: SM >= 50, RR >= 0.10
        # DIAMOND gate: bio >= 80 and sm >= 50 and rr >= 0.25 and ba <= 20
        
        # risk_level from ADMET: 4-tier (DIAMOND, GOLD, SILVER, FAILURE); legacy HIGH_RISK_TRASH = reject
        risk_level = stages.get("admet", {}).get("toxicity", {}).get("risk_level", "UNKNOWN")
        # Diamond, Gold, Silver are viable; only FAILURE (and legacy HIGH_RISK_TRASH) = hard rejection
        toxicity_ok = risk_level in ("TOXICITY_DIAMOND", "TOXICITY_GOLD", "TOXICITY_SILVER")
        if risk_level == "HIGH_RISK_TRASH":
            toxicity_ok = False  # Legacy label
        is_gold = sm >= 50 and rr >= 0.10 and toxicity_ok

        if is_gold:
            print(f"  ❌ FALSE POSITIVE: {smiles} ({control['reason']}) passed as GOLD")
            false_positives += 1
        else:
            print(f"  ✅ Correctly rejected: {smiles} ({control['reason']})")
            
    fpr = false_positives / total
    print(f"  ✅ False Positive Rate: {fpr*100:.1f}%")
    
    results = {
        "total_negative_controls": total,
        "false_positives": false_positives,
        "false_positive_rate": fpr,
        "accuracy_passed": fpr < 0.1
    }
    
    return results

if __name__ == "__main__":
    accuracy_results = test_accuracy()
    repo_root = Path(__file__).resolve().parents[2]
    out_dir = repo_root / "PX_Warehouse" / "Operations" / "reports"
    out_dir.mkdir(parents=True, exist_ok=True)
    output_path = out_dir / "accuracy_metrics.json"
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(accuracy_results, f, indent=2)
    print(f"Accuracy metrics saved to {output_path}")
