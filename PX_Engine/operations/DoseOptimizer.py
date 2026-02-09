"""
DoseOptimizer.py
Dose Optimization Module (Deterministic v1)

Provides grid search and optimization methods for finding optimal doses
based on exposure targets.

Constitutional Note: This is a deterministic implementation for exposure-based
dose optimization. Clinical efficacy endpoints require clinical trial data.
"""

from typing import Dict, Any, List, TypedDict
from PX_Engine.operations import TrialEngine


class ToxicityInfo(TypedDict):
    toxicity_index: float


class AdmetInfo(TypedDict, total=False):
    toxicity: ToxicityInfo
    safety_margins: Dict[str, float]


def grid_search_dose(
    smiles: str,
    admet: AdmetInfo,
    dose_grid_mg: List[float],
    target_auc_mg_h_per_L: float,
    n_patients: int = 20,
) -> Dict[str, Any]:
    """
    Simple exposure-based dose optimization via grid search.
    
    Searches through a grid of doses to find the one that produces
    an AUC closest to the target.
    
    Args:
        smiles: SMILES string (for reference)
        admet: ADMET analysis dict
        dose_grid_mg: List of doses to test (mg)
        target_auc_mg_h_per_L: Target AUC (mg·h/L)
        n_patients: Number of virtual patients per dose
    
    Returns:
        Dict with optimal dose, achieved AUC, and error
    """
    
    engine = TrialEngine(time_step_h=1.0)
    best = None
    
    # CONSTITUTIONAL HARD-LOCK (Law L4/L11/U27): Harm is impossible
    if "toxicity" not in admet or "toxicity_index" not in admet["toxicity"]:
        raise ValueError("TraceabilityError: toxicity_index missing from ADMET")
    tox = float(admet["toxicity"]["toxicity_index"])
    if tox >= 0.0210:
        return {
            "status": "COLLAPSED",
            "reason": "LAW L4 HARD-LOCK: Toxicity index >= 0.0210",
            "toxicity_index": tox
        }
    
    for dose in dose_grid_mg:
        if tox >= 0.0210:
            return {
                "status": "COLLAPSED",
                "reason": "LAW L4 HARD-LOCK: Toxicity index >= 0.0210",
                "toxicity_index": tox,
            }
        # CONSTITUTIONAL HARD-LOCK: Toxicity Filter for target AUC
        # If the dose produces an AUC that would require exceeding the 0.0210 toxicity threshold, reject.
        # Since toxicity is currently calculated from SMILES/descriptors, it's mostly dose-independent in this simplified model,
        # but we must ensure the safety margin holds.
        
        protocol = {
            "trial_id": f"TRIAL-DOSE-{dose}",
            "duration_days": 3.0,
            "arms": [
                {
                    "arm_id": "A1",
                    "label": f"QD {dose} mg",
                    "dose_mg": dose,
                    "dosing_interval_h": 24.0,
                    "n_patients": n_patients,
                }
            ],
        }
        
        trial = engine.run_trial(protocol, admet)
        auc = trial["arms"][0]["exposure_summary"]["auc_mg_h_per_L"]["mean"]
        
        # Verify Law L11 for this dose level
        # If AUC exceeds max_tolerated_auc, it's a violation
        max_auc = admet.get("safety_margins", {}).get("max_tolerated_auc_mg_h_per_L", 5000.0)
        if auc > max_auc:
            continue  # This dose is too high, skip it
            
        error = abs(auc - target_auc_mg_h_per_L)
        
        candidate = {
            "dose_mg": dose,
            "auc_mg_h_per_L": auc,
            "error": error,
        }
        
        if best is None or error < best["error"]:
            best = candidate
    
    if best is None:
        return {
            "status": "NO_VALID_DOSE",
            "reason": "No dose satisfied constitutional AUC limits",
            "toxicity_index": tox,
        }

    # Add constitutional metadata
    best["constitutional"] = {
        "status": "OPTIMIZED",
        "engine": "DOSE_OPTIMIZER_GRID_SEARCH_V1",
        "notes": "Exposure-based optimization only. Clinical efficacy not considered.",
        "target_auc": target_auc_mg_h_per_L,
        "doses_tested": len(dose_grid_mg),
    }
    
    return best


__all__ = ["grid_search_dose"]
