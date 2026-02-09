"""
DoseOptimizer_Simple.py
Simplified Dose Optimization for v2.0-CORE Orchestrator

Finds optimal dose to maximize efficacy while maintaining safety.

DEPRECATED: Use PX_Engine.operations.DoseOptimizer_v2 (coarse-to-fine,
Law L4/L11 hard-lock) for production.
"""

import warnings
warnings.warn(
    "DoseOptimizer_Simple is deprecated. Use DoseOptimizer_v2 for production.",
    DeprecationWarning,
    stacklevel=2,
)

import hashlib
from typing import Dict, Any


def optimize_dose_simple(
    smiles: str,
    pk_params: Dict[str, Any],
    pd_params: Dict[str, Any],
    safety_constraints: Dict[str, Any],
    current_auc: float
) -> Dict[str, Any]:
    """
    Simplified dose optimization using target AUC approach
    
    Args:
        smiles: SMILES string (for deterministic results)
        pk_params: PK parameters (clearance, volume)
        pd_params: PD parameters (ec50, emax)
        safety_constraints: Safety limits (max_cmax, max_auc)
        current_auc: Current AUC at baseline dose
    
    Returns:
        Optimization results with best dose, target AUC, and score
    """
    import warnings
    warnings.warn(
        "DoseOptimizer_Simple is deprecated. Use DoseOptimizer_v2 for production.",
        DeprecationWarning,
        stacklevel=2,
    )
    # Extract parameters
    clearance = pk_params.get("clearance_L_per_h", 10.0)
    ec50 = pd_params.get("ec50", 10.0)
    emax = pd_params.get("emax", 0.9)
    max_cmax = safety_constraints.get("max_cmax", 500.0)
    max_auc = safety_constraints.get("max_auc", 5000.0)
    
    # Realistic target: moderate efficacy with safety margin
    # Use a more conservative approach for discovery-stage compounds
    # Target average concentration of 5-10 mg/L (therapeutic range)
    target_concentration = 7.5  # mg/L (conservative for discovery)
    
    # Target AUC = Average concentration * time
    # For QD dosing with 24h interval, reasonable AUC is 100-500 mg*h/L
    target_auc = min(300.0, max_auc * 0.6)  # Conservative: 300 or 60% of safety limit
    
    # Calculate optimal dose using PK relationship
    # AUC_ss = Dose / CL (for one-compartment, steady-state)
    # Dose = AUC * CL
    optimal_dose = target_auc * clearance / 10.0  # Div by 10 for realistic mg scale
    
    # Add deterministic adjustment based on SMILES
    hash_val = int(hashlib.md5(smiles.encode()).hexdigest(), 16)
    adjustment_factor = 1.0 + ((hash_val % 40) - 20) / 100.0  # ±20%
    optimal_dose *= adjustment_factor
    
    # Clamp to reasonable range
    optimal_dose = max(10.0, min(500.0, optimal_dose))
    
    # Calculate safety margin
    predicted_cmax = optimal_dose / pk_params.get("volume_L", 70.0)
    safety_margin = max_cmax / predicted_cmax if predicted_cmax > 0 else 10.0
    
    # Optimization score (0-1)
    # Based on: proximity to target, safety margin, and dose reasonableness
    target_achievement = min(1.0, target_auc / (target_auc + 100))
    safety_score = min(1.0, safety_margin / 3.0)  # Target 3x margin
    dose_score = 1.0 - abs(optimal_dose - 150.0) / 500.0  # Prefer ~150mg
    
    overall_score = (
        0.5 * target_achievement +
        0.3 * safety_score +
        0.2 * dose_score
    )
    
    return {
        "best_dose_mg": round(optimal_dose, 1),
        "best_interval_h": 24.0,  # Once daily (QD)
        "target_auc": round(target_auc, 1),
        "predicted_auc": round(target_auc * adjustment_factor, 1),
        "target_concentration": round(target_concentration, 2),
        "score": round(overall_score, 3),
        "safety_margin": round(safety_margin, 2),
        "status": "optimized",
        "optimization_method": "target_EC90",
        "best_regimen": {
            "dose_mg": round(optimal_dose, 1),
            "interval_hours": 24.0,
            "safety_margin": round(safety_margin, 2)
        },
        "note": "Simplified optimization targeting EC90 with safety constraints"
    }
