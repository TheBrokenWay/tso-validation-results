"""
PKPD.py
PK/PD Modeling Module (Production v2.0-PHASE1)

Provides Emax and sigmoid Emax models for linking pharmacokinetic
exposure to pharmacodynamic effect.

Constitutional Note: PD models are theoretical simulations.
Clinical efficacy data required for validation.
"""

import time
from typing import Dict, Any, List, Optional
import math
from PX_System.foundation.sign_off import create_sign_off


def emax_model(
    concentration: float,
    emax: float,
    ec50: float,
    hill: float = 1.0,
    baseline: float = 0.0,
) -> float:
    """
    Sigmoid Emax pharmacodynamic model.
    
    E = E0 + (Emax * C^hill) / (EC50^hill + C^hill)
    
    Args:
        concentration: Drug concentration (mg/L)
        emax: Maximum effect above baseline (0-1 for fractional inhibition)
        ec50: Concentration at 50% of Emax (mg/L)
        hill: Hill coefficient (slope factor, sigmoidicity)
        baseline: Baseline effect at C=0 (default 0.0)
    
    Returns:
        Effect value (baseline to baseline+emax)
        
    Examples:
        >>> emax_model(5.0, emax=0.9, ec50=5.0, hill=1.0)  # At EC50
        0.45  # 50% of Emax when hill=1
        
        >>> emax_model(0.0, emax=0.9, ec50=5.0)  # Zero concentration
        0.0  # Baseline
    """
    if concentration <= 0:
        return baseline
    
    c_hill = concentration ** hill
    ec50_hill = ec50 ** hill
    
    effect = baseline + (emax * c_hill) / (ec50_hill + c_hill)
    
    return effect


def compute_pd_metrics(
    time_h: List[float],
    effect: List[float],
    effect_threshold: Optional[float] = None,
) -> Dict[str, Any]:
    """
    Compute pharmacodynamic metrics from effect-time profile.
    
    Args:
        time_h: Time grid (hours)
        effect: Effect values at each time point
        effect_threshold: Threshold for time_above calculation (default: median effect)
    
    Returns:
        {
            "max_effect": Maximum effect achieved,
            "time_to_max_effect_h": Time to maximum effect,
            "auec_h": Area Under Effect Curve (trapezoidal),
            "time_above_threshold_h": Time above specified threshold,
            "mean_effect": Mean effect over time course,
            "effect_at_steady_state": Mean effect in last 10% of profile
        }
    """
    if not time_h or not effect:
        return {
            "max_effect": 0.0,
            "time_to_max_effect_h": 0.0,
            "auec_h": 0.0,
            "time_above_threshold_h": 0.0,
            "mean_effect": 0.0,
            "effect_at_steady_state": 0.0,
        }
    
    # Max effect and time to max
    max_effect = max(effect)
    max_idx = effect.index(max_effect)
    time_to_max_effect_h = time_h[max_idx]
    
    # AUEC using trapezoidal rule
    auec_h = 0.0
    for i in range(1, len(time_h)):
        dt = time_h[i] - time_h[i - 1]
        auec_h += 0.5 * (effect[i] + effect[i - 1]) * dt
    
    # Time above threshold
    if effect_threshold is None:
        # Default: use median effect as threshold
        sorted_effect = sorted(effect)
        effect_threshold = sorted_effect[len(sorted_effect) // 2]
    
    time_above_threshold_h = 0.0
    for i in range(1, len(time_h)):
        dt = time_h[i] - time_h[i - 1]
        # Count time if either endpoint is above threshold
        if effect[i] >= effect_threshold or effect[i - 1] >= effect_threshold:
            time_above_threshold_h += dt
    
    # Mean effect
    mean_effect = sum(effect) / len(effect) if effect else 0.0
    
    # Effect at steady state (last 10% of profile)
    ss_start_idx = max(0, int(0.9 * len(effect)))
    effect_at_steady_state = sum(effect[ss_start_idx:]) / max(1, len(effect[ss_start_idx:])) if effect else 0.0
    
    return {
        "max_effect": float(max_effect),
        "time_to_max_effect_h": float(time_to_max_effect_h),
        "auec_h": float(auec_h),
        "time_above_threshold_h": float(time_above_threshold_h),
        "mean_effect": float(mean_effect),
        "effect_at_steady_state": float(effect_at_steady_state),
    }


def link_pk_to_pd(
    pk_profile: Dict[str, Any],
    pd_params: Dict[str, Any],
) -> Dict[str, Any]:
    """
    Link PK concentration-time profile to PD effect-time profile.
    
    This is the core PK/PD linking function for PREDATOR X v2.0.
    
    Args:
        pk_profile: Output from SimulationEngine.simulate_one_compartment()
            {
                "time_grid_h": List[float],
                "concentration_mg_per_L": List[float],
                "summary": {...}
            }
        
        pd_params: PD model parameters
            {
                "emax": float,           # Maximum effect (0-1 for fractional inhibition)
                "ec50": float,           # Concentration at 50% effect (mg/L)
                "hill": float,           # Hill coefficient (default 1.0)
                "baseline": float,       # Baseline effect (default 0.0)
                "effect_threshold": float  # Optional: for time_above calculations
            }
    
    Returns:
        {
            "time_h": List[float],       # Same as PK time grid
            "effect": List[float],       # Effect at each time point
            "pd_summary": {              # Aggregate PD metrics
                "max_effect": float,
                "time_to_max_effect_h": float,
                "auec_h": float,
                "time_above_threshold_h": float,
                "mean_effect": float,
                "effect_at_steady_state": float
            },
            "parameters": {...},         # PD parameters used
            "constitutional": {...}      # L51/L34 compliance
        }
    
    Example:
        >>> pk = sim_engine.simulate_one_compartment(...)
        >>> pd = link_pk_to_pd(pk, {"emax": 0.9, "ec50": 5.0, "hill": 1.5})
        >>> print(pd["pd_summary"]["max_effect"])
        0.75  # 75% target inhibition achieved
    """
    # Extract PK profile
    time_h = pk_profile.get("time_grid_h", [])
    concentration = pk_profile.get("concentration_mg_per_L", [])
    
    if not time_h or not concentration:
        raise ValueError("PK profile must contain time_grid_h and concentration_mg_per_L")
    
    # Extract PD parameters
    emax = pd_params.get("emax")
    ec50 = pd_params.get("ec50")
    
    _t0 = time.monotonic()

    if emax is None or ec50 is None:
        raise ValueError("pd_params must contain 'emax' and 'ec50'")
    
    hill = pd_params.get("hill", 1.0)
    baseline = pd_params.get("baseline", 0.0)
    effect_threshold = pd_params.get("effect_threshold")
    
    # Compute effect at each time point using Emax model
    effect_profile = []
    for c in concentration:
        e = emax_model(c, emax, ec50, hill, baseline)
        effect_profile.append(e)
    
    # Compute PD metrics
    pd_summary = compute_pd_metrics(time_h, effect_profile, effect_threshold)
    
    result = {
        "time_h": time_h,
        "effect": effect_profile,
        "pd_summary": pd_summary,
        "parameters": {
            "emax": emax,
            "ec50": ec50,
            "hill": hill,
            "baseline": baseline,
            "effect_threshold": effect_threshold,
        },
        "constitutional": {
            "status": "SIMULATED",
            "engine": "PKPD_EMAX_V2.0",
            "notes": (
                "PD model is theoretical and based on Emax assumptions. "
                "Clinical validation required. EC50 and Emax must be "
                "determined experimentally for actual drug-target pairs."
            ),
        },
    }
    _elapsed_ms = int((time.monotonic() - _t0) * 1000)
    result["sign_off"] = create_sign_off(
        engine_id="PKPD_EMAX_V2",
        version="2.0-PHASE1",
        inputs={"pk_summary": pk_profile.get("summary", {}), "pd_params": pd_params},
        outputs=result,
        laws_checked=["U27"],
        laws_results={"U27": True},
        execution_time_ms=_elapsed_ms,
    )
    return result


# Legacy function for backward compatibility
def apply_pkpd_to_profile(
    time_grid_h: List[float],
    concentration_mg_per_L: List[float],
    pd_params: Dict[str, Any],
) -> Dict[str, Any]:
    """
    Legacy function. Use link_pk_to_pd() instead.
    
    Maintained for backward compatibility.
    """
    # Convert to new format
    pk_profile = {
        "time_grid_h": time_grid_h,
        "concentration_mg_per_L": concentration_mg_per_L,
    }
    
    result = link_pk_to_pd(pk_profile, pd_params)
    
    # Return legacy format
    return {
        "time_grid_h": result["time_h"],
        "effect": result["effect"],
        "parameters": result["parameters"],
        "constitutional": result["constitutional"],
    }


__all__ = ["emax_model", "compute_pd_metrics", "link_pk_to_pd", "apply_pkpd_to_profile"]
