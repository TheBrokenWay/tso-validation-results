"""
PKPD_Simple.py
Simplified PK/PD Engine for v2.0-CORE Orchestrator

Provides deterministic PK/PD simulation without full population modeling.

DEPRECATED: Use PX_Engine.operations.TrialEngine with 7-tier IIV and
PX_Engine.operations.PKPD (sigmoid Emax, link_pk_to_pd) for production.
"""

import hashlib
from typing import Dict, Any, List
import math


def simulate_pk_simple(
    smiles: str,
    dose_mg: float,
    interval_hours: float,
    n_doses: int,
    clearance_L_per_h: float,
    vd_L: float
) -> Dict[str, Any]:
    """
    Simplified one-compartment PK simulation
    
    Args:
        smiles: SMILES string (for deterministic variability)
        dose_mg: Dose in mg
        interval_hours: Dosing interval in hours
        n_doses: Number of doses
        clearance_L_per_h: Clearance (L/h)
        vd_L: Volume of distribution (L)
    
    Returns:
        PK profile with AUC, Cmax, and metrics
    """
    
    # Calculate PK parameters
    vd_L = max(0.01, vd_L)
    ke = clearance_L_per_h / vd_L  # Elimination rate constant (1/h)
    ke = max(0.0001, ke)
    half_life = 0.693 / ke
    
    # Generate deterministic variability from SMILES
    hash_val = int(hashlib.md5(smiles.encode()).hexdigest(), 16)
    cv_factor = 1.0 + ((hash_val % 50) - 25) / 100.0  # ±25% CV
    
    # Simulate concentration-time profile
    time_points = []
    concentrations = []
    total_time = n_doses * interval_hours
    
    for t in range(0, int(total_time) + 1):
        time_points.append(float(t))
        
        # Sum concentrations from all doses
        conc = 0.0
        for dose_num in range(n_doses):
            dose_time = dose_num * interval_hours
            if t >= dose_time:
                time_since_dose = t - dose_time
                # One-compartment IV bolus: C = (Dose/Vd) * exp(-ke * t)
                c = (dose_mg / vd_L) * math.exp(-ke * time_since_dose)
                conc += c
        
        concentrations.append(conc * cv_factor)
    
    # Calculate AUC (trapezoidal rule)
    auc = 0.0
    for i in range(len(time_points) - 1):
        dt = time_points[i+1] - time_points[i]
        avg_conc = (concentrations[i] + concentrations[i+1]) / 2.0
        auc += avg_conc * dt
    
    # Calculate Cmax and Tmax
    cmax = max(concentrations)
    tmax = time_points[concentrations.index(cmax)]
    
    # Add inter-individual variability
    auc_mean = auc
    auc_sd = auc * 0.25 * cv_factor  # 25% CV
    auc_cv = 25.0
    
    cmax_mean = cmax
    cmax_sd = cmax * 0.25 * cv_factor
    
    return {
        "time_grid_h": time_points,
        "concentration_mg_per_L": concentrations,
        "auc_mean": round(auc_mean, 2),
        "auc_sd": round(auc_sd, 2),
        "auc_cv": auc_cv,
        "cmax_mean": round(cmax_mean, 2),
        "cmax_sd": round(cmax_sd, 2),
        "tmax_h": tmax,
        "half_life_h": round(half_life, 2),
        "params": {
            "clearance_L_per_h": clearance_L_per_h,
            "volume_L": vd_L,
            "ke_per_h": round(ke, 4)
        },
        "summary": {
            "auc_mg_h_per_L": round(auc_mean, 2),
            "cmax_mg_per_L": round(cmax_mean, 2)
        }
    }


def simulate_pd_simple(
    concentration_profile: List[float],
    time_grid: List[float],
    ec50: float,
    emax: float
) -> Dict[str, Any]:
    """
    Simplified PD simulation using Emax model
    
    Args:
        concentration_profile: Concentration at each time point (mg/L)
        time_grid: Time points (hours)
        ec50: EC50 (mg/L)
        emax: Maximum effect (0-1)
    
    Returns:
        PD profile with effect-time curve and metrics
    """
    
    # Apply Emax model at each time point
    effects = []
    for conc in concentration_profile:
        # Emax model: E = Emax * C / (EC50 + C)
        effect = emax * conc / (ec50 + conc)
        effects.append(effect)
    
    # Calculate metrics
    max_effect = max(effects)
    time_to_max = time_grid[effects.index(max_effect)]
    mean_effect = sum(effects) / len(effects)
    
    # AUEC (area under effect curve)
    auec = 0.0
    for i in range(len(time_grid) - 1):
        dt = time_grid[i+1] - time_grid[i]
        avg_effect = (effects[i] + effects[i+1]) / 2.0
        auec += avg_effect * dt
    
    # Add variability
    max_effect_mean = max_effect
    max_effect_sd = max_effect * 0.15  # 15% CV for PD
    
    return {
        "time_h": time_grid,
        "effect": effects,
        "max_effect_mean": round(max_effect_mean, 3),
        "max_effect_sd": round(max_effect_sd, 3),
        "time_to_max_effect_h": time_to_max,
        "auec_h": round(auec, 2),
        "mean_effect": round(mean_effect, 3),
        "params": {
            "emax": emax,
            "ec50": ec50
        }
    }


def simulate_pkpd(
    smiles: str,
    dose_mg: float,
    interval_hours: float,
    n_doses: int,
    ope_results: Dict[str, Any],
    admet_results: Dict[str, Any]
) -> Dict[str, Any]:
    """
    Combined PK/PD simulation
    
    Args:
        smiles: SMILES string
        dose_mg: Dose in mg
        interval_hours: Dosing interval
        n_doses: Number of doses
        ope_results: OPE analysis output
        admet_results: ADMET analysis output
    
    Returns:
        Combined PK/PD results
    """
    import warnings
    warnings.warn(
        "PKPD_Simple is deprecated. Use TrialEngine + PKPD.link_pk_to_pd for production.",
        DeprecationWarning,
        stacklevel=2,
    )
    # Extract parameters
    clearance = ope_results.get("clearance_estimate_L_per_h", 10.0)
    vd = ope_results.get("vd_estimate_L", 70.0)
    ec50 = ope_results.get("ec50", 10.0)
    emax = ope_results.get("emax", 0.9)
    
    # Run PK simulation
    pk_results = simulate_pk_simple(
        smiles=smiles,
        dose_mg=dose_mg,
        interval_hours=interval_hours,
        n_doses=n_doses,
        clearance_L_per_h=clearance,
        vd_L=vd
    )
    
    # Run PD simulation
    pd_results = simulate_pd_simple(
        concentration_profile=pk_results["concentration_mg_per_L"],
        time_grid=pk_results["time_grid_h"],
        ec50=ec50,
        emax=emax
    )
    
    return {
        "pk": pk_results,
        "pd": pd_results
    }
