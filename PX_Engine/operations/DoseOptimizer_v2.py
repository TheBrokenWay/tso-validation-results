"""
DoseOptimizer_v2.py
Dose Optimization Module v2.1 - Holistic Assessment

Production implementation with Whole Product view:
- No single-point toxicity kill; high-tox candidates are optimized and judged by
  safety margin and PRV eligibility (Smart Offset, strategic pass).
- scoring_function accepts optional admet_data and penalizes toxicity heavily unless
  safety margin is massive (>= 50) or regimen is otherwise acceptable.

Constitutional: All optimization based on virtual trial simulations. L51/L34 compliant.
"""

import time
from typing import Dict, Any, List, Tuple, Optional
import statistics
from PX_Engine.operations.TrialEngine import TrialEngine
from PX_System.foundation.sign_off import create_sign_off


def optimize_dose(
    smiles: str,
    admet: Dict[str, Any],
    protocol_template: Dict[str, Any],
    target_pk_range: Optional[Dict[str, Tuple[float, float]]] = None,
    target_pd_range: Optional[Dict[str, Tuple[float, float]]] = None,
    dose_bounds: Tuple[float, float] = (10.0, 500.0),
    interval_options: List[float] = None,
    variability: Optional[Dict[str, Any]] = None,
    pd_params: Optional[Dict[str, Any]] = None,
    search_strategy: str = "coarse_to_fine",
    n_eval_patients: int = 10,
) -> Dict[str, Any]:
    """
    Optimize dose and dosing interval for target PK/PD ranges.
    
    v2.0-PHASE4: Production dose optimization with multi-dimensional search.
    
    Args:
        smiles: SMILES string (for reference)
        admet: ADMET analysis dict
        protocol_template: Base trial protocol (will modify dose/interval)
        target_pk_range: Dict of PK metric → (min, max) ranges
            Example: {"auc_mg_h_per_L": (200.0, 400.0)}
        target_pd_range: Dict of PD metric → (min, max) ranges
            Example: {"max_effect": (0.6, 0.8)}
        dose_bounds: (min_dose, max_dose) in mg
        interval_options: List of dosing intervals to test (hours)
            Default: [8.0, 12.0, 24.0] (TID, BID, QD)
        variability: IIV parameters for realistic populations
        pd_params: PD parameters if PD optimization desired
        search_strategy: "coarse_to_fine" or "binary_search"
        n_eval_patients: Number of patients per mini-trial
    
    Returns:
        Dict with:
            - best_regimen: {dose_mg, interval_h, score}
            - search_history: List of all evaluated regimens
            - target_achievement: How close to target ranges
            - constitutional: Metadata and notes
    
    Example:
        >>> result = optimize_dose(
        ...     smiles="...",
        ...     admet={...},
        ...     protocol_template={...},
        ...     target_pk_range={"auc_mg_h_per_L": (200.0, 400.0)},
        ...     dose_bounds=(50.0, 300.0),
        ...     search_strategy="coarse_to_fine"
        ... )
        >>> print(f"Optimal dose: {result['best_regimen']['dose_mg']}mg")
    """
    
    _t0 = time.monotonic()

    if interval_options is None:
        interval_options = [24.0, 12.0, 8.0]  # QD, BID, TID

    # v2.1 HOLISTIC: No hard lock. High-tox candidates proceed to optimization; scoring and
    # whole-product assessment (e.g. in PX_Refinery) judge risk/benefit (Smart Offset, PRV).

    # Determine search strategy
    if search_strategy == "binary_search":
        # Check if metric is monotonic for binary search
        if target_pk_range and len(target_pk_range) == 1:
            metric = list(target_pk_range.keys())[0]
            if is_monotonic_metric(metric):
                return _attach_sign_off(binary_search_dose(
                    admet=admet,
                    protocol_template=protocol_template,
                    target_pk_range=target_pk_range,
                    dose_bounds=dose_bounds,
                    interval_options=interval_options,
                    variability=variability,
                    pd_params=pd_params,
                    n_eval_patients=n_eval_patients,
                ), smiles, _t0)
        # Fall back to coarse-to-fine if not suitable for binary search
        search_strategy = "coarse_to_fine"

    if search_strategy == "coarse_to_fine":
        return _attach_sign_off(coarse_to_fine_search(
            admet=admet,
            protocol_template=protocol_template,
            target_pk_range=target_pk_range,
            target_pd_range=target_pd_range,
            dose_bounds=dose_bounds,
            interval_options=interval_options,
            variability=variability,
            pd_params=pd_params,
            n_eval_patients=n_eval_patients,
        ), smiles, _t0)

    raise ValueError(f"Unknown search strategy: {search_strategy}")


def _attach_sign_off(result: Dict[str, Any], smiles: str, _t0: float) -> Dict[str, Any]:
    """Attach sign-off block to optimize_dose result."""
    _elapsed_ms = int((time.monotonic() - _t0) * 1000)
    result["sign_off"] = create_sign_off(
        engine_id="DOSE_OPT_V2",
        version="2.1-HOLISTIC",
        inputs={"smiles": smiles},
        outputs=result,
        laws_checked=["L11"],
        laws_results={"L11": result.get("best_regimen") is not None},
        execution_time_ms=_elapsed_ms,
    )
    return result


def evaluate_regimen(
    dose_mg: float,
    interval_h: float,
    admet: Dict[str, Any],
    protocol_template: Dict[str, Any],
    variability: Optional[Dict[str, Any]] = None,
    pd_params: Optional[Dict[str, Any]] = None,
    n_patients: int = 10,
    duration_days: float = 7.0,
) -> Dict[str, Any]:
    """
    Evaluate a specific dose/interval regimen via mini-trial.
    v2.1: No toxicity hard lock; holistic scoring and whole-product assessment handle risk.
    Runs a small virtual trial (5-10 patients) to assess PK/PD performance.
    
    Args:
        dose_mg: Dose in mg
        interval_h: Dosing interval in hours
        admet: ADMET parameters
        protocol_template: Base protocol structure
        variability: IIV parameters
        pd_params: PD parameters
        n_patients: Number of virtual patients
        duration_days: Trial duration
    
    Returns:
        Dict with:
            - dose_mg: Dose tested
            - interval_h: Interval tested
            - pk_summary: PK metrics (mean, std, CV%)
            - pd_summary: PD metrics (if pd_params provided)
            - trial_result: Full TrialEngine output
    """
    trial_engine = TrialEngine(time_step_h=1.0)
    
    # Build protocol for this regimen
    protocol = {
        "trial_id": f"OPT-{dose_mg}mg-Q{interval_h}h",
        "duration_days": duration_days,
        "arms": [{
            "arm_id": "OPT1",
            "label": f"{dose_mg}mg Q{interval_h}h",
            "dose_mg": dose_mg,
            "dosing_interval_h": interval_h,
            "n_patients": n_patients,
        }]
    }
    
    # Run mini-trial
    trial_result = trial_engine.run_trial(
        protocol=protocol,
        admet=admet,
        pd_params=pd_params,
        variability=variability,
    )
    
    arm = trial_result["arms"][0]
    
    # Extract PK summary with CV%
    pk_summary = {}
    for metric, values in arm["exposure_summary"].items():
        pk_summary[metric] = {
            "mean": values["mean"],
            "std": values["std"],
            "cv_percent": (values["std"] / values["mean"] * 100) if values["mean"] > 0 else 0,
            "min": values["min"],
            "max": values["max"],
        }
    
    # Extract PD summary if present
    pd_summary = None
    if "pd_summary" in arm:
        pd_summary = {}
        for metric, values in arm["pd_summary"].items():
            pd_summary[metric] = {
                "mean": values["mean"],
                "std": values["std"],
                "cv_percent": (values["std"] / values["mean"] * 100) if values["mean"] > 0 else 0,
                "min": values["min"],
                "max": values["max"],
            }
    
    return {
        "dose_mg": dose_mg,
        "interval_h": interval_h,
        "pk_summary": pk_summary,
        "pd_summary": pd_summary,
        "trial_result": trial_result,
    }


def scoring_function(
    pk_summary: Dict[str, Dict[str, float]],
    pd_summary: Optional[Dict[str, Dict[str, float]]],
    target_pk_range: Optional[Dict[str, Tuple[float, float]]],
    target_pd_range: Optional[Dict[str, Tuple[float, float]]],
    admet_data: Optional[Dict[str, Any]] = None,
) -> float:
    """
    Score a regimen based on how well it meets target ranges (mission) and safety (cost).
    v2.1 HOLISTIC: When admet_data is provided, applies dynamic toxicity/safety-margin
    penalty: high tox without massive safety margin is heavily penalized (Smart Offset).
    
    Lower score is better; 0.0 is perfect.
    """
    score = 0.0

    # 1. Target attainment (PK/PD ranges)
    if target_pk_range:
        for metric, (min_target, max_target) in target_pk_range.items():
            if metric not in pk_summary:
                continue
            
            mean_value = pk_summary[metric]["mean"]
            cv_percent = pk_summary[metric]["cv_percent"]
            
            # Distance from target range
            if mean_value < min_target:
                # Below target
                distance = (min_target - mean_value) / min_target
                score += distance ** 2
            elif mean_value > max_target:
                # Above target
                distance = (mean_value - max_target) / max_target
                score += distance ** 2
            # else: within range, no penalty
            
            # Variability penalty (prefer low CV%)
            if cv_percent > 30.0:  # High variability
                score += (cv_percent - 30.0) / 100.0
    
    # Score PD metrics
    if target_pd_range and pd_summary:
        for metric, (min_target, max_target) in target_pd_range.items():
            if metric not in pd_summary:
                continue
            
            mean_value = pd_summary[metric]["mean"]
            cv_percent = pd_summary[metric]["cv_percent"]
            
            # Distance from target range
            if mean_value < min_target:
                distance = (min_target - mean_value) / min_target
                score += distance ** 2
            elif mean_value > max_target:
                distance = (mean_value - max_target) / max_target
                score += distance ** 2
            
            # Variability penalty
            if cv_percent > 25.0:  # PD typically less variable
                score += (cv_percent - 25.0) / 100.0

    # 2. Safety/toxicity penalty (Whole Product cost). No hard lock; penalize unless margin is huge.
    if admet_data and pk_summary:
        tox = (admet_data.get("toxicity") or {}).get("toxicity_index")
        safety_margins = admet_data.get("safety_margins") or {}
        max_cmax = safety_margins.get("max_tolerated_cmax_mg_per_L")
        cmax_dist = pk_summary.get("cmax_mg_per_L", {})
        cmax_mean = cmax_dist.get("mean", 0.0) if isinstance(cmax_dist, dict) else 0.0
        if tox is not None and max_cmax is not None and cmax_mean > 0:
            safety_margin = float(max_cmax) / cmax_mean
            if safety_margin < 1.0:
                score += 1000.0  # Unsafe dose (critical failure)
            elif safety_margin < 10.0:
                score += (10.0 - safety_margin) * 0.5  # Progressive penalty
            # Whole Product: high tox allowed only if safety margin is massive
            if float(tox) >= 0.0210 and safety_margin < 50.0:
                score += 50.0  # High penalty for risky molecules without massive buffer

    return score


def is_monotonic_metric(metric: str) -> bool:
    """
    Check if a PK/PD metric is monotonic with dose.
    
    Monotonic metrics increase predictably with dose, enabling binary search.
    
    Args:
        metric: Metric name (e.g., "auc_mg_h_per_L")
    
    Returns:
        True if monotonic, False otherwise
    """
    
    # Known monotonic PK metrics
    monotonic_pk = {
        "auc_mg_h_per_L",
        "cmax_mg_per_L",
        "cmin_steady_state_mg_per_L",
    }
    
    # Known monotonic PD metrics
    monotonic_pd = {
        "max_effect",
        "auec_h",
        "mean_effect",
    }
    
    return metric in monotonic_pk or metric in monotonic_pd


def binary_search_dose(
    admet: Dict[str, Any],
    protocol_template: Dict[str, Any],
    target_pk_range: Dict[str, Tuple[float, float]],
    dose_bounds: Tuple[float, float],
    interval_options: List[float],
    variability: Optional[Dict[str, Any]],
    pd_params: Optional[Dict[str, Any]],
    n_eval_patients: int,
    tolerance: float = 0.05,
) -> Dict[str, Any]:
    """
    Binary search for optimal dose (monotonic metrics only).
    
    Efficiently finds dose within target range for monotonic metrics.
    
    Args:
        (same as optimize_dose)
        tolerance: Stop when within this fraction of target midpoint
    
    Returns:
        Same structure as optimize_dose()
    """
    
    metric = list(target_pk_range.keys())[0]
    min_target, max_target = target_pk_range[metric]
    target_midpoint = (min_target + max_target) / 2.0
    
    best_regimen = None
    best_score = float('inf')
    search_history = []
    
    # Test each interval
    for interval_h in interval_options:
        low_dose, high_dose = dose_bounds
        
        # Binary search for this interval
        while (high_dose - low_dose) / target_midpoint > tolerance:
            mid_dose = (low_dose + high_dose) / 2.0
            
            # Evaluate regimen
            result = evaluate_regimen(
                dose_mg=mid_dose,
                interval_h=interval_h,
                admet=admet,
                protocol_template=protocol_template,
                variability=variability,
                pd_params=pd_params,
                n_patients=n_eval_patients,
            )
            
            mean_value = result["pk_summary"][metric]["mean"]
            score = scoring_function(
                pk_summary=result["pk_summary"],
                pd_summary=result["pd_summary"],
                target_pk_range=target_pk_range,
                target_pd_range=None,
                admet_data=admet,
            )
            
            search_history.append({
                "dose_mg": mid_dose,
                "interval_h": interval_h,
                "metric_value": mean_value,
                "score": score,
            })
            
            if score < best_score:
                best_score = score
                best_regimen = {
                    "dose_mg": mid_dose,
                    "interval_h": interval_h,
                    "score": score,
                    "pk_summary": result["pk_summary"],
                    "pd_summary": result["pd_summary"],
                }
            
            # Adjust search bounds
            if mean_value < target_midpoint:
                low_dose = mid_dose  # Need higher dose
            else:
                high_dose = mid_dose  # Need lower dose
    
    return {
        "best_regimen": best_regimen,
        "search_history": search_history,
        "search_strategy": "binary_search",
        "evaluations": len(search_history),
        "target_achievement": {
            "metric": metric,
            "target_range": target_pk_range[metric],
            "achieved_mean": best_regimen["pk_summary"][metric]["mean"] if best_regimen else None,
        },
        "constitutional": {
            "status": "OPTIMIZED",
            "engine": "DOSE_OPTIMIZER_V2_BINARY_SEARCH",
            "notes": (
                "Binary search for monotonic metric. "
                "All evaluations based on virtual trials. "
                "L51: No fabricated values. "
                "L34: Optimization explicitly labeled VIRTUAL."
            ),
        },
    }


def coarse_to_fine_search(
    admet: Dict[str, Any],
    protocol_template: Dict[str, Any],
    target_pk_range: Optional[Dict[str, Tuple[float, float]]],
    target_pd_range: Optional[Dict[str, Tuple[float, float]]],
    dose_bounds: Tuple[float, float],
    interval_options: List[float],
    variability: Optional[Dict[str, Any]],
    pd_params: Optional[Dict[str, Any]],
    n_eval_patients: int,
) -> Dict[str, Any]:
    """
    Coarse-to-fine search for optimal dose/interval.
    
    Strategy:
    1. Coarse grid: Test 5 doses × N intervals
    2. Select top 3 candidates
    3. Fine grid: Refine around each candidate
    4. Return best overall
    
    Args:
        (same as optimize_dose)
    
    Returns:
        Same structure as optimize_dose()
    """
    
    min_dose, max_dose = dose_bounds
    
    # PHASE 1: Coarse grid
    coarse_doses = [
        min_dose,
        min_dose + (max_dose - min_dose) * 0.25,
        min_dose + (max_dose - min_dose) * 0.50,
        min_dose + (max_dose - min_dose) * 0.75,
        max_dose,
    ]
    
    coarse_results = []
    
    for dose in coarse_doses:
        for interval in interval_options:
            result = evaluate_regimen(
                dose_mg=dose,
                interval_h=interval,
                admet=admet,
                protocol_template=protocol_template,
                variability=variability,
                pd_params=pd_params,
                n_patients=n_eval_patients,
            )
            
            score = scoring_function(
                pk_summary=result["pk_summary"],
                pd_summary=result["pd_summary"],
                target_pk_range=target_pk_range,
                target_pd_range=target_pd_range,
                admet_data=admet,
            )
            
            coarse_results.append({
                "dose_mg": dose,
                "interval_h": interval,
                "score": score,
                "pk_summary": result["pk_summary"],
                "pd_summary": result["pd_summary"],
            })
    
    # Sort by score (lower is better)
    coarse_results.sort(key=lambda x: x["score"])
    
    # PHASE 2: Fine grid around top 3
    top_candidates = coarse_results[:3]
    fine_results = []
    
    for candidate in top_candidates:
        center_dose = candidate["dose_mg"]
        interval = candidate["interval_h"]
        
        # Refine around this dose (±20%)
        fine_doses = [
            center_dose * 0.85,
            center_dose * 0.925,
            center_dose,
            center_dose * 1.075,
            center_dose * 1.15,
        ]
        
        for dose in fine_doses:
            # Skip if out of bounds
            if dose < min_dose or dose > max_dose:
                continue
            
            result = evaluate_regimen(
                dose_mg=dose,
                interval_h=interval,
                admet=admet,
                protocol_template=protocol_template,
                variability=variability,
                pd_params=pd_params,
                n_patients=n_eval_patients,
            )
            
            score = scoring_function(
                pk_summary=result["pk_summary"],
                pd_summary=result["pd_summary"],
                target_pk_range=target_pk_range,
                target_pd_range=target_pd_range,
                admet_data=admet,
            )
            
            fine_results.append({
                "dose_mg": dose,
                "interval_h": interval,
                "score": score,
                "pk_summary": result["pk_summary"],
                "pd_summary": result["pd_summary"],
            })
    
    # Combine all results
    all_results = coarse_results + fine_results
    all_results.sort(key=lambda x: x["score"])
    
    best_regimen = all_results[0]
    
    # Calculate target achievement
    target_achievement = {}
    if target_pk_range:
        for metric, (min_t, max_t) in target_pk_range.items():
            achieved = best_regimen["pk_summary"][metric]["mean"]
            within_range = min_t <= achieved <= max_t
            target_achievement[f"pk_{metric}"] = {
                "target_range": (min_t, max_t),
                "achieved": achieved,
                "within_range": within_range,
            }
    
    if target_pd_range and best_regimen["pd_summary"]:
        for metric, (min_t, max_t) in target_pd_range.items():
            achieved = best_regimen["pd_summary"][metric]["mean"]
            within_range = min_t <= achieved <= max_t
            target_achievement[f"pd_{metric}"] = {
                "target_range": (min_t, max_t),
                "achieved": achieved,
                "within_range": within_range,
            }
    
    return {
        "best_regimen": best_regimen,
        "search_history": all_results,
        "search_strategy": "coarse_to_fine",
        "evaluations": len(all_results),
        "coarse_evaluations": len(coarse_results),
        "fine_evaluations": len(fine_results),
        "target_achievement": target_achievement,
        "constitutional": {
            "status": "OPTIMIZED",
            "engine": "DOSE_OPTIMIZER_V2_COARSE_TO_FINE",
            "notes": (
                "Coarse-to-fine grid search. "
                "All evaluations based on virtual trials with IIV. "
                "L51: No fabricated values. "
                "L34: Optimization explicitly labeled VIRTUAL."
            ),
        },
    }


__all__ = [
    "optimize_dose",
    "evaluate_regimen",
    "scoring_function",
    "is_monotonic_metric",
    "binary_search_dose",
    "coarse_to_fine_search",
]
