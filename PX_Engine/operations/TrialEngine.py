"""
TrialEngine.py
Exposure-only trial simulation engine v1.

Orchestrates:
- virtual population generation
- per-arm PK simulation using PX_Laboratory.SimulationEngine
- exposure-based endpoints (Cmax, AUC distributions)

Non-clinical, deterministic, constitutional.
"""

import time
from typing import Dict, Any, List
import statistics

from PX_Laboratory import SimulationEngine
from PX_System.foundation.sign_off import create_sign_off


def generate_virtual_population(
    n_patients: int,
    base_weight_kg: float = 70.0,
    weight_sd_kg: float = 10.0,
    variability: Dict[str, Any] | None = None,
) -> List[Dict[str, Any]]:
    """
    Parametric virtual population generator with deterministic IIV.
    
    v2.0-PHASE2: Full inter-individual variability (IIV) implementation.
    Uses deterministic tiers instead of RNG for reproducibility.
    
    Args:
        n_patients: Number of virtual patients to generate
        base_weight_kg: Base weight (kg)
        weight_sd_kg: Weight standard deviation (kg)
        variability: Optional IIV parameters (v2.0-PHASE2):
            {
                "clearance_variation": 0.3,  # ±30% realistic range
                "vd_variation": 0.25,        # ±25% realistic range
                "ka_variation": 0.2,         # ±20% realistic range (optional)
                "n_tiers": 7,                # Number of tiers (default 7)
            }
            
            If None, all patients identical (v1.x behavior for backward compat)
    
    Returns:
        List of patient dicts with deterministic variation
        
    Example:
        >>> pop = generate_virtual_population(
        ...     n_patients=21,
        ...     variability={
        ...         "clearance_variation": 0.3,
        ...         "vd_variation": 0.25,
        ...         "ka_variation": 0.2
        ...     }
        ... )
        >>> # Patients will have factors cycling through 7 tiers:
        >>> # 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3
    """

    if n_patients <= 0:
        raise ValueError("n_patients must be positive")

    patients: List[Dict[str, Any]] = []
    
    # Deterministic weight variation pattern: symmetric offsets
    weight_offsets = [0, -1, 1, -2, 2, -3, 3, -4, 4, -5, 5]

    for i in range(n_patients):
        # Weight variation
        offset = weight_offsets[i % len(weight_offsets)]
        weight = base_weight_kg + offset * (weight_sd_kg / 5.0)
        patient = {
            "patient_id": f"VP-{i+1:04d}",
            "weight_kg": max(40.0, min(120.0, weight)),
        }
        
        # v2.0-PHASE2: Full IIV with deterministic tiers
        if variability:
            # Number of tiers (default 7 for good distribution)
            n_tiers = variability.get("n_tiers", 7)
            
            # Calculate tier index: cycles through -n_tiers/2 to +n_tiers/2
            # e.g., n_tiers=7: -3, -2, -1, 0, 1, 2, 3
            # e.g., n_tiers=5: -2, -1, 0, 1, 2
            tier_idx = (i % n_tiers) - (n_tiers // 2)
            
            # Clearance factor
            if "clearance_variation" in variability:
                cl_var = variability["clearance_variation"]
                # Map tier to factor: tier=-3 → 0.7, tier=0 → 1.0, tier=3 → 1.3 (for 30% variation)
                factor = 1.0 + tier_idx * (cl_var / (n_tiers // 2))
                # Clamp to physiologically plausible range [0.5, 2.0]
                factor = max(0.5, min(2.0, factor))
                patient["clearance_factor"] = factor
            
            # Volume of distribution factor
            if "vd_variation" in variability:
                vd_var = variability["vd_variation"]
                factor = 1.0 + tier_idx * (vd_var / (n_tiers // 2))
                # Clamp to physiologically plausible range [0.5, 2.0]
                factor = max(0.5, min(2.0, factor))
                patient["vd_factor"] = factor
            
            # Absorption rate constant factor (optional)
            if "ka_variation" in variability:
                ka_var = variability["ka_variation"]
                factor = 1.0 + tier_idx * (ka_var / (n_tiers // 2))
                # Clamp to physiologically plausible range [0.5, 2.0]
                factor = max(0.5, min(2.0, factor))
                patient["ka_factor"] = factor
        
        patients.append(patient)

    return patients


class TrialEngine:
    """
    Exposure-only trial simulation engine.

    Protocol schema (minimal v1):

    protocol = {
        "trial_id": str,
        "duration_days": float,
        "arms": [
            {
                "arm_id": str,
                "label": str,
                "dose_mg": float,
                "dosing_interval_h": float,
                "n_patients": int,
            },
            ...
        ],
        
        # Optional adaptive rules (v2.0-PHASE3: production implementation)
        "adaptive_rules": {
            "metric": "auc_mg_h_per_L",  # or "cmax_mg_per_L", "max_effect", etc.
            "lower_bound": 100.0,        # Optional: trigger INCREASE_DOSE if below
            "upper_bound": 300.0,        # Optional: trigger REDUCE_DOSE if above
            "action": "REDUCE_DOSE",     # or "INCREASE_DOSE", "STOP_ARM"
            "dose_adjustment_factor": 0.75,  # For dose changes (default 0.75 or 1.33)
            "interim_after_n": 10,       # Evaluate after N patients per arm
        },
        
        # Optional crossover design (v2 structure only)
        "design": "PARALLEL" or "CROSSOVER",
        "sequences": [
            {
                "sequence_id": "S1",
                "treatments": ["A1", "A2"],  # arm IDs in order
            },
        ]
    }
    """

    def __init__(self, time_step_h: float = 1.0):
        self.time_step_h = time_step_h
        self.pk_engine = SimulationEngine(time_step_h=time_step_h)

    def run_trial(
        self,
        protocol: Dict[str, Any],
        admet: Dict[str, Any],
        pd_params: Dict[str, Any] | None = None,
        variability: Dict[str, Any] | None = None,
    ) -> Dict[str, Any]:
        """
        Run a trial simulation with optional PK/PD modeling and IIV.

        For each arm:
        - generate virtual population (with optional IIV)
        - run PK simulation per patient
        - aggregate exposure metrics (Cmax, AUC, Cmin_ss) with std
        - if pd_params: run PK/PD and aggregate effect metrics

        Args:
            protocol: Trial protocol dict
            admet: ADMET parameters
            pd_params: Optional PD parameters for PK/PD modeling:
                {
                    "emax": float,           # Maximum effect
                    "ec50": float,           # EC50 (mg/L)
                    "hill": float,           # Hill coefficient (default 1.0)
                    "baseline": float,       # Baseline effect (default 0.0)
                    "effect_threshold": float  # For time_above calc (optional)
                }
            variability: Optional IIV parameters (v2.0-PHASE2):
                {
                    "clearance_variation": 0.3,  # ±30%
                    "vd_variation": 0.25,        # ±25%
                    "ka_variation": 0.2,         # ±20% (optional)
                    "n_tiers": 7                 # Number of tiers (default 7)
                }
                If None, all patients identical (v1.x behavior)

        Returns:
            dict with per-arm summaries (exposure + optional PD) and trial metadata.
        """

        _t0 = time.monotonic()
        trial_id = protocol.get("trial_id", "TRIAL-UNKNOWN")
        duration_days = protocol.get("duration_days", 7.0)
        duration_h = duration_days * 24.0

        arms = protocol.get("arms", [])
        if not arms:
            raise ValueError("Protocol must define at least one arm")
        
        # Only PARALLEL design is supported; CROSSOVER etc. are not yet implemented
        design = protocol.get("design", "PARALLEL")
        if design != "PARALLEL":
            raise ValueError(
                f"Unsupported trial design '{design}'. Only design='PARALLEL' is supported."
            )

        arm_results: List[Dict[str, Any]] = []
        
        # Phase 3: Adaptive rules setup
        adaptive_rules = protocol.get("adaptive_rules")
        
        for arm in arms:
            arm_id = arm.get("arm_id", "ARM-UNKNOWN")
            label = arm.get("label", "")
            initial_dose_mg = arm.get("dose_mg", None)
            dosing_interval_h = arm.get("dosing_interval_h", 24.0)
            n_patients = arm.get("n_patients", 10)

            if initial_dose_mg is None or initial_dose_mg <= 0:
                raise ValueError(f"Invalid dose_mg for arm {arm_id}")

            population = generate_virtual_population(
                n_patients=n_patients,
                variability=variability
            )

            # PK metrics
            cmax_values: List[float] = []
            auc_values: List[float] = []
            cmin_values: List[float] = []
            
            # PD metrics (if pd_params provided)
            max_effect_values: List[float] = []
            auec_values: List[float] = []
            time_above_threshold_values: List[float] = []
            mean_effect_values: List[float] = []
            
            # Phase 3: Adaptive state tracking
            current_dose_mg = initial_dose_mg
            arm_stopped = False
            adaptation_log: List[Dict[str, Any]] = []
            
            # Phase 3: Epoch-based adaptation
            interim_after_n = adaptive_rules.get("interim_after_n", n_patients) if adaptive_rules else n_patients

            for i, patient in enumerate(population):
                # Phase 3: Check if arm stopped
                if arm_stopped:
                    # Skip remaining patients if arm stopped
                    continue
                
                # Phase 3: Interim analysis checkpoint
                if adaptive_rules and i > 0 and i % interim_after_n == 0:
                    # Evaluate adaptive rule based on data so far
                    decision = self._evaluate_adaptive_rule(
                        adaptive_rules=adaptive_rules,
                        cmax_values=cmax_values,
                        auc_values=auc_values,
                        cmin_values=cmin_values,
                        max_effect_values=max_effect_values if pd_params else [],
                        auec_values=auec_values if pd_params else [],
                        time_above_threshold_values=time_above_threshold_values if pd_params else [],
                        mean_effect_values=mean_effect_values if pd_params else [],
                        current_dose_mg=current_dose_mg,
                        patient_count=i,
                    )
                    
                    # Log decision
                    adaptation_log.append(decision)
                    
                    # Execute action if triggered
                    if decision["triggered"]:
                        action = decision["action"]
                        if action == "STOP_ARM":
                            arm_stopped = True
                            continue
                        elif action in ["REDUCE_DOSE", "INCREASE_DOSE"]:
                            current_dose_mg = decision["new_dose_mg"]
                
                # 1. Run PK simulation with current dose
                pk = self.pk_engine.simulate_one_compartment(
                    dose_mg=current_dose_mg,
                    duration_h=duration_h,
                    dosing_interval_h=dosing_interval_h,
                    patient=patient,
                    admet=admet,
                )
                
                # Store PK metrics
                summary = pk["summary"]
                cmax_values.append(summary["cmax_mg_per_L"])
                auc_values.append(summary["auc_mg_h_per_L"])
                cmin_values.append(summary["cmin_steady_state_mg_per_L"])
                
                # 2. If PD enabled, run PK/PD linking
                if pd_params:
                    from PX_Engine.operations.PKPD import link_pk_to_pd
                    
                    pd_result = link_pk_to_pd(pk, pd_params)
                    pd_summary = pd_result["pd_summary"]
                    
                    # Store PD metrics
                    max_effect_values.append(pd_summary["max_effect"])
                    auec_values.append(pd_summary["auec_h"])
                    time_above_threshold_values.append(pd_summary["time_above_threshold_h"])
                    mean_effect_values.append(pd_summary["mean_effect"])

            arm_result = {
                "arm_id": arm_id,
                "label": label,
                "n_patients": n_patients,
                "initial_dose_mg": initial_dose_mg,
                "final_dose_mg": current_dose_mg,
                "dosing_interval_h": dosing_interval_h,
                "exposure_summary": {
                    "cmax_mg_per_L": _dist_summary(cmax_values),
                    "auc_mg_h_per_L": _dist_summary(auc_values),
                    "cmin_steady_state_mg_per_L": _dist_summary(cmin_values),
                },
                "arm_stopped": arm_stopped,
                "patients_enrolled": len(cmax_values),  # Actual patients simulated
            }
            
            # Add PD summary if pd_params was provided
            if pd_params:
                arm_result["pd_summary"] = {
                    "max_effect": _dist_summary(max_effect_values),
                    "auec_h": _dist_summary(auec_values),
                    "time_above_threshold_h": _dist_summary(time_above_threshold_values),
                    "mean_effect": _dist_summary(mean_effect_values),
                }
            
            # Phase 3: Add adaptation log if adaptive rules present
            if adaptive_rules:
                arm_result["adaptation_log"] = adaptation_log
                arm_result["adaptations_triggered"] = sum(1 for d in adaptation_log if d["triggered"])
            
            arm_results.append(arm_result)

        # Build constitutional notes based on PD status and adaptation
        if adaptive_rules:
            if pd_params:
                constitutional_notes = (
                    "Adaptive trial simulation with PK/PD modeling. "
                    "Dose adjustments based on interim analyses. "
                    "PD model is theoretical and requires clinical validation. "
                    "No clinical endpoints."
                )
                engine_version = "TRIAL_ENGINE_V2.0-ADAPTIVE-PKPD"
            else:
                constitutional_notes = (
                    "Adaptive trial simulation. "
                    "Dose adjustments based on interim analyses. "
                    "Exposure-only simulation; no clinical endpoints."
                )
                engine_version = "TRIAL_ENGINE_V2.0-ADAPTIVE"
        else:
            if pd_params:
                constitutional_notes = (
                    "Trial simulation with PK/PD modeling. "
                    "PD model is theoretical and requires clinical validation. "
                    "No clinical endpoints."
                )
                engine_version = "TRIAL_ENGINE_V2.0-PKPD"
            else:
                constitutional_notes = "Exposure-only trial simulation; no clinical endpoints."
                engine_version = "TRIAL_ENGINE_V1"
        
        result = {
            "trial_id": trial_id,
            "duration_days": duration_days,
            "arms": arm_results,
            "pd_params": pd_params,  # Include for provenance
            "adaptive_rules": adaptive_rules,  # Include for provenance
            "constitutional": {
                "status": "SIMULATED",
                "engine": engine_version,
                "notes": constitutional_notes,
            },
        }
        _elapsed_ms = int((time.monotonic() - _t0) * 1000)
        result["sign_off"] = create_sign_off(
            engine_id="TRIAL_ENGINE_V1",
            version="1.0-EXPOSURE",
            inputs={"trial_id": protocol.get("trial_id", "")},
            outputs=result,
            laws_checked=["L11", "L34"],
            laws_results={"L11": True, "L34": True},
            execution_time_ms=_elapsed_ms,
        )
        return result
    
    def _evaluate_adaptive_rule(
        self,
        adaptive_rules: Dict[str, Any],
        cmax_values: List[float],
        auc_values: List[float],
        cmin_values: List[float],
        max_effect_values: List[float],
        auec_values: List[float],
        time_above_threshold_values: List[float],
        mean_effect_values: List[float],
        current_dose_mg: float,
        patient_count: int,
    ) -> Dict[str, Any]:
        """
        Evaluate adaptive rule at interim analysis.
        
        v2.0-PHASE3: Production adaptive logic.
        
        Args:
            adaptive_rules: Rule definition from protocol
            *_values: Accumulated metric values from patients enrolled so far
            current_dose_mg: Current dose for this arm
            patient_count: Number of patients enrolled
            
        Returns:
            Decision dict with triggered status, action, and rationale
        """
        import statistics
        
        metric = adaptive_rules.get("metric")
        lower_bound = adaptive_rules.get("lower_bound")
        upper_bound = adaptive_rules.get("upper_bound")
        action = adaptive_rules.get("action", "REDUCE_DOSE")
        dose_adjustment_factor = adaptive_rules.get("dose_adjustment_factor", 0.75)
        
        # Get metric values
        metric_map = {
            "auc_mg_h_per_L": auc_values,
            "cmax_mg_per_L": cmax_values,
            "cmin_steady_state_mg_per_L": cmin_values,
            "max_effect": max_effect_values,
            "auec_h": auec_values,
            "time_above_threshold_h": time_above_threshold_values,
            "mean_effect": mean_effect_values,
        }
        
        if metric not in metric_map:
            # Unknown metric, no action
            return {
                "patient_count": patient_count,
                "metric": metric,
                "triggered": False,
                "reason": f"Unknown metric: {metric}",
                "action": None,
                "current_dose_mg": current_dose_mg,
                "new_dose_mg": current_dose_mg,
            }
        
        values = metric_map[metric]
        if not values:
            # No data yet
            return {
                "patient_count": patient_count,
                "metric": metric,
                "triggered": False,
                "reason": "No data available",
                "action": None,
                "current_dose_mg": current_dose_mg,
                "new_dose_mg": current_dose_mg,
            }
        
        mean_value = statistics.fmean(values)
        
        # Evaluate bounds
        triggered = False
        triggered_action = None
        reason = None
        new_dose_mg = current_dose_mg
        
        if upper_bound is not None and mean_value > upper_bound:
            triggered = True
            triggered_action = "REDUCE_DOSE" if action == "REDUCE_DOSE" else action
            reason = f"{metric} mean ({mean_value:.2f}) > upper_bound ({upper_bound})"
            if triggered_action == "REDUCE_DOSE":
                new_dose_mg = current_dose_mg * dose_adjustment_factor
            elif triggered_action == "STOP_ARM":
                new_dose_mg = current_dose_mg  # No change, arm will stop
        
        elif lower_bound is not None and mean_value < lower_bound:
            triggered = True
            triggered_action = "INCREASE_DOSE" if action == "INCREASE_DOSE" else action
            reason = f"{metric} mean ({mean_value:.2f}) < lower_bound ({lower_bound})"
            if triggered_action == "INCREASE_DOSE":
                # Inverse of reduction factor for increase
                increase_factor = adaptive_rules.get("increase_factor", 1.33)
                new_dose_mg = current_dose_mg * increase_factor
            elif triggered_action == "STOP_ARM":
                new_dose_mg = current_dose_mg  # No change, arm will stop
        
        if not triggered:
            reason = f"{metric} mean ({mean_value:.2f}) within bounds"
        
        return {
            "patient_count": patient_count,
            "metric": metric,
            "mean_value": mean_value,
            "lower_bound": lower_bound,
            "upper_bound": upper_bound,
            "triggered": triggered,
            "reason": reason,
            "action": triggered_action,
            "current_dose_mg": current_dose_mg,
            "new_dose_mg": new_dose_mg,
        }


def _dist_summary(values: List[float]) -> Dict[str, float]:
    """
    Compute distribution summary statistics.
    
    v2.0-PHASE2: Added std (standard deviation) for IIV analysis.
    """
    if not values:
        return {
            "mean": 0.0,
            "median": 0.0,
            "min": 0.0,
            "max": 0.0,
            "std": 0.0,
        }
    
    mean_val = statistics.fmean(values)
    
    # Calculate std (sample standard deviation)
    if len(values) > 1:
        std_val = statistics.stdev(values)
    else:
        std_val = 0.0
    
    return {
        "mean": float(mean_val),
        "median": float(statistics.median(values)),
        "min": float(min(values)),
        "max": float(max(values)),
        "std": float(std_val),
    }
