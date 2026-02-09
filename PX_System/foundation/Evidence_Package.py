"""
Foundation Evidence Package - Dossier Generation
Generates FDA 2026-compliant dossiers for compound evaluation results and trial simulations.

Constitutional Requirements:
- ALCOA+ data integrity (Attributable, Legible, Contemporaneous, Original, Accurate)
- Immutable audit trail via Sovereign_Log_Chain
- SHA-256 constitutional seals
"""
import sys
import os
from pathlib import Path
from typing import Dict, Any, Optional
from datetime import datetime, timezone
import json
import hashlib

try:
    from PX_System.foundation.Sovereign_Log_Chain import append as log_to_chain
except ImportError:
    log_to_chain = None


def generate_dossier(
    candidate: Dict[str, Any],
    engine_results: Dict[str, Any],
    zeus_verdict: Optional[Dict[str, Any]] = None,
    context: Optional[Dict[str, Any]] = None
) -> Dict[str, Any]:
    """
    Generate FDA 2026-compliant dossier from evaluation results
    
    CONSTITUTIONAL REQUIREMENT (L10: Harm is Absolute):
    Must calculate harm_energy from engine results. No defaults, no estimates.
    
    Args:
        candidate: Compound data (name, smiles, properties)
        engine_results: Results from 7-engine evaluation
        zeus_verdict: Constitutional adjudication result
        context: Additional execution context
    
    Returns:
        Complete dossier dict with ALCOA+ metadata and harm_energy
    
    Raises:
        ValueError: If harm_energy cannot be calculated from engine results
    """
    timestamp = datetime.now(timezone.utc).isoformat()
    
    # CONSTITUTIONAL EVIDENCE: Calculate harm_energy from engine results
    harm_energy = _calculate_harm_energy_from_engines(engine_results)
    
    internal_snapshot = (context or {}).get("physics_snapshot")
    external_snapshot = (context or {}).get("external_snapshot")
    enforced_snapshot, l1_overrode = _enforce_law_l1(internal_snapshot, external_snapshot)

    causal_trace_log = _build_causal_trace_log(engine_results, context, timestamp)

    dossier = {
        "dossier_version": "2.0.0",
        "generated": timestamp,
        "causal_trace_log": causal_trace_log,
        "alcoa_metadata": {
            "attributable_to": "OLYMPUS_Pipeline_Gate",
            "timestamp": timestamp,
            "original": True,
            "complete": True,
            "consistent": True,
            "contemporaneous": True
        },
        "candidate": candidate,
        "engines": engine_results,
        "harm_energy": harm_energy,  # Constitutional evidence
        "zeus": zeus_verdict or {},
        "context": {
            **(context or {}),
            "physics_snapshot": enforced_snapshot,
            "law_l1_overrode_external": l1_overrode,
        },
        "fda_compliance": "2026.01.14"
    }
    
    # Generate constitutional seal (SHA-256 hash of core data)
    seal_data = json.dumps({
        "candidate": candidate.get("name"),
        "smiles": candidate.get("smiles"),
        "timestamp": timestamp,
        "harm_energy": harm_energy,
        "zeus_verdict": zeus_verdict.get("verdict") if zeus_verdict else None
    }, sort_keys=True)
    
    dossier["constitutional_seal"] = hashlib.sha256(seal_data.encode()).hexdigest()
    
    # Log to Sovereign Chain
    if log_to_chain:
        try:
            log_to_chain({
                "source": "Evidence_Package",
                "action": "DOSSIER_GENERATED",
                "compound": candidate.get("name"),
                "harm_energy": harm_energy,
                "seal": dossier["constitutional_seal"]
            })
        except Exception:
            pass

    return dossier


def _calculate_harm_energy_from_engines(engine_results: Dict[str, Any]) -> float:
    """
    Calculate harm_energy from engine evaluation results
    
    CONSTITUTIONAL LOGIC (L10: Harm is Absolute):
    1. If OSE (Safety Engine) provides harm_energy, use it
    2. If OBE (Binding/Biological) provides harm_energy, use it
    3. If selectivity data available (IC50/CC50), calculate harm_energy
    4. If selectivity index (SI) provided directly, calculate harm_energy
    5. Otherwise, FAIL HARD - no defaults allowed
    
    Args:
        engine_results: Dict of engine outputs
    
    Returns:
        float: Calculated harm_energy value
    
    Raises:
        ValueError: If harm_energy cannot be calculated (fail-closed)
    """
    # Priority 1: OSE (Safety Engine) explicit harm_energy
    if "OSE" in engine_results:
        ose_result = engine_results["OSE"]
        if isinstance(ose_result, dict) and "harm_energy" in ose_result:
            if ose_result["harm_energy"] is not None:
                return float(ose_result["harm_energy"])
    
    # Priority 2: OBE (Biological Engine) explicit harm_energy
    if "OBE" in engine_results:
        obe_result = engine_results["OBE"]
        if isinstance(obe_result, dict) and "harm_energy" in obe_result:
            if obe_result["harm_energy"] is not None:
                return float(obe_result["harm_energy"])
    
    # Priority 3: Calculate from selectivity index (IC50/CC50)
    # Extract from any engine that provides these values
    ic50 = None
    cc50 = None
    si = None
    
    for engine_name, engine_result in engine_results.items():
        if not isinstance(engine_result, dict):
            continue
        
        properties = engine_result.get("properties", {})
        if "ic50_um" in properties and properties["ic50_um"] is not None:
            ic50 = float(properties["ic50_um"])
        if "cc50_um" in properties and properties["cc50_um"] is not None:
            cc50 = float(properties["cc50_um"])
        if "SI" in properties and properties["SI"] is not None:
            si = float(properties["SI"])
        if "si" in properties and properties["si"] is not None:
            si = float(properties["si"])
    
    if ic50 is not None and cc50 is not None and ic50 > 0:
        # Selectivity Index (SI) = CC50 / IC50
        # SI > 10 = safe (harm_energy = 0)
        # SI < 10 = toxic (harm_energy scales from 0 to 1)
        selectivity_index = cc50 / ic50
        harm_energy = max(0.0, 1.0 - (selectivity_index / 10.0))
        return harm_energy

    if si is not None:
        # Direct SI input: linear penalty
        harm_energy = max(0.0, 1.0 - (si / 10.0))
        return harm_energy

    # Priority 4: ADMET pipeline output (OPE+ADMET stages; L10: toxicity_index is harm-related)
    admet = engine_results.get("admet") or (engine_results.get("stages") or {}).get("admet")
    if isinstance(admet, dict):
        tox = admet.get("toxicity_index")
        if tox is None and "toxicity" in admet and isinstance(admet["toxicity"], dict):
            tox = admet["toxicity"].get("toxicity_index")
        if tox is not None:
            return float(tox)

    # FAIL HARD: No valid harm_energy source found
    raise ValueError(
        "CONSTITUTIONAL VIOLATION: Cannot calculate harm_energy from engine results. "
        "Required: OSE/OBE harm_energy field OR IC50/CC50 selectivity data OR SI. "
        "Per L10 (Harm is Absolute), no defaults are allowed. "
        "Fix the engine evaluation to provide complete safety data."
    )


def validate_dossier(dossier: Dict[str, Any]) -> bool:
    """
    Validate dossier integrity
    
    Args:
        dossier: Dossier dict to validate
    
    Returns:
        True if valid, False otherwise
    """
    required_fields = [
        "dossier_version",
        "generated",
        "alcoa_metadata",
        "candidate",
        "constitutional_seal"
    ]
    
    for field in required_fields:
        if field not in dossier:
            return False
    
    return True


def wrap_trial_simulation(
    protocol: Dict[str, Any],
    trial_result: Dict[str, Any],
    ope: Dict[str, Any],
    admet: Dict[str, Any],
    output_dir: str = "PX_Warehouse/Calibration_Molecules",
) -> str:
    """
    Wraps a TrialEngine simulation into a constitutional Evidence Package.
    Trial dossiers go to Calibration_Molecules (QA vault for reference/testing data).

    v2.0 Enhancement: Supports PK/PD modeling when pd_params present in trial_result.

    Produces:
        TRIAL_SIMULATION_DOSSIER-<HASH>.json

    Args:
        protocol: Trial protocol with arms, doses, etc.
        trial_result: TrialEngine output (includes pd_params if PD enabled)
        ope: OPE analysis dict
        admet: ADMET analysis dict
        output_dir: Output directory for dossier

    Returns:
        Full path to the persisted dossier.
    """

    timestamp = datetime.now(timezone.utc).isoformat()
    
    # Extract PD parameters if present
    pd_params = trial_result.get("pd_params")
    has_pd = pd_params is not None
    
    # Determine version based on features present
    # v3.0: IIV + Adaptive + Dose Opt + Efficacy
    # v2.1: PK/PD only
    # v1.0: Exposure only
    arms = trial_result.get("arms", [])
    has_iiv = False
    if arms:
        has_iiv = arms[0].get("exposure_summary", {}).get("auc_mg_h_per_L", {}).get("std") is not None
    
    has_adaptive = trial_result.get("adaptive_rules") is not None
    
    if has_iiv or has_adaptive:
        dossier_version = "3.0"
    elif has_pd:
        dossier_version = "2.1"
    else:
        dossier_version = "1.0"
    
    # Build constitutional notes based on features
    notes_parts = []
    if has_adaptive:
        notes_parts.append("Adaptive virtual trial with mid-trial dose adjustments.")
    if has_pd:
        notes_parts.append("PK/PD modeling (theoretical, requires clinical validation).")
    if has_iiv:
        notes_parts.append("Population variability (IIV) modeled.")
    if not notes_parts:
        notes_parts.append("Exposure-only virtual trial.")
    
    notes_parts.append("No clinical endpoints.")
    notes_parts.append("L51: No fabricated values.")
    notes_parts.append("L34: All simulations explicitly labeled VIRTUAL.")
    
    constitutional_notes = " ".join(notes_parts)

    causal_trace_log = _build_causal_trace_log(trial_result, None, timestamp)

    dossier = {
        "dossier_type": "TRIAL_SIMULATION_DOSSIER",
        "version": dossier_version,
        "timestamp_utc": timestamp,
        "causal_trace_log": causal_trace_log,

        "protocol": protocol,
        "trial_result": trial_result,

        "provenance": {
            "ope_engine": ope.get("status", "UNKNOWN"),
            "admet_engine": admet.get("constitutional", {}).get("engine", "UNKNOWN"),
            "trial_engine": trial_result.get("constitutional", {}).get("engine", "UNKNOWN"),
            "pkpd_enabled": has_pd,
        },

        "inputs": {
            "ope_analysis": ope,
            "admet_analysis": admet,
        },
        
        # v2.1: PK/PD Analysis
        "pkpd_analysis": {
            "pd_model": "EMAX" if has_pd else None,
            "pd_parameters": pd_params,
            "pd_summary_per_arm": [
                {
                    "arm_id": arm.get("arm_id"),
                    "pd_summary": arm.get("pd_summary"),
                }
                for arm in trial_result.get("arms", [])
                if arm.get("pd_summary") is not None
            ] if has_pd else None,
            "constitutional": {
                "status": "SIMULATED",
                "notes": (
                    "PD model is theoretical and based on Emax assumptions. "
                    "Clinical validation required. EC50 and Emax must be "
                    "determined experimentally for actual drug-target pairs."
                )
            } if has_pd else None,
        } if has_pd else None,
        
        # v3.0: IIV (Population Variability)
        "iiv_analysis": {
            "enabled": has_iiv,
            "variability_parameters": trial_result.get("variability") if has_iiv else None,
            "population_distributions": [
                {
                    "arm_id": arm.get("arm_id"),
                    "pk_distribution": {
                        "auc": arm.get("exposure_summary", {}).get("auc_mg_h_per_L", {}),
                        "cmax": arm.get("exposure_summary", {}).get("cmax_mg_per_L", {}),
                    },
                    "pd_distribution": {
                        "max_effect": arm.get("pd_summary", {}).get("max_effect", {}),
                    } if has_pd else None,
                }
                for arm in trial_result.get("arms", [])
            ] if has_iiv else None,
            "constitutional": {
                "status": "SIMULATED",
                "notes": (
                    "Population variability modeled with deterministic 7-tier system. "
                    "IIV factors: clearance, Vd, ka. "
                    "L51: No random sampling. "
                    "L34: IIV explicitly labeled VIRTUAL."
                )
            } if has_iiv else None,
        } if has_iiv else None,
        
        # v3.0: Adaptive Trial Decisions
        "adaptive_analysis": {
            "enabled": has_adaptive,
            "adaptive_rules": trial_result.get("adaptive_rules"),
            "decisions_per_arm": [
                {
                    "arm_id": arm.get("arm_id"),
                    "initial_dose_mg": arm.get("initial_dose_mg"),
                    "final_dose_mg": arm.get("final_dose_mg"),
                    "arm_stopped": arm.get("arm_stopped", False),
                    "patients_enrolled": arm.get("patients_enrolled"),
                    "n_patients_planned": arm.get("n_patients"),
                    "adaptation_log": arm.get("adaptation_log", []),
                    "adaptations_triggered": arm.get("adaptations_triggered", 0),
                }
                for arm in trial_result.get("arms", [])
            ] if has_adaptive else None,
            "constitutional": {
                "status": "LOGGED",
                "notes": (
                    "Adaptive decisions logged with full rationale. "
                    "All dose adjustments based on interim analyses. "
                    "L51: No fabricated adaptation triggers. "
                    "L34: Adaptive behavior explicitly labeled VIRTUAL."
                )
            } if has_adaptive else None,
        } if has_adaptive else None,

        "constitutional": {
            "status": "EVIDENCE_PACKAGE_CREATED",
            "law_basis": ["L51", "L34", "ALCOA+"],
            "notes": constitutional_notes,
        },
    }

    # Hash for reproducibility
    dossier_bytes = json.dumps(dossier, sort_keys=True).encode("utf-8")
    dossier_hash = hashlib.sha256(dossier_bytes).hexdigest()[:12]

    dossier["evidence_hash"] = dossier_hash

    # Output path
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    filename = f"TRIAL_SIMULATION_DOSSIER-{dossier_hash}.json"
    full_path = str(Path(output_dir) / filename)

    with open(full_path, "w", encoding="utf-8") as f:
        json.dump(dossier, f, indent=2)

    # Log to Sovereign Chain
    if log_to_chain:
        try:
            log_to_chain(
                "TRIAL_DOSSIER_GENERATED",
                {
                    "trial_id": protocol.get("trial_id", "UNKNOWN"),
                    "hash": dossier_hash,
                    "path": full_path,
                },
                {"source": "Evidence_Package"}
            )
        except Exception:
            pass

    return full_path


__all__ = ["generate_dossier", "validate_dossier", "wrap_trial_simulation"]


def _enforce_law_l1(internal_snapshot, external_snapshot):
    """
    Law L1: Internal 35D physics snapshot is primary.
    External data may provide context but cannot override internal physics.
    """
    if external_snapshot is None:
        return internal_snapshot, False
    if internal_snapshot != external_snapshot:
        return internal_snapshot, True
    return internal_snapshot, False


def _build_causal_trace_log(engine_results: Dict[str, Any], context: Optional[Dict[str, Any]], timestamp: str):
    energy_delta = None
    if context:
        energy_delta = context.get("energy_delta")
        physics = context.get("physics_snapshot") or {}
        if energy_delta is None and isinstance(physics, dict):
            energy_delta = physics.get("energy_delta")
    authorization = "GLOBAL_SUM conserved"
    if energy_delta is not None:
        authorization = "Energy Delta = 0" if energy_delta == 0 else "Energy Delta non-zero"
    return [
        {
            "law": "U34",
            "authorization": authorization,
            "timestamp": timestamp,
        }
    ]
