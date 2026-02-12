"""
Finalization Pipeline — When a dossier is "Finalized" in Predator X.

Governed by PX_Warehouse/Finalization_Spec.FINALIZATION_SPEC. A dossier is only
finalized when it has all finalized_requires (promotion, MoA, disease context,
PRV scoring, lineage, ancestral_trace, physics_map, trial_simulation binding,
worldline, discovery grading, constitutional overrides, Zeus verdict, seals,
tier, ALCOA, causal_trace_log, fda_compliance). All engines are used.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional, Tuple

_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in __import__("sys").path:
    __import__("sys").path.insert(0, str(_REPO_ROOT))

from PX_System.foundation.quint.converter import emit
from PX_System.foundation.quint.engine_adapter import (
    q_run_ole, q_run_grading, q_run_dose_optimizer, q_run_virtual_efficacy,
    q_run_trial,
)

# Required keys that must be present on a dossier to be finalizable (base Evidence Package)
_REQUIRED_BASE = ["candidate", "engines", "harm_energy", "constitutional_seal", "dossier_version"]


def _get_spec():
    try:
        from PX_Warehouse.Finalization_Spec import FINALIZATION_SPEC
        return FINALIZATION_SPEC
    except Exception:
        return {}


def _ensure_repo_path():
    if str(_REPO_ROOT) not in __import__("sys").path:
        import sys
        sys.path.insert(0, str(_REPO_ROOT))


def _get_engines(dossier: Dict[str, Any]) -> Tuple[Dict, Dict]:
    """Return (ope, admet) from dossier."""
    engines = dossier.get("engines") or {}
    ope = engines.get("ope") or {}
    admet = engines.get("admet") or {}
    return ope, admet


def _add_moa_hypothesis(dossier: Dict[str, Any], repo_root: Path) -> str:
    """Generate MoA hypothesis from OPE/ADMET and optional Disease_Constraint_Model."""
    ope, admet = _get_engines(dossier)
    parts = []
    if ope.get("ec50") is not None:
        parts.append(f"OPE EC50 estimate {ope.get('ec50'):.4f}; binding_affinity_nM={ope.get('binding_affinity_nM')}")
    tox = (admet.get("toxicity") or {}) if isinstance(admet.get("toxicity"), dict) else {}
    if tox:
        parts.append(f"Toxicity index {tox.get('toxicity_index', 'N/A')}; risk_level={tox.get('risk_level', 'N/A')}")
    try:
        from PX_System.foundation.Disease_Constraint_Model import get_disease_constraints
        d = get_disease_constraints("NiV_Malaysia")
        if d:
            parts.append(f"Disease constraints: ic50_max_um={d.get('ic50_max_um')}, SI_min={d.get('selectivity_index_min')}")
    except Exception:
        pass
    return " | ".join(parts) if parts else "Mechanism inferred from OPE/ADMET; disease constraints not applied"


def _add_disease_space_anchoring(dossier: Dict[str, Any], repo_root: Path) -> list:
    """Anchor to PRV disease space from manifest."""
    manifest_path = repo_root / "PX_Domain" / "PRV_Diseases" / "manifest.json"
    if not manifest_path.exists():
        return ["PRV disease manifest not found; anchoring deferred"]
    try:
        data = json.loads(manifest_path.read_text(encoding="utf-8"))
        diseases = data.get("prv_diseases") or []
        return [d.get("name") for d in diseases if d.get("prv_eligible")][:15]
    except Exception:
        return ["Anchoring failed"]


def _add_disease_context_prv_category_rationale(
    dossier: Dict[str, Any], repo_root: Path, ole_result: Dict[str, Any]
) -> Tuple[List[str], str, str]:
    """Return (disease_context list, prv_category, prv_rationale)."""
    manifest_path = repo_root / "PX_Domain" / "PRV_Diseases" / "manifest.json"
    disease_context = []
    prv_category = "STANDARD_FDA"
    prv_rationale = "Indication not in PRV-eligible set"
    if manifest_path.exists():
        try:
            data = json.loads(manifest_path.read_text(encoding="utf-8"))
            diseases = data.get("prv_diseases") or []
            disease_context = [d.get("name") for d in diseases if d.get("prv_eligible")][:20]
        except Exception:
            pass
    if ole_result.get("prv_eligible"):
        prv_category = ole_result.get("regulatory_pathway", "PRV_ACCELERATED")
        prv_rationale = "Indication in PRV-eligible disease space (tropical/rare_pediatric/mcm)"
    return disease_context, prv_category, prv_rationale


def _add_physics_map(dossier: Dict[str, Any]) -> Tuple[str, str]:
    """Compute physics_map_id and physics_map_hash from engines (OPE+ADMET)."""
    engines = dossier.get("engines") or {}
    payload = {
        "ope": engines.get("ope"),
        "admet": engines.get("admet"),
    }
    raw = json.dumps(payload, sort_keys=True)
    h = hashlib.sha256(raw.encode()).hexdigest()
    return f"PM_{h[:12]}", h


def _run_virtual_trial_for_dossier(dossier: Dict[str, Any], item_id: str) -> Tuple[Optional[Dict[str, Any]], Optional[Dict[str, Any]], Optional[Dict[str, Any]]]:
    """
    Run full in-silico trial when dossier has OPE/ADMET but no trial_result/virtual_efficacy.
    Returns (trial_result, virtual_efficacy, dose_optimization) or (None, None, None) on skip/error.
    """
    _ensure_repo_path()
    engines = dossier.get("engines") or {}
    if engines.get("trial_result") and engines.get("virtual_efficacy"):
        return None, None, None
    ope = engines.get("ope") or {}
    admet = engines.get("admet") or {}
    if not admet:
        return None, None, None
    try:
        ec50 = float(ope.get("ec50", 1.0) or 1.0)
        emax = float(ope.get("emax", 0.8) or 0.8)
        pd_params = {"emax": emax, "ec50": ec50, "hill": 1.0, "baseline": 0.0, "effect_threshold": 0.5}
        variability = {"clearance_variation": 0.3, "vd_variation": 0.25, "ka_variation": 0.2, "n_tiers": 7}
        protocol = {
            "trial_id": f"FINAL-{item_id}",
            "duration_days": 7.0,
            "arms": [{"arm_id": "A1", "label": "100mg Q24h", "dose_mg": 100.0, "dosing_interval_h": 24.0, "n_patients": 21}],
        }
        trial_qf = q_run_trial(protocol, admet, pd_params, variability, time_step_h=1.0)
        trial_result = emit(trial_qf, target_label="finalization_trial")
        dose_optimization = None
        smiles = (dossier.get("candidate") or {}).get("smiles") or ""
        if smiles:
            try:
                dose_qf = q_run_dose_optimizer(
                    smiles, admet, protocol, pd_params, n_eval_patients=10,
                )
                dose_optimization = emit(dose_qf, target_label="finalization_dose")
                best = (dose_optimization or {}).get("best_regimen") or {}
                d_mg = best.get("dose_mg", 100.0)
                d_int = best.get("dosing_interval_h", 24.0)
                protocol_opt = {
                    "trial_id": f"PTA-{item_id}",
                    "duration_days": 7.0,
                    "arms": [{"arm_id": "P1", "label": f"{d_mg}mg Q{d_int}h", "dose_mg": d_mg, "dosing_interval_h": d_int, "n_patients": 21}],
                }
                trial_qf2 = q_run_trial(protocol_opt, admet, pd_params, variability, time_step_h=1.0)
                trial_result = emit(trial_qf2, target_label="finalization_trial_optimized")
            except Exception:
                pass
        ve_qf = q_run_virtual_efficacy(trial_result, metric="auc_mg_h_per_L", target_threshold=200.0)
        ve_raw = emit(ve_qf, target_label="finalization_virtual_efficacy")
        pta_auc = ve_raw.get("pta", {})
        responder = ve_raw.get("responder_rate", {})
        ve = {
            "pta": pta_auc,
            "pta_auc": pta_auc,
            "pk_pta": {"auc_mg_h_per_L": pta_auc} if isinstance(pta_auc, dict) else {},
            "responder_rate": responder,
            "variability": ve_raw.get("variability"),
        }
        return trial_result, ve, dose_optimization
    except Exception:
        return None, None, None


def _add_ancestral_trace(item_id: str, worldline_path: Optional[str], first_seen: str) -> Dict[str, Any]:
    """Build ancestral_trace: full_chain, first_seen, worldline_ids."""
    full_chain = [item_id]
    worldline_ids = [worldline_path] if worldline_path else []
    return {
        "full_chain": full_chain,
        "first_seen": first_seen,
        "worldline_ids": worldline_ids,
    }


def _add_trial_simulation_binding(dossier: Dict[str, Any], item_id: str) -> Dict[str, Any]:
    """
    Trial simulation binding: id, hash, full trial_outcome_summary per TRIAL_INTEGRATED_DOSSIER_SPEC.
    When trial/virtual_efficacy data present: real id/hash and status COMPLETED with responder_rate, pta,
    toxicity_rate, effect_size, uncertainty_bounds, cohort_statistics. Otherwise NOT_RUN with same structure.
    """
    engines = dossier.get("engines") or {}
    trial_result = dossier.get("trial_result") or engines.get("trial_result")
    virtual_efficacy = engines.get("virtual_efficacy") or dossier.get("virtual_efficacy")
    pkpd = engines.get("pkpd") or dossier.get("pkpd")
    dose_opt = engines.get("dose_optimization") or dossier.get("dose_optimization")
    admet = engines.get("admet") or dossier.get("admet") or {}

    # Build full outcome summary structure (spec 2)
    def _empty_outcome() -> Dict[str, Any]:
        return {
            "status": "VIRTUAL_TRIAL_NOT_RUN",
            "note": "Trial simulation optional at finalization",
            "responder_rate": None,
            "pta": None,
            "toxicity_rate": None,
            "effect_size": None,
            "uncertainty_bounds": None,
            "cohort_statistics": None,
        }

    if trial_result and trial_result.get("arms"):
        raw = json.dumps(trial_result, sort_keys=True, default=str)
        h = hashlib.sha256(raw.encode()).hexdigest()[:16]
        trial_hash = hashlib.sha256(raw.encode()).hexdigest()
        arms = trial_result.get("arms", [])
        # Extract from arms/virtual_efficacy/pkpd/admet
        responder_rate = None
        pta_val = None
        if virtual_efficacy:
            rr = virtual_efficacy.get("responder_rate") or virtual_efficacy.get("pd_responders") or {}
            if isinstance(rr, dict) and "effect_ge_threshold" in rr:
                responder_rate = round(float(rr["effect_ge_threshold"]) * 100.0, 4)
            elif isinstance(rr, dict) and rr.get("max_effect", {}).get("responder_rate") is not None:
                responder_rate = round(float(rr["max_effect"]["responder_rate"]) * 100.0, 4)
            pk_pta = virtual_efficacy.get("pk_pta") or virtual_efficacy.get("pk_pta_auc")
            if isinstance(pk_pta, dict) and "auc_mg_h_per_L" in pk_pta:
                pta_val = round(float(pk_pta["auc_mg_h_per_L"].get("pta", 0)), 4)
            elif isinstance(pk_pta, (int, float)):
                pta_val = round(float(pk_pta), 4)
        if pta_val is None and pkpd and isinstance(pkpd.get("auc"), dict):
            pta_val = round(float(pkpd["auc"].get("pta", 0)), 4)
        tox_dict = admet.get("toxicity") if isinstance(admet.get("toxicity"), dict) else {}
        toxicity_rate = None
        if tox_dict and tox_dict.get("toxicity_index") is not None:
            toxicity_rate = round(float(tox_dict["toxicity_index"]), 6)
        if toxicity_rate is None and dossier.get("harm_energy") is not None:
            toxicity_rate = round(float(dossier["harm_energy"]), 6)
        n_total = sum(arm.get("patients_enrolled", arm.get("n_patients", 0)) for arm in arms)
        cohort_statistics = {"n_arms": len(arms), "n_total": n_total}
        outcome_summary = {
            "status": "COMPLETED",
            "responder_rate": responder_rate,
            "pta": pta_val,
            "toxicity_rate": toxicity_rate,
            "effect_size": virtual_efficacy.get("effect_size") if isinstance(virtual_efficacy, dict) else None,
            "uncertainty_bounds": virtual_efficacy.get("uncertainty_bounds") if isinstance(virtual_efficacy, dict) else None,
            "cohort_statistics": cohort_statistics,
            "arms": len(arms),
        }
        return {
            "trial_simulation_id": f"TRIAL_{item_id}_{h}",
            "trial_simulation_hash": trial_hash,
            "trial_outcome_summary": outcome_summary,
        }

    # Virtual efficacy only (no trial_result arms) — still treat as simulated if we have key metrics
    if virtual_efficacy or (pkpd and admet):
        payload = {"ve": virtual_efficacy, "pkpd": pkpd, "admet": admet, "item_id": item_id}
        raw = json.dumps(payload, sort_keys=True, default=str)
        h = hashlib.sha256(raw.encode()).hexdigest()[:16]
        trial_hash = hashlib.sha256(raw.encode()).hexdigest()
        rr = (virtual_efficacy or {}).get("responder_rate") or (virtual_efficacy or {}).get("pd_responders") or {}
        responder_rate = None
        if isinstance(rr, dict) and "effect_ge_threshold" in rr:
            responder_rate = round(float(rr["effect_ge_threshold"]) * 100.0, 4)
        pta_val = None
        pk_pta = (virtual_efficacy or {}).get("pk_pta") or (virtual_efficacy or {}).get("pk_pta_auc")
        if isinstance(pk_pta, dict) and pk_pta.get("auc_mg_h_per_L"):
            pta_val = round(float(pk_pta["auc_mg_h_per_L"].get("pta", 0)), 4)
        tox_dict = admet.get("toxicity") if isinstance(admet.get("toxicity"), dict) else {}
        toxicity_rate = round(float(tox_dict.get("toxicity_index", dossier.get("harm_energy") or 0.5)), 6) if tox_dict or dossier.get("harm_energy") is not None else None
        outcome_summary = {
            "status": "COMPLETED",
            "responder_rate": responder_rate,
            "pta": pta_val,
            "toxicity_rate": toxicity_rate,
            "effect_size": (virtual_efficacy or {}).get("effect_size"),
            "uncertainty_bounds": (virtual_efficacy or {}).get("uncertainty_bounds"),
            "cohort_statistics": {"n_arms": 0, "n_total": 0},
            "arms": 0,
        }
        return {
            "trial_simulation_id": f"TRIAL_{item_id}_{h}",
            "trial_simulation_hash": trial_hash,
            "trial_outcome_summary": outcome_summary,
        }

    outcome = _empty_outcome()
    return {
        "trial_simulation_id": f"TRIAL_{item_id}_NOT_RUN",
        "trial_simulation_hash": hashlib.sha256(f"{item_id}_VIRTUAL_NOT_RUN".encode()).hexdigest(),
        "trial_outcome_summary": outcome,
    }


def _build_trial_results_block(dossier: Dict[str, Any], trial_binding: Dict[str, Any]) -> Dict[str, Any]:
    """Spec 10: trial_results — cohort summary, dose–response table, PTA curve, toxicity curve, effect size, uncertainty intervals, recommended_dose, therapeutic_window."""
    summary = trial_binding.get("trial_outcome_summary") or {}
    if summary.get("status") != "COMPLETED":
        return {"status": "VIRTUAL_TRIAL_NOT_RUN", "cohort_summary": None, "dose_response_table": None, "pta_curve": None, "toxicity_curve": None, "effect_size": None, "uncertainty_intervals": None, "recommended_dose": None, "therapeutic_window": None}
    engines = dossier.get("engines") or {}
    ve = engines.get("virtual_efficacy") or dossier.get("virtual_efficacy") or {}
    trial_result = dossier.get("trial_result") or engines.get("trial_result") or {}
    dose_opt = engines.get("dose_optimization") or dossier.get("dose_optimization") or {}
    arms = trial_result.get("arms", [])
    cohort_summary = {
        "n_arms": len(arms),
        "n_total": summary.get("cohort_statistics") or {},
        "responder_rate": summary.get("responder_rate"),
        "pta": summary.get("pta"),
        "toxicity_rate": summary.get("toxicity_rate"),
    }
    dose_response_table = None
    if arms:
        dose_response_table = [{"arm_id": a.get("arm_id"), "dose": a.get("final_dose_mg") or a.get("initial_dose_mg"), "exposure_summary": a.get("exposure_summary"), "pd_summary": a.get("pd_summary")} for a in arms]
    pta_curve = ve.get("pk_pta") or ve.get("pta_curve")
    toxicity_curve = None
    if isinstance(engines.get("admet"), dict) and isinstance(engines["admet"].get("toxicity"), dict):
        toxicity_curve = {"toxicity_index": engines["admet"]["toxicity"].get("toxicity_index")}
    best_reg = dose_opt.get("best_regimen") or {}
    recommended_dose = {"dose_mg": best_reg.get("dose_mg"), "interval_h": best_reg.get("dosing_interval_h"), "score": best_reg.get("score")} if best_reg else None
    therapeutic_window = None
    if ve.get("pk_pta") or ve.get("pd_responders"):
        therapeutic_window = {"pk_target_auc": 200.0, "pd_target_max_effect": 0.5, "safety_threshold": 0.95}
    return {
        "status": "COMPLETED",
        "cohort_summary": cohort_summary,
        "dose_response_table": dose_response_table,
        "pta_curve": pta_curve,
        "toxicity_curve": toxicity_curve,
        "effect_size": summary.get("effect_size") or ve.get("effect_size"),
        "uncertainty_intervals": summary.get("uncertainty_bounds") or ve.get("uncertainty_bounds"),
        "recommended_dose": recommended_dose,
        "therapeutic_window": therapeutic_window,
    }


def _get_effect_size_vs_soc(dossier: Dict[str, Any], soc_benchmarking: Dict[str, Any]) -> Dict[str, Any]:
    """Spec 7: effect_size_vs_soc — positive/neutral/negative, magnitude, confidence interval, superiority/non-inferiority classification."""
    engines = dossier.get("engines") or {}
    ve = engines.get("virtual_efficacy") or dossier.get("virtual_efficacy") or {}
    effect_size = ve.get("effect_size")
    better = soc_benchmarking.get("better_than_soc")
    if effect_size is None:
        return {"direction": "neutral", "magnitude": None, "confidence_interval": None, "superiority": None, "non_inferiority": None, "note": "No trial effect size", "better_than_soc": better}
    try:
        es = float(effect_size)
        direction = "positive" if es > 0 else ("negative" if es < 0 else "neutral")
        superiority = "superior" if es > 0.1 and better else ("non_superior" if es <= 0.1 else None)
        non_inferiority = "non_inferior" if es >= 0 and not (es < -0.05) else ("inferior" if es < -0.05 else None)
    except (TypeError, ValueError):
        direction = "neutral"
        es = None
        superiority = None
        non_inferiority = None
    return {"direction": direction, "magnitude": es, "confidence_interval": ve.get("uncertainty_bounds"), "superiority": superiority, "non_inferiority": non_inferiority, "better_than_soc": better}


def _get_trial_adjusted_prv_score(prv_eligibility: Dict[str, Any], trial_outcome_summary: Dict[str, Any]) -> Dict[str, Any]:
    """Spec 6: trial_adjusted_prv_score — efficacy_weight, safety_weight, novelty_weight, soc_weight."""
    if (trial_outcome_summary or {}).get("status") != "COMPLETED":
        return {"efficacy_weight": None, "safety_weight": None, "novelty_weight": None, "soc_weight": None, "note": "Trial not run"}
    return {
        "efficacy_weight": 0.35,
        "safety_weight": 0.35,
        "novelty_weight": 0.15,
        "soc_weight": 0.15,
        "prv_eligible": prv_eligibility.get("prv_eligible", False),
    }


def _get_worldline_trial_event(trial_binding: Dict[str, Any], finalized_at: str, tier_after_trial: str) -> Dict[str, Any]:
    """Spec 5: worldline_trial_event — timestamp, simulation_id, trial_hash, outcome_summary, tier_after_trial."""
    summary = trial_binding.get("trial_outcome_summary") or {}
    if summary.get("status") != "COMPLETED":
        return {"timestamp": finalized_at, "simulation_id": trial_binding.get("trial_simulation_id"), "trial_hash": trial_binding.get("trial_simulation_hash"), "outcome_summary": summary, "status": "NOT_RUN", "tier_after_trial": None}
    return {
        "timestamp": finalized_at,
        "simulation_id": trial_binding.get("trial_simulation_id"),
        "trial_hash": trial_binding.get("trial_simulation_hash"),
        "outcome_summary": summary,
        "status": "COMPLETED",
        "tier_after_trial": tier_after_trial,
    }


def _get_trial_adjusted_sa_score(synthetic_accessibility_score: Dict[str, Any], trial_outcome_summary: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    """Spec 8 (optional): trial_adjusted_sa_score."""
    if (trial_outcome_summary or {}).get("status") != "COMPLETED":
        return None
    base = (synthetic_accessibility_score or {}).get("score")
    if base is None:
        return None
    return {"score": base, "trial_adjusted": False, "note": "Optional; high dose requirements could lower SA in future"}


def _get_tier_from_trial_outcome(trial_outcome_summary: Dict[str, Any], discovery_tier: str) -> str:
    """Spec 4: when trial status COMPLETED, optionally refine tier from trial outcomes."""
    if (trial_outcome_summary or {}).get("status") != "COMPLETED":
        return discovery_tier
    pta = trial_outcome_summary.get("pta")
    rr = trial_outcome_summary.get("responder_rate")
    tox = trial_outcome_summary.get("toxicity_rate")
    if pta is not None and rr is not None and tox is not None:
        try:
            pta_f, rr_f, tox_f = float(pta), float(rr), float(tox)
            if tox_f >= 0.0210:
                return "Bronze"
            if pta_f >= 30 and rr_f >= 30 and tox_f < 0.021:
                return "Diamond" if discovery_tier == "Diamond" else "Gold"
            if pta_f >= 20 and rr_f >= 20:
                return "Gold" if discovery_tier in ("Diamond", "Gold") else "Silver"
            if pta_f >= 10 and rr_f >= 10:
                return "Silver" if discovery_tier in ("Diamond", "Gold", "Silver") else "Bronze"
        except (TypeError, ValueError):
            pass
    return discovery_tier


def _add_prv_eligibility(dossier: Dict[str, Any], repo_root: Path) -> Dict[str, Any]:
    """Run OLE for PRV eligibility scoring (via QUINT)."""
    _ensure_repo_path()
    try:
        candidate = dossier.get("candidate") or {}
        payload = {
            "compound_id": candidate.get("name") or "UNKNOWN",
            "indication": "PRV_eligible_disease_space",
        }
        ole_qf = q_run_ole(payload)
        return emit(ole_qf, target_label="finalization_ole")
    except Exception as e:
        return {"prv_eligible": False, "status": "OLE_ERROR", "error": str(e)}


def _add_discovery_grading(dossier: Dict[str, Any], repo_root: Path) -> Dict[str, Any]:
    """Run GradingEngine for discovery-stage grade (via QUINT)."""
    _ensure_repo_path()
    try:
        engines = dossier.get("engines") or {}
        d = {**dossier, "ope": engines.get("ope"), "admet": engines.get("admet")}
        grade_qf = q_run_grading(d)
        return emit(grade_qf, target_label="finalization_grading")
    except Exception as e:
        return {
            "grade": "NEEDS_REVIEW",
            "reason": str(e),
            "metrics": {
                "pta": 0.0,
                "responder_rate": 0.0,
                "toxicity": 0.5,
                "dose_score": 0.5,
                "variability_cv": 0.0,
                "binding_affinity": 0.5,
            },
            "reasoning": {"grade": "NEEDS_REVIEW", "reason": str(e), "criteria_met_count": 0},
        }


# ---------------------------------------------------------------------------
# Grading metrics: SoC comparison, novelty fingerprint, synthetic accessibility
# ---------------------------------------------------------------------------

def _compute_soc_benchmarking_score(
    dossier: Dict[str, Any], disease_context: List[str]
) -> Dict[str, Any]:
    """
    Standard-of-Care comparison: deterministic comparison against best-in-class
    for the disease space (binding affinity vs reference). Score 0-1, 1 = beats SoC.
    """
    engines = dossier.get("engines") or {}
    ope = engines.get("ope") or {}
    candidate_nM = ope.get("binding_affinity_nM")
    if candidate_nM is None:
        candidate_nM = ope.get("ec50", 100.0) * 1000.0  # ec50 in uM -> nM proxy
    try:
        candidate_nM = float(candidate_nM)
    except (TypeError, ValueError):
        candidate_nM = 100.0
    # Reference: best-in-class proxy (nM). Lower affinity = better. Typical SoC 1-50 nM.
    soc_reference_nM = 10.0
    if disease_context:
        # Optional: load disease-specific SoC from manifest or table later
        pass
    better_than_soc = candidate_nM <= soc_reference_nM
    if candidate_nM <= 0:
        score = 1.0
    else:
        ratio = candidate_nM / max(1e-6, soc_reference_nM)
        score = max(0.0, min(1.0, 1.0 - (ratio - 1.0) * 0.2))  # 1 at SoC, decays above
    return {
        "score": round(score, 6),
        "candidate_binding_nM": round(candidate_nM, 4),
        "soc_reference_nM": soc_reference_nM,
        "better_than_soc": better_than_soc,
        "note": "Deterministic comparison vs best-in-class reference for disease space",
    }


def _compute_novelty_fingerprint_score(
    dossier: Dict[str, Any], repo_root: Path
) -> Dict[str, Any]:
    """
    IP Readiness: uniqueness vs known SMILES (NoveltyBudgetEngine + optional ref set).
    Score 0-1, higher = more novel / patent-eligible basis.
    """
    _ensure_repo_path()
    candidate = dossier.get("candidate") or {}
    reference_set: Optional[List[Dict[str, Any]]] = None
    known_path = repo_root / "PX_Data" / "known_smiles.json"
    if known_path.exists():
        try:
            data = json.loads(known_path.read_text(encoding="utf-8"))
            reference_set = data if isinstance(data, list) else data.get("smiles_list", [])
        except Exception:
            pass
    try:
        from PX_System.foundation.Novelty_Budget_Engine import NoveltyBudgetEngine
        engine = NoveltyBudgetEngine(max_evaluations=1)
        score = engine.calculate_novelty_score(candidate, reference_set=reference_set)
    except Exception as e:
        score = 0.5
        reference_set = None
    return {
        "score": round(float(score), 6),
        "reference_set_used": reference_set is not None and len(reference_set) > 0,
        "note": "Uniqueness vs known SMILES; basis for patent eligibility",
    }


def _compute_synthetic_accessibility_score(dossier: Dict[str, Any]) -> Dict[str, Any]:
    """
    Economic Scalability: ensure molecule is synthesizable (RDKit SA_Score or proxy).
    SA raw 1-10 (lower = easier). synthetic_accessibility_score 0-1 (1 = easy).
    """
    smiles = (dossier.get("candidate") or {}).get("smiles") or ""
    sa_raw: Optional[float] = None
    source = "proxy"
    try:
        from rdkit import Chem
        from rdkit.Chem.rdMolDescriptors import CalcNumRings, CalcNumAtomStereoCenters
        mol = Chem.MolFromSmiles(smiles) if smiles else None
        if mol is not None:
            n_rings = CalcNumRings(mol)
            n_stereo = CalcNumAtomStereoCenters(mol)
            complexity = n_rings * 0.5 + n_stereo * 1.0 + (mol.GetNumHeavyAtoms() / 50.0)
            sa_raw = min(10.0, max(1.0, 1.0 + complexity))
        else:
            sa_raw = 5.0
    except Exception:
        try:
            from rdkit import Chem
            mol = Chem.MolFromSmiles(smiles) if smiles else None
            sa_raw = min(10.0, max(1.0, 1.0 + (mol.GetNumHeavyAtoms() / 20.0 if mol else 5.0)))
        except Exception:
            sa_raw = 5.0
    if sa_raw is None:
        sa_raw = 5.0
    # Invert: 1 (easy) -> score 1.0; 10 (hard) -> score 0.0
    accessibility = max(0.0, min(1.0, 1.0 - (sa_raw - 1.0) / 9.0))
    return {
        "score": round(accessibility, 6),
        "sa_raw": round(sa_raw, 4),
        "source": source,
        "note": "Synthetic accessibility; 1 = lab-feasible",
    }


def _add_constitutional_overrides(dossier: Dict[str, Any], zeus_verdict: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """Log constitutional risk overrides (L1, Zeus)."""
    ctx = dossier.get("context") or {}
    zeus = zeus_verdict or dossier.get("zeus") or {}
    return {
        "law_l1_overrode_external": ctx.get("law_l1_overrode_external", False),
        "zeus_verdict": zeus.get("rationale") if isinstance(zeus.get("rationale"), str) else (zeus.get("authorized") and "PASS" or "FAIL"),
        "zeus_laws_results": zeus.get("laws_results") if isinstance(zeus.get("laws_results"), dict) else None,
        "harm_energy": dossier.get("harm_energy"),
    }


def _add_rep_seal(finalized_payload: Dict[str, Any]) -> str:
    """REP seal: hash of critical finalized fields."""
    seal_data = json.dumps({
        "candidate": finalized_payload.get("candidate", {}).get("name"),
        "harm_energy": finalized_payload.get("harm_energy"),
        "promotion_stage": finalized_payload.get("finalization", {}).get("promotion_stage"),
        "discovery_grade": finalized_payload.get("finalization", {}).get("discovery_grading", {}).get("grade"),
        "constitutional_seal": finalized_payload.get("constitutional_seal"),
        "finalized_at": finalized_payload.get("finalization", {}).get("finalized_at"),
    }, sort_keys=True)
    return hashlib.sha256(seal_data.encode()).hexdigest()


def _ensure_worldline(
    item_id: str,
    dossier: Dict[str, Any],
    tier: str,
    is_novel: bool,
    repo_root: Path,
) -> Optional[str]:
    """Ensure a worldline exists for this dossier; create if missing. Return path or None."""
    _ensure_repo_path()
    from PX_Warehouse.warehouse_layout import get_worldline_dir, TIERS
    wl_dir = get_worldline_dir(tier if tier in TIERS else "Bronze", repo_root)
    safe_id = "".join(c if c.isalnum() or c in "-_" else "_" for c in str(item_id))[:64]
    path = wl_dir / f"{safe_id}.worldline"
    if path.exists():
        return str(path)
    try:
        from PX_Warehouse.WorldLine_Database import WorldLineDatabase
        wl_db = WorldLineDatabase()
        return wl_db.record_discovery(item_id, dossier, tier, is_novel, repo_root)
    except Exception:
        return None


def _validate_finalized_dossier(
    finalized: Dict[str, Any], spec: Dict[str, Any]
) -> Tuple[bool, List[str], Dict[str, Any]]:
    """
    Check finalized dossier has all finalized_requires.
    Returns (ok, missing_or_errors, spec_compliance_report).
    """
    required = spec.get("finalized_requires") or []
    missing: List[str] = []
    passed_required: List[str] = []
    fin = finalized.get("finalization") or {}
    for key in required:
        if "." in key:
            part, sub = key.split(".", 1)
            obj = fin.get(part) if part == "ancestral_trace" else finalized.get(part)
            if not isinstance(obj, dict) or sub not in obj:
                missing.append(key)
            else:
                passed_required.append(key)
        else:
            if key in ("engines", "alcoa_metadata", "causal_trace_log", "fda_compliance", "constitutional_seal"):
                if key not in finalized:
                    missing.append(key)
                else:
                    passed_required.append(key)
            elif key not in fin and key not in finalized:
                missing.append(key)
            else:
                passed_required.append(key)
    checklist = spec.get("finalization_checklist") or {}
    must_pass = checklist.get("must_pass") or []
    required_missing_count = len(missing)
    failed_checks: List[str] = []
    passed_checks: List[str] = list(passed_required)
    for c in must_pass:
        fail = False
        if c == "all_required_fields_present":
            fail = required_missing_count > 0
        elif c == "zeus_gate_passed":
            if not fin.get("zeus_verdict", {}).get("authorized"):
                missing.append("zeus_gate_passed")
                fail = True
        elif c == "tier_computed":
            if not fin.get("tier"):
                missing.append("tier_computed")
                fail = True
        elif c == "rep_seal_generated":
            if not fin.get("rep_seal"):
                missing.append("rep_seal_generated")
                fail = True
        elif c == "worldline_written":
            if not fin.get("worldline_path"):
                missing.append("worldline_written")
                fail = True
        elif c == "grading_metrics_complete":
            for k in ("soc_benchmarking_score", "novelty_fingerprint_score", "synthetic_accessibility_score"):
                if k not in fin:
                    missing.append(k)
            fail = any(k not in fin for k in ("soc_benchmarking_score", "novelty_fingerprint_score", "synthetic_accessibility_score"))
        if fail:
            failed_checks.append(c)
        else:
            passed_checks.append(c)
    # Always record Zeus gate in report for audit (write is gated by zeus_verdict.authorized)
    if not fin.get("zeus_verdict", {}).get("authorized"):
        failed_checks.append("zeus_gate_passed")
    else:
        passed_checks.append("zeus_gate_passed")
    report = {
        "missing_fields": list(missing),
        "failed_checks": failed_checks,
        "passed_checks": passed_checks,
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }
    return (len(missing) == 0, missing, report)


def run_finalization(
    dossier: Dict[str, Any],
    item_id: str,
    is_novel: bool,
    repo_root: Optional[Path] = None,
) -> Tuple[Optional[Dict[str, Any]], Optional[str], Optional[str]]:
    """
    Run the FULL finalization pipeline for every dossier (NO_SINGLE_METRIC_REJECTION).
    Does not write to disk. Zeus gate does not cause early return; all metrics,
    worldline, and spec_compliance_report are computed. Only the caller gates
    WRITE to Finalized_Dossiers on zeus_verdict.authorized.

    Args:
        dossier: Evidence Package dossier (from Prv_Dossiers or Novel_Dossiers).
        item_id: Candidate id (e.g. PRV_NOV_xxx, PRV_REP_xxx).
        is_novel: True if novel, False if repurposed.
        repo_root: Repo root; default from this file.

    Returns:
        (finalized_dossier_dict, tier, error_message)
        On full pipeline success: (dict, tier, None). Caller checks zeus_verdict.authorized to gate write.
        On invalid input only: (None, None, reason).
    """
    root = Path(repo_root) if repo_root else _REPO_ROOT
    _ensure_repo_path()
    spec = _get_spec()
    zeus_cfg = (spec.get("zeus_gate") or {})
    tox_threshold = zeus_cfg.get("toxicity_threshold", 0.0210)
    harm_threshold = zeus_cfg.get("harm_energy_threshold", 0.0210)

    for key in _REQUIRED_BASE:
        if key not in dossier:
            return None, None, f"Missing required field: {key}"

    # Upgrade: run full in-silico trial when dossier has OPE/ADMET but no trial/virtual_efficacy
    tr, ve, do = _run_virtual_trial_for_dossier(dossier, item_id)
    if tr is not None or ve is not None or do is not None:
        eng = dict(dossier.get("engines") or {})
        if tr is not None:
            eng["trial_result"] = tr
        if ve is not None:
            eng["virtual_efficacy"] = ve
        if do is not None:
            eng["dose_optimization"] = do
        dossier = {**dossier, "engines": eng}

    # IMPLEMENTATION_ROADMAP: No early exits. Zeus runs at the END; all stages run first.
    # Tier (placement)
    try:
        from PX_Warehouse.warehouse_layout import get_tier
        tier = get_tier(dossier)
    except Exception:
        tier = "Bronze"

    # Promotion stage: NOV → DISC → REP
    promotion_stage = "NOV" if is_novel else "REP"
    finalized_at = datetime.now(timezone.utc).isoformat()

    # MoA hypothesis
    moa_hypothesis = _add_moa_hypothesis(dossier, root)

    # Disease-space anchoring (legacy key) + disease_context, prv_category, prv_rationale
    disease_anchoring = _add_disease_space_anchoring(dossier, root)
    prv_eligibility = _add_prv_eligibility(dossier, root)
    disease_context, prv_category, prv_rationale = _add_disease_context_prv_category_rationale(dossier, root, prv_eligibility)

    # Discovery grading (GradingEngine)
    discovery_grading = _add_discovery_grading(dossier, root)

    # Grading metrics: SoC comparison, novelty fingerprint, synthetic accessibility
    soc_benchmarking_score = _compute_soc_benchmarking_score(dossier, disease_context)
    novelty_fingerprint_score = _compute_novelty_fingerprint_score(dossier, root)
    synthetic_accessibility_score = _compute_synthetic_accessibility_score(dossier)

    # Lineage + ancestral_trace
    lineage_assignment = {
        "source": "prv_24h",
        "item_id": item_id,
        "stage": promotion_stage,
        "is_novel": is_novel,
    }
    worldline_path = _ensure_worldline(item_id, dossier, tier, is_novel, root)
    ancestral_trace = _add_ancestral_trace(item_id, worldline_path, finalized_at)

    # Physics map binding
    physics_map_id, physics_map_hash = _add_physics_map(dossier)

    # Trial simulation binding (spec 2: full trial_outcome_summary)
    trial_binding = _add_trial_simulation_binding(dossier, item_id)
    # Spec 4: recompute tier from trial outcomes when COMPLETED
    tier = _get_tier_from_trial_outcome(trial_binding.get("trial_outcome_summary") or {}, tier)

    # Grade/tier consistency gate: NEEDS_REVIEW and REJECTED grades cannot go to Diamond or Gold
    _grade = (discovery_grading or {}).get("grade", "")
    if _grade in ("NEEDS_REVIEW", "REJECTED"):
        if tier in ("Diamond", "Gold"):
            tier = "Silver" if _grade == "NEEDS_REVIEW" else "Bronze"

    # Governance stage: Zeus gate at the END (IMPLEMENTATION_ROADMAP). Spec 9: pass trial_outcome_summary for Zeus trial review.
    try:
        from PX_System.foundation.ZeusLaws import run_zeus_gate
        dossier_for_zeus = {**dossier, "trial_outcome_summary": trial_binding.get("trial_outcome_summary")}
        zeus_verdict = run_zeus_gate(dossier_for_zeus, toxicity_threshold=tox_threshold, harm_energy_threshold=harm_threshold)
    except Exception as e:
        zeus_verdict = {"authorized": False, "rationale": str(e), "laws_required": [], "laws_results": {}}

    # Constitutional overrides (include full zeus_verdict)
    constitutional_overrides = _add_constitutional_overrides(dossier, zeus_verdict)

    # Ensure ALCOA/traceability on base dossier
    if "alcoa_metadata" not in dossier:
        dossier = {**dossier, "alcoa_metadata": dossier.get("alcoa_metadata") or {"attributable_to": "finalization", "timestamp": finalized_at}}
    if "causal_trace_log" not in dossier:
        dossier = {**dossier, "causal_trace_log": dossier.get("causal_trace_log") or []}
    if "fda_compliance" not in dossier:
        dossier = {**dossier, "fda_compliance": dossier.get("fda_compliance") or "2026.01.14"}

    # Spec 10: trial_results block; spec 7: effect_size_vs_soc; spec 6: trial_adjusted_prv_score; spec 5: worldline_trial_event; spec 8: trial_adjusted_sa_score
    trial_results_block = _build_trial_results_block(dossier, trial_binding)
    effect_size_vs_soc = _get_effect_size_vs_soc(dossier, soc_benchmarking_score)
    trial_adjusted_prv = _get_trial_adjusted_prv_score(prv_eligibility, trial_binding.get("trial_outcome_summary") or {})
    worldline_trial_event = _get_worldline_trial_event(trial_binding, finalized_at, tier)
    trial_adjusted_sa = _get_trial_adjusted_sa_score(synthetic_accessibility_score, trial_binding.get("trial_outcome_summary") or {})

    finalization = {
        "promotion_stage": promotion_stage,
        "moa_hypothesis": moa_hypothesis,
        "disease_space_anchoring": disease_anchoring,
        "disease_context": disease_context,
        "prv_category": prv_category,
        "prv_rationale": prv_rationale,
        "prv_eligibility_scoring": prv_eligibility,
        "lineage_assignment": lineage_assignment,
        "ancestral_trace": ancestral_trace,
        "physics_map_id": physics_map_id,
        "physics_map_hash": physics_map_hash,
        "trial_simulation_id": trial_binding["trial_simulation_id"],
        "trial_simulation_hash": trial_binding["trial_simulation_hash"],
        "trial_outcome_summary": trial_binding["trial_outcome_summary"],
        "trial_results": trial_results_block,
        "effect_size_vs_soc": effect_size_vs_soc,
        "trial_adjusted_prv_score": trial_adjusted_prv,
        "worldline_trial_event": worldline_trial_event,
        "trial_adjusted_sa_score": trial_adjusted_sa,
        "discovery_grading": discovery_grading,
        "soc_benchmarking_score": soc_benchmarking_score,
        "novelty_fingerprint_score": novelty_fingerprint_score,
        "synthetic_accessibility_score": synthetic_accessibility_score,
        "constitutional_risk_overrides": constitutional_overrides,
        "zeus_verdict": zeus_verdict,
        "worldline_generated": worldline_path is not None,
        "worldline_path": worldline_path,
        "tier": tier,
        "finalized_at": finalized_at,
    }

    finalized_dossier = {**dossier, "finalization": finalization}
    rep_seal = _add_rep_seal(finalized_dossier)
    finalized_dossier["finalization"]["rep_seal"] = rep_seal

    # Validate and build spec compliance report; store the report for final state (after version + report added)
    if spec:
        finalization_version = spec.get("finalization_version", "1.0.0")
        finalized_dossier["finalization"]["finalization_version"] = finalization_version
        ok1, missing1, report1 = _validate_finalized_dossier(finalized_dossier, spec)
        finalized_dossier["finalization"]["spec_compliance_report"] = report1
        ok2, missing2, report2 = _validate_finalized_dossier(finalized_dossier, spec)
        finalized_dossier["finalization"]["spec_compliance_report"] = report2
        if not ok2:
            return None, None, f"Finalization checklist missing: {missing2}"
    return finalized_dossier, tier, None


def write_finalized_dossier(
    finalized_dossier: Dict[str, Any],
    item_id: str,
    tier: str,
    repo_root: Optional[Path] = None,
) -> Path:
    """
    Write a finalized dossier to PX_Warehouse/Finalized_Dossiers/<tier>/<safe_id>.json.
    Returns the path written.
    """
    root = Path(repo_root) if repo_root else _REPO_ROOT
    from PX_Warehouse.warehouse_layout import get_finalized_dossier_dir, TIERS
    t = tier if tier in TIERS else "Bronze"
    out_dir = get_finalized_dossier_dir(t, root)
    out_dir.mkdir(parents=True, exist_ok=True)
    safe_id = "".join(c if c.isalnum() or c in "-_" else "_" for c in str(item_id))[:64]
    out_path = out_dir / f"{safe_id}.json"
    out_path.write_text(json.dumps(finalized_dossier, indent=2, default=str), encoding="utf-8")
    return out_path


def passes_discovery_bar(finalized_dossier: Dict[str, Any]) -> bool:
    """
    True if whole-product toxicity <= DISCOVERY_ACCEPTED_TOX_THRESHOLD (5%).
    Used for PATH B: Zeus fail but tox <= 0.05 -> Discovery_Accepted.
    """
    from PX_Warehouse.warehouse_layout import DISCOVERY_ACCEPTED_TOX_THRESHOLD
    zv = (finalized_dossier.get("finalization") or {}).get("zeus_verdict") or {}
    tox = zv.get("toxicity_index")
    if tox is None:
        tox = finalized_dossier.get("harm_energy")
    try:
        t = float(tox) if tox is not None else 1.0
        return t <= DISCOVERY_ACCEPTED_TOX_THRESHOLD
    except (TypeError, ValueError):
        return False


def write_discovery_accepted_dossier(
    finalized_dossier: Dict[str, Any],
    item_id: str,
    repo_root: Optional[Path] = None,
) -> Path:
    """
    Write to Discovery_Accepted (The Workshop) instead of Finalized (The Vault).
    Tags placement_status and placement_note so we know why it's here.
    """
    from PX_Warehouse.warehouse_layout import get_discovery_accepted_dir
    root = Path(repo_root) if repo_root else _REPO_ROOT
    out_dir = get_discovery_accepted_dir(root)
    out_dir.mkdir(parents=True, exist_ok=True)
    safe_id = "".join(c if c.isalnum() or c in "-_" else "_" for c in str(item_id))[:64]
    fin = (finalized_dossier.get("finalization") or {}).copy()
    fin["placement_status"] = "DISCOVERY_ACCEPTED"
    fin["placement_note"] = "Failed Zeus (0.021) but passed Discovery (0.05)"
    fin["discovery_accepted"] = True
    fin["zeus_authorized"] = False
    out_dossier = {**finalized_dossier, "finalization": {**(finalized_dossier.get("finalization") or {}), **fin}}
    out_path = out_dir / f"{safe_id}.json"
    out_path.write_text(json.dumps(out_dossier, indent=2, default=str), encoding="utf-8")
    return out_path


def finalize_and_place(
    dossier: Dict[str, Any],
    item_id: str,
    is_novel: bool,
    repo_root: Optional[Path] = None,
) -> Optional[Path]:
    """
    Run finalization.
    - If Zeus pass -> Write to Finalized_Dossiers/<tier>
    - If Zeus fail but tox <= 5% -> Write to Finalized_Dossiers/Discovery_Accepted
    - If > 5% tox -> Drop (Garbage Collection)
    """
    finalized, tier, err = run_finalization(dossier, item_id, is_novel, repo_root)
    if err or finalized is None:
        return None
    zeus_verdict = (finalized.get("finalization") or {}).get("zeus_verdict") or {}
    # PATH A: The Gold Standard (Filing Ready)
    if zeus_verdict.get("authorized"):
        return write_finalized_dossier(finalized, item_id, tier, repo_root)
    # PATH B: The Discovery Funnel (Refinement Ready)
    if passes_discovery_bar(finalized):
        return write_discovery_accepted_dossier(finalized, item_id, repo_root)
    # PATH C: The Hard Floor (Garbage)
    return None


__all__ = [
    "run_finalization",
    "write_finalized_dossier",
    "write_discovery_accepted_dossier",
    "passes_discovery_bar",
    "finalize_and_place",
    "validate_finalized_dossier",
]


def validate_finalized_dossier(finalized: Dict[str, Any]) -> Tuple[bool, List[str]]:
    """Validate a finalized dossier against FINALIZATION_SPEC. Returns (ok, missing_or_errors)."""
    spec = _get_spec()
    ok, missing, _ = _validate_finalized_dossier(finalized, spec)
    return (ok, missing)
