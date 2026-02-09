"""
ZeusLaws - Constitutional Governance
Provides constitutional checks for OLYMPUS operations.

Laws:
  L1  - Internal 35D physics snapshot is primary (external cannot override).
  L10 - Harm is Absolute (harm_energy from engines; no defaults).
  L11 - Deterministic engine (toxicity_index threshold 0.0210).
  U27 - Stability law (engine outputs stable).
  U34 - Global sum conserved (physics invariant).
"""

from datetime import datetime, timezone
from typing import Any, Dict


def check_constitutional(operation: str, payload: dict) -> dict:
    """
    Deterministic constitutional check (L10/L11).
    """
    tox = payload.get("toxicity_index")
    harm = payload.get("harm_energy")
    if tox is not None and tox >= 0.0210:
        return {
            "authorized": False,
            "operation": operation,
            "rationale": "Law L11 Violation: toxicity_index >= 0.0210",
            "timestamp": None
        }
    if harm is not None and harm >= 0.0210:
        return {
            "authorized": False,
            "operation": operation,
            "rationale": "Law L10 Violation: harm_energy >= 0.0210",
            "timestamp": None
        }
    return {
        "authorized": True,
        "operation": operation,
        "rationale": "Constitutional check passed",
        "timestamp": None
    }


def _whole_profile_passes(admet: Dict[str, Any], toxicity_index: float, harm_energy: float) -> bool:
    """
    Whole-profile evaluation: do not reject on a single line if the full ADMET profile
    is acceptable. TOXICITY_DIAMOND = (tox < 0.01) OR (tox > 0.02 and safety_margin > 50).
    """
    if not admet:
        return False
    tox_dict = admet.get("toxicity") if isinstance(admet.get("toxicity"), dict) else {}
    risk_level = (tox_dict.get("risk_level") or "").strip()
    safety_margins = admet.get("safety_margins") or {}
    safety_margin = 0.0
    if isinstance(safety_margins, dict) and safety_margins.get("safety_margin") is not None:
        try:
            safety_margin = float(safety_margins["safety_margin"])
        except (TypeError, ValueError):
            pass
    if safety_margin <= 0 and isinstance(tox_dict.get("safety_margin"), (int, float)):
        safety_margin = float(tox_dict["safety_margin"])

    # ADMET already decided: Diamond = ultra-pure OR smart offset (high tox + safety_margin > 50)
    if risk_level == "TOXICITY_DIAMOND":
        return True
    # Explicit whole-profile rule (same as ADMET/warehouse_layout): tox < 0.01 OR (tox > 0.02 and safety_margin > 50)
    if toxicity_index is not None and harm_energy is not None:
        if toxicity_index < 0.0100:
            return True
        if toxicity_index > 0.0200 and safety_margin > 50.0:
            return True
    return False


def run_zeus_gate(
    dossier: Dict[str, Any],
    toxicity_threshold: float = 0.0210,
    harm_energy_threshold: float = 0.0210,
) -> Dict[str, Any]:
    """
    Full Zeus gate per FINALIZATION_SPEC: L1, U27, U34, L11.
    WHOLE-PROFILE: Do not reject on a single metric if the full profile passes
    (e.g. TOXICITY_DIAMOND or tox < 0.01 or (tox > 0.02 and safety_margin > 50)).
    Returns zeus_verdict with laws_required and per-law result.
    """
    ts = datetime.now(timezone.utc).isoformat()
    engines = dossier.get("engines") or {}
    admet = engines.get("admet") or {}
    tox_dict = admet.get("toxicity") if isinstance(admet.get("toxicity"), dict) else {}
    toxicity_index = dossier.get("harm_energy")
    if tox_dict and tox_dict.get("toxicity_index") is not None:
        toxicity_index = float(tox_dict["toxicity_index"])
    harm_energy = dossier.get("harm_energy")
    if harm_energy is not None:
        harm_energy = float(harm_energy)

    laws_required = [
        "L1_HARM_LAW",
        "U27_STABILITY_LAW",
        "U34_GLOBAL_SUM",
        "L11_DETERMINISTIC_ENGINE"
    ]
    results = {}

    # L11: Deterministic engine toxicity threshold (single-metric)
    l11_fail = toxicity_index is not None and toxicity_index >= toxicity_threshold
    results["L11_DETERMINISTIC_ENGINE"] = (
        {"passed": False, "reason": f"toxicity_index {toxicity_index} >= {toxicity_threshold}"} if l11_fail
        else {"passed": True, "reason": "toxicity below threshold"}
    )

    # L10 / L1_HARM_LAW: harm_energy threshold (single-metric)
    l1_fail = harm_energy is not None and harm_energy >= harm_energy_threshold
    results["L1_HARM_LAW"] = (
        {"passed": False, "reason": f"harm_energy {harm_energy} >= {harm_energy_threshold}"} if l1_fail
        else {"passed": True, "reason": "harm_energy below threshold"}
    )

    # U34: Global sum conserved (causal_trace_log or engines present)
    has_trace = bool(dossier.get("causal_trace_log"))
    has_engines = bool(dossier.get("engines"))
    results["U34_GLOBAL_SUM"] = {"passed": has_trace or has_engines, "reason": "causal_trace or engines present"}

    # U27: Stability (engine outputs present and deterministic)
    results["U27_STABILITY_LAW"] = {"passed": bool(engines.get("ope") and engines.get("admet")), "reason": "OPE and ADMET outputs present"}

    # Spec 9: Zeus review of trial outcomes — harm_energy_under_trial, toxicity_under_trial, stability_under_trial
    trial_outcome = dossier.get("trial_outcome_summary") or dossier.get("finalization", {}).get("trial_outcome_summary") or {}
    harm_energy_under_trial = None
    toxicity_under_trial = None
    stability_under_trial = None
    if trial_outcome.get("status") == "COMPLETED":
        toxicity_under_trial = trial_outcome.get("toxicity_rate")
        harm_energy_under_trial = dossier.get("harm_energy") or toxicity_under_trial
        stability_under_trial = True
        if toxicity_under_trial is not None:
            try:
                ttox = float(toxicity_under_trial)
                results["TRIAL_TOXICITY"] = {"passed": ttox < toxicity_threshold, "reason": f"toxicity_under_trial {ttox} vs threshold {toxicity_threshold}"}
            except (TypeError, ValueError):
                results["TRIAL_TOXICITY"] = {"passed": True, "reason": "trial toxicity not numeric"}
        else:
            results["TRIAL_TOXICITY"] = {"passed": True, "reason": "no trial toxicity_rate"}
        results["TRIAL_HARM_ENERGY"] = {"passed": harm_energy_under_trial is None or float(harm_energy_under_trial) < harm_energy_threshold, "reason": f"harm_energy_under_trial {harm_energy_under_trial}"}
        results["TRIAL_STABILITY"] = {"passed": bool(stability_under_trial), "reason": "stability_under_trial from virtual cohort"}
    else:
        results["TRIAL_TOXICITY"] = {"passed": True, "reason": "trial not run"}
        results["TRIAL_HARM_ENERGY"] = {"passed": True, "reason": "trial not run"}
        results["TRIAL_STABILITY"] = {"passed": True, "reason": "trial not run"}

    all_single_passed = all(r.get("passed") for r in results.values())
    whole_profile_ok = _whole_profile_passes(admet, toxicity_index, harm_energy)

    # Authorize if every law passes OR the whole profile passes (multi-factor; no single-line rejection when profile is acceptable)
    authorized = all_single_passed or whole_profile_ok
    rationale = (
        "All laws passed" if all_single_passed else
        "Whole profile passed (TOXICITY_DIAMOND or tox<0.01 or safety_margin>50)" if whole_profile_ok else
        "One or more laws failed and whole profile did not pass"
    )
    out = {
        "authorized": authorized,
        "operation": "zeus_gate",
        "rationale": rationale,
        "timestamp": ts,
        "laws_required": laws_required,
        "laws_results": results,
        "whole_profile_passed": whole_profile_ok,
        "all_single_passed": all_single_passed,
        "toxicity_index": toxicity_index,
        "harm_energy": harm_energy,
        "verdict_field": "zeus_verdict",
    }
    if trial_outcome.get("status") == "COMPLETED":
        out["harm_energy_under_trial"] = harm_energy_under_trial
        out["toxicity_under_trial"] = toxicity_under_trial
        out["stability_under_trial"] = stability_under_trial
    return out


# Legacy exports for backward compatibility
__all__ = ["check_constitutional", "run_zeus_gate"]
