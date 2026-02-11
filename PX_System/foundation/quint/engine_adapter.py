"""
QUINT Engine Adapter — Thin wrappers for all 12 PRV pipeline engines.

Each q_run_* function accepts EITHER a QFrame or a plain dict (backward compatible),
calls the existing engine function unchanged, and wraps the result in a QFrame
(QType.QRESULT) with lineage tracking.

Design:
  - Backward compatible: plain dicts work identically to before
  - QFrame-aware: extracts payload via converter.emit() when given a QFrame
  - Preserves sign_off blocks inside the QFrame payload
  - Python stdlib only for the adapter itself (constitutional module)

The 12 engines in PRV pipeline order:
  OPE -> OBE -> OCE -> OLE -> OME -> OSE -> ADMET -> PKPD ->
  DoseOptimizer_v2 -> VirtualEfficacyAnalytics -> GradingEngine -> ZeusLaws
"""

from __future__ import annotations

from typing import Any, Dict, Optional

from PX_System.foundation.quint.kernel import QFrame, QType, create_qframe
from PX_System.foundation.quint.converter import emit, ingest


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def _unwrap(data):
    """Extract plain dict from QFrame or pass through dict unchanged."""
    if isinstance(data, QFrame):
        return emit(data, target_label="engine_unwrap")
    if isinstance(data, dict):
        return data
    raise TypeError(f"Expected QFrame or dict, got {type(data).__name__}")


def _wrap_result(result, engine_id: str, qid_prefix: str, lineage_tag: str) -> QFrame:
    """Wrap an engine result dict into a sealed QRESULT QFrame."""
    result_payload = dict(result) if isinstance(result, dict) else {"raw": result}
    result_payload["engine_id"] = engine_id
    result_payload.setdefault("status", "PASSED")

    frame = create_qframe(
        qtype=QType.QRESULT,
        qid=f"QRESULT-{qid_prefix}",
        payload=result_payload,
    )
    frame.lineage.append(f"engine:{lineage_tag}")
    frame.seal()
    return frame


# ---------------------------------------------------------------------------
# 1. OPE — run_ope(smiles) -> dict
# ---------------------------------------------------------------------------

def q_run_ope(input_data):
    """QUINT-wrapped OPE engine."""
    from PX_Engine.operations.OPE import run_ope

    payload = _unwrap(input_data)
    smiles = payload.get("smiles", "CCO")
    result = run_ope(smiles)
    return _wrap_result(result, "OPE_V3_DETERMINISTIC", f"OPE-{smiles[:16]}", "OPE")


# ---------------------------------------------------------------------------
# 2. OBE — execute({"smiles": smiles}) -> dict
# ---------------------------------------------------------------------------

def q_run_obe(input_data):
    """QUINT-wrapped OBE engine."""
    from PX_Engine.operations.OBE import execute as obe_execute

    payload = _unwrap(input_data)
    result = obe_execute(payload)
    smiles = payload.get("smiles", "UNK")
    return _wrap_result(result, "OBE_V3_DETERMINISTIC", f"OBE-{smiles[:16]}", "OBE")


# ---------------------------------------------------------------------------
# 3. OCE — execute({"p_vector": [...], "csa_scores": [...]}) -> dict
# ---------------------------------------------------------------------------

def q_run_oce(input_data):
    """QUINT-wrapped OCE engine."""
    from PX_Engine.operations.OCE import execute as oce_execute

    payload = _unwrap(input_data)
    result = oce_execute(payload)
    return _wrap_result(result, "OCE_V3_DETERMINISTIC", "OCE", "OCE")


# ---------------------------------------------------------------------------
# 4. OLE — execute({"compound_id": ..., "indication": ..., ...}) -> dict
# ---------------------------------------------------------------------------

def q_run_ole(input_data):
    """QUINT-wrapped OLE engine."""
    from PX_Engine.operations.OLE import execute as ole_execute

    payload = _unwrap(input_data)
    result = ole_execute(payload)
    cid = payload.get("compound_id", "UNK")
    return _wrap_result(result, "OLE_V3_DETERMINISTIC", f"OLE-{cid[:16]}", "OLE")


# ---------------------------------------------------------------------------
# 5. OME — execute({"smiles": smiles}) -> dict
# ---------------------------------------------------------------------------

def q_run_ome(input_data):
    """QUINT-wrapped OME engine."""
    from PX_Engine.operations.OME import execute as ome_execute

    payload = _unwrap(input_data)
    result = ome_execute(payload)
    smiles = payload.get("smiles", "UNK")
    return _wrap_result(result, "OME_V3_DETERMINISTIC", f"OME-{smiles[:16]}", "OME")


# ---------------------------------------------------------------------------
# 6. OSE — execute({"smiles": smiles}) -> dict
# ---------------------------------------------------------------------------

def q_run_ose(input_data):
    """QUINT-wrapped OSE engine."""
    from PX_Engine.operations.OSE import execute as ose_execute

    payload = _unwrap(input_data)
    result = ose_execute(payload)
    smiles = payload.get("smiles", "UNK")
    return _wrap_result(result, "OSE_V3_DETERMINISTIC", f"OSE-{smiles[:16]}", "OSE")


# ---------------------------------------------------------------------------
# 7. ADMET — run_admet(smiles, ope_result, ome_result=None, ose_result=None)
# ---------------------------------------------------------------------------

def q_run_admet(input_data, ope_result=None, ome_result=None, ose_result=None):
    """
    QUINT-wrapped ADMET engine.

    Args:
        input_data: QFrame or dict containing at least {"smiles": ...}.
        ope_result: QFrame or dict with OPE engine output (required by ADMET).
        ome_result: QFrame or dict with OME engine output (optional).
        ose_result: QFrame or dict with OSE engine output (optional).
    """
    from PX_Engine.operations.ADMET import run_admet

    payload = _unwrap(input_data)
    smiles = payload.get("smiles", "CCO")

    ope_dict = _unwrap(ope_result) if ope_result is not None else payload
    ome_dict = _unwrap(ome_result) if ome_result is not None else None
    ose_dict = _unwrap(ose_result) if ose_result is not None else None

    result = run_admet(smiles, ope_dict, ome_result=ome_dict, ose_result=ose_dict)
    return _wrap_result(result, "ADMET_V3_DETERMINISTIC", f"ADMET-{smiles[:16]}", "ADMET")


# ---------------------------------------------------------------------------
# 8. PKPD — link_pk_to_pd(pk_profile, pd_params) -> dict
# ---------------------------------------------------------------------------

def q_run_pkpd(pk_profile, pd_params):
    """
    QUINT-wrapped PKPD engine.

    Args:
        pk_profile: QFrame or dict with PK concentration-time profile.
        pd_params:  QFrame or dict with PD model parameters.
    """
    from PX_Engine.operations.PKPD import link_pk_to_pd

    pk_dict = _unwrap(pk_profile)
    pd_dict = _unwrap(pd_params)

    result = link_pk_to_pd(pk_dict, pd_dict)
    return _wrap_result(result, "PKPD_V2_DETERMINISTIC", "PKPD", "PKPD")


# ---------------------------------------------------------------------------
# 9. DoseOptimizer_v2 — optimize_dose(smiles, admet, protocol_template, ...)
# ---------------------------------------------------------------------------

def q_run_dose_optimizer(
    smiles_or_frame,
    admet,
    protocol_template,
    pd_params=None,
    n_eval_patients: int = 10,
):
    """
    QUINT-wrapped DoseOptimizer_v2 engine.

    Args:
        smiles_or_frame: QFrame or dict containing {"smiles": ...}, or plain SMILES str.
        admet:            QFrame or dict with ADMET output.
        protocol_template: QFrame or dict with trial protocol template.
        pd_params:        QFrame or dict with PD parameters (optional).
        n_eval_patients:  Number of evaluation patients per mini-trial.
    """
    from PX_Engine.operations.DoseOptimizer_v2 import optimize_dose

    if isinstance(smiles_or_frame, str) and not isinstance(smiles_or_frame, dict):
        # Caller passed a plain SMILES string directly
        smiles = smiles_or_frame
    else:
        smiles_payload = _unwrap(smiles_or_frame)
        smiles = smiles_payload.get("smiles", "CCO")

    admet_dict = _unwrap(admet)
    protocol_dict = _unwrap(protocol_template)
    pd_dict = _unwrap(pd_params) if pd_params is not None else None

    result = optimize_dose(
        smiles=smiles,
        admet=admet_dict,
        protocol_template=protocol_dict,
        pd_params=pd_dict,
        n_eval_patients=n_eval_patients,
    )
    return _wrap_result(result, "DOSE_OPTIMIZER_V2", f"DOSE-{smiles[:16]}", "DoseOptimizer_v2")


# ---------------------------------------------------------------------------
# 10. VirtualEfficacyAnalytics — compute_pta, virtual_responder_rate,
#                                  effect_variability_risk
# ---------------------------------------------------------------------------

def q_run_virtual_efficacy(
    trial_result_or_frame,
    metric: str = "auc_mg_h_per_L",
    target_threshold: float = 250.0,
    response_metric: str = "max_effect",
    responder_threshold: float = 0.5,
    effect_metric: str = "max_effect",
):
    """
    QUINT-wrapped VirtualEfficacyAnalytics — runs all three analytics
    (PTA, responder rate, variability risk) and returns a combined QFrame.

    Args:
        trial_result_or_frame: QFrame or dict with TrialEngine output.
        metric:              PK/PD metric for PTA calculation.
        target_threshold:    PTA threshold.
        response_metric:     Metric for responder rate.
        responder_threshold: Responder rate threshold.
        effect_metric:       Metric for variability risk.
    """
    from PX_Engine.operations.VirtualEfficacyAnalytics import (
        compute_pta,
        virtual_responder_rate,
        effect_variability_risk,
    )

    trial_dict = _unwrap(trial_result_or_frame)

    pta_result = compute_pta(trial_dict, metric, target_threshold)
    responder_result = virtual_responder_rate(trial_dict, response_metric, responder_threshold)
    variability_result = effect_variability_risk(trial_dict, effect_metric)

    combined = {
        "pta": pta_result,
        "responder_rate": responder_result,
        "variability": variability_result,
        "engine_id": "VIRTUAL_EFFICACY_V2",
        "status": "PASSED",
    }

    frame = create_qframe(
        qtype=QType.QRESULT,
        qid="QRESULT-VirtualEfficacy",
        payload=combined,
    )
    frame.lineage.append("engine:VirtualEfficacyAnalytics")
    frame.seal()
    return frame


# ---------------------------------------------------------------------------
# 11. GradingEngine — GradingEngine(verbose=False).grade_dossier(dossier)
# ---------------------------------------------------------------------------

def q_run_grading(dossier_or_frame):
    """
    QUINT-wrapped GradingEngine.

    Args:
        dossier_or_frame: QFrame or dict containing a complete dossier.
    """
    from PX_Engine.operations.GradingEngine import GradingEngine

    dossier = _unwrap(dossier_or_frame)
    ge = GradingEngine(verbose=False)
    result = ge.grade_dossier(dossier)
    return _wrap_result(result, "GRADING_ENGINE_V2", "Grading", "GradingEngine")


# ---------------------------------------------------------------------------
# 12. ZeusLaws — run_zeus_gate(dossier) -> dict
# ---------------------------------------------------------------------------

def q_run_zeus(dossier_or_frame):
    """
    QUINT-wrapped ZeusLaws governance gate.

    Args:
        dossier_or_frame: QFrame or dict containing a complete dossier.
    """
    from PX_System.foundation.ZeusLaws import run_zeus_gate

    dossier = _unwrap(dossier_or_frame)
    result = run_zeus_gate(dossier)
    return _wrap_result(result, "ZEUS_GATE_V1", "Zeus", "ZeusLaws")


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

__all__ = [
    "_unwrap",
    "q_run_ope",
    "q_run_obe",
    "q_run_oce",
    "q_run_ole",
    "q_run_ome",
    "q_run_ose",
    "q_run_admet",
    "q_run_pkpd",
    "q_run_dose_optimizer",
    "q_run_virtual_efficacy",
    "q_run_grading",
    "q_run_zeus",
]
