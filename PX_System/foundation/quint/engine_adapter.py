"""
QUINT Engine Adapter — Thin wrappers for ALL engines in the Predator X system.

Each q_run_* function accepts EITHER a QFrame or a plain dict (backward compatible),
calls the existing engine function unchanged, and wraps the result in a QFrame
(QType.QRESULT) with lineage tracking.

Design:
  - Backward compatible: plain dicts work identically to before
  - QFrame-aware: extracts payload via converter.emit() when given a QFrame
  - Preserves sign_off blocks inside the QFrame payload
  - Python stdlib only for the adapter itself (constitutional module)

PRV pipeline engines (1-12):
  OPE -> OBE -> OCE -> OLE -> OME -> OSE -> ADMET -> PKPD ->
  DoseOptimizer_v2 -> VirtualEfficacyAnalytics -> GradingEngine -> ZeusLaws

Auxiliary engines (13-20):
  VectorCore, SimulationEngine, TrialEngine, CheckConstitutional,
  EvidencePackage, WrapTrialSimulation, Metabolism, BlockOrchestrator
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


# ===========================================================================
# AUXILIARY ENGINES (13-20) — Physics, Simulation, Trial, Governance, etc.
# ===========================================================================


# ---------------------------------------------------------------------------
# 13. VectorCore — VectorCore().execute(p_array, physical_descriptors)
# ---------------------------------------------------------------------------

def q_run_vector_core(input_data, physical_descriptors=None):
    """
    QUINT-wrapped VectorCore physics engine.

    Args:
        input_data: QFrame or dict containing {"p_vector": [...]}.
                    Also accepts a plain list as the p_vector directly.
        physical_descriptors: Optional dict of physical descriptors.
    """
    from PX_Engine.Vector_Core import VectorCore

    if isinstance(input_data, list):
        p_vector = input_data
    else:
        payload = _unwrap(input_data)
        p_vector = payload.get("p_vector", [0.1, 0.0, 35.0, 1.0])
        if physical_descriptors is None:
            physical_descriptors = payload.get("physical_descriptors")

    core = VectorCore()
    result = core.execute(p_vector, physical_descriptors=physical_descriptors)
    return _wrap_result(result, "VECTOR_CORE_V1", "VectorCore", "VectorCore")


# ---------------------------------------------------------------------------
# 14. SimulationEngine — SimulationEngine().materialize_candidate(wl_id, coh)
# ---------------------------------------------------------------------------

def q_run_simulation(input_data):
    """
    QUINT-wrapped SimulationEngine materialization.

    Args:
        input_data: QFrame or dict containing {"task_id": ..., "coherence": ...}.
    """
    from PX_Laboratory.Simulation_Engine import SimulationEngine

    payload = _unwrap(input_data)
    task_id = payload.get("task_id", payload.get("worldline_id", "SIM-UNK"))
    coherence = payload.get("coherence", payload.get("amplitude", 0.0))
    engine = SimulationEngine()
    result = engine.materialize_candidate(task_id, coherence)
    if not isinstance(result, dict):
        result = {"materialization": result, "task_id": task_id}
    return _wrap_result(result, "SIMULATION_ENGINE_V1", f"SIM-{str(task_id)[:16]}", "SimulationEngine")


# ---------------------------------------------------------------------------
# 14b. PK Simulation — SimulationEngine.simulate_one_compartment(...)
# ---------------------------------------------------------------------------

def q_run_pk_simulation(dose_mg, duration_h, dosing_interval_h, patient, admet, time_step_h=0.5):
    """
    QUINT-wrapped PK one-compartment simulation.

    Args:
        dose_mg:            Dose in mg.
        duration_h:         Duration in hours.
        dosing_interval_h:  Dosing interval in hours.
        patient:            Dict with patient parameters (e.g. weight_kg).
        admet:              QFrame or dict with ADMET output.
        time_step_h:        Time step in hours for simulation.
    """
    from PX_Laboratory.Simulation_Engine import SimulationEngine

    admet_dict = _unwrap(admet)
    sim = SimulationEngine(time_step_h=time_step_h)
    result = sim.simulate_one_compartment(
        dose_mg=dose_mg,
        duration_h=duration_h,
        dosing_interval_h=dosing_interval_h,
        patient=patient,
        admet=admet_dict,
    )
    return _wrap_result(result, "PK_SIMULATION_V1", "PKSim", "PKSimulation")


# ---------------------------------------------------------------------------
# 15. TrialEngine — TrialEngine(time_step_h).run_trial(protocol, admet, ...)
# ---------------------------------------------------------------------------

def q_run_trial(protocol_or_frame, admet=None, pd_params=None, variability=None, time_step_h=1.0):
    """
    QUINT-wrapped TrialEngine clinical trial simulation.

    Args:
        protocol_or_frame: QFrame or dict with trial protocol.
        admet:       QFrame or dict with ADMET output.
        pd_params:   QFrame or dict with PD parameters (optional).
        variability: QFrame or dict with IIV parameters (optional).
        time_step_h: Time step in hours for PK simulation.
    """
    from PX_Engine.operations.TrialEngine import TrialEngine

    protocol = _unwrap(protocol_or_frame)
    admet_dict = _unwrap(admet) if admet is not None else None
    pd_dict = _unwrap(pd_params) if pd_params is not None else None
    var_dict = _unwrap(variability) if variability is not None else None

    engine = TrialEngine(time_step_h=time_step_h)
    result = engine.run_trial(protocol=protocol, admet=admet_dict, pd_params=pd_dict, variability=var_dict)
    return _wrap_result(result, "TRIAL_ENGINE_V1", "Trial", "TrialEngine")


# ---------------------------------------------------------------------------
# 16. CheckConstitutional — check_constitutional(operation, payload)
# ---------------------------------------------------------------------------

def q_run_check_constitutional(operation, data_or_frame):
    """
    QUINT-wrapped ZeusLaws constitutional check (L10/L11).

    Args:
        operation: String identifying the operation being checked.
        data_or_frame: QFrame or dict with toxicity/harm data.
    """
    from PX_System.foundation.ZeusLaws import check_constitutional

    data = _unwrap(data_or_frame)
    result = check_constitutional(operation, data)
    return _wrap_result(result, "CONSTITUTIONAL_CHECK_V1", f"CONST-{operation[:16]}", "CheckConstitutional")


# ---------------------------------------------------------------------------
# 17. EvidencePackage — generate_dossier(candidate, engine_results, ...)
# ---------------------------------------------------------------------------

def q_run_evidence_dossier(candidate_or_frame, engine_results=None, zeus_verdict=None, context=None):
    """
    QUINT-wrapped Evidence_Package dossier generation.

    Args:
        candidate_or_frame: QFrame or dict with candidate data.
        engine_results: QFrame or dict with engine outputs.
        zeus_verdict:   QFrame or dict with Zeus gate result.
        context:        QFrame or dict with additional context.
    """
    from PX_System.foundation.Evidence_Package import generate_dossier

    candidate = _unwrap(candidate_or_frame)
    eng = _unwrap(engine_results) if engine_results is not None else None
    zeus = _unwrap(zeus_verdict) if zeus_verdict is not None else None
    ctx = _unwrap(context) if context is not None else None

    result = generate_dossier(candidate, eng, zeus_verdict=zeus, context=ctx)
    return _wrap_result(result, "EVIDENCE_PACKAGE_V1", "Evidence", "EvidencePackage")


# ---------------------------------------------------------------------------
# 18. WrapTrialSimulation — wrap_trial_simulation(protocol, trial, ope, admet)
# ---------------------------------------------------------------------------

def q_run_wrap_trial(protocol_or_frame, trial_result=None, ope_result=None, admet_result=None, output_dir=None):
    """
    QUINT-wrapped wrap_trial_simulation for Evidence_Package.

    Args:
        protocol_or_frame: QFrame or dict with trial protocol.
        trial_result: QFrame or dict with TrialEngine output.
        ope_result:   QFrame or dict with OPE output.
        admet_result: QFrame or dict with ADMET output.
        output_dir:   Output directory for trial simulation files.
    """
    from PX_System.foundation.Evidence_Package import wrap_trial_simulation

    protocol = _unwrap(protocol_or_frame)
    trial = _unwrap(trial_result) if trial_result is not None else None
    ope = _unwrap(ope_result) if ope_result is not None else None
    admet = _unwrap(admet_result) if admet_result is not None else None

    result = wrap_trial_simulation(protocol, trial, ope, admet, output_dir=output_dir)
    return _wrap_result(result, "WRAP_TRIAL_V1", "WrapTrial", "WrapTrialSimulation")


# ---------------------------------------------------------------------------
# 19. Metabolism — Metabolism().pulse(task_id, p_vec, csa_s)
# ---------------------------------------------------------------------------

def q_run_metabolism_pulse(input_data):
    """
    QUINT-wrapped Metabolism pulse (heartbeat cycle).

    Args:
        input_data: QFrame or dict containing {"task_id": ..., "p_vec": [...], "csa_s": [...]}.
                    Also accepts a plain string as task_id.
    """
    from PX_Engine.Metabolism import Metabolism

    if isinstance(input_data, str):
        payload = {"task_id": input_data}
    else:
        payload = _unwrap(input_data)

    task_id = payload.get("task_id", payload.get("event", "pulse"))
    p_vec = payload.get("p_vec")
    csa_s = payload.get("csa_s")

    m = Metabolism()
    result = m.pulse(task_id, p_vec, csa_s)
    if not isinstance(result, dict):
        if isinstance(result, tuple) and len(result) == 2:
            result = {"cycle_age": result[0], "status": result[1], "task_id": task_id}
        else:
            result = {"pulse_result": result, "task_id": task_id}
    return _wrap_result(result, "METABOLISM_V1", f"METAB-{str(task_id)[:16]}", "Metabolism")


# ---------------------------------------------------------------------------
# 20. BlockOrchestrator — execute_living_pulse(task_id, ctx, p_vec, csa_s)
# ---------------------------------------------------------------------------

def q_run_block_pulse(input_data):
    """
    QUINT-wrapped Block_Orchestrator living pulse.

    Args:
        input_data: QFrame or dict containing
                    {"task_id": ..., "ctx": {...}, "p_vec": [...], "csa_s": [...]}.
    """
    from PX_Engine.Block_Orchestrator import execute_living_pulse

    payload = _unwrap(input_data)
    task_id = payload.get("task_id", "UNK")
    ctx = payload.get("ctx", payload.get("source_ctx", {}))
    p_vec = payload.get("p_vec", payload.get("p_vector", [0.1, 0.0, 35.0, 1.0]))
    csa_s = payload.get("csa_s", payload.get("csa_scores", []))

    result = execute_living_pulse(task_id, ctx, p_vec, csa_s)
    if not isinstance(result, dict):
        result = {"pulse_result": result, "task_id": task_id}
    return _wrap_result(result, "BLOCK_ORCHESTRATOR_V1", f"BLOCK-{str(task_id)[:16]}", "BlockOrchestrator")


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

__all__ = [
    "_unwrap",
    # PRV pipeline (1-12)
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
    # Auxiliary (13-20)
    "q_run_vector_core",
    "q_run_simulation",
    "q_run_pk_simulation",
    "q_run_trial",
    "q_run_check_constitutional",
    "q_run_evidence_dossier",
    "q_run_wrap_trial",
    "q_run_metabolism_pulse",
    "q_run_block_pulse",
]
