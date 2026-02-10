"""
ADMET - Absorption, Distribution, Metabolism, Excretion, Toxicity
Predicts pharmacokinetic and safety properties using deterministic physical laws.

TAXONOMY: 4-TIER SYSTEM (Diamond, Gold, Silver, Failure).
  TOXICITY_DIAMOND — Ultra-pure (tox < 0.01) OR Smart Offset (high tox + safety_margin > 50)
  TOXICITY_GOLD    — Standard safe (tox < 0.0200)
  TOXICITY_SILVER  — Borderline / warning (0.0200 <= tox < 0.0210)
  TOXICITY_FAILURE — Hard stop (tox >= 0.0210)
"""

import time
from typing import Dict, Any
from PX_System.foundation.sign_off import create_sign_off


def run_admet(
    smiles: str,
    ope_analysis: Dict[str, Any],
    ome_result: Dict[str, Any] | None = None,
    ose_result: Dict[str, Any] | None = None,
) -> Dict[str, Any]:
    """
    Run ADMET analysis on a SMILES string using OPE molecular descriptors.
    When ome_result is provided, metabolism clearance and half-life are taken from OME
    (Operational Metabolism Engine) instead of being derived from OPE.
    When ose_result is provided, Selectivity Index (SI) is taken from OSE and included
    in the toxicity/risk level label (mandatory field for final risk level).
    Returns full structure for pipeline (absorption, distribution, metabolism, excretion, toxicity, etc.).
    """
    _t0 = time.monotonic()
    logp_provided = "logp" in ope_analysis
    logp = ope_analysis.get("logp", 2.0)
    mw = ope_analysis.get("molecular_weight", 350.0)
    hbd = ope_analysis.get("hbd", 2)
    hba = ope_analysis.get("hba", 4)
    tpsa = ope_analysis.get("tpsa", 70.0)

    violations = 0
    if mw > 500: violations += 1
    if logp > 5: violations += 1
    if hbd > 5: violations += 1
    if hba > 10: violations += 1
    base_abs = 90.0 - (violations * 20.0)
    if tpsa > 140: base_abs -= 15.0
    absorption = max(5.0, min(98.0, base_abs))

    ppb = min(99.9, 40.0 + (logp * 10.0))
    if ppb < 0: ppb = 5.0

    # Ingest OME metabolic stability and half-life when provided (bridge OME → ADMET)
    if ome_result and "clearance_L_per_h" in ome_result and "half_life_h" in ome_result:
        clearance = float(ome_result["clearance_L_per_h"])
        half_life_from_ome = float(ome_result["half_life_h"])
    else:
        clearance = ope_analysis.get("clearance_estimate_L_per_h", 10.0)
        half_life_from_ome = None
    hepatic_flow = 90.0
    extraction_ratio = min(0.95, clearance / hepatic_flow)

    renal_fraction = 0.1 + (tpsa / 200.0) - (logp / 20.0)
    renal_fraction = max(0.05, min(0.9, renal_fraction))

    herg_risk = 0.1 + (logp * 0.1)
    if tpsa < 40: herg_risk += 0.2
    herg_risk = min(0.9, herg_risk)

    hepatotox_risk = 0.1 + (logp * 0.05)
    if mw > 500: hepatotox_risk += 0.15
    hepatotox_risk = min(0.9, hepatotox_risk)

    ames_risk = 0.05
    if any(x in smiles for x in ["N=N", "N(=O)=O", "ClC(Cl)Cl"]):
        ames_risk += 0.4

    toxicity_index = (0.4 * herg_risk + 0.3 * hepatotox_risk + 0.3 * ames_risk)

    max_cmax = 100.0 / (toxicity_index + 0.1)
    max_auc = max_cmax * 24.0

    ec50 = ope_analysis.get("ec50", 1.0)
    if not isinstance(ec50, (int, float)) or ec50 <= 0:
        ec50 = 1.0
    safety_margin = (float(ec50) * 10.0) / max_cmax if max_cmax > 0 else 0.0

    # Selectivity Index (SI): mandatory field from OSE when provided; otherwise deterministic proxy
    if ose_result is not None and "selectivity_index" in ose_result:
        selectivity_index = float(ose_result["selectivity_index"])
    else:
        selectivity_index = max(1.0, 20.0 - (logp * 2.0))

    if toxicity_index < 0.0100 or (toxicity_index > 0.0200 and safety_margin > 50.0):
        risk_label = "TOXICITY_DIAMOND"
    elif toxicity_index < 0.0200:
        risk_label = "TOXICITY_GOLD"
    elif toxicity_index < 0.0210:
        risk_label = "TOXICITY_SILVER"
    else:
        risk_label = "TOXICITY_FAILURE"
    risk_level_display = f"{risk_label} (SI={selectivity_index:.2f})"

    half_life_h = (
        half_life_from_ome
        if half_life_from_ome is not None
        else (0.693 * ope_analysis.get("vd_estimate_L", 70.0) / max(0.1, clearance))
    )

    result = {
        "absorption": {
            "oral_bioavailability_percent": float(f"{absorption:.8f}"),
            "intestinal_permeability": "high" if absorption > 70 else "moderate",
            "predicted_bioavailability": float(f"{absorption/100.0:.8f}"),
            "bcs_class": "I" if absorption > 80 and logp < 3 else "II" if logp >= 3 else "III",
            "note": "Calculated from Rule of Five and TPSA"
        },
        "distribution": {
            "plasma_protein_binding_percent": float(f"{ppb:.8f}"),
            "vd_L": ope_analysis.get("vd_estimate_L", 70.0),
            "blood_brain_barrier": "permeable" if tpsa < 90 and logp > 1 else "limited",
            "predicted_vd_L_per_kg": float(f"{ope_analysis.get('vd_estimate_L', 70.0)/70.0:.8f}"),
            "bbb_penetration": "HIGH" if tpsa < 60 else "LOW"
        },
        "metabolism": {
            "hepatic_extraction_ratio": float(f"{extraction_ratio:.8f}"),
            "clearance_L_per_h": float(f"{clearance:.8f}"),
            "half_life_h": float(f"{half_life_h:.8f}"),
            "source": "OME" if ome_result and "clearance_L_per_h" in ome_result else "OPE"
        },
        "excretion": {
            "renal_clearance_fraction": float(f"{renal_fraction:.8f}"),
            "hepatic_clearance_fraction": float(f"{1.0 - renal_fraction:.8f}")
        },
        "toxicity": {
            "toxicity_index": float(f"{toxicity_index:.8f}"),
            "herg_risk": float(f"{herg_risk:.8f}"),
            "hepatotoxicity_risk": float(f"{hepatotox_risk:.8f}"),
            "ames_mutagenicity": float(f"{ames_risk:.8f}"),
            "risk_level": risk_label,
            "risk_level_display": risk_level_display,
            "selectivity_index": float(f"{selectivity_index:.2f}"),
            "safety_margin": float(f"{safety_margin:.8f}")
        },
        "toxicity_flags": {
            "hepatotoxicity_risk": (
                "UNKNOWN" if not logp_provided else
                "HIGH" if logp > 4.5 else
                "MEDIUM" if logp >= 3.5 else
                "LOW"
            ),
            "cardiotoxicity_risk": "HIGH" if herg_risk > 0.5 else "LOW"
        },
        "constitutional": {
            "status": "VERIFIED",
            "engine": "OPE_ADMET_V3_DETERMINISTIC",
            "law_basis": ["L11", "U27"]
        },
        "safety_margins": {
            "max_tolerated_cmax_mg_per_L": float(f"{max_cmax:.8f}"),
            "max_tolerated_auc_mg_h_per_L": float(f"{max_auc:.8f}"),
            "safety_margin": float(f"{safety_margin:.8f}")
        },
        "toxicity_threshold": float(f"{max_cmax:.8f}"),
        "note": "ADMET engine using deterministic physical laws (v3.0-CORE-DETERMINISTIC)",
        "version": "3.0-CORE-DETERMINISTIC"
    }
    _elapsed_ms = int((time.monotonic() - _t0) * 1000)
    result["sign_off"] = create_sign_off(
        engine_id="OPE_ADMET_V3_DETERMINISTIC",
        version="3.0-CORE-DETERMINISTIC",
        inputs={"smiles": smiles},
        outputs=result,
        laws_checked=["L11"],
        laws_results={"L11": toxicity_index < 0.0210},
        execution_time_ms=_elapsed_ms,
    )
    return result
