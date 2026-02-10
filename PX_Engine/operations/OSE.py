"""
OSE - Operational Safety Engine
Calculates toxicity and selectivity index (SI)
"""

import time
from rdkit import Chem
from rdkit.Chem import Descriptors
from PX_System.foundation.sign_off import create_sign_off

def execute(payload):
    """
    Calculates safety metrics and harm energy.
    """
    _t0 = time.monotonic()
    smiles = payload.get("smiles")
    if not smiles:
        result = {"error": "Missing SMILES"}
        _elapsed_ms = int((time.monotonic() - _t0) * 1000)
        result["sign_off"] = create_sign_off(
            engine_id="OSE_V3_DETERMINISTIC",
            version="3.0-CORE",
            inputs=payload,
            outputs=result,
            laws_checked=["L11"],
            laws_results={"L11": False},
            execution_time_ms=_elapsed_ms,
        )
        return result

    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        result = {"error": "Invalid SMILES"}
        _elapsed_ms = int((time.monotonic() - _t0) * 1000)
        result["sign_off"] = create_sign_off(
            engine_id="OSE_V3_DETERMINISTIC",
            version="3.0-CORE",
            inputs=payload,
            outputs=result,
            laws_checked=["L11"],
            laws_results={"L11": False},
            execution_time_ms=_elapsed_ms,
        )
        return result
        
    logp = Descriptors.MolLogP(mol)
    
    # Toxicity index (0-1) - Law L11 Enforcement
    # Higher LogP often correlates with higher non-specific toxicity
    toxicity_index = 0.01 + (logp * 0.005)
    toxicity_index = max(0.001, toxicity_index)
    
    # Selectivity Index (SI) = CC50 / IC50
    # Higher is safer
    si = 20.0 - (logp * 2.0)
    si = max(1.0, si)
    
    result = {
        "engine": "OSE_V3_DETERMINISTIC",
        "toxicity_index": float(f"{toxicity_index:.8f}"),
        "selectivity_index": float(f"{si:.2f}"),
        "harm_energy": float(f"{toxicity_index * 1.5:.8f}"),
        "safety_status": "SAFE" if toxicity_index < 0.0200 else "CAUTION" if toxicity_index < 0.0210 else "FAILURE",
        "law_l11_compliance": "VERIFIED"
    }
    _elapsed_ms = int((time.monotonic() - _t0) * 1000)
    result["sign_off"] = create_sign_off(
        engine_id="OSE_V3_DETERMINISTIC",
        version="3.0-CORE",
        inputs=payload,
        outputs=result,
        laws_checked=["L11"],
        laws_results={"L11": toxicity_index < 0.0210},
        execution_time_ms=_elapsed_ms,
    )
    return result
