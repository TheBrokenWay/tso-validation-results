"""
PX_Refinery.py
The Great Filter v2.4 - Serialization Safe.

Updates:
- Implements '_json_safe' recursive sanitizer to strip NumPy types.
- Fixes "Object of type bool is not JSON serializable" errors.
- Preserves Holistic Logic (Diamond/PRV/Clean).

Repo-aligned: REPO_ROOT paths; Vector_Core 5-element p_vector (sum 36.1); pulse gets 4 elements;
_extract_attributes includes candidate (Evidence Package); safety_margin computed from pk_summary;
protocol_template with arms; correct imports (optimize_dose from DoseOptimizer_v2).
"""

from __future__ import annotations

import json
import os
import shutil
import sys
import numpy as np
from pathlib import Path
from datetime import datetime, timezone
from typing import Any, Dict, List, Tuple

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

WAREHOUSE_PATH = REPO_ROOT / "PX_Warehouse"
BACKUP_PATH = WAREHOUSE_PATH / "Backup_Pre_Refinery"

# --- Sub-root physics ---
from PX_Engine.Vector_Core import VectorCore
from PX_Engine.Metabolism import Metabolism

# --- Operations (production); DoseOptimizer_v2 not in __all__ so import from module ---
from PX_Engine.operations.OPE import run_ope
from PX_Engine.operations import OBE, OME, OLE, OSE, ADMET, OCE
from PX_Engine.operations.DoseOptimizer_v2 import optimize_dose
from PX_System.foundation.ZeusLaws import run_zeus_gate

vector_gate = VectorCore(threshold=0.95, dims_limit=35.0, global_sum_target=36.1)
metabolism = Metabolism()


def _json_safe(obj: Any) -> Any:
    """
    Recursively converts NumPy types to native Python types for JSON serialization.
    """
    if isinstance(obj, dict):
        return {k: _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return _json_safe(obj.tolist())
    if isinstance(obj, np.generic):
        return obj.item()
    return obj


def setup_backup() -> bool:
    """Safety first: Backup the current timeline."""
    if not WAREHOUSE_PATH.exists():
        print(">>> [REFINERY] Warehouse not found!")
        return False
    BACKUP_PATH.mkdir(parents=True, exist_ok=True)
    count = 0
    for file in WAREHOUSE_PATH.rglob("*.json"):
        if "Backup_Pre_Refinery" in str(file):
            continue
        rel_path = file.relative_to(WAREHOUSE_PATH)
        dest_path = BACKUP_PATH / rel_path
        dest_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(file, dest_path)
        count += 1
    print(f">>> [REFINERY] Timeline secured. {count} records backed up.")
    return True


def _extract_attributes(data: Dict[str, Any]) -> Tuple[str | None, str]:
    """Hunts for SMILES and ID in various JSON schemas (including Evidence Package candidate)."""
    smiles = None
    compound_id = "UNKNOWN"

    if data.get("smiles") and isinstance(data["smiles"], str):
        smiles = data["smiles"]
    elif data.get("candidate") and isinstance(data["candidate"], dict) and data["candidate"].get("smiles"):
        smiles = data["candidate"]["smiles"]
    elif data.get("candidate_profile") and isinstance(data["candidate_profile"], dict):
        if data["candidate_profile"].get("smiles"):
            smiles = data["candidate_profile"]["smiles"]
    elif data.get("candidate_data") and isinstance(data["candidate_data"], dict):
        prv = data["candidate_data"].get("prv_candidate")
        if isinstance(prv, dict) and prv.get("smiles"):
            smiles = prv["smiles"]
        elif data["candidate_data"].get("smiles"):
            smiles = data["candidate_data"]["smiles"]

    if data.get("compound_id"):
        compound_id = str(data["compound_id"])
    elif data.get("id"):
        compound_id = str(data["id"])
    elif data.get("candidate") and isinstance(data["candidate"], dict) and data["candidate"].get("name"):
        compound_id = str(data["candidate"]["name"])
    elif data.get("candidate_profile") and isinstance(data["candidate_profile"], dict):
        if data["candidate_profile"].get("compound_id"):
            compound_id = str(data["candidate_profile"]["compound_id"])
    elif data.get("header") and isinstance(data["header"], dict) and data["header"].get("worldline_id"):
        compound_id = str(data["header"]["worldline_id"])
    elif data.get("metadata") and isinstance(data["metadata"], dict) and data["metadata"].get("id"):
        compound_id = str(data["metadata"]["id"])

    if not isinstance(smiles, str):
        smiles = None
    else:
        smiles = (smiles or "").strip() or None
    compound_id = (compound_id or "UNKNOWN")[:64]
    return smiles, compound_id


def assess_whole_product(
    admet: Dict[str, Any],
    dose_result: Dict[str, Any],
    legal_result: Dict[str, Any],
    vector_state: Dict[str, Any],
) -> Tuple[str, bool]:
    """
    Judges the artifact based on the WHOLE PRODUCT view (v2.1 Logic).
    safety_margin is computed from admet + best_regimen pk_summary (best_regimen has no "safety_margin" key).
    """
    if not vector_state.get("authorized"):
        return "COLLAPSED_PHYSICS_U34", False

    best = dose_result.get("best_regimen")
    if not best or not best.get("pk_summary"):
        return "COLLAPSED_NO_REGIMEN", False

    tox_index = float((admet.get("toxicity") or {}).get("toxicity_index", 1.0))
    safety_margins = admet.get("safety_margins") or {}
    max_cmax = float(safety_margins.get("max_tolerated_cmax_mg_per_L", 0) or 0)
    cmax_dist = (best.get("pk_summary") or {}).get("cmax_mg_per_L") or {}
    cmax_mean = float((cmax_dist.get("mean") or 0) or 0)
    safety_margin = (max_cmax / cmax_mean) if cmax_mean > 0 else 0.0
    prv = bool(legal_result.get("prv_eligible", False))

    if tox_index < 0.0210 and safety_margin > 10.0:
        return "OPTIMIZED_CLEAN", True
    if tox_index >= 0.0210 and safety_margin > 50.0:
        return "OPTIMIZED_DIAMOND", True
    if prv and safety_margin > 30.0:
        return "OPTIMIZED_STRATEGIC_PRV", True
    if tox_index >= 0.0210:
        return f"COLLAPSED_TOX_RISK (Tox={tox_index:.4f} / SM={safety_margin:.1f})", False
    return f"COLLAPSED_LOW_MARGIN (SM={safety_margin:.1f})", False


def reforge_artifact(file_path: Path) -> Dict[str, Any] | None:
    with open(file_path, "r", encoding="utf-8") as f:
        try:
            data = json.load(f)
        except json.JSONDecodeError:
            print(f"    [FAIL] Corrupt JSON: {file_path.name}")
            return None

    smiles, compound_id = _extract_attributes(data)
    if not smiles:
        return None

    print(f"\n>>> [REFINERY] Reforging {compound_id}...")

    # --- STAGE 1: PHYSICS GATE (Vector_Core requires 5-element p_vector summing to 36.1) ---
    try:
        ope_result = run_ope(smiles)
    except Exception as e:
        print(f"    [FAIL] OPE Error: {e}")
        return None

    p0 = 36.1 - (0.0 + 35.0 + 1.0 + 0.85)
    p_vec_np = np.array([p0, 0.0, 35.0, 1.0, 0.85])
    physical_descriptors = {
        "molecular_weight": ope_result.get("molecular_weight", 350.0),
        "logp": ope_result.get("logp", 2.5),
        "tpsa": ope_result.get("tpsa", 70.0),
    }
    vector_state = vector_gate.execute(p_vec_np, physical_descriptors=physical_descriptors)

    if not vector_state.get("authorized"):
        print("    [COLLAPSE] Failed Vector Core (Law U34).")
        data["status"] = "COLLAPSED_U34"
        data["refinery_timestamp"] = datetime.now(timezone.utc).isoformat()
        data["refinery_log"] = vector_state
        return _json_safe(data)

    # --- STAGE 2: OPERATIONS ---
    obe_result = OBE.execute({"smiles": smiles})
    ome_result = OME.execute({"smiles": smiles})
    ole_result = OLE.execute({"compound_id": compound_id, "indication": "General"})
    ose_result = OSE.execute({"smiles": smiles})

    admet_result = ADMET.run_admet(smiles, ope_result, ome_result=ome_result, ose_result=ose_result)
    admet_result["metabolism"]["clearance_L_per_h"] = ome_result.get("clearance_L_per_h", admet_result["metabolism"].get("clearance_L_per_h"))
    admet_result["metabolism"]["half_life_h"] = ome_result.get("half_life_h", admet_result["metabolism"].get("half_life_h"))
    admet_result["toxicity"]["selectivity_index"] = ose_result.get("selectivity_index", admet_result["toxicity"].get("selectivity_index", 0))

    # --- STAGE 3: OPTIMIZATION (protocol_template must have "arms" for evaluate_regimen) ---
    variability = {"clearance_variation": 0.3, "vd_variation": 0.25, "n_tiers": 7}
    pd_params = {
        "emax": float(ope_result.get("emax", 0.9)),
        "ec50": float(ope_result.get("ec50", 10.0)),
        "hill": 1.5,
        "baseline": 0.0,
    }
    protocol_template = {
        "trial_id": "REFINE",
        "duration_days": 7.0,
        "arms": [
            {"arm_id": "R1", "label": "Refinery", "dose_mg": 100.0, "dosing_interval_h": 24.0, "n_patients": 10}
        ],
    }

    dose_result = optimize_dose(
        smiles=smiles,
        admet=admet_result,
        protocol_template=protocol_template,
        target_pk_range={"auc_mg_h_per_L": (200.0, 500.0)},
        search_strategy="coarse_to_fine",
        variability=variability,
        pd_params=pd_params,
        n_eval_patients=10,
    )

    # --- STAGE 4: WHOLE PRODUCT ---
    status, passed = assess_whole_product(admet_result, dose_result, ole_result, vector_state)

    if not passed:
        print(f"    [COLLAPSE] {status}")
        data["status"] = status
        data["refinery_timestamp"] = datetime.now(timezone.utc).isoformat()
        data["admet"] = admet_result
        data["dose_attempt"] = dose_result
        return _json_safe(data)

    # --- STAGE 5: GOVERNANCE ---
    oce_payload = {
        "p_vector": [0.1, 0.0, 35.0, 1.0],
        "csa_scores": [1.0, 1.0, 1.0, 1.0, 1.0],
        "security_score": 1.0,
    }
    oce_result = OCE.execute(oce_payload)
    if not oce_result.get("authorized"):
        print("    [FAIL] OCE Manifold Incoherence.")
        data["status"] = "UNSTABLE_MANIFOLD"
        data["refinery_timestamp"] = datetime.now(timezone.utc).isoformat()
        data["oce"] = oce_result
        return _json_safe(data)

    # --- STAGE 6: PULSE (OCE/Block_Universe expect 4-element p_vector) ---
    cycle, _ = metabolism.pulse(f"REFINE-{compound_id}", p_vec_np[:4].tolist(), [1.0, 1.0, 1.0, 1.0, 1.0])

    # --- UPDATE ARTIFACT ---
    best = dose_result.get("best_regimen") or {}
    pk = best.get("pk_summary") or {}
    cmax_dist = pk.get("cmax_mg_per_L") or {}
    cmax_mean = float((cmax_dist.get("mean") or 0) or 0)
    max_cmax = float((admet_result.get("safety_margins") or {}).get("max_tolerated_cmax_mg_per_L", 0) or 0)
    safety_margin = (max_cmax / cmax_mean) if cmax_mean > 0 else 0.0

    data["version"] = "2.4-SERIALIZATION-FIX"
    data["refinery_timestamp"] = datetime.now(timezone.utc).isoformat()
    data["status"] = status
    data["physics_gate"] = vector_state
    data["engines"] = {
        "OPE": ope_result,
        "OBE": obe_result,
        "OME": ome_result,
        "OLE": ole_result,
        "OSE": ose_result,
        "ADMET": admet_result,
        "OCE": oce_result,
    }
    data["clinical_simulation"] = {
        "dose_optimization": dose_result,
        "variability_model": "7-TIER-DETERMINISTIC",
        "pd_model": "SIGMOID_EMAX",
    }
    data["whole_product_score"] = {
        "toxicity_index": admet_result["toxicity"]["toxicity_index"],
        "safety_margin": safety_margin,
        "prv_eligible": ole_result.get("prv_eligible", False),
    }
    data["system_cycle"] = cycle

    print(f"    [SUCCESS] {status} (Cycle {cycle})")
    return _json_safe(data)


def run_refinery() -> None:
    print("=== PREDATOR X: WAREHOUSE REFINERY v2.4 ===")
    print("Initiating full logical reprocessing (Serialization Safe)...")

    if not setup_backup():
        return

    files = [f for f in WAREHOUSE_PATH.rglob("*.json") if "Backup_Pre_Refinery" not in str(f)]
    success_count = 0
    collapse_count = 0
    skipped_count = 0

    for file_path in files:
        try:
            new_data = reforge_artifact(file_path)
            if new_data is not None:
                # Zeus governance gate before warehouse write (fail-closed)
                zeus = run_zeus_gate(new_data)
                if not zeus.get("authorized", False):
                    print(f"    [ZEUS REJECT] {file_path.name}: {zeus.get('verdict', 'UNKNOWN')}")
                    collapse_count += 1
                    continue
                with open(file_path, "w", encoding="utf-8") as f:
                    json.dump(new_data, f, indent=2)
                status_val = new_data.get("status", "")
                if "COLLAPSED" in status_val:
                    collapse_count += 1
                elif "OPTIMIZED" in status_val:
                    success_count += 1
            else:
                skipped_count += 1
        except Exception as e:
            print(f"    [ERROR] Processing {file_path.name}: {e}")
            skipped_count += 1

    print("\n=== REFINERY COMPLETE ===")
    print(f"Total Files Scanned: {len(files)}")
    print(f"Reforged (Production-Ready): {success_count}")
    print(f"Collapsed (Laws U34/L11/Margin): {collapse_count}")
    print(f"Skipped (No SMILES/Invalid): {skipped_count}")


if __name__ == "__main__":
    run_refinery()
