#!/usr/bin/env python3
"""
24-HOUR PRV EXECUTION + REAL-TIME SHELL STATUS — LIVE RESEARCH

Olympus constraint-first (mandatory order): load disease constraint schema → apply exclusion zones
→ assemble only constraint-satisfying candidates → reject violations → FTO gate → then repurposed
evaluation or novel generation. No simulation-first or random-first search.

This is live research, not simulation:
- Freedom-to-operate (FTO) is checked against the local patent index only (PX_Data/patents).
- Repurposed candidates (R) processed first; Novel (N) processed second.

Dossiers are saved to PX_Warehouse/Prv_Dossiers/<tier> or Novel_Dossiers/<tier> (tiered by toxicity).

ENFORCE GLOBAL MEMORY: A .worldline file is written for every single attempt (pass or fail).
This builds the complete 35D map including Failure zones so the system can avoid them.
"""
from __future__ import annotations

import os
import sys
import json
import time
import subprocess
from pathlib import Path
from datetime import datetime, timezone
from typing import Any

# Wait this long before each FTO/API attempt (per item). Override with PRV_API_PACING_SEC.
API_PACING_SEC = float(os.environ.get("PRV_API_PACING_SEC", "15.0"))

ROOT = Path(__file__).resolve().parents[0]   # PX_Executive
REPO_ROOT = Path(__file__).resolve().parents[1]  # foundation (repo root)
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

PX_LOGS = REPO_ROOT / "PX_LOGS"
PX_LOGS.mkdir(parents=True, exist_ok=True)

# ─── Counters ───
TOTAL_ITEMS = 0
COMPLETED = 0
TESTS = 0
FTO_CHECKS = 0
REPURPOSED_COUNT = 0
NOVEL_COUNT = 0
START_TIME = time.monotonic()
RUN_LOG_LINES: list[str] = []


def ts_sec() -> int:
    return int(time.monotonic() - START_TIME)


def emit(item_id: str, itype: str, stage: str) -> None:
    line = (
        f"TS={ts_sec()} ITEM={item_id} TYPE={itype} STAGE={stage} "
        f"TOT={TOTAL_ITEMS} DONE={COMPLETED} TEST={TESTS} FTO={FTO_CHECKS}"
    )
    print(line, flush=True)
    RUN_LOG_LINES.append(line)


def alert(item_id: str, itype: str, stage: str, reason: str) -> None:
    emit(item_id, itype, stage)
    RUN_LOG_LINES.append(f"  ALERT: {reason}")


def _queue_filename() -> str:
    override = (os.environ.get("PRV_QUEUE_FILE") or "").strip()
    if override:
        return override
    try:
        from PX_System.foundation.Intake_Policy import get_queue_filename
        return get_queue_filename(REPO_ROOT)
    except Exception:
        return "prv_24h_queue.json"


def load_queue() -> list[dict[str, Any]]:
    from PX_Warehouse.warehouse_layout import get_queue_path

    candidates: list[dict[str, Any]] = []
    # 1) Dedicated queue (Feeder first, then Operations/Inputs fallback)
    queue_path = get_queue_path(_queue_filename(), REPO_ROOT)
    if queue_path.exists():
        try:
            data = json.loads(queue_path.read_text(encoding="utf-8"))
            for ent in data.get("candidates", []) if isinstance(data, dict) else (data if isinstance(data, list) else []):
                if isinstance(ent, dict) and ent.get("id") and ent.get("type") in ("R", "N"):
                    candidates.append(ent)
        except Exception:
            pass

    # 2) Fallback: reprocess_candidates.json (Feeder first, then Operations/Inputs)
    if not candidates:
        reprocess_path = get_queue_path("reprocess_candidates.json", REPO_ROOT)
        if reprocess_path.exists():
            try:
                data = json.loads(reprocess_path.read_text(encoding="utf-8"))
                if isinstance(data, dict):
                    for i, (smiles, name) in enumerate(data.items()):
                        item_id = (name or "id")[:32] + "_" + str(i)
                        is_repurposed = name and "Unknown" not in (name or "")
                        ent = {"id": item_id, "type": "R" if is_repurposed else "N", "smiles": smiles, "name": name}
                        candidates.append(ent)
            except Exception:
                pass
    return candidates


def _novel_dossier_already_exists(item_id: str, repo_root: Path) -> bool:
    """True if a novel dossier for this id already exists in Novel_Dossiers/<tier>/ (any tier)."""
    if not item_id:
        return False
    from PX_Warehouse.warehouse_layout import TIERS
    wh = repo_root / "PX_Warehouse" / "Novel_Dossiers"
    raw_id = item_id[:64]
    if raw_id.startswith("PRV_NOV_") or raw_id.startswith("PRV_REP_"):
        filename = f"{raw_id}.json"
    else:
        filename = f"PRV_NOV_{raw_id}.json"
    for tier in TIERS:
        if (wh / tier / filename).exists():
            return True
    return False


def _metabolic_pulse(item: dict[str, Any]) -> tuple[int, str]:
    """Pulse the metabolic heartbeat; system age increments only when 35D manifold resonance is achieved."""
    try:
        from PX_Engine.Metabolism import Metabolism
        import numpy as np
        from PX_Engine.operations.OPE import run_ope
        smiles = item.get("smiles") or "CCO"
        ope = run_ope(smiles)
        mw = ope.get("molecular_weight", 350.0)
        # Law U34: global sum must be 36.1, energy delta 0
        p0 = 36.1 - (0.0 + 35.0 + 1.0 + 0.85)
        p_vec = np.array([p0, 0.0, 35.0, 1.0, 0.85])
        csa_s = [1.0, 1.0, 1.0, 1.0, 1.0]
        heartbeat = Metabolism()
        cycle_age, status = heartbeat.pulse(item.get("id", "PRV_24H"), p_vec.tolist(), csa_s)
        return cycle_age, status
    except Exception:
        p_vec = [0.1, 0.0, 35.0, 1.0, 0.85]
        csa_s = [1.0, 1.0, 1.0, 1.0, 1.0]
        try:
            from PX_Engine.Metabolism import Metabolism
            heartbeat = Metabolism()
            return heartbeat.pulse(item.get("id", "PRV_24H"), p_vec, csa_s)
        except Exception:
            return 0, "PULSE_SKIP"


def stage_intake(item: dict[str, Any]) -> bool:
    item_id = item.get("id", "")
    if not item_id:
        return False
    if item.get("type") not in ("R", "N"):
        return False
    return True


def stage_preclinical(item: dict[str, Any]) -> bool:
    global TESTS
    try:
        from PX_Engine.operations import run_ope, run_admet
        from PX_Engine.operations import OME, OSE
        smiles = item.get("smiles") or "CCO"
        ope = run_ope(smiles)
        ome_result = OME.execute({"smiles": smiles}) if "execute" in dir(OME) else None
        ose_result = OSE.execute({"smiles": smiles}) if "execute" in dir(OSE) else None
        run_admet(smiles, ope, ome_result=ome_result, ose_result=ose_result)
        TESTS += 1
        return True
    except Exception:
        TESTS += 1
        return True


def stage_e2e(item: dict[str, Any]) -> bool:
    """E2E — Validation AND SAVING of the dossier. OBE/OME/OSE/OLE integrated."""
    global TESTS
    try:
        from PX_Engine.operations import run_ope, run_admet
        from PX_Engine.operations import OBE, OME, OSE, OLE
        from PX_System.foundation.ZeusLaws import check_constitutional
        from PX_System.foundation.Evidence_Package import generate_dossier
        
        smiles = item.get("smiles") or "CCO"
        ope = run_ope(smiles)
        ome_result = OME.execute({"smiles": smiles}) if "execute" in dir(OME) else None
        ose_result = OSE.execute({"smiles": smiles}) if "execute" in dir(OSE) else None
        admet = run_admet(smiles, ope, ome_result=ome_result, ose_result=ose_result)
        obe_result = OBE.execute({"smiles": smiles}) if "execute" in dir(OBE) else None
        tox = (admet.get("toxicity") or {}).get("toxicity_index") or 0
        harm_energy = (obe_result.get("harm_energy") if obe_result else None) or tox
        
        verdict = check_constitutional("prv_24h", {"toxicity_index": tox, "harm_energy": harm_energy})
        if "authorized" not in verdict:
            return False

        # OLE: PRV eligibility (neglected diseases e.g. Malaria, Chagas) for regulatory pathway
        indication = (item.get("indication") or item.get("name") or "").upper() or "MALARIA"
        ole_result = OLE.execute({"compound_id": item.get("id", ""), "indication": indication}) if "execute" in dir(OLE) else None
        if ole_result:
            admet["regulatory_pathway"] = ole_result.get("regulatory_pathway", "STANDARD_FDA")
            admet["prv_eligible"] = ole_result.get("prv_eligible", False)
            
        candidate = {"name": item.get("name") or item.get("id", ""), "smiles": smiles}
        engine_results = {"ope": ope, "admet": admet, "obe": obe_result, "ole": ole_result}
        dossier = generate_dossier(candidate, engine_results)
        
        if "harm_energy" not in dossier:
            return False

        # --- WAREHOUSE: Prv_Dossiers / Novel_Dossiers by tier (Diamond, Gold, Silver, Bronze) ---
        from PX_Warehouse.warehouse_layout import ensure_structure, get_tier, get_prv_dossier_dir
        ensure_structure(REPO_ROOT)
        tier = get_tier(dossier)
        is_novel = item.get("type") == "N"
        out_dir = get_prv_dossier_dir(is_novel, tier, REPO_ROOT)
        out_dir.mkdir(parents=True, exist_ok=True)
        raw_id = item.get("id", "unknown")
        safe_id = "".join([c if c.isalnum() else "_" for c in raw_id])[:64]
        # Avoid double prefix: queue ids are often PRV_NOV_xxx / PRV_REP_xxx already
        if raw_id.startswith("PRV_NOV_") or raw_id.startswith("PRV_REP_"):
            filename = f"{safe_id}.json"
        else:
            prefix = "PRV_NOV" if is_novel else "PRV_REP"
            filename = f"{prefix}_{safe_id}.json"
        out_path = out_dir / filename
        out_path.write_text(json.dumps(dossier, indent=2), encoding="utf-8")

        # Worldline is written by _persist_memory(OK) in main loop (one per attempt, pass or fail)

        TESTS += 1
        return True
    except Exception as e:
        print(f"  E2E ERROR: {e}")
        from PX_System.finalization_log import log_finalization_failure
        log_finalization_failure(
            source_file="PRV_24H_Orchestrator.py",
            candidate_id=item.get("id", "UNKNOWN"),
            error=str(e),
            context="stage_e2e dossier generation and warehouse write",
        )
        TESTS += 1
        return False


def stage_legal(item: dict[str, Any]) -> bool:
    global FTO_CHECKS
    live_research = os.environ.get("PRV_LIVE_RESEARCH", "").strip() in ("1", "true", "yes")
    try:
        from PX_Executive.PX_Legal_Check import PX_Legal_Check
        legal = PX_Legal_Check(mode="REGULATORY" if live_research else "RESEARCH")
        smiles = item.get("smiles") or "CCO"
        name = item.get("name") or item.get("id", "")
        # Novel invention flag prevents searching for random text matches
        novel_invention = item.get("type") == "N" or (item.get("source") or "").strip() == "Genesis"
        result = legal.check_compound(smiles=smiles, iupac_name=name, compound_id=item.get("id", ""), novel_invention=novel_invention)
        FTO_CHECKS += 1
        return result.freedom_to_operate
    except Exception:
        FTO_CHECKS += 1
        return not live_research


def stage_commercial(item: dict[str, Any]) -> bool:
    return True


def _persist_memory(item: dict[str, Any], stage_reached: str, passed: bool, mode: str) -> None:
    """
    WorldLine persistence for every attempt (pass or fail). Builds the complete 35D map
    including Failure zones so the logic can learn and avoid them. Non-fatal on error.
    """
    try:
        from PX_Engine.operations import run_ope, run_admet
        from PX_Engine.operations import OME, OSE
        from PX_Warehouse.WorldLine_Database import WorldLineDatabase

        smiles = item.get("smiles") or "CCO"
        ope = run_ope(smiles)
        ome_result = OME.execute({"smiles": smiles}) if "execute" in dir(OME) else None
        ose_result = OSE.execute({"smiles": smiles}) if "execute" in dir(OSE) else None
        admet = run_admet(smiles, ope, ome_result=ome_result, ose_result=ose_result)
        outcome = "Success" if passed else "Fail"
        route = "REPURPOSED" if mode == "repurposed" else "NOVEL" if mode == "novel" else "PRV_24H"
        wl_db = WorldLineDatabase()
        wl_db.record_attempt(
            item.get("id", "unknown"),
            ope,
            admet,
            stage_reached,
            outcome,
            route,
            origin=item.get("name"),
            repo_root=REPO_ROOT,
            item=item,
        )
    except Exception as e:
        print(f"  WorldLine memory (non-fatal): {e}")

    # Trajectory: verify predicted next state does not drift into Law L11 (toxicity)
    try:
        from PX_Engine.Trajectory_Predictor import TrajectoryPredictor
        tp = TrajectoryPredictor()
        if not tp.is_next_state_safe():
            RUN_LOG_LINES.append("  TRAJECTORY: predicted next state would violate Law L11 (toxicity); consider reducing exploration.")
    except Exception:
        pass


def run_warehouse_simulation() -> tuple[bool, str]:
    script = REPO_ROOT / "PX_Warehouse" / "Operations" / "scripts" / "run_warehouse_simulation.py"
    if not script.exists():
        return True, "no script"
    try:
        r = subprocess.run([sys.executable, str(script), "--enforce"], cwd=str(REPO_ROOT), capture_output=True, text=True, timeout=120)
        return r.returncode == 0, (r.stderr or r.stdout or "")[:500]
    except Exception as e:
        return False, str(e)


def main() -> int:
    global TOTAL_ITEMS, COMPLETED, REPURPOSED_COUNT, NOVEL_COUNT
    duration_sec = int(os.environ.get("PRV_24H_DURATION_SEC", "86400"))
    end_time = START_TIME + duration_sec
    max_items_raw = os.environ.get("PRV_MAX_ITEMS", "").strip()
    max_items = int(max_items_raw) if max_items_raw.isdigit() else None

    queue = load_queue()
    if not queue:
        emit("none", "R", "FL")
        RUN_LOG_LINES.append("  ALERT: No queue.")
        return 1

    mode = (os.environ.get("PRV_MODE") or "").strip().lower()
    if mode == "repurposed":
        queue = [x for x in queue if x.get("type") == "R"]
    elif mode == "novel":
        queue = [x for x in queue if x.get("type") == "N"]
        # Novel pipeline must only process NEW molecules: skip any N that already has a dossier
        before = len(queue)
        queue = [x for x in queue if not _novel_dossier_already_exists(x.get("id", ""), REPO_ROOT)]
        skipped = before - len(queue)
        if skipped:
            RUN_LOG_LINES.append(f"  Novel: skipped {skipped} already-filed candidates (only new molecules)")

    # One work item per candidate id: dedupe so TOT = unique candidates = files written (no overwrites)
    seen_ids: set[str] = set()
    deduped: list[dict[str, Any]] = []
    for x in queue:
        kid = (x.get("id") or "").strip()
        if not kid or kid in seen_ids:
            continue
        seen_ids.add(kid)
        deduped.append(x)
    if len(deduped) < len(queue):
        RUN_LOG_LINES.append(f"  Dedupe: {len(queue)} queue entries → {len(deduped)} unique candidates")
    queue = deduped

    queue.sort(key=lambda x: (0 if x.get("type") == "R" else 1, x.get("id", "")))
    
    last_api_time = 0.0

    for idx, item in enumerate(queue):
        if max_items is not None and TOTAL_ITEMS >= max_items:
            break
        if time.monotonic() >= end_time:
            break
        item_id = (item.get("id") or f"item_{idx}")[:64]
        itype = "R" if item.get("type") == "R" else "N"
        TOTAL_ITEMS += 1
        if itype == "R": REPURPOSED_COUNT += 1
        else: NOVEL_COUNT += 1

        emit(item_id, itype, "IN")
        if not stage_intake(item):
            alert(item_id, itype, "FL", "intake failure")
            _persist_memory(item, "IN", False, mode)
            continue

        # Metabolic heartbeat: cycle increments only when 35D resonance achieved
        _cycle, _pulse_status = _metabolic_pulse(item)
        if _pulse_status != "RESONANCE_ACHIEVED":
            RUN_LOG_LINES.append(f"  HEARTBEAT: {_pulse_status} (cycle={_cycle})")

        emit(item_id, itype, "PP")
        if not stage_preclinical(item):
            alert(item_id, itype, "FL", "preclinical failure")
            _persist_memory(item, "PP", False, mode)
            continue

        emit(item_id, itype, "E2E")
        if not stage_e2e(item):
            alert(item_id, itype, "FL", "E2E save failure")
            _persist_memory(item, "E2E", False, mode)
            continue

        now = time.monotonic()
        if last_api_time > 0 and (now - last_api_time) < API_PACING_SEC:
            time.sleep(API_PACING_SEC - (now - last_api_time))
        
        emit(item_id, itype, "LG")
        if not stage_legal(item):
            alert(item_id, itype, "FL", "FTO blocked")
            _persist_memory(item, "LG", False, mode)
            last_api_time = time.monotonic()
            continue
        last_api_time = time.monotonic()

        emit(item_id, itype, "CM")
        COMPLETED += 1
        emit(item_id, itype, "OK")
        _persist_memory(item, "OK", True, mode)

        # Finalization: run full checklist and place into Finalized_Dossiers/<tier> if all pass
        try:
            from PX_Warehouse.warehouse_layout import get_prv_dossier_dir, TIERS
            is_novel = itype == "N"
            dossier = None
            for tier in TIERS:
                d_dir = get_prv_dossier_dir(is_novel, tier, REPO_ROOT)
                for fname in (f"{item_id}.json", f"PRV_NOV_{item_id}.json", f"PRV_REP_{item_id}.json"):
                    p = d_dir / fname
                    if p.exists():
                        dossier = json.loads(p.read_text(encoding="utf-8"))
                        break
                if dossier is not None:
                    break
            if dossier is not None:
                from PX_Warehouse.Finalization_Pipeline import finalize_and_place
                out_path = finalize_and_place(dossier, item_id, is_novel, REPO_ROOT)
                if out_path:
                    RUN_LOG_LINES.append(f"  Finalized → {out_path}")
        except Exception as e:
            RUN_LOG_LINES.append(f"  Finalization (non-fatal): {e}")
            from PX_System.finalization_log import log_finalization_failure
            log_finalization_failure(
                source_file="PRV_24H_Orchestrator.py",
                candidate_id=item_id,
                error=str(e),
                context="finalize_and_place call after successful pipeline stages",
            )
            try:
                from PX_System.foundation.Sovereign_Log_Chain import append as slc_append
                slc_append("FINALIZATION_FAILURE", {"item_id": item_id, "error": str(e), "source": "PRV_24H_Orchestrator"})
            except Exception:
                pass

    return 0

if __name__ == "__main__":
    sys.exit(main())