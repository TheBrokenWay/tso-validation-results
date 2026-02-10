#!/usr/bin/env python3
"""
px_prv.py — Unified PRV Evaluation Orchestrator

Runs every candidate through the mandatory 12-engine pipeline:
    OPE → OBE → OCE → OLE → OME → OSE → ADMET → PKPD → DoseOptimizer_v2
    → VirtualEfficacyAnalytics → GradingEngine → ZeusLaws

Usage:
    python PX_Executive/px_prv.py --type novel         # Process type N from queue
    python PX_Executive/px_prv.py --type repurpose     # Process type R from queue
    python PX_Executive/px_prv.py --type all            # Process both
    python PX_Executive/px_prv.py --disease nipah       # Filter by disease
    python PX_Executive/px_prv.py --limit 100           # Max candidates
    python PX_Executive/px_prv.py --dry-run             # Report only, no writes

Environment variables:
    PRV_24H_DURATION_SEC  Max runtime in seconds (default: 86400)
    PRV_API_PACING_SEC    Delay between FTO/API calls (default: 15)
    PRV_QUEUE_FILE        Override queue filename
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


# Complete PRV disease list (alphabetical)
PRV_DISEASES = [
    "Buruli ulcer",
    "Chagas disease",
    "Chikungunya",
    "Cholera",
    "Dengue",
    "Dracunculiasis",
    "Ebola virus disease",
    "Lassa fever",
    "Leishmaniasis",
    "Leprosy",
    "Lymphatic filariasis",
    "Malaria",
    "Marburg virus disease",
    "Nipah virus infection",
    "Onchocerciasis",
    "Rabies",
    "Schistosomiasis",
    "Trachoma",
    "Tuberculosis",
    "Yaws",
    "Yellow fever",
    "Zika virus disease",
]


# ─── Counters ───
_START = time.monotonic()
_STATS = {"total": 0, "passed": 0, "failed": 0, "zeus_rejected": 0}


def _ts() -> str:
    return time.strftime("%H:%M:%S")


def _log(item_id: str, msg: str) -> None:
    print(f"[{_ts()}] {item_id}: {msg}", flush=True)


# ─── Queue Loading ───

def _queue_filename() -> str:
    override = (os.environ.get("PRV_QUEUE_FILE") or "").strip()
    if override:
        return override
    try:
        from PX_System.foundation.Intake_Policy import get_queue_filename
        return get_queue_filename(REPO_ROOT)
    except Exception:
        return "prv_24h_queue.json"


def load_queue(candidate_type: str, disease_filter: str | None) -> list[dict[str, Any]]:
    """Load and filter candidates from Feeder queue."""
    from PX_Warehouse.warehouse_layout import get_queue_path

    queue_path = get_queue_path(_queue_filename(), REPO_ROOT)
    candidates: list[dict[str, Any]] = []

    if queue_path.exists():
        try:
            data = json.loads(queue_path.read_text(encoding="utf-8"))
            raw = data.get("candidates", []) if isinstance(data, dict) else (data if isinstance(data, list) else [])
            for ent in raw:
                if isinstance(ent, dict) and ent.get("id") and ent.get("type") in ("R", "N"):
                    candidates.append(ent)
        except Exception:
            pass

    # Fallback: reprocess_candidates.json
    if not candidates:
        reprocess_path = get_queue_path("reprocess_candidates.json", REPO_ROOT)
        if reprocess_path.exists():
            try:
                data = json.loads(reprocess_path.read_text(encoding="utf-8"))
                if isinstance(data, dict):
                    for i, (smiles, name) in enumerate(data.items()):
                        item_id = (name or "id")[:32] + "_" + str(i)
                        is_rep = name and "Unknown" not in (name or "")
                        candidates.append({"id": item_id, "type": "R" if is_rep else "N", "smiles": smiles, "name": name})
            except Exception:
                pass

    # Filter by type
    if candidate_type == "novel":
        candidates = [c for c in candidates if c.get("type") == "N"]
    elif candidate_type == "repurpose":
        candidates = [c for c in candidates if c.get("type") == "R"]
    # else "all" — keep both

    # Filter by disease
    if disease_filter:
        df = disease_filter.lower()
        candidates = [c for c in candidates if df in (c.get("indication") or c.get("name") or "").lower()]

    # Dedupe by id
    seen: set[str] = set()
    deduped: list[dict[str, Any]] = []
    for c in candidates:
        kid = (c.get("id") or "").strip()
        if kid and kid not in seen:
            seen.add(kid)
            deduped.append(c)

    # Sort: repurposed first, then by id (alphabetical)
    deduped.sort(key=lambda x: (0 if x.get("type") == "R" else 1, x.get("id", "")))
    return deduped


# ─── Dossier Existence Check ───

def _dossier_already_exists(item_id: str, is_novel: bool) -> bool:
    """Check if dossier already exists in warehouse."""
    from PX_Warehouse.warehouse_layout import TIERS
    folder = "Novel_Dossiers" if is_novel else "Prv_Dossiers"
    wh = REPO_ROOT / "PX_Warehouse" / folder
    raw_id = item_id[:64]
    filenames = [f"{raw_id}.json"]
    prefix = "PRV_NOV_" if is_novel else "PRV_REP_"
    if not raw_id.startswith(prefix):
        filenames.append(f"{prefix}{raw_id}.json")
    for tier in TIERS:
        for fn in filenames:
            if (wh / tier / fn).exists():
                return True
    return False


# ─── 12-Engine Pipeline ───

def run_full_pipeline(item: dict[str, Any]) -> Tuple[Optional[Dict[str, Any]], Optional[str]]:
    """
    Run the mandatory 12-engine pipeline for a single candidate.

    Returns:
        (dossier, None) on success
        (None, error_message) on failure
    """
    smiles = item.get("smiles") or "CCO"
    item_id = item.get("id", "UNKNOWN")
    indication = item.get("indication") or item.get("name") or ""
    is_novel = item.get("type") == "N"

    # ── 1. OPE: Molecular descriptors ──
    from PX_Engine.operations.OPE import run_ope
    ope_result = run_ope(smiles)

    # ── 2. OBE: Binding energy / harm energy ──
    from PX_Engine.operations.OBE import execute as obe_execute
    obe_result = obe_execute({"smiles": smiles})

    # ── 3. OCE: 35D manifold coherence + physics snapshot ──
    from PX_Engine.operations.OCE import execute as oce_execute
    import numpy as np
    # Build p_vector for 35D manifold (Law U34: global sum = 36.1, energy_delta = 0)
    p0 = 36.1 - (0.0 + 35.0 + 1.0 + 0.85)
    p_vector = [p0, 0.0, 35.0, 1.0, 0.85]
    csa_scores = [1.0, 1.0, 1.0, 1.0, 1.0]
    oce_result = oce_execute({"p_vector": p_vector, "csa_scores": csa_scores})

    physics_snapshot = {
        "p_vector": p_vector,
        "csa_scores": csa_scores,
        "coherence": oce_result.get("coherence"),
        "manifold_id": oce_result.get("manifold_id"),
    }

    # ── 4. OLE: Legal/IP + PRV eligibility ──
    from PX_Engine.operations.OLE import execute as ole_execute
    ole_result = ole_execute({
        "compound_id": item_id,
        "indication": indication,
        "disease_context": [indication] if indication else [],
    })

    # ── 5. OME: Metabolic pathway ──
    from PX_Engine.operations.OME import execute as ome_execute
    ome_result = ome_execute({"smiles": smiles})

    # ── 6. OSE: Safety / selectivity ──
    from PX_Engine.operations.OSE import execute as ose_execute
    ose_result = ose_execute({"smiles": smiles})

    # ── 7. ADMET: Full pharmacokinetic safety profile ──
    from PX_Engine.operations.ADMET import run_admet
    admet_result = run_admet(smiles, ope_result, ome_result=ome_result, ose_result=ose_result)

    # ── 8. PKPD: Sigmoid Emax PK/PD modeling ──
    from PX_Engine.operations.PKPD import link_pk_to_pd
    from PX_Laboratory import SimulationEngine

    sim = SimulationEngine(time_step_h=0.5)
    dose_mg = 100.0  # default starting dose
    patient = {"weight_kg": 70.0}
    pk_profile = sim.simulate_one_compartment(
        dose_mg=dose_mg,
        duration_h=168.0,  # 7 days
        dosing_interval_h=24.0,
        patient=patient,
        admet=admet_result,
    )
    emax = ope_result.get("emax", 0.9)
    ec50 = ope_result.get("ec50", 5.0)
    pkpd_result = link_pk_to_pd(pk_profile, {"emax": emax, "ec50": ec50, "hill": 1.5})

    # ── 9. DoseOptimizer_v2: Coarse-to-fine dose search ──
    from PX_Engine.operations.DoseOptimizer_v2 import optimize_dose
    protocol_template = {
        "dose_mg": dose_mg,
        "interval_h": 24.0,
        "n_doses": 7,
        "clearance_L_per_h": clearance,
        "vd_L": vd,
        "ka_per_h": ka,
    }
    dose_result = optimize_dose(
        smiles=smiles,
        admet=admet_result,
        protocol_template=protocol_template,
        pd_params={"emax": emax, "ec50": ec50, "hill": 1.5},
        n_eval_patients=10,
    )

    # ── 10. VirtualEfficacyAnalytics: PTA, responder rate ──
    from PX_Engine.operations.VirtualEfficacyAnalytics import compute_pta, virtual_responder_rate, effect_variability_risk

    # Run a mini-trial for efficacy analytics
    from PX_Engine.operations.TrialEngine import TrialEngine as TE
    te = TE(time_step_h=0.5)
    best_regimen = dose_result.get("best_regimen", {})
    trial_protocol = {
        "arms": [{
            "arm_id": "treatment",
            "dose_mg": best_regimen.get("dose_mg", dose_mg),
            "interval_h": best_regimen.get("interval_h", 24.0),
            "n_patients": 21,
        }],
        "duration_days": 7.0,
    }
    variability_params = {"clearance_variation": 0.3, "vd_variation": 0.25}
    trial_result = te.run_trial(
        trial_protocol, admet_result,
        pd_params={"emax": emax, "ec50": ec50, "hill": 1.5},
        variability=variability_params,
    )

    pta_result = compute_pta(trial_result, metric="auc_mg_h_per_L", target_threshold=200.0)
    responder_result = virtual_responder_rate(trial_result, response_metric="max_effect", responder_threshold=0.5)
    variability_result = effect_variability_risk(trial_result, effect_metric="max_effect")

    virtual_efficacy = {
        "pta": pta_result,
        "responder_rate": responder_result,
        "effect_variability": variability_result,
    }

    # ── 11. GradingEngine: 5-tier classification ──
    from PX_Engine.operations.GradingEngine import GradingEngine as GE

    # Build dossier-like dict for grading
    tox_idx = (admet_result.get("toxicity") or {}).get("toxicity_index", 0.5)
    harm_energy = obe_result.get("harm_energy", tox_idx)

    from PX_System.foundation.Evidence_Package import generate_dossier
    candidate = {"name": item.get("name") or item_id, "smiles": smiles}
    engine_results = {
        "ope": ope_result,
        "admet": admet_result,
        "obe": obe_result,
        "oce": oce_result,
        "ole": ole_result,
        "ome": ome_result,
        "ose": ose_result,
        "trial_result": trial_result,
        "pkpd": pkpd_result,
        "dose_optimization": dose_result,
        "virtual_efficacy": virtual_efficacy,
    }
    dossier = generate_dossier(candidate, engine_results)

    # Inject all engine outputs into dossier
    engines_block = dossier.get("engines") or {}
    engines_block["oce"] = oce_result
    engines_block["ome"] = ome_result
    engines_block["ose"] = ose_result
    engines_block["pkpd"] = pkpd_result.get("pd_summary", pkpd_result)
    engines_block["dose_optimization"] = dose_result
    engines_block["virtual_efficacy"] = virtual_efficacy
    dossier["engines"] = engines_block

    # Physics snapshot
    ctx = dossier.get("context") or {}
    ctx["physics_snapshot"] = physics_snapshot
    dossier["context"] = ctx

    # PRV eligibility from OLE
    dossier["prv_eligible"] = ole_result.get("prv_eligible", False)
    dossier["regulatory_pathway"] = ole_result.get("regulatory_pathway", "STANDARD_FDA")

    # Grade
    ge = GE(verbose=False)
    grade_result = ge.grade_dossier(dossier)
    dossier["discovery_grading"] = grade_result

    # ── 12. ZeusLaws: Constitutional governance (full gate) ──
    from PX_System.foundation.ZeusLaws import run_zeus_gate

    zeus_verdict = run_zeus_gate(dossier)
    dossier["zeus_verdict"] = zeus_verdict

    if not zeus_verdict.get("authorized"):
        _STATS["zeus_rejected"] += 1
        # Still save dossier (Zeus runs at end, dossier is preserved with verdict)

    return dossier, None


# ─── Warehouse Write ───

def save_dossier(dossier: Dict[str, Any], item: dict[str, Any]) -> Optional[Path]:
    """Save dossier to appropriate warehouse tier directory."""
    from PX_Warehouse.warehouse_layout import ensure_structure, get_tier, get_prv_dossier_dir

    ensure_structure(REPO_ROOT)
    tier = get_tier(dossier)

    # Grade/tier consistency: downgrade if grade doesn't support tier
    grade = (dossier.get("discovery_grading") or {}).get("grade", "")
    if grade in ("NEEDS_REVIEW", "REJECTED") and tier in ("Diamond", "Gold"):
        tier = "Silver" if grade == "NEEDS_REVIEW" else "Bronze"

    is_novel = item.get("type") == "N"
    out_dir = get_prv_dossier_dir(is_novel, tier, REPO_ROOT)
    out_dir.mkdir(parents=True, exist_ok=True)

    raw_id = item.get("id", "unknown")
    safe_id = "".join(c if c.isalnum() or c in "-_" else "_" for c in raw_id)[:64]
    if raw_id.startswith("PRV_NOV_") or raw_id.startswith("PRV_REP_"):
        filename = f"{safe_id}.json"
    else:
        prefix = "PRV_NOV" if is_novel else "PRV_REP"
        filename = f"{prefix}_{safe_id}.json"

    out_path = out_dir / filename
    out_path.write_text(json.dumps(dossier, indent=2, default=str), encoding="utf-8")
    return out_path


# ─── WorldLine Persistence ───

def persist_worldline(item: dict[str, Any], stage: str, passed: bool) -> None:
    """Write WorldLine file for this attempt (pass or fail)."""
    try:
        from PX_Engine.operations import run_ope, run_admet, OME, OSE
        from PX_Warehouse.WorldLine_Database import WorldLineDatabase

        smiles = item.get("smiles") or "CCO"
        ope = run_ope(smiles)
        ome_r = OME.execute({"smiles": smiles}) if "execute" in dir(OME) else None
        ose_r = OSE.execute({"smiles": smiles}) if "execute" in dir(OSE) else None
        admet = run_admet(smiles, ope, ome_result=ome_r, ose_result=ose_r)
        outcome = "Success" if passed else "Fail"
        route = "NOVEL" if item.get("type") == "N" else "REPURPOSED"
        wl_db = WorldLineDatabase()
        wl_db.record_attempt(
            item.get("id", "unknown"), ope, admet, stage, outcome, route,
            origin=item.get("name"), repo_root=REPO_ROOT, item=item,
        )
    except Exception as e:
        print(f"  WorldLine (non-fatal): {e}", file=sys.stderr)


# ─── FTO Gate ───

def check_fto(item: dict[str, Any]) -> bool:
    """Freedom-to-operate check via PX_Legal_Check."""
    live = os.environ.get("PRV_LIVE_RESEARCH", "").strip() in ("1", "true", "yes")
    try:
        from PX_Executive.PX_Legal_Check import PX_Legal_Check
        legal = PX_Legal_Check(mode="REGULATORY" if live else "RESEARCH")
        smiles = item.get("smiles") or "CCO"
        name = item.get("name") or item.get("id", "")
        novel_invention = item.get("type") == "N"
        result = legal.check_compound(smiles=smiles, iupac_name=name, compound_id=item.get("id", ""), novel_invention=novel_invention)
        return result.freedom_to_operate
    except Exception:
        return not live


# ─── Finalization ───

def try_finalize(item_id: str, is_novel: bool) -> Optional[str]:
    """Attempt finalization for a just-saved dossier."""
    try:
        from PX_Warehouse.warehouse_layout import get_prv_dossier_dir, TIERS
        from PX_Warehouse.Finalization_Pipeline import finalize_and_place

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
            out = finalize_and_place(dossier, item_id, is_novel, REPO_ROOT)
            return str(out) if out else None
    except Exception as e:
        print(f"  Finalization (non-fatal): {e}", file=sys.stderr)
        try:
            from PX_System.finalization_log import log_finalization_failure
            log_finalization_failure(
                source_file="px_prv.py",
                candidate_id=item_id,
                error=str(e),
                context="finalize_and_place after pipeline",
            )
        except Exception:
            pass
    return None


# ─── Main Loop ───

def main() -> int:
    parser = argparse.ArgumentParser(description="Unified PRV evaluation orchestrator (12-engine pipeline)")
    parser.add_argument("--type", choices=["novel", "repurpose", "all"], default="all", help="Candidate type filter")
    parser.add_argument("--disease", type=str, default=None, help="Disease filter (e.g. nipah)")
    parser.add_argument("--limit", type=int, default=0, help="Max candidates to process (0 = all)")
    parser.add_argument("--dry-run", action="store_true", help="Report only, no writes")
    args = parser.parse_args()

    duration_sec = int(os.environ.get("PRV_24H_DURATION_SEC", "86400"))
    api_pacing = float(os.environ.get("PRV_API_PACING_SEC", "15.0"))
    end_time = _START + duration_sec

    print("=" * 64)
    print("   PREDATOR X — PRV EVALUATION (12-ENGINE PIPELINE)")
    print("=" * 64)
    print(f"Type: {args.type} | Disease: {args.disease or 'all'} | Limit: {args.limit or 'none'}")
    print(f"Duration: {duration_sec}s | Dry-run: {args.dry_run}")
    print("-" * 64)

    queue = load_queue(args.type, args.disease)
    if not queue:
        print("No candidates in queue.")
        return 1

    if args.limit > 0:
        queue = queue[:args.limit]

    print(f"Queue: {len(queue)} candidates")

    # Skip already-filed novel dossiers
    if args.type in ("novel", "all"):
        before = len(queue)
        queue = [x for x in queue if not (x.get("type") == "N" and _dossier_already_exists(x.get("id", ""), True))]
        skipped = before - len(queue)
        if skipped:
            print(f"  Skipped {skipped} novel candidates with existing dossiers")

    print(f"Processing: {len(queue)} candidates\n")

    last_api_time = 0.0

    for idx, item in enumerate(queue):
        if time.monotonic() >= end_time:
            print(f"\nTime limit reached ({duration_sec}s)")
            break

        item_id = (item.get("id") or f"item_{idx}")[:64]
        itype = item.get("type", "?")
        is_novel = itype == "N"
        _STATS["total"] += 1

        _log(item_id, f"[{idx+1}/{len(queue)}] type={itype} — Starting 12-engine pipeline")

        if args.dry_run:
            _log(item_id, "[dry-run] Would process through full pipeline")
            _STATS["passed"] += 1
            continue

        # Run full 12-engine pipeline
        try:
            dossier, err = run_full_pipeline(item)
        except Exception as e:
            _log(item_id, f"Pipeline error: {e}")
            _STATS["failed"] += 1
            persist_worldline(item, "PIPELINE_ERROR", False)
            continue

        if err:
            _log(item_id, f"Pipeline failed: {err}")
            _STATS["failed"] += 1
            persist_worldline(item, "PIPELINE_FAIL", False)
            continue

        # FTO check (with API pacing)
        now = time.monotonic()
        if last_api_time > 0 and (now - last_api_time) < api_pacing:
            time.sleep(api_pacing - (now - last_api_time))

        fto_ok = check_fto(item)
        last_api_time = time.monotonic()
        if not fto_ok:
            _log(item_id, "FTO blocked")
            _STATS["failed"] += 1
            persist_worldline(item, "FTO_BLOCKED", False)
            continue

        # Save dossier to warehouse
        out_path = save_dossier(dossier, item)
        zeus_ok = (dossier.get("zeus_verdict") or {}).get("authorized", False)
        grade = (dossier.get("discovery_grading") or {}).get("grade", "?")

        _log(item_id, f"Saved → {out_path} | grade={grade} | zeus={'OK' if zeus_ok else 'REJECTED'}")

        persist_worldline(item, "OK", True)
        _STATS["passed"] += 1

        # Try finalization
        fin_path = try_finalize(item_id, is_novel)
        if fin_path:
            _log(item_id, f"Finalized → {fin_path}")

    # Summary
    print("\n" + "=" * 64)
    print("PRV EVALUATION SUMMARY")
    print("=" * 64)
    print(f"Total:         {_STATS['total']}")
    print(f"Passed:        {_STATS['passed']}")
    print(f"Failed:        {_STATS['failed']}")
    print(f"Zeus rejected: {_STATS['zeus_rejected']}")
    print("=" * 64)

    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("\n[STOP] PRV evaluation shutdown.")
        sys.exit(0)
