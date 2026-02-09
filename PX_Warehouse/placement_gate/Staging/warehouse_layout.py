"""
Canonical warehouse layout:
  Feeder              — Loading dock for raw queues (prv_24h_queue.json, reprocess_candidates.json, etc.)
  Calibration_Molecules — QA vault for reference standards (TRIAL_SIMULATION, Ibuprofen, Ethanol, etc.)
  Prv_Dossiers / Novel_Dossiers — Tiered output (Diamond, Gold, Silver, Bronze)
  Finalized_Dossiers  — Tiered finalized dossiers (Bronze, Silver, Diamond, Gold)
  Learning_Material   — General logs and non-calibration trial output.

All placement must go through PX_Warehouse.placement_gate before being sorted (place_* functions + audit).
This module provides path resolution and tier logic; placement_gate is the single entry point for writes.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict

# Staging(0) -> placement_gate(1) -> PX_Warehouse(2) -> repo(3)
_REPO_ROOT = Path(__file__).resolve().parents[3]
WAREHOUSE_ROOT = _REPO_ROOT / "PX_Warehouse"

# Phase 2 zones
FEEDER_ROOT = WAREHOUSE_ROOT / "Feeder"
CALIBRATION_MOLECULES_ROOT = WAREHOUSE_ROOT / "Calibration_Molecules"

# Dossier and learning
PRV_DOSSIERS_ROOT = WAREHOUSE_ROOT / "Prv_Dossiers"
NOVEL_DOSSIERS_ROOT = WAREHOUSE_ROOT / "Novel_Dossiers"
FINALIZED_DOSSIERS_ROOT = WAREHOUSE_ROOT / "Finalized_Dossiers"
LEARNING_MATERIAL_ROOT = WAREHOUSE_ROOT / "Learning_Material"
WORLD_LINES_ROOT = WAREHOUSE_ROOT / "WorldLines"

# Legacy fallback (read-only fallback for queue files)
OPERATIONS_INPUTS_ROOT = WAREHOUSE_ROOT / "Operations" / "Inputs"

TIERS = ("Diamond", "Gold", "Silver", "Bronze")

# Two-tier bar: Zeus at 0.021 → Finalized_Dossiers/<tier>; discovery bar → Discovery_Accepted
ZEUS_TOXICITY_THRESHOLD = 0.0210  # Governance lock; do not relax for regulatory path
DISCOVERY_ACCEPTED_TOX_THRESHOLD = 0.05  # 5% Toxicity Limit for Discovery Review


def get_tier(data: Dict[str, Any]) -> str:
    """
    Tier from dossier/result data (toxicity_index + optional safety_margin).
    Diamond: tox < 0.01 OR (tox > 0.02 and safety_margin > 50)
    Gold: tox < 0.021
    Silver: tox < 0.10
    Bronze: else
    """
    try:
        tox = 1.0
        if isinstance(data.get("engines"), dict) and isinstance(data["engines"].get("admet"), dict):
            admet_tox = (data["engines"]["admet"].get("toxicity") or {}).get("toxicity_index")
            if admet_tox is not None:
                tox = float(admet_tox)
        if tox == 1.0 and "candidate_profile" in data and data["candidate_profile"]:
            cp_tox = data["candidate_profile"].get("toxicity_index")
            if cp_tox is not None:
                tox = float(cp_tox)
        if tox == 1.0 and "admet_analysis" in data and data["admet_analysis"]:
            admet = data["admet_analysis"]
            if isinstance(admet.get("toxicity"), dict) and "toxicity_index" in admet["toxicity"]:
                tox = float(admet["toxicity"]["toxicity_index"])
        if tox == 1.0:
            physics = data.get("physics", {}) or data.get("engines", {}).get("ope", {}) or {}
            tox = physics.get("toxicity_index", 1.0)
        if data.get("toxicity_index") is not None:
            tox = float(data["toxicity_index"])

        safety_margin = float(data.get("safety_margin", 0))

        if tox < 0.01 or (tox > 0.02 and safety_margin > 50.0):
            return "Diamond"
        if tox < 0.021:
            return "Gold"
        if tox < 0.10:
            return "Silver"
        return "Bronze"
    except Exception:
        return "Bronze"


def get_prv_dossier_dir(is_novel: bool, tier: str, repo_root: Path | None = None) -> Path:
    """Root directory for a PRV dossier: Prv_Dossiers/<tier> or Novel_Dossiers/<tier>."""
    root = repo_root or _REPO_ROOT
    wh = root / "PX_Warehouse"
    if is_novel:
        return wh / "Novel_Dossiers" / tier
    return wh / "Prv_Dossiers" / tier


def get_finalized_dossier_dir(tier: str, repo_root: Path | None = None) -> Path:
    """Root directory for a finalized dossier: Finalized_Dossiers/<tier>."""
    root = repo_root or _REPO_ROOT
    wh = root / "PX_Warehouse" / "Finalized_Dossiers"
    return wh / tier if tier in TIERS else wh


def get_discovery_accepted_dir(repo_root: Path | None = None) -> Path:
    """Root for dossiers that failed Zeus (0.021) but passed Discovery Bar (0.05)."""
    root = repo_root or _REPO_ROOT
    wh = root / "PX_Warehouse" / "Finalized_Dossiers"
    return wh / "Discovery_Accepted"


def get_learning_material_dir(repo_root: Path | None = None) -> Path:
    """Root for trial simulations and general learning material."""
    root = repo_root or _REPO_ROOT
    return root / "PX_Warehouse" / "Learning_Material"


def get_feeder_dir(repo_root: Path | None = None) -> Path:
    """Loading dock for queue files (prv_24h_queue.json, reprocess_candidates.json, etc.)."""
    root = repo_root or _REPO_ROOT
    return root / "PX_Warehouse" / "Feeder"


def get_calibration_molecules_dir(repo_root: Path | None = None) -> Path:
    """QA vault for reference standards (TRIAL_SIMULATION, Ibuprofen, Ethanol, etc.)."""
    root = repo_root or _REPO_ROOT
    return root / "PX_Warehouse" / "Calibration_Molecules"


def get_worldline_dir(tier: str, repo_root: Path | None = None) -> Path:
    """WorldLines folder by tier (Diamond, Gold, Silver, Bronze). Every discovery gets a worldline here for 35D memory."""
    root = repo_root or _REPO_ROOT
    wh = root / "PX_Warehouse" / "WorldLines"
    if tier in TIERS:
        return wh / tier
    return wh


def get_queue_path(queue_filename: str, repo_root: Path | None = None) -> Path:
    """
    Resolve queue file path: Feeder first (canonical), then Operations/Inputs (fallback).
    Use this for reading. For writing, use get_feeder_dir() / filename so Feeder is canonical.
    """
    root = repo_root or _REPO_ROOT
    wh = root / "PX_Warehouse"
    feeder_file = wh / "Feeder" / queue_filename
    inputs_file = wh / "Operations" / "Inputs" / queue_filename
    if feeder_file.exists():
        return feeder_file
    return inputs_file


def ensure_structure(repo_root: Path | None = None) -> None:
    """Create all canonical zones and tier subdirs (Feeder, Calibration_Molecules, Prv/Novel, Finalized_Dossiers, Learning_Material, WorldLines)."""
    root = repo_root or _REPO_ROOT
    wh = root / "PX_Warehouse"
    (wh / "Feeder").mkdir(parents=True, exist_ok=True)
    (wh / "Calibration_Molecules").mkdir(parents=True, exist_ok=True)
    for folder in ("Prv_Dossiers", "Novel_Dossiers"):
        (wh / folder).mkdir(parents=True, exist_ok=True)
        for tier in TIERS:
            (wh / folder / tier).mkdir(exist_ok=True)
    (wh / "Finalized_Dossiers").mkdir(parents=True, exist_ok=True)
    for tier in TIERS:
        (wh / "Finalized_Dossiers" / tier).mkdir(exist_ok=True)
    (get_finalized_dossier_dir("Bronze", root).parent / "Discovery_Accepted").mkdir(parents=True, exist_ok=True)
    (wh / "Learning_Material").mkdir(parents=True, exist_ok=True)
    (wh / "WorldLines").mkdir(parents=True, exist_ok=True)
    for tier in TIERS:
        (wh / "WorldLines" / tier).mkdir(exist_ok=True)


__all__ = [
    "WAREHOUSE_ROOT",
    "FEEDER_ROOT",
    "CALIBRATION_MOLECULES_ROOT",
    "PRV_DOSSIERS_ROOT",
    "NOVEL_DOSSIERS_ROOT",
    "FINALIZED_DOSSIERS_ROOT",
    "LEARNING_MATERIAL_ROOT",
    "WORLD_LINES_ROOT",
    "OPERATIONS_INPUTS_ROOT",
    "TIERS",
    "ZEUS_TOXICITY_THRESHOLD",
    "DISCOVERY_ACCEPTED_TOX_THRESHOLD",
    "get_tier",
    "get_prv_dossier_dir",
    "get_finalized_dossier_dir",
    "get_discovery_accepted_dir",
    "get_learning_material_dir",
    "get_worldline_dir",
    "get_feeder_dir",
    "get_calibration_molecules_dir",
    "get_queue_path",
    "ensure_structure",
]
