"""
Placement Gatekeeper — routes dossiers to canonical warehouse folders (Prv_Dossiers / Novel_Dossiers by tier).
Uses PX_Warehouse.warehouse_layout for single source of truth.
"""
from pathlib import Path
from typing import Any, Dict, Tuple

def _get_tier(dossier: Dict[str, Any]) -> str:
    """Derive tier from dossier (toxicity_index + safety_margin for Diamond)."""
    try:
        from PX_Warehouse.warehouse_layout import get_tier
        return get_tier(dossier)
    except Exception:
        pass
    tox = 1.0
    engines = dossier.get("engines") or {}
    admet = engines.get("admet") or {}
    if isinstance(admet, dict):
        t = (admet.get("toxicity") or {}).get("toxicity_index")
        if t is not None:
            tox = float(t)
        elif admet.get("toxicity_index") is not None:
            tox = float(admet["toxicity_index"])
    safety = float(dossier.get("safety_margin") or 0)
    if tox < 0.0100 or (tox > 0.0200 and safety > 50.0):
        return "Diamond"
    if tox < 0.0200:
        return "Gold"
    if tox < 0.0210:
        return "Silver"
    return "Bronze"


class PlacementGatekeeper:
    """Routes dossiers to Prv_Dossiers/<tier> or Novel_Dossiers/<tier>."""

    @staticmethod
    def resolve_placement(dossier: Dict[str, Any], candidate_type: str) -> Tuple[str, str]:
        """
        Returns (root_folder, tier_name).
        candidate_type "R" = repurposed -> Prv_Dossiers; "N" = novel -> Novel_Dossiers.
        """
        tier = _get_tier(dossier)
        root = "Novel_Dossiers" if candidate_type == "N" else "Prv_Dossiers"
        return root, tier

    @staticmethod
    def get_secure_path(repo_root: Path, filename: str, dossier: Dict[str, Any], candidate_type: str) -> Path:
        """Return Path: PX_Warehouse/<root_folder>/<tier>/<filename>."""
        root_folder, tier_name = PlacementGatekeeper.resolve_placement(dossier, candidate_type)
        return Path(repo_root) / "PX_Warehouse" / root_folder / tier_name / filename
