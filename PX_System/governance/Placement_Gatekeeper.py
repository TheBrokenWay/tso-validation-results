"""
Placement Gatekeeper — routes dossiers to canonical warehouse folders (Prv_Dossiers / Novel_Dossiers by tier).
Uses PX_Warehouse.warehouse_layout for single source of truth.
"""
from pathlib import Path
from typing import Any, Dict, Tuple

def _get_tier(dossier: Dict[str, Any]) -> str:
    """Delegates to canonical PX_Warehouse.warehouse_layout.get_tier."""
    from PX_Warehouse.warehouse_layout import get_tier
    return get_tier(dossier)


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
