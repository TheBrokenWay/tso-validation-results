"""
WorldLine Database - persistence for 35D manifold materializations.
Provides path to WorldLines directory and record_materialization for audit/orchestrator integration.
Every discovery creates a worldline in WorldLines/<tier> so the logic remembers 35D work and improves accuracy.
"""
import os
import json
from pathlib import Path
from typing import Any, Dict
from datetime import datetime, timezone

# Default path: PX_Warehouse/WorldLines (repo root = Staging -> placement_gate -> PX_Warehouse -> repo)
_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
DEFAULT_WORLDLINES_PATH = os.path.join(_REPO_ROOT, "PX_Warehouse", "WorldLines")


class WorldLineDatabase:
    """
    WorldLine persistence: .path for scan/listdir, record_materialization to write .worldline files.
    """

    def __init__(self, path: str | None = None):
        self._path = path or DEFAULT_WORLDLINES_PATH
        os.makedirs(self._path, exist_ok=True)

    @property
    def path(self) -> str:
        return self._path

    def record_materialization(
        self,
        task_id: str,
        block,
        coherence: float,
        lab_results: dict,
        route: str = "BLUE",
        cycle: int | None = None,
        origin: str | None = None,
        parent_id: str | None = None,
    ) -> str:
        """
        Write a WorldLine artifact. Returns path to the written file.
        block: list or list-like of 35 floats (coordinate_35d).
        """
        safe_id = "".join(c if c.isalnum() or c in "-_" else "_" for c in str(task_id))[:64]
        fname = f"{safe_id}.worldline"
        filepath = os.path.join(self._path, fname)
        payload = {
            "header": {
                "worldline_id": task_id,
                "system_version": "1.2.0-GAIP",
                "route": route,
            },
            "coordinate_35d": list(block)[:35] if block else [0.0] * 35,
            "coherence_amplitude": coherence,
            "route": route,
            "physical_realization": lab_results,
        }
        if cycle is not None:
            payload["metabolic_cycle"] = cycle
        if origin:
            payload["origin"] = origin
        if parent_id is not None:
            payload["parent_id"] = parent_id
        with open(filepath, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
        return filepath

    def record_attempt(
        self,
        item_id: str,
        ope: Dict[str, Any],
        admet: Dict[str, Any],
        stage_reached: str,
        outcome: str,
        route: str,
        origin: str | None = None,
        repo_root: Path | str | None = None,
        item: Dict[str, Any] | None = None,
    ) -> str:
        """
        Write a worldline for every attempt (pass or fail). Builds the complete 35D map
        including Failure zones so the system can avoid them. Called from orchestrator _persist_memory.
        """
        from PX_Warehouse.warehouse_layout import get_worldline_dir, TIERS

        tox = (admet.get("toxicity") or {})
        risk_level = tox.get("risk_level", "TOXICITY_GOLD")
        toxicity_index = float(tox.get("toxicity_index", 0.02))
        safety_margins = admet.get("safety_margins") or {}
        safety_margin = float(safety_margins.get("safety_margin", 0.0))
        # Map risk_level to tier for folder placement
        tier = "Diamond" if risk_level == "TOXICITY_DIAMOND" else "Gold" if risk_level == "TOXICITY_GOLD" else "Silver" if risk_level == "TOXICITY_SILVER" else "Bronze"
        if tier not in TIERS:
            tier = "Bronze"

        root = Path(repo_root) if repo_root else Path(_REPO_ROOT)
        wl_dir = get_worldline_dir(tier, root)
        wl_dir.mkdir(parents=True, exist_ok=True)

        safe_id = "".join(c if c.isalnum() or c in "-_" else "_" for c in str(item_id))[:64]
        task_id = f"WL-{safe_id}"
        fname = f"{safe_id}.worldline"
        filepath = os.path.join(str(wl_dir), fname)

        coord = [
            float(ope.get("logp", 2.0)),
            float(ope.get("molecular_weight", 350.0)),
            float(ope.get("tpsa", 70.0)),
            float(ope.get("ec50", 1.0)),
            toxicity_index,
            float(tox.get("herg_risk", 0.2)),
            float(tox.get("hepatotoxicity_risk", 0.2)),
            float(tox.get("ames_mutagenicity", 0.05)),
            float(admet.get("absorption", {}).get("oral_bioavailability_percent", 70) if isinstance(admet.get("absorption"), dict) else 70),
            float(admet.get("distribution", {}).get("plasma_protein_binding_percent", 50) if isinstance(admet.get("distribution"), dict) else 50),
            safety_margin,
            toxicity_index,
        ]
        while len(coord) < 35:
            coord.append(0.0)
        coordinate_35d = coord[:35]
        coherence = max(0.0, 1.0 - toxicity_index)

        payload = {
            "header": {
                "worldline_id": task_id,
                "timestamp": datetime.now(timezone.utc).isoformat(),
                "system_version": "1.2.0-GAIP",
                "route": route,
                "compliance": "GAIP-2026",
            },
            "coordinate_35d": coordinate_35d,
            "coherence_amplitude": coherence,
            "route": route,
            "physical_realization": {
                "status": risk_level,
                "toxicity_index": toxicity_index,
                "safety_margin": safety_margin,
                "outcome": outcome,
                "stage": stage_reached,
            },
            "physics_snapshot": {
                "coordinate_35d": coordinate_35d,
                "toxicity_index": toxicity_index,
                "safety_margin": safety_margin,
            },
        }
        if origin:
            payload["origin"] = origin
        if item:
            payload["candidate_data"] = {
                "prv_candidate": {"smiles": item.get("smiles", ""), "name": item.get("name", "") or str(item_id)},
                "task_id": task_id,
                "coherence": coherence,
                "toxicity": toxicity_index,
                "affinity": float(ope.get("ec50", 1.0)),
            }
        with open(filepath, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
        return filepath

    def record_discovery(
        self,
        item_id: str,
        dossier: Dict[str, Any],
        tier: str,
        is_novel: bool,
        repo_root: Path | str | None = None,
    ) -> str:
        """
        Create a worldline for every discovery. Writes to WorldLines/<tier>/ so the 35D logic
        remembers what work was done and can improve research accuracy over time.
        """
        from PX_Warehouse.warehouse_layout import get_worldline_dir, TIERS

        root = Path(repo_root) if repo_root else Path(_REPO_ROOT)
        tier = tier if tier in TIERS else "Bronze"
        wl_dir = get_worldline_dir(tier, root)
        wl_dir.mkdir(parents=True, exist_ok=True)

        safe_id = "".join(c if c.isalnum() or c in "-_" else "_" for c in str(item_id))[:64]
        task_id = f"WL-{'NOV' if is_novel else 'REP'}_{safe_id}"
        fname = f"{safe_id}.worldline"
        filepath = os.path.join(str(wl_dir), fname)

        engines = dossier.get("engines") or {}
        ope = engines.get("ope") or {}
        admet = engines.get("admet") or {}
        tox = (admet.get("toxicity") or {})
        toxicity_index = float(tox.get("toxicity_index", 0.02))
        risk_level = tox.get("risk_level", "TOXICITY_GOLD")
        safety_margins = admet.get("safety_margins") or {}
        safety_margin = float(safety_margins.get("safety_margin", 0.0))
        harm_energy = float(dossier.get("harm_energy", toxicity_index))

        # 35D coordinate from OPE + ADMET for manifold memory (efficacy, safety, geometry, etc.)
        coord = [
            float(ope.get("logp", 2.0)),
            float(ope.get("molecular_weight", 350.0)),
            float(ope.get("tpsa", 70.0)),
            float(ope.get("ec50", 1.0)),
            toxicity_index,
            float(tox.get("herg_risk", 0.2)),
            float(tox.get("hepatotoxicity_risk", 0.2)),
            float(tox.get("ames_mutagenicity", 0.05)),
            float(admet.get("absorption", {}).get("oral_bioavailability_percent", 70) if isinstance(admet.get("absorption"), dict) else 70),
            float(admet.get("distribution", {}).get("plasma_protein_binding_percent", 50) if isinstance(admet.get("distribution"), dict) else 50),
            safety_margin,
            harm_energy,
        ]
        # Pad to 35 dimensions
        while len(coord) < 35:
            coord.append(0.0)
        coordinate_35d = coord[:35]

        candidate = dossier.get("candidate") or {}
        smiles = candidate.get("smiles", "")
        name = candidate.get("name", "") or item_id
        coherence = max(0.0, 1.0 - toxicity_index)

        payload = {
            "header": {
                "worldline_id": task_id,
                "timestamp": datetime.now(timezone.utc).isoformat(),
                "system_version": "1.2.0-GAIP",
                "route": "DISCOVERY",
                "compliance": "GAIP-2026",
                "reality_type": "REAL",
            },
            "coordinate_35d": coordinate_35d,
            "coherence_amplitude": coherence,
            "route": "DISCOVERY",
            "physical_realization": {
                "toxicity_index": toxicity_index,
                "risk_level": risk_level,
                "safety_margin": safety_margin,
                "harm_energy": harm_energy,
                "source_dossier": "prv_24h",
            },
            "physics_snapshot": {
                "coordinate_35d": coordinate_35d,
                "toxicity_index": toxicity_index,
                "safety_margin": safety_margin,
            },
            "candidate_data": {
                "prv_candidate": {"smiles": smiles, "name": name},
                "task_id": task_id,
                "coherence": coherence,
                "toxicity": toxicity_index,
                "affinity": float(ope.get("ec50", 1.0)),
            },
        }
        with open(filepath, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
        return filepath

    def get_all_candidates(self):
        """
        Load all worldline files (root and tier subdirs) for mining/trajectory.
        Returns list of dicts with physics, name, smiles, worldline_path.
        """
        from PX_Warehouse.warehouse_layout import TIERS

        candidates = []
        root = Path(self._path)
        if not root.exists():
            return candidates
        for wl_path in root.rglob("*.worldline"):
            try:
                with open(wl_path, "r", encoding="utf-8") as f:
                    data = json.load(f)
            except Exception:
                continue
            header = data.get("header") or {}
            phys = data.get("physical_realization") or data.get("physics_snapshot") or {}
            cand = data.get("candidate_data") or {}
            prv = cand.get("prv_candidate") or {}
            candidates.append({
                "name": prv.get("name", header.get("worldline_id", "Unknown")),
                "smiles": prv.get("smiles", ""),
                "physics": {
                    "toxicity_index": phys.get("toxicity_index", 0.02),
                    "safety_margin": phys.get("safety_margin", 0.0),
                },
                "worldline_path": str(wl_path),
            })
        return candidates


__all__ = ["WorldLineDatabase", "DEFAULT_WORLDLINES_PATH"]
