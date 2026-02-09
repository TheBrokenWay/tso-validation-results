"""
Warehouse Placement Gate — single entry point for all file placement.

EVERYTHING must go through this module before being sorted into the warehouse.
No code should write directly to PX_Warehouse paths; use the place_* functions below.
This prevents misplaced files and keeps Feeder / Calibration_Molecules / Dossiers / Learning_Material correct.

Canonical zones (only these receive content):
  - Feeder                 — queue files (prv_24h_queue.json, reprocess_candidates.json, etc.)
  - Calibration_Molecules  — TRIAL_SIMULATION dossiers, LiveRuns, BatchRuns (reference/QA)
  - Prv_Dossiers/<tier>    — PRV_REP dossiers by tier (Diamond, Gold, Silver, Bronze)
  - Novel_Dossiers/<tier>  — PRV_NOV dossiers by tier
  - Learning_Material      — general logs (non-calibration)

Monitoring: run audit_placement() periodically to detect files outside these zones.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Tuple

from PX_Warehouse.warehouse_layout import (
    WAREHOUSE_ROOT,
    FEEDER_ROOT,
    CALIBRATION_MOLECULES_ROOT,
    PRV_DOSSIERS_ROOT,
    NOVEL_DOSSIERS_ROOT,
    LEARNING_MATERIAL_ROOT,
    TIERS,
    ensure_structure,
    get_tier,
    get_feeder_dir,
    get_calibration_molecules_dir,
    get_prv_dossier_dir,
    get_learning_material_dir,
    get_queue_path,
)


# --- Sanctioned zones (paths that are valid placement targets) ---
def _canonical_zone_roots(repo_root: Path | None = None) -> List[Path]:
    """Roots of all zones where files are allowed. Used by audit."""
    root = repo_root or WAREHOUSE_ROOT.parent
    wh = root / "PX_Warehouse"
    zones = [
        wh / "Feeder",
        wh / "Calibration_Molecules",
        wh / "Learning_Material",
        wh / "Prv_Dossiers",
        wh / "Novel_Dossiers",
        wh / "Finalized_Dossiers",
    ]
    for t in TIERS:
        zones.append(wh / "Prv_Dossiers" / t)
        zones.append(wh / "Novel_Dossiers" / t)
        zones.append(wh / "Finalized_Dossiers" / t)
    zones.append(wh / "Finalized_Dossiers" / "Discovery_Accepted")
    return zones


def is_path_sanctioned(path: Path, repo_root: Path | None = None) -> bool:
    """True if path is under a canonical zone (Feeder, Calibration_Molecules, Prv/Novel/<tier>, Learning_Material)."""
    path = Path(path).resolve()
    for zone in _canonical_zone_roots(repo_root):
        zone = zone.resolve()
        if path == zone or path.is_relative_to(zone):
            return True
    return False


def resolve_zone(path: Path, repo_root: Path | None = None) -> str | None:
    """Return zone name for path ('Feeder', 'Calibration_Molecules', 'Prv_Dossiers', 'Novel_Dossiers', 'Learning_Material') or None if not sanctioned."""
    path = Path(path).resolve()
    root = repo_root or WAREHOUSE_ROOT.parent
    wh = root / "PX_Warehouse"
    if path == (wh / "Feeder") or (wh / "Feeder") in path.parents:
        return "Feeder"
    if path == (wh / "Calibration_Molecules") or (wh / "Calibration_Molecules") in path.parents:
        return "Calibration_Molecules"
    if path == (wh / "Learning_Material") or (wh / "Learning_Material") in path.parents:
        return "Learning_Material"
    if (wh / "Prv_Dossiers") in path.parents:
        return "Prv_Dossiers"
    if (wh / "Novel_Dossiers") in path.parents:
        return "Novel_Dossiers"
    if (wh / "Finalized_Dossiers" / "Discovery_Accepted") in path.parents or path.parent == (wh / "Finalized_Dossiers" / "Discovery_Accepted"):
        return "Discovery_Accepted"
    if (wh / "Finalized_Dossiers") in path.parents:
        return "Finalized_Dossiers"
    return None


# --- Placement API (use these instead of writing to PX_Warehouse directly) ---

def place_queue_file(
    filename: str,
    payload: Dict[str, Any],
    repo_root: Path | None = None,
    *,
    atomic: bool = True,
) -> Path:
    """
    Write a queue file into Feeder. Use for prv_24h_queue.json, reprocess_candidates.json, e2e_test_queue.json, etc.
    Returns the path written.
    """
    ensure_structure(repo_root)
    root = repo_root or WAREHOUSE_ROOT.parent
    feeder = get_feeder_dir(root)
    feeder.mkdir(parents=True, exist_ok=True)
    out = feeder / filename
    if atomic and out.suffix.lower() == ".json":
        tmp = out.with_suffix(out.suffix + ".tmp")
        tmp.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
        tmp.replace(out)
    else:
        out.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    return out


def place_calibration_dossier(
    filename: str,
    payload: Dict[str, Any],
    subdir: str | None = None,
    repo_root: Path | None = None,
) -> Path:
    """
    Write a file into Calibration_Molecules (optionally under subdir e.g. LiveRuns/run_ts or BatchRuns).
    Use for TRIAL_SIMULATION_DOSSIER-*.json and other reference/QA outputs.
    """
    ensure_structure(repo_root)
    root = repo_root or WAREHOUSE_ROOT.parent
    cal = get_calibration_molecules_dir(root)
    if subdir:
        cal = cal / subdir
    cal.mkdir(parents=True, exist_ok=True)
    out = cal / filename
    out.write_text(json.dumps(payload, indent=2, ensure_ascii=False, default=str), encoding="utf-8")
    return out


def place_prv_dossier(
    dossier: Dict[str, Any],
    item: Dict[str, Any],
    repo_root: Path | None = None,
) -> Path:
    """
    Tier and place a PRV dossier into Prv_Dossiers/<tier> or Novel_Dossiers/<tier>.
    item must have 'type' ('R' or 'N') and 'id'. Returns path written.
    """
    ensure_structure(repo_root)
    tier = get_tier(dossier)
    is_novel = item.get("type") == "N"
    out_dir = get_prv_dossier_dir(is_novel, tier, repo_root)
    out_dir.mkdir(parents=True, exist_ok=True)
    safe_id = "".join(c if c.isalnum() else "_" for c in item.get("id", "unknown"))
    prefix = "PRV_NOV" if is_novel else "PRV_REP"
    filename = f"{prefix}_{safe_id}.json"
    out = out_dir / filename
    out.write_text(json.dumps(dossier, indent=2), encoding="utf-8")
    return out


def place_learning_material_file(
    filename: str,
    content: str | bytes,
    repo_root: Path | None = None,
    *,
    binary: bool = False,
) -> Path:
    """Write a file into Learning_Material (general logs). content is str or bytes."""
    ensure_structure(repo_root)
    root = repo_root or WAREHOUSE_ROOT.parent
    lm = get_learning_material_dir(root)
    lm.mkdir(parents=True, exist_ok=True)
    out = lm / filename
    if binary:
        out.write_bytes(content if isinstance(content, bytes) else content.encode("utf-8"))
    else:
        out.write_text(content if isinstance(content, str) else content.decode("utf-8"), encoding="utf-8")
    return out


def get_queue_path_for_reading(filename: str, repo_root: Path | None = None) -> Path:
    """Resolve queue file for reading (Feeder first, then Operations/Inputs). Use for load_queue logic."""
    return get_queue_path(filename, repo_root)


# --- Monitoring: audit placement so nothing is written outside sanctioned zones ---

def audit_placement(
    repo_root: Path | None = None,
    *,
    exclude_dirs: Tuple[str, ...] = ("__pycache__", ".git", "Operations", "WorldLines"),
) -> Tuple[List[Path], List[Path], Dict[str, int]]:
    """
    Scan PX_Warehouse and report:
      - misplaced: files/dirs not under any canonical zone (should be moved or removed)
      - sanctioned: sample of paths under canonical zones (for sanity)
      - counts: { "Feeder": n, "Calibration_Molecules": n, ... } by zone

    exclude_dirs: top-level dir names to skip (e.g. Operations, WorldLines are special-use).
    Returns (misplaced_list, sanctioned_sample, zone_counts).
    """
    root = repo_root or WAREHOUSE_ROOT.parent
    wh = root / "PX_Warehouse"
    if not wh.exists():
        return [], [], {}

    zone_roots = {z.resolve(): z for z in _canonical_zone_roots(repo_root)}
    misplaced: List[Path] = []
    sanctioned_sample: List[Path] = []
    counts: Dict[str, int] = {}

    def _count_zone(p: Path) -> None:
        z = resolve_zone(p, repo_root)
        if z:
            counts[z] = counts.get(z, 0) + 1

    for item in wh.rglob("*"):
        if not item.exists():
            continue
        rel = item.relative_to(wh) if wh in item.parents or item == wh else None
        if rel and len(rel.parts) >= 1 and rel.parts[0] in exclude_dirs:
            continue
        # Warehouse root infra (layout, gate, README, .py at root) — not "placed" content
        if rel and len(rel.parts) == 1:
            continue
        if item.is_dir():
            continue  # count files only
        if not is_path_sanctioned(item, repo_root):
            misplaced.append(item)
        else:
            _count_zone(item)
            if len(sanctioned_sample) < 50:  # cap sample
                sanctioned_sample.append(item)

    return misplaced, sanctioned_sample, counts


def run_placement_audit(repo_root: Path | None = None, verbose: bool = True) -> int:
    """
    Run audit_placement and print a short report. Returns number of misplaced items (0 = OK).
    Call this from CI or a cron job to ensure no drift.
    """
    misplaced, _sample, counts = audit_placement(repo_root)
    if verbose:
        print("Warehouse placement audit")
        print("  Zone counts:", counts)
        if misplaced:
            print(f"  MISPLACED ({len(misplaced)}):")
            for p in misplaced[:30]:
                print(f"    - {p}")
            if len(misplaced) > 30:
                print(f"    ... and {len(misplaced) - 30} more")
        else:
            print("  All scanned files are in sanctioned zones.")
    return len(misplaced)


__all__ = [
    "is_path_sanctioned",
    "resolve_zone",
    "place_queue_file",
    "place_calibration_dossier",
    "place_prv_dossier",
    "place_learning_material_file",
    "get_queue_path_for_reading",
    "audit_placement",
    "run_placement_audit",
]
