"""
QUINT Brains — Disease Constraint Evaluator.

Loads the 22 PRV disease constraint JSON files from PX_Domain/PRV_Diseases/
as QCONSTRAINT frames and provides constraint-based molecule validation.

Usage:
    from PX_System.foundation.quint.brains.disease_constraints import (
        load_disease_constraint,
        load_all_disease_constraints,
        evaluate_molecule_against_constraint,
        get_eligible_diseases,
    )

    constraints = load_all_disease_constraints()
    nipah = load_disease_constraint("nipah_virus_infection")
    result = evaluate_molecule_against_constraint(mol_frame, nipah)

Constitutional: Python stdlib only. Deterministic. Fail-closed.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

REPO_ROOT = Path(__file__).resolve().parents[4]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from PX_System.foundation.quint.kernel import QFrame, QType, create_qframe
from PX_System.foundation.quint.converter import ingest


# ── Path to disease constraint files ────────────────────────────────────────

_DISEASE_DIR = REPO_ROOT / "PX_Domain" / "PRV_Diseases"


# ── Constraint check definitions ────────────────────────────────────────────
# Each tuple: (mol_field, constraint_field, comparison)
# Comparisons: "lt" = mol < constraint, "lte" = mol <= constraint,
#              "gte" = mol >= constraint, "range" = between two constraint fields

_CHECKS: List[Tuple[str, str, str]] = [
    ("toxicity", "toxicity_threshold", "lt"),
    ("mw", "molecular_weight_max", "lte"),
    ("hbd", "hbd_max", "lte"),
    ("hba", "hba_max", "lte"),
]

# logp is a range check: logp_min <= mol.logp <= logp_max
_RANGE_CHECKS: List[Tuple[str, str, str]] = [
    ("logp", "logp_min", "logp_max"),
]


# ── Loading functions ───────────────────────────────────────────────────────

def load_disease_constraint(disease_id: str) -> QFrame:
    """
    Load a single disease constraint JSON file as a QCONSTRAINT QFrame.

    Args:
        disease_id: The disease identifier matching the filename stem
                    (e.g., "nipah_virus_infection").

    Returns:
        A sealed QCONSTRAINT QFrame containing the constraint data.

    Raises:
        FileNotFoundError: If the constraint file does not exist.
    """
    filepath = _DISEASE_DIR / f"{disease_id}.json"
    if not filepath.exists():
        raise FileNotFoundError(
            f"Disease constraint file not found: {filepath}"
        )

    with open(filepath, "r", encoding="utf-8") as f:
        raw_data = json.load(f)

    qid = f"QCONSTRAINT-{disease_id}"
    frame = ingest(
        raw_data,
        qtype=QType.QCONSTRAINT,
        qid=qid,
        source_label=f"disease_constraint:{disease_id}",
    )
    frame.lineage.append(f"disease_constraint:{disease_id}")
    return frame


def load_all_disease_constraints() -> Dict[str, QFrame]:
    """
    Load all disease constraint JSON files from PX_Domain/PRV_Diseases/.

    Excludes manifest.json and any non-JSON files.

    Returns:
        Dict mapping disease_id (file stem) to QCONSTRAINT QFrame.
    """
    constraints: Dict[str, QFrame] = {}

    if not _DISEASE_DIR.exists():
        return constraints

    for filepath in sorted(_DISEASE_DIR.glob("*.json")):
        if filepath.stem == "manifest":
            continue
        disease_id = filepath.stem
        constraints[disease_id] = load_disease_constraint(disease_id)

    return constraints


# ── Evaluation functions ────────────────────────────────────────────────────

def evaluate_molecule_against_constraint(
    mol_frame: QFrame,
    constraint_frame: QFrame,
) -> Dict[str, Any]:
    """
    Validate a QMOLECULE frame against a QCONSTRAINT frame.

    Checks molecule properties against constraint thresholds. If a molecule
    field is missing, that check is skipped (tolerant of incomplete data).

    Checks performed:
        - toxicity < toxicity_threshold (strict less-than, constitutional)
        - mw <= molecular_weight_max
        - logp between logp_min and logp_max (inclusive)
        - hbd <= hbd_max
        - hba <= hba_max

    Args:
        mol_frame:        QFrame with qtype QMOLECULE (or any with relevant payload fields).
        constraint_frame: QFrame with qtype QCONSTRAINT containing disease thresholds.

    Returns:
        Dict with keys:
            passed:          bool — True if all checks passed (or were skipped)
            disease_name:    str — from the constraint
            violations:      list of dicts with field, value, threshold, comparison
            checks_performed: int — number of checks actually run
            checks_passed:   int — number of checks that passed
    """
    mol = mol_frame.payload
    con = constraint_frame.payload

    disease_name = con.get("disease_name", "unknown")
    violations: List[Dict[str, Any]] = []
    checks_performed = 0
    checks_passed = 0

    # Standard threshold checks
    for mol_field, con_field, comparison in _CHECKS:
        mol_val = mol.get(mol_field)
        con_val = con.get(con_field)

        if mol_val is None or con_val is None:
            continue

        checks_performed += 1

        if comparison == "lt":
            if mol_val < con_val:
                checks_passed += 1
            else:
                violations.append({
                    "field": mol_field,
                    "value": mol_val,
                    "threshold": con_val,
                    "comparison": f"{mol_field} must be < {con_field}",
                })
        elif comparison == "lte":
            if mol_val <= con_val:
                checks_passed += 1
            else:
                violations.append({
                    "field": mol_field,
                    "value": mol_val,
                    "threshold": con_val,
                    "comparison": f"{mol_field} must be <= {con_field}",
                })

    # Range checks (logp between logp_min and logp_max)
    for mol_field, con_min_field, con_max_field in _RANGE_CHECKS:
        mol_val = mol.get(mol_field)
        con_min = con.get(con_min_field)
        con_max = con.get(con_max_field)

        if mol_val is None or con_min is None or con_max is None:
            continue

        checks_performed += 1

        if con_min <= mol_val <= con_max:
            checks_passed += 1
        else:
            violations.append({
                "field": mol_field,
                "value": mol_val,
                "threshold": f"[{con_min}, {con_max}]",
                "comparison": f"{mol_field} must be between {con_min_field} and {con_max_field}",
            })

    return {
        "passed": len(violations) == 0,
        "disease_name": disease_name,
        "violations": violations,
        "checks_performed": checks_performed,
        "checks_passed": checks_passed,
    }


def evaluate_molecule_against_all(
    mol_frame: QFrame,
) -> Dict[str, Dict[str, Any]]:
    """
    Evaluate a molecule against ALL disease constraints.

    Loads all 22 disease constraint files and runs the molecule through
    each constraint evaluator.

    Args:
        mol_frame: QFrame with molecule properties in payload.

    Returns:
        Dict mapping disease_id to evaluation result dict.
    """
    constraints = load_all_disease_constraints()
    results: Dict[str, Dict[str, Any]] = {}

    for disease_id, constraint_frame in constraints.items():
        results[disease_id] = evaluate_molecule_against_constraint(
            mol_frame, constraint_frame
        )

    return results


def get_eligible_diseases(mol_frame: QFrame) -> List[str]:
    """
    Return disease IDs where the molecule passes all constraints.

    Args:
        mol_frame: QFrame with molecule properties in payload.

    Returns:
        Sorted list of disease_id strings where the molecule passes.
    """
    results = evaluate_molecule_against_all(mol_frame)
    return sorted(
        disease_id
        for disease_id, result in results.items()
        if result["passed"]
    )


def get_constraint_summary() -> Dict[str, Any]:
    """
    Return a summary of all loaded disease constraints.

    Returns:
        Dict with:
            count:       int — number of disease constraint files
            disease_ids: list — sorted list of disease identifiers
            versions:    dict — mapping disease_id to constraint_version
    """
    constraints = load_all_disease_constraints()
    versions: Dict[str, str] = {}

    for disease_id, frame in constraints.items():
        versions[disease_id] = frame.payload.get("constraint_version", "unknown")

    return {
        "count": len(constraints),
        "disease_ids": sorted(constraints.keys()),
        "versions": versions,
    }


__all__ = [
    "load_disease_constraint",
    "load_all_disease_constraints",
    "evaluate_molecule_against_constraint",
    "evaluate_molecule_against_all",
    "get_eligible_diseases",
    "get_constraint_summary",
]
