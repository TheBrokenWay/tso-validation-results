"""
Disease Registry — single entry point for all PRV disease data.

Thin adapter wrapping the QUINT disease_constraints.py evaluator (already
tested with 17 tests) and the PRV_Diseases manifest.json.

Usage:
    from PX_Domain.PRV_Diseases.disease_registry import (
        get_prv_diseases,
        get_disease_constraint,
        evaluate_candidate,
        get_all_constraints,
    )

    diseases = get_prv_diseases()
    result = evaluate_candidate({"mw": 350, "logp": 2.5, "hbd": 2, "hba": 5}, "nipah_virus_infection")

Constitutional: Python stdlib + QUINT only. Deterministic. Fail-closed.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from PX_System.foundation.quint.brains.disease_constraints import (
    load_disease_constraint as _quint_load,
    load_all_disease_constraints as _quint_load_all,
    evaluate_molecule_against_constraint as _quint_evaluate,
    get_eligible_diseases as _quint_eligible,
    get_constraint_summary as _quint_summary,
)
from PX_System.foundation.quint.kernel import QType
from PX_System.foundation.quint.converter import ingest as _quint_ingest


# ── Path to manifest ─────────────────────────────────────────────────────────

_DISEASE_DIR = _REPO_ROOT / "PX_Domain" / "PRV_Diseases"
_MANIFEST_PATH = _DISEASE_DIR / "manifest.json"


# ── Manifest access ──────────────────────────────────────────────────────────

def get_prv_diseases() -> List[Dict[str, Any]]:
    """
    Load the PRV disease manifest and return all disease entries.

    Returns:
        List of dicts, each with: name, category, pathogen_type, prv_eligible.
        Empty list if manifest not found.
    """
    if not _MANIFEST_PATH.exists():
        return []
    with open(_MANIFEST_PATH, "r", encoding="utf-8") as f:
        manifest = json.load(f)
    return manifest.get("prv_diseases", [])


def get_prv_disease_names() -> List[str]:
    """Return sorted list of all PRV disease names from manifest."""
    return sorted(d["name"] for d in get_prv_diseases() if "name" in d)


def get_prv_disease_ids() -> List[str]:
    """
    Return disease_ids matching constraint file stems.

    Derives ID from disease name: lowercase, spaces → underscores.
    """
    return sorted(
        d["name"].lower().replace(" ", "_")
        for d in get_prv_diseases()
        if "name" in d
    )


# ── Constraint access (delegates to QUINT) ───────────────────────────────────

def get_disease_constraint(disease_id: str):
    """
    Load a single disease constraint as a QCONSTRAINT QFrame.

    Args:
        disease_id: File stem, e.g. "nipah_virus_infection".

    Returns:
        QCONSTRAINT QFrame, or None if file not found.
    """
    try:
        return _quint_load(disease_id)
    except FileNotFoundError:
        return None


def get_all_constraints():
    """
    Load all disease constraint files as QCONSTRAINT QFrames.

    Returns:
        Dict mapping disease_id to QFrame.
    """
    return _quint_load_all()


def get_constraint_summary():
    """Summary of all loaded constraints (count, ids, versions)."""
    return _quint_summary()


# ── Candidate evaluation (the key pipeline integration point) ────────────────

def evaluate_candidate(
    mol_properties: Dict[str, Any],
    disease_id: str,
) -> Optional[Dict[str, Any]]:
    """
    Evaluate a molecule's properties against a disease constraint.

    This is the main integration point for px_prv.py. It accepts plain dicts
    (from OPE results) and handles the QUINT QFrame conversion internally.

    Args:
        mol_properties: Dict with keys like mw, logp, hbd, hba, toxicity.
                        Missing keys are tolerated (checks are skipped).
        disease_id:     File stem, e.g. "nipah_virus_infection".

    Returns:
        Dict with: passed (bool), disease_name, violations (list),
        checks_performed, checks_passed.
        None if the disease constraint file does not exist.
    """
    constraint_frame = get_disease_constraint(disease_id)
    if constraint_frame is None:
        return None

    # Convert plain dict to QMOLECULE QFrame via QUINT gateway
    # QMOLECULE requires 'smiles' field — inject placeholder if not present
    props = dict(mol_properties)
    if "smiles" not in props:
        props["smiles"] = props.get("SMILES", "")
    mol_frame = _quint_ingest(
        props,
        qtype=QType.QMOLECULE,
        qid=f"candidate-eval-{disease_id}",
        source_label=f"registry:evaluate:{disease_id}",
    )

    return _quint_evaluate(mol_frame, constraint_frame)


def evaluate_candidate_all(
    mol_properties: Dict[str, Any],
) -> Dict[str, Dict[str, Any]]:
    """
    Evaluate a molecule against ALL disease constraints.

    Returns:
        Dict mapping disease_id to evaluation result.
    """
    constraints = get_all_constraints()
    if not constraints:
        return {}

    props = dict(mol_properties)
    if "smiles" not in props:
        props["smiles"] = props.get("SMILES", "")
    mol_frame = _quint_ingest(
        props,
        qtype=QType.QMOLECULE,
        qid="candidate-eval-all",
        source_label="registry:evaluate:all",
    )

    results: Dict[str, Dict[str, Any]] = {}
    for disease_id, constraint_frame in constraints.items():
        results[disease_id] = _quint_evaluate(mol_frame, constraint_frame)
    return results


def get_eligible_diseases_for(mol_properties: Dict[str, Any]) -> List[str]:
    """
    Return disease IDs where the molecule passes all constraints.
    """
    results = evaluate_candidate_all(mol_properties)
    return sorted(
        did for did, r in results.items() if r.get("passed", False)
    )


__all__ = [
    "get_prv_diseases",
    "get_prv_disease_names",
    "get_prv_disease_ids",
    "get_disease_constraint",
    "get_all_constraints",
    "get_constraint_summary",
    "evaluate_candidate",
    "evaluate_candidate_all",
    "get_eligible_diseases_for",
]
