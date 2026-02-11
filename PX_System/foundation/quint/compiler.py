"""
QUINT Compiler — Normalizes external data into QFrame payload format.

The compiler is an INTERNAL module used by the Converter gateway.
It does NOT touch external I/O directly — it only transforms data structures.

Responsibilities:
  - Detect the QType from raw data shape
  - Normalize field names to QUINT canonical form
  - Strip unknown/dangerous fields
  - Validate required fields per QType
  - Produce a clean payload dict ready for QFrame wrapping

Constitutional: Python stdlib only. Deterministic. No side effects.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

from PX_System.foundation.quint.kernel import QType


# ── Required fields per QType ────────────────────────────────────────────────

REQUIRED_FIELDS: Dict[QType, List[str]] = {
    QType.QMOLECULE: ["smiles"],
    QType.QRESULT: ["engine_id", "status"],
    QType.QDOSSIER: ["compound_id"],
    QType.QPHYSICS: ["dimensions"],
    QType.QSIGNAL: ["signal_type", "source"],
    QType.QTRIAL: ["protocol"],
    QType.QCONSTRAINT: ["disease_name"],
    QType.QRAW: [],  # no requirements
}

# ── Field aliases → canonical name ───────────────────────────────────────────

_ALIASES: Dict[str, str] = {
    "SMILES": "smiles",
    "Smiles": "smiles",
    "molecule": "smiles",
    "mol": "smiles",
    "molecular_weight": "mw",
    "MW": "mw",
    "mol_weight": "mw",
    "compound": "compound_id",
    "compound_name": "compound_id",
    "name": "compound_id",
    "disease": "disease_name",
    "indication": "disease_name",
    "engine": "engine_id",
    "engine_name": "engine_id",
    "dim": "dimensions",
    "dims": "dimensions",
    "tox": "toxicity",
    "toxicity_index": "toxicity",
    "harm": "harm_energy",
    "harm_energy_score": "harm_energy",
}

# ── Blocked fields (never pass through) ──────────────────────────────────────

_BLOCKED: frozenset = frozenset({
    "password", "secret", "token", "api_key", "credential",
    "private_key", "access_token", "refresh_token",
})


class CompilationError(Exception):
    """Raised when raw data cannot be compiled into a valid QFrame payload."""
    pass


def detect_qtype(data: Dict[str, Any]) -> QType:
    """
    Infer QType from the shape of raw data.

    Returns QType.QRAW if no specific type can be determined.
    """
    keys = set(data.keys())

    # Check for explicit type hint
    explicit = data.get("_qtype") or data.get("qtype")
    if explicit:
        try:
            return QType(explicit)
        except ValueError:
            pass

    # Heuristic detection by field presence
    if "smiles" in keys or "SMILES" in keys or "molecule" in keys:
        return QType.QMOLECULE
    if "engine_id" in keys or "engine" in keys or "engine_name" in keys:
        return QType.QRESULT
    if "protocol" in keys and ("arm" in str(data.get("protocol", "")) or "trial" in str(keys)):
        return QType.QTRIAL
    if "dimensions" in keys or "physics_vector" in keys:
        return QType.QPHYSICS
    if "signal_type" in keys:
        return QType.QSIGNAL
    if "disease_name" in keys or "disease" in keys or "indication" in keys:
        return QType.QCONSTRAINT
    if "compound_id" in keys or "compound" in keys:
        return QType.QDOSSIER

    return QType.QRAW


def normalize(data: Dict[str, Any], qtype: Optional[QType] = None) -> Tuple[QType, Dict[str, Any]]:
    """
    Normalize raw data into a canonical QUINT payload.

    Steps:
      1. Strip blocked fields
      2. Apply field aliases
      3. Detect or confirm QType
      4. Validate required fields
      5. Return (qtype, clean_payload)

    Raises CompilationError if required fields are missing.
    """
    # Step 1: Strip blocked fields
    clean = {k: v for k, v in data.items() if k.lower() not in _BLOCKED}

    # Step 2: Apply aliases
    normalized: Dict[str, Any] = {}
    for key, value in clean.items():
        canonical = _ALIASES.get(key, key)
        # Don't overwrite if canonical already set by a prior key
        if canonical not in normalized:
            normalized[canonical] = value

    # Step 3: Detect or confirm QType
    resolved_type = qtype if qtype is not None else detect_qtype(normalized)

    # Step 4: Validate required fields
    required = REQUIRED_FIELDS.get(resolved_type, [])
    missing = [f for f in required if f not in normalized]
    if missing:
        raise CompilationError(
            f"QType {resolved_type.value} requires fields {missing} — "
            f"got keys: {sorted(normalized.keys())}"
        )

    return resolved_type, normalized


def compile_batch(items: List[Dict[str, Any]], qtype: Optional[QType] = None) -> List[Tuple[QType, Dict[str, Any]]]:
    """Normalize a list of raw dicts. Returns list of (qtype, payload) tuples."""
    return [normalize(item, qtype) for item in items]


__all__ = [
    "CompilationError",
    "REQUIRED_FIELDS",
    "compile_batch",
    "detect_qtype",
    "normalize",
]
