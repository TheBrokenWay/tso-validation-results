"""
QUINT Resonance — Round-trip fidelity measurement.

Resonance measures how faithfully data survives the journey:
    External → Converter.ingest() → QFrame → Converter.emit() → External

A resonance of 100% means every field and value is perfectly preserved.
A resonance of 99% means 99% of fields survived unchanged.

The platform requirement is >= 99% resonance on all standard data types.

Constitutional: Python stdlib only. Deterministic.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Any, Dict, List, Optional

from PX_System.foundation.quint.converter import ingest, emit, reset_stats
from PX_System.foundation.quint.kernel import QType


@dataclass
class ResonanceResult:
    """Result of a resonance measurement."""
    total_fields: int
    preserved_fields: int
    lost_fields: List[str]
    altered_fields: List[str]
    added_fields: List[str]
    resonance_pct: float
    round_trip_us: float  # microseconds for the round-trip

    @property
    def passed(self) -> bool:
        return self.resonance_pct >= 99.0

    def __repr__(self) -> str:
        return (
            f"Resonance({self.resonance_pct:.2f}% | "
            f"{self.preserved_fields}/{self.total_fields} fields | "
            f"{self.round_trip_us:.1f}us)"
        )


def measure(
    data: Dict[str, Any],
    qtype: Optional[QType] = None,
) -> ResonanceResult:
    """
    Measure round-trip resonance for a single data dict.

    Performs: ingest(data) → QFrame → emit(frame) → compare with original.
    """
    t0 = time.perf_counter()
    frame = ingest(data, qtype=qtype, source_label="resonance_in")
    result = emit(frame, target_label="resonance_out")
    elapsed_us = (time.perf_counter() - t0) * 1_000_000

    # Compare original keys vs result keys
    orig_keys = set(data.keys())
    result_keys = set(result.keys())

    preserved = []
    altered = []
    lost = []
    added = []

    for key in orig_keys:
        if key in result_keys:
            if result[key] == data[key]:
                preserved.append(key)
            else:
                altered.append(key)
        else:
            lost.append(key)

    for key in result_keys:
        if key not in orig_keys:
            added.append(key)

    total = len(orig_keys)
    pct = (len(preserved) / total * 100) if total > 0 else 100.0

    return ResonanceResult(
        total_fields=total,
        preserved_fields=len(preserved),
        lost_fields=sorted(lost),
        altered_fields=sorted(altered),
        added_fields=sorted(added),
        resonance_pct=pct,
        round_trip_us=elapsed_us,
    )


def measure_batch(
    items: List[Dict[str, Any]],
    qtype: Optional[QType] = None,
) -> List[ResonanceResult]:
    """Measure resonance for a batch of data dicts."""
    return [measure(item, qtype) for item in items]


def measure_aggregate(
    items: List[Dict[str, Any]],
    qtype: Optional[QType] = None,
) -> Dict[str, Any]:
    """
    Measure resonance across a batch and return aggregate statistics.

    Returns:
        count: number of items measured
        min_resonance: worst case
        max_resonance: best case
        mean_resonance: average
        all_passed: True if every item >= 99%
        mean_round_trip_us: average round-trip time in microseconds
        total_fields: total fields across all items
        total_preserved: total preserved fields across all items
    """
    results = measure_batch(items, qtype)
    if not results:
        return {
            "count": 0,
            "min_resonance": 0.0,
            "max_resonance": 0.0,
            "mean_resonance": 0.0,
            "all_passed": True,
            "mean_round_trip_us": 0.0,
            "total_fields": 0,
            "total_preserved": 0,
        }

    resonances = [r.resonance_pct for r in results]
    times = [r.round_trip_us for r in results]
    total_f = sum(r.total_fields for r in results)
    total_p = sum(r.preserved_fields for r in results)

    return {
        "count": len(results),
        "min_resonance": min(resonances),
        "max_resonance": max(resonances),
        "mean_resonance": sum(resonances) / len(resonances),
        "all_passed": all(r.passed for r in results),
        "mean_round_trip_us": sum(times) / len(times),
        "total_fields": total_f,
        "total_preserved": total_p,
    }


# ── Standard test vectors ────────────────────────────────────────────────────

STANDARD_VECTORS: List[Dict[str, Any]] = [
    # QMOLECULE
    {"smiles": "CCO", "mw": 46.07, "logp": -0.31, "hbd": 1, "hba": 1, "tpsa": 20.23},
    {"smiles": "CC(=O)O", "mw": 60.05, "logp": -0.17, "toxicity": 0.005},
    {"smiles": "c1ccccc1", "mw": 78.11, "logp": 1.56, "hba": 0, "hbd": 0},
    # QRESULT
    {"engine_id": "OPE_V3_DETERMINISTIC", "status": "PASSED", "score": 0.95, "authorized": True},
    {"engine_id": "ADMET", "status": "PASSED", "score": 0.88, "toxicity": 0.019},
    # QDOSSIER
    {"compound_id": "PRV_NOV_abc123", "smiles": "CCO", "tier": "GOLD", "grade": "GOLD_TIER"},
    # QCONSTRAINT
    {"disease_name": "Nipah Virus Infection", "pathogen_type": "virus", "ic50_max_um": 1.0},
    # QPHYSICS
    {"dimensions": 35, "energy": 0.42, "global_sum": 1.0},
    # QRAW
    {"arbitrary_key": "arbitrary_value", "count": 42, "nested": {"a": 1, "b": 2}},
    # Edge cases
    {"smiles": "", "mw": 0.0},  # empty smiles
    {"smiles": "C" * 1000, "mw": 14000.0},  # very long SMILES
]


def run_standard_resonance_test() -> Dict[str, Any]:
    """
    Run the standard resonance test suite.

    Returns aggregate results. All standard vectors must achieve >= 99% resonance.
    """
    return measure_aggregate(STANDARD_VECTORS)


__all__ = [
    "ResonanceResult",
    "STANDARD_VECTORS",
    "measure",
    "measure_aggregate",
    "measure_batch",
    "run_standard_resonance_test",
]
