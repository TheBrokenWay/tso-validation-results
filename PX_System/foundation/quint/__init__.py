"""
PX_System.foundation.quint — QUINT Internal Mathematical Language

QUINT is the platform's internal data language. All external data enters
through the Converter (single gateway in), all output exits through the
Converter (single gateway out). No other module touches raw external data.

Architecture:
    External World
        | Converter.ingest()   [one way in]
        v
    QFrame (typed, hashed, validated)
        | Runtime / Engines operate on QFrames
        v
    QFrame (transformed, lineage tracked)
        | Converter.emit()     [one way out]
        v
    External World

Modules:
    kernel     — QFrame, QType core types
    compiler   — Data normalization into QFrame payloads
    converter  — THE single gateway (ingest / emit)
    runtime    — QFrame operations, validation, integrity
    registry   — Schema registry for QTypes
    resonance  — Round-trip fidelity measurement

Constitutional: Python stdlib only. No external dependencies.
"""

from PX_System.foundation.quint.kernel import (
    QFrame,
    QType,
    SCHEMA_VERSION,
    create_qframe,
)

from PX_System.foundation.quint.converter import (
    ingest,
    emit,
    emit_json,
    ingest_file,
    ingest_batch,
    emit_batch,
    round_trip,
    ConversionError,
)

from PX_System.foundation.quint.runtime import (
    validate,
    transform,
    merge,
    extract,
    annotate,
    chain,
)

from PX_System.foundation.quint.resonance import (
    measure,
    measure_aggregate,
    run_standard_resonance_test,
    ResonanceResult,
)

__all__ = [
    # Kernel
    "QFrame", "QType", "SCHEMA_VERSION", "create_qframe",
    # Converter (THE gateway)
    "ingest", "emit", "emit_json", "ingest_file", "ingest_batch",
    "emit_batch", "round_trip", "ConversionError",
    # Runtime
    "validate", "transform", "merge", "extract", "annotate", "chain",
    # Resonance
    "measure", "measure_aggregate", "run_standard_resonance_test", "ResonanceResult",
]
