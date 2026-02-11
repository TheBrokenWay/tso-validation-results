"""
QUINT Converter — THE Single Gateway.

ONE WAY IN:  ingest()  — external data → QFrame
ONE WAY OUT: emit()    — QFrame → external data

No other module in the platform should parse, deserialize, or construct
data from external sources. Everything enters through ingest(), everything
exits through emit(). The converter is the membrane between the outside
world and the QUINT internal representation.

Constitutional: Python stdlib only. Deterministic. Fail-closed.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from PX_System.foundation.quint.compiler import CompilationError, normalize
from PX_System.foundation.quint.kernel import (
    QFrame,
    QType,
    _canonical_json,
    _hash_payload,
    create_qframe,
)


class ConversionError(Exception):
    """Raised when external data cannot be converted to/from QUINT."""
    pass


# ── Counters for diagnostics ─────────────────────────────────────────────────

_stats = {
    "ingested": 0,
    "emitted": 0,
    "rejected": 0,
}


def get_stats() -> Dict[str, int]:
    """Return converter throughput counters."""
    return dict(_stats)


def reset_stats() -> None:
    """Reset converter counters (for testing)."""
    _stats["ingested"] = 0
    _stats["emitted"] = 0
    _stats["rejected"] = 0


# ── ONE WAY IN ───────────────────────────────────────────────────────────────

def ingest(
    data: Union[Dict[str, Any], str, bytes],
    qtype: Optional[QType] = None,
    qid: Optional[str] = None,
    source_label: str = "external",
) -> QFrame:
    """
    THE single entry point for external data into QUINT.

    Accepts:
      - dict: used directly
      - str:  parsed as JSON
      - bytes: decoded UTF-8, then parsed as JSON

    Returns a sealed QFrame with:
      - source_hash: SHA-256 of the raw input
      - frame_hash:  integrity seal of the normalized payload
      - lineage:     ["ingest:<source_label>"]

    Raises ConversionError on invalid input.
    """
    # Step 1: Parse raw input to dict
    raw_bytes: bytes
    if isinstance(data, bytes):
        raw_bytes = data
        try:
            parsed = json.loads(data.decode("utf-8"))
        except (UnicodeDecodeError, json.JSONDecodeError) as e:
            _stats["rejected"] += 1
            raise ConversionError(f"Cannot parse bytes as JSON: {e}") from e
    elif isinstance(data, str):
        raw_bytes = data.encode("utf-8")
        try:
            parsed = json.loads(data)
        except json.JSONDecodeError as e:
            _stats["rejected"] += 1
            raise ConversionError(f"Cannot parse string as JSON: {e}") from e
    elif isinstance(data, dict):
        raw_bytes = _canonical_json(data).encode("utf-8")
        parsed = data
    else:
        _stats["rejected"] += 1
        raise ConversionError(f"Unsupported input type: {type(data).__name__}")

    if not isinstance(parsed, dict):
        _stats["rejected"] += 1
        raise ConversionError(f"Expected JSON object, got {type(parsed).__name__}")

    # Step 2: Compute source hash (fingerprint of raw external data)
    source_hash = hashlib.sha256(raw_bytes).hexdigest()

    # Step 3: Compile (normalize + type-detect + validate)
    try:
        resolved_type, payload = normalize(parsed, qtype)
    except CompilationError as e:
        _stats["rejected"] += 1
        raise ConversionError(str(e)) from e

    # Step 4: Generate QID if not provided
    if qid is None:
        qid = f"Q-{source_hash[:12]}"

    # Step 5: Create sealed QFrame
    now = datetime.now(timezone.utc).isoformat()
    frame = create_qframe(
        qtype=resolved_type,
        qid=qid,
        payload=payload,
        source_hash=source_hash,
        ingested_at=now,
    )
    frame.lineage.append(f"ingest:{source_label}")

    _stats["ingested"] += 1
    return frame


def ingest_file(
    path: Union[str, Path],
    qtype: Optional[QType] = None,
    qid: Optional[str] = None,
) -> QFrame:
    """
    Convenience: ingest a JSON file through the gateway.

    The file path is recorded in lineage.
    """
    p = Path(path)
    if not p.exists():
        raise ConversionError(f"File not found: {p}")
    raw = p.read_bytes()
    return ingest(raw, qtype=qtype, qid=qid, source_label=f"file:{p.name}")


def ingest_batch(
    items: List[Union[Dict[str, Any], str, bytes]],
    qtype: Optional[QType] = None,
    source_label: str = "external_batch",
) -> List[QFrame]:
    """Ingest a list of items. Each gets its own QFrame."""
    return [
        ingest(item, qtype=qtype, source_label=f"{source_label}[{i}]")
        for i, item in enumerate(items)
    ]


# ── ONE WAY OUT ──────────────────────────────────────────────────────────────

def emit(
    frame: QFrame,
    include_meta: bool = False,
    target_label: str = "external",
) -> Dict[str, Any]:
    """
    THE single exit point for QUINT data to the external world.

    By default, emits only the payload (clean external representation).
    With include_meta=True, includes QUINT envelope metadata prefixed with _.

    Updates frame lineage with the emit operation.
    """
    if not isinstance(frame, QFrame):
        raise ConversionError(f"Can only emit QFrame, got {type(frame).__name__}")

    frame.lineage.append(f"emit:{target_label}")

    if include_meta:
        result = frame.to_dict()
    else:
        result = dict(frame.payload)

    _stats["emitted"] += 1
    return result


def emit_json(
    frame: QFrame,
    include_meta: bool = False,
    target_label: str = "external",
) -> str:
    """Emit as a JSON string."""
    d = emit(frame, include_meta=include_meta, target_label=target_label)
    return json.dumps(d, indent=2, default=str)


def emit_batch(
    frames: List[QFrame],
    include_meta: bool = False,
    target_label: str = "external_batch",
) -> List[Dict[str, Any]]:
    """Emit a list of QFrames."""
    return [
        emit(f, include_meta=include_meta, target_label=f"{target_label}[{i}]")
        for i, f in enumerate(frames)
    ]


# ── ROUND-TRIP VERIFICATION ─────────────────────────────────────────────────

def round_trip(data: Dict[str, Any], qtype: Optional[QType] = None) -> Dict[str, Any]:
    """
    Full round-trip: ingest external dict → QFrame → emit back to dict.

    Used for resonance testing: compare input vs output.
    """
    frame = ingest(data, qtype=qtype, source_label="round_trip_in")
    return emit(frame, target_label="round_trip_out")


__all__ = [
    "ConversionError",
    "emit",
    "emit_batch",
    "emit_json",
    "get_stats",
    "ingest",
    "ingest_batch",
    "ingest_file",
    "reset_stats",
    "round_trip",
]
