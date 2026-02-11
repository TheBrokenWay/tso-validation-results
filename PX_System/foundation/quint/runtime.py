"""
QUINT Runtime — Operations on QFrames.

Provides validated operations that engines use to work with QUINT data:
  - validate:   check a QFrame's integrity and type compliance
  - transform:  apply a function to a QFrame's payload, producing a new QFrame
  - merge:      combine multiple QFrames into one
  - extract:    pull specific fields from a QFrame
  - annotate:   add metadata to a QFrame without altering payload

All operations preserve lineage tracking and integrity seals.

Constitutional: Python stdlib only. Deterministic.
"""

from __future__ import annotations

from typing import Any, Callable, Dict, List, Optional, Set

from PX_System.foundation.quint.compiler import REQUIRED_FIELDS
from PX_System.foundation.quint.kernel import QFrame, QType, create_qframe


class RuntimeError(Exception):
    """Raised when a runtime operation fails."""
    pass


def validate(frame: QFrame) -> bool:
    """
    Full validation of a QFrame:
      1. Integrity seal (frame_hash matches)
      2. Required fields present for its QType
      3. QID is non-empty

    Returns True if valid, raises RuntimeError if not.
    """
    if not frame.qid:
        raise RuntimeError("QFrame has empty qid")

    if not frame.verify():
        raise RuntimeError(
            f"QFrame {frame.qid} integrity check failed — "
            f"frame_hash does not match current state"
        )

    required = REQUIRED_FIELDS.get(frame.qtype, [])
    missing = [f for f in required if f not in frame.payload]
    if missing:
        raise RuntimeError(
            f"QFrame {frame.qid} ({frame.qtype.value}) missing required fields: {missing}"
        )

    return True


def transform(
    frame: QFrame,
    fn: Callable[[Dict[str, Any]], Dict[str, Any]],
    label: str = "transform",
    new_qtype: Optional[QType] = None,
) -> QFrame:
    """
    Apply a pure function to a QFrame's payload, producing a new sealed QFrame.

    The original frame is NOT modified. The new frame inherits the source_hash
    and extends the lineage.
    """
    new_payload = fn(dict(frame.payload))
    if not isinstance(new_payload, dict):
        raise RuntimeError(f"Transform function must return dict, got {type(new_payload).__name__}")

    new_frame = create_qframe(
        qtype=new_qtype or frame.qtype,
        qid=frame.qid,
        payload=new_payload,
        source_hash=frame.source_hash,
        ingested_at=frame.ingested_at,
    )
    new_frame.lineage = list(frame.lineage) + [label]
    new_frame.seal()  # reseal after lineage update
    return new_frame


def merge(
    frames: List[QFrame],
    qtype: QType = QType.QDOSSIER,
    qid: Optional[str] = None,
    label: str = "merge",
) -> QFrame:
    """
    Merge multiple QFrames into a single QFrame.

    Payloads are combined (later frames override earlier on key conflict).
    Lineage tracks all source QIDs.
    """
    if not frames:
        raise RuntimeError("Cannot merge empty frame list")

    combined: Dict[str, Any] = {}
    source_qids: List[str] = []
    for f in frames:
        combined.update(f.payload)
        source_qids.append(f.qid)

    merged_qid = qid or f"MERGE-{'-'.join(source_qids[:4])}"
    merged = create_qframe(
        qtype=qtype,
        qid=merged_qid,
        payload=combined,
        source_hash=frames[0].source_hash,
        ingested_at=frames[0].ingested_at,
    )
    merged.lineage = [f"{label}({','.join(source_qids)})"]
    merged.seal()
    return merged


def extract(
    frame: QFrame,
    fields: Set[str],
    label: str = "extract",
) -> QFrame:
    """
    Extract a subset of fields from a QFrame into a new QRAW frame.
    """
    subset = {k: v for k, v in frame.payload.items() if k in fields}
    extracted = create_qframe(
        qtype=QType.QRAW,
        qid=f"{frame.qid}:extract",
        payload=subset,
        source_hash=frame.source_hash,
        ingested_at=frame.ingested_at,
    )
    extracted.lineage = list(frame.lineage) + [label]
    extracted.seal()
    return extracted


def annotate(
    frame: QFrame,
    annotations: Dict[str, Any],
    label: str = "annotate",
) -> QFrame:
    """
    Add metadata fields to a QFrame's payload without removing existing data.

    Returns a new sealed QFrame. Original is not modified.
    """
    new_payload = dict(frame.payload)
    new_payload.update(annotations)
    annotated = create_qframe(
        qtype=frame.qtype,
        qid=frame.qid,
        payload=new_payload,
        source_hash=frame.source_hash,
        ingested_at=frame.ingested_at,
    )
    annotated.lineage = list(frame.lineage) + [label]
    annotated.seal()
    return annotated


def chain(
    frame: QFrame,
    operations: List[Callable[[QFrame], QFrame]],
) -> QFrame:
    """
    Apply a sequence of QFrame → QFrame operations.
    Each operation receives the output of the previous one.
    """
    current = frame
    for op in operations:
        current = op(current)
    return current


__all__ = [
    "RuntimeError",
    "annotate",
    "chain",
    "extract",
    "merge",
    "transform",
    "validate",
]
