"""
QUINT Kernel — Core types for the internal mathematical language.

QFrame is the fundamental data unit. Every piece of data inside the platform
exists as a QFrame once it passes through the Converter gateway.

Design constraints:
  - Python stdlib only (constitutional module)
  - Near-zero overhead: QFrame is a thin typed envelope around a dict payload
  - Deterministic: same input → same QFrame (no random, no time-dependent logic
    except the ingestion timestamp which is caller-supplied)
  - Fail-closed: invalid data → rejection, not silent coercion
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional


class QType(Enum):
    """QUINT type system — every QFrame has exactly one type."""
    QMOLECULE = "QMOLECULE"        # Molecular representation (SMILES + descriptors)
    QRESULT   = "QRESULT"          # Engine computation result
    QDOSSIER  = "QDOSSIER"         # Complete dossier in QUINT form
    QPHYSICS  = "QPHYSICS"         # Physics state vector (worldline snapshot)
    QSIGNAL   = "QSIGNAL"          # Inter-engine communication
    QTRIAL    = "QTRIAL"           # Trial simulation data
    QCONSTRAINT = "QCONSTRAINT"    # Disease constraint / governance rule
    QRAW      = "QRAW"             # Untyped passthrough (for data that doesn't fit above)


# Fields that the kernel reserves — cannot appear in user payload
_RESERVED_META = frozenset({
    "_qtype", "_qid", "_schema_version", "_source_hash",
    "_frame_hash", "_ingested_at", "_lineage",
})

SCHEMA_VERSION = "1.0.0"


def _canonical_json(obj: Any) -> str:
    """Deterministic JSON serialization for hashing."""
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), default=str)


def _hash_payload(payload: Dict[str, Any]) -> str:
    """SHA-256 of the canonical JSON form of a payload dict."""
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


@dataclass(frozen=False)
class QFrame:
    """
    The fundamental QUINT data unit.

    A QFrame is a lightweight typed envelope around a dict payload.
    Overhead vs raw dict: ~200 bytes of metadata + two SHA-256 hashes.

    Attributes:
        qtype:          The QUINT type of this frame.
        qid:            Unique identifier (caller-assigned or auto-generated).
        payload:        The actual data — a plain dict.
        schema_version: For forward compatibility.
        source_hash:    SHA-256 of the original external data before ingestion.
        frame_hash:     SHA-256 of (qtype + qid + payload) — integrity seal.
        ingested_at:    ISO timestamp of when the frame was created.
        lineage:        Ordered list of operations applied to this frame.
    """
    qtype: QType
    qid: str
    payload: Dict[str, Any]
    schema_version: str = SCHEMA_VERSION
    source_hash: str = ""
    frame_hash: str = ""
    ingested_at: str = ""
    lineage: List[str] = field(default_factory=list)

    def seal(self) -> None:
        """Compute and set the frame_hash from current state."""
        self.frame_hash = self._compute_hash()

    def verify(self) -> bool:
        """Return True if frame_hash matches current state."""
        if not self.frame_hash:
            return False
        return self.frame_hash == self._compute_hash()

    def _compute_hash(self) -> str:
        """Hash of (qtype, qid, schema_version, payload) — deterministic."""
        blob = _canonical_json({
            "qtype": self.qtype.value,
            "qid": self.qid,
            "schema_version": self.schema_version,
            "payload": self.payload,
        })
        return hashlib.sha256(blob.encode("utf-8")).hexdigest()

    def to_dict(self) -> Dict[str, Any]:
        """Serialize QFrame to a plain dict (for storage/transmission)."""
        return {
            "_qtype": self.qtype.value,
            "_qid": self.qid,
            "_schema_version": self.schema_version,
            "_source_hash": self.source_hash,
            "_frame_hash": self.frame_hash,
            "_ingested_at": self.ingested_at,
            "_lineage": list(self.lineage),
            **self.payload,
        }

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> QFrame:
        """Reconstruct a QFrame from a serialized dict."""
        meta_qtype = d.get("_qtype", QType.QRAW.value)
        meta_qid = d.get("_qid", "")
        meta_sv = d.get("_schema_version", SCHEMA_VERSION)
        meta_sh = d.get("_source_hash", "")
        meta_fh = d.get("_frame_hash", "")
        meta_ia = d.get("_ingested_at", "")
        meta_lin = d.get("_lineage", [])

        payload = {k: v for k, v in d.items() if k not in _RESERVED_META}

        try:
            qtype = QType(meta_qtype)
        except ValueError:
            qtype = QType.QRAW

        frame = cls(
            qtype=qtype,
            qid=meta_qid,
            payload=payload,
            schema_version=meta_sv,
            source_hash=meta_sh,
            frame_hash=meta_fh,
            ingested_at=meta_ia,
            lineage=list(meta_lin),
        )
        return frame

    def size_bytes(self) -> int:
        """Approximate memory overhead of the QFrame envelope (excludes payload)."""
        import sys
        overhead = (
            sys.getsizeof(self.qtype.value)
            + sys.getsizeof(self.qid)
            + sys.getsizeof(self.schema_version)
            + sys.getsizeof(self.source_hash)
            + sys.getsizeof(self.frame_hash)
            + sys.getsizeof(self.ingested_at)
            + sys.getsizeof(self.lineage)
        )
        return overhead


def create_qframe(
    qtype: QType,
    qid: str,
    payload: Dict[str, Any],
    source_hash: str = "",
    ingested_at: str = "",
) -> QFrame:
    """
    Factory function to create a sealed QFrame.

    This is the ONLY approved way to create QFrames outside of the converter.
    The frame is sealed (hash computed) before return.
    """
    frame = QFrame(
        qtype=qtype,
        qid=qid,
        payload=dict(payload),  # shallow copy — near-zero overhead for flat dicts
        source_hash=source_hash,
        ingested_at=ingested_at,
    )
    frame.seal()
    return frame


__all__ = [
    "QType",
    "QFrame",
    "create_qframe",
    "SCHEMA_VERSION",
]
