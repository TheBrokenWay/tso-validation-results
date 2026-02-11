"""
Sovereign_Log_Chain - Audit Trail
Provides immutable audit logging for OLYMPUS operations.
"""
import json
import hashlib
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Optional

_REPO_ROOT = Path(__file__).resolve().parents[2]
LOG_PATH = _REPO_ROOT / "PX_Audit" / "sovereign_log_chain.jsonl"

def append(event_type: str, data: Dict[str, Any], context: Optional[Dict] = None, qframe=None) -> bool:
    """
    Append a record to the sovereign log chain.
    """
    LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    prev_hash = get_chain_hash()
    log_entry = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "event_type": event_type,
        "data": data,
        "context": context or {},
        "prev_hash": prev_hash,
    }
    # QUINT integration: add frame metadata when available
    if qframe is not None:
        try:
            quint_meta = {
                "frame_hash": getattr(qframe, "frame_hash", ""),
                "qtype": getattr(qframe, "qtype", None),
                "qid": getattr(qframe, "qid", ""),
                "source_hash": getattr(qframe, "source_hash", ""),
                "lineage": list(getattr(qframe, "lineage", [])),
                "schema_version": getattr(qframe, "schema_version", ""),
            }
            # Convert QType enum to string if present
            if hasattr(quint_meta["qtype"], "value"):
                quint_meta["qtype"] = quint_meta["qtype"].value
            log_entry["quint_metadata"] = quint_meta
        except Exception:
            pass  # Fail-open for audit — never block logging
    entry_bytes = json.dumps(log_entry, sort_keys=True).encode("utf-8")
    log_entry["record_hash"] = hashlib.sha256(entry_bytes).hexdigest()

    with LOG_PATH.open("a", encoding="utf-8") as f:
        f.write(json.dumps(log_entry) + "\n")
    return True

def get_chain_hash() -> str:
    """Get the current chain hash."""
    if not LOG_PATH.exists():
        return "0" * 64
    try:
        with LOG_PATH.open("r", encoding="utf-8") as f:
            lines = f.readlines()
            if not lines:
                return "0" * 64
            last = json.loads(lines[-1])
            return last.get("record_hash", "0" * 64)
    except (OSError, json.JSONDecodeError):
        return "0" * 64

# Legacy exports for backward compatibility
__all__ = ['append', 'get_chain_hash']
