"""
Finalization failure audit logger.

Appends structured JSONL entries to PX_LOGS/finalization_failures.jsonl.
Each entry records: timestamp (ISO 8601 UTC), source file, candidate_id,
error message, and context string.

Dependencies: json, datetime, pathlib (stdlib only).
"""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[1]
_LOG_PATH = _REPO_ROOT / "PX_LOGS" / "finalization_failures.jsonl"


def log_finalization_failure(
    *,
    source_file: str,
    candidate_id: str,
    error: str,
    context: str,
) -> None:
    """Append a single finalization-failure record to the JSONL audit trail.

    Parameters
    ----------
    source_file : str
        Which orchestrator/pipeline file the failure occurred in.
    candidate_id : str
        The candidate or dossier ID being finalized (empty string if unknown).
    error : str
        The exception message or error description.
    context : str
        What was being done when the failure happened
        (e.g. "finalize_and_place call", "dossier read").
    """
    entry = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "file": source_file,
        "candidate_id": candidate_id,
        "error": error,
        "context": context,
    }
    try:
        _LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
        with open(_LOG_PATH, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(entry, ensure_ascii=False) + "\n")
    except Exception:
        # Last resort: if even the log write fails, there is nothing safe to do.
        # The caller already prints to stderr as a fallback.
        pass
