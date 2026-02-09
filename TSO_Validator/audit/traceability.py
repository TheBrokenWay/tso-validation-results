"""Traceability — end-to-end audit trail linking inputs to validation results."""

import hashlib
import json
from datetime import datetime, timezone


def build_trace_chain(records: list, provenance: dict, run_id: str) -> dict:
    """Build an audit trail linking input records to the validation run."""
    entries = []
    for rec in records:
        rec_id = rec.get("record_id", "UNKNOWN")
        entry = {
            "record_id": rec_id,
            "input_hash": hashlib.sha256(json.dumps(rec, sort_keys=True).encode()).hexdigest(),
            "run_id": run_id,
        }
        entries.append(entry)

    chain_payload = json.dumps(entries, sort_keys=True).encode()
    chain_hash = hashlib.sha256(chain_payload).hexdigest()

    return {
        "run_id": run_id,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "record_count": len(records),
        "provenance_file_count": provenance.get("file_count", 0),
        "entries": entries,
        "chain_hash": chain_hash,
    }
