"""Reproducibility verification — determinism checks."""

import json
from pathlib import Path


def check_determinism(record: dict) -> list:
    """Check that a record contains no non-deterministic markers. Returns failures."""
    failures = []
    # Flag records that claim random generation without a seed
    if record.get("generation_method") == "random" and "rng_seed" not in record:
        failures.append("Non-deterministic: generation_method=random but no rng_seed")
    # Flag NaN/Inf values (non-reproducible computation artifacts)
    for key, val in record.items():
        if isinstance(val, float):
            if val != val:  # NaN check
                failures.append(f"NaN detected in field: {key}")
            elif val == float("inf") or val == float("-inf"):
                failures.append(f"Infinity detected in field: {key}")
    return failures


def check_records_deterministic(records: list) -> dict:
    """Check all records for reproducibility. Returns summary."""
    all_failures = []
    for rec in records:
        rec_failures = check_determinism(rec)
        if rec_failures:
            all_failures.append({
                "record_id": rec.get("record_id", "UNKNOWN"),
                "failures": rec_failures,
            })
    return {
        "records_checked": len(records),
        "records_failed": len(all_failures),
        "failures": all_failures,
        "passed": len(all_failures) == 0,
    }
