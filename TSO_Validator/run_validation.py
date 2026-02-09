#!/usr/bin/env python
"""
TSO_Validator — Standalone physics validation safety gate.

Entry point: discovers JSON files in data/raw/, validates against schema
and constitutional constraints, writes results to data/results/run_<timestamp>/.

Exit codes:
  0 — TSO_VALIDATED or TSO_INDETERMINATE (no data)
  1 — TSO_FAILED (at least one hard failure)

Architectural rule: ZERO imports from PX_System, PX_Engine, PX_Executive,
PX_Warehouse, or any PX_* module. Uses only Python stdlib.
"""

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

TSO_ROOT = Path(__file__).resolve().parent
RAW_DIR = TSO_ROOT / "data" / "raw"
RESULTS_DIR = TSO_ROOT / "data" / "results"
CONFIG_DIR = TSO_ROOT / "config"

# Import audit modules (sibling packages — no PX_* imports)
sys.path.insert(0, str(TSO_ROOT))
from audit.ethical_guard import run_ethical_guard
from audit.provenance import build_provenance
from audit.reproducibility import check_records_deterministic
from audit.drift import run_drift_analysis
from audit.traceability import build_trace_chain


def _load_raw_records(raw_dir: Path) -> list:
    """Load all JSON files from raw directory. Returns list of (path, record) tuples."""
    records = []
    for fp in sorted(raw_dir.glob("*.json")):
        try:
            with open(fp, "r", encoding="utf-8") as f:
                data = json.load(f)
            if isinstance(data, dict):
                records.append((fp, data))
            elif isinstance(data, list):
                for item in data:
                    if isinstance(item, dict):
                        records.append((fp, item))
        except (json.JSONDecodeError, IOError) as e:
            records.append((fp, {"_parse_error": str(e), "record_id": fp.name}))
    return records


def run_validation() -> dict:
    """Run full TSO validation pipeline. Returns run_summary dict."""
    timestamp = datetime.now(timezone.utc)
    run_id = f"run_{timestamp.strftime('%Y%m%d_%H%M%S')}"

    # 1. Discover raw files
    if not RAW_DIR.exists():
        RAW_DIR.mkdir(parents=True, exist_ok=True)

    raw_pairs = _load_raw_records(RAW_DIR)
    records = [rec for _, rec in raw_pairs]

    if len(records) == 0:
        return {
            "status": "TSO_VALIDATED",
            "run_id": run_id,
            "timestamp": timestamp.isoformat(),
            "files_checked": 0,
            "failures": [],
            "message": "No input files in data/raw/; validation trivially passes.",
        }

    # 2. Provenance (checksums for all raw files)
    provenance = build_provenance(RAW_DIR)

    # 3. Ethical guard (no-harm, no-mock per record)
    all_failures = []
    ethical_results = []
    for rec in records:
        if "_parse_error" in rec:
            all_failures.append(f"Parse error in {rec.get('record_id', '?')}: {rec['_parse_error']}")
            continue
        result = run_ethical_guard(rec)
        ethical_results.append(result)
        if not result["passed"]:
            for fail in result["failures"]:
                all_failures.append(f"[{result['record_id']}] {fail}")

    # 4. Reproducibility
    clean_records = [r for r in records if "_parse_error" not in r]
    repro = check_records_deterministic(clean_records)
    if not repro["passed"]:
        for item in repro["failures"]:
            for fail in item["failures"]:
                all_failures.append(f"[{item['record_id']}] {fail}")

    # 5. Schema drift
    drift = run_drift_analysis(clean_records, CONFIG_DIR)
    if not drift.get("passed", True):
        for item in drift.get("issues", []):
            for issue in item["issues"]:
                all_failures.append(f"[{item['record_id']}] DRIFT: {issue}")

    # 6. Traceability
    trace = build_trace_chain(clean_records, provenance, run_id)

    # 7. Determine status
    if len(all_failures) > 0:
        status = "TSO_FAILED"
    else:
        status = "TSO_VALIDATED"

    run_summary = {
        "status": status,
        "run_id": run_id,
        "timestamp": timestamp.isoformat(),
        "schema_version": drift.get("schema_version", "unknown"),
        "files_checked": len(raw_pairs),
        "records_checked": len(records),
        "failures": all_failures,
        "ethical_guard": {
            "records_checked": len(ethical_results),
            "records_passed": sum(1 for r in ethical_results if r["passed"]),
        },
        "reproducibility": {
            "passed": repro["passed"],
            "records_failed": repro["records_failed"],
        },
        "drift": {
            "passed": drift.get("passed", True),
            "records_with_drift": drift.get("records_with_drift", 0),
        },
    }

    # 8. Write results
    run_dir = RESULTS_DIR / run_id
    run_dir.mkdir(parents=True, exist_ok=True)

    with open(run_dir / "run_summary.json", "w", encoding="utf-8") as f:
        json.dump(run_summary, f, indent=2)

    with open(run_dir / "provenance.json", "w", encoding="utf-8") as f:
        json.dump(provenance, f, indent=2)

    with open(run_dir / "historical_context.json", "w", encoding="utf-8") as f:
        json.dump({
            "drift_analysis": drift,
            "trace_chain": trace,
        }, f, indent=2, default=str)

    return run_summary


def main() -> int:
    """Entry point. Returns 0 on pass, 1 on TSO_FAILED."""
    summary = run_validation()
    status = summary["status"]
    files = summary["files_checked"]
    failures = summary["failures"]

    print(f"TSO_Validator: {status} ({files} files, {len(failures)} failures)")
    if failures:
        for f in failures[:20]:
            print(f"  - {f}")
        if len(failures) > 20:
            print(f"  ... and {len(failures) - 20} more")

    if status == "TSO_FAILED":
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
