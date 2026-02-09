"""Schema drift detection — compare input data against validation schema."""

import json
from pathlib import Path


def load_schema(config_dir: Path) -> dict:
    """Load the canonical validation schema."""
    schema_path = config_dir / "validation_schema.json"
    if not schema_path.exists():
        return {}
    with open(schema_path, "r", encoding="utf-8") as f:
        return json.load(f)


def check_schema_drift(record: dict, schema: dict) -> list:
    """Check a single record against the schema. Returns list of drift issues."""
    failures = []
    required = schema.get("required_fields", [])
    for field in required:
        if field not in record:
            failures.append(f"Missing required field: {field}")

    physics = schema.get("physics_fields", {})
    for field_name, spec in physics.items():
        val = record.get(field_name)
        if val is None:
            continue  # Optional physics fields
        expected_type = spec.get("type")
        if expected_type == "number" and not isinstance(val, (int, float)):
            failures.append(f"Type drift: {field_name} expected number, got {type(val).__name__}")
            continue
        if isinstance(val, (int, float)):
            fmin = spec.get("min")
            fmax = spec.get("max")
            if fmin is not None and val < fmin:
                failures.append(f"Range drift: {field_name}={val} < min={fmin}")
            if fmax is not None and val > fmax:
                failures.append(f"Range drift: {field_name}={val} > max={fmax}")

    constraints = schema.get("constraints", {})
    tox_limit = constraints.get("toxicity_hard_limit")
    if tox_limit is not None:
        tox = record.get("toxicity_index")
        if isinstance(tox, (int, float)) and tox >= tox_limit:
            failures.append(f"Constraint violation: toxicity_index {tox} >= hard limit {tox_limit}")

    return failures


def run_drift_analysis(records: list, config_dir: Path) -> dict:
    """Run drift analysis on all records. Returns summary."""
    schema = load_schema(config_dir)
    if not schema:
        return {"status": "NO_SCHEMA", "message": "No validation_schema.json found", "passed": True}

    all_issues = []
    for rec in records:
        issues = check_schema_drift(rec, schema)
        if issues:
            all_issues.append({
                "record_id": rec.get("record_id", "UNKNOWN"),
                "issues": issues,
            })
    return {
        "schema_version": schema.get("schema_version", "unknown"),
        "records_checked": len(records),
        "records_with_drift": len(all_issues),
        "issues": all_issues,
        "passed": len(all_issues) == 0,
    }
