"""Ethical guard — constitutional enforcement for validated data.

Rules enforced (mirroring ZeusLaws L10/L11 without importing PX_*):
  - No-harm: toxicity_index < 0.0210, harm_energy < 0.0210
  - No-mock: records must contain real numeric values, not placeholders
  - No-sycophancy: thresholds are absolute, no rounding or adjustment
"""

TOXICITY_HARD_LIMIT = 0.0210
HARM_ENERGY_HARD_LIMIT = 0.0210
HARMONIC_OVERDRIVE = 1.02  # Fixed constant, never adaptive


def check_no_harm(record: dict) -> list:
    """Check toxicity and harm energy hard limits. Returns list of failures."""
    failures = []
    tox = record.get("toxicity_index")
    if tox is not None:
        if not isinstance(tox, (int, float)):
            failures.append(f"toxicity_index is not numeric: {type(tox).__name__}")
        elif tox >= TOXICITY_HARD_LIMIT:
            failures.append(f"L11 violation: toxicity_index {tox} >= {TOXICITY_HARD_LIMIT}")

    harm = record.get("harm_energy")
    if harm is not None:
        if not isinstance(harm, (int, float)):
            failures.append(f"harm_energy is not numeric: {type(harm).__name__}")
        elif harm >= HARM_ENERGY_HARD_LIMIT:
            failures.append(f"L10 violation: harm_energy {harm} >= {HARM_ENERGY_HARD_LIMIT}")

    return failures


def check_no_mock(record: dict) -> list:
    """Detect placeholder or mock data. Returns list of failures."""
    failures = []
    for key in ("toxicity_index", "binding_affinity_kj", "coherence", "harm_energy"):
        val = record.get(key)
        if val is None:
            continue
        if isinstance(val, str):
            failures.append(f"Mock data: {key} is string '{val}', expected numeric")
        if isinstance(val, (int, float)) and val == -1:
            failures.append(f"Placeholder detected: {key} = -1")
    record_id = record.get("record_id", "")
    if isinstance(record_id, str) and any(tag in record_id.upper() for tag in ("MOCK", "DUMMY", "TEST", "FAKE")):
        failures.append(f"Mock record_id detected: {record_id}")
    return failures


def run_ethical_guard(record: dict) -> dict:
    """Run all ethical checks on a single record. Returns result dict."""
    failures = []
    failures.extend(check_no_harm(record))
    failures.extend(check_no_mock(record))
    return {
        "record_id": record.get("record_id", "UNKNOWN"),
        "passed": len(failures) == 0,
        "failures": failures,
    }
