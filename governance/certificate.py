"""
Governance Readiness Certificate — auto-generated from coverage + invariants.

Purely informational, but powerful for audits: combines governance coverage,
constitutional invariants, and stamp status into a single certificate.
"""
import json
from pathlib import Path
from datetime import datetime, timezone

_GOVERNANCE_DIR = Path(__file__).resolve().parent
_ROOT_DIR = _GOVERNANCE_DIR.parent
COVERAGE_PATH = _GOVERNANCE_DIR / "governance_coverage.json"
STAMP_PATH = _GOVERNANCE_DIR / ".poison_pill_gate_passed"
CERTIFICATE_PATH = _GOVERNANCE_DIR / "governance_readiness_certificate.json"
CERTIFICATE_MD_PATH = _GOVERNANCE_DIR / "governance_readiness_certificate.md"


def generate_governance_readiness_certificate() -> dict:
    """
    Generate Governance Readiness Certificate from coverage + invariants + stamp.
    Writes JSON and Markdown; returns the certificate dict.
    """
    # Coverage
    coverage = {}
    if COVERAGE_PATH.exists():
        with open(COVERAGE_PATH, encoding="utf-8") as f:
            coverage = json.load(f)

    # Invariants
    invariants = {}
    try:
        import sys
        sys.path.insert(0, str(_ROOT_DIR))
        from PX_Constitution import mandatory_pre_execution_tests, LAYER_MONOTONICITY_INVARIANT
        invariants = {
            "mandatory_pre_execution_tests": mandatory_pre_execution_tests,
            "layer_monotonicity_invariant": LAYER_MONOTONICITY_INVARIANT,
        }
    except Exception:
        invariants = {"mandatory_pre_execution_tests": [], "layer_monotonicity_invariant": ""}

    # Stamp
    stamp_passed = STAMP_PATH.exists()
    stamp_mtime = None
    if stamp_passed and STAMP_PATH.exists():
        stamp_mtime = datetime.fromtimestamp(STAMP_PATH.stat().st_mtime, tz=timezone.utc).isoformat()

    cert = {
        "generated": datetime.now(timezone.utc).isoformat(),
        "governance_coverage": coverage,
        "invariants": invariants,
        "stamp": {
            "passed": stamp_passed,
            "stamp_file": ".poison_pill_gate_passed",
            "last_passed_at": stamp_mtime,
        },
        "total_pills": coverage.get("total_pills", 0),
        "active_pills": coverage.get("active", 0),
        "coverage_ratio": coverage.get("coverage", 0.0),
    }

    with open(CERTIFICATE_PATH, "w", encoding="utf-8") as f:
        json.dump(cert, f, indent=2)

    # Markdown summary
    md_lines = [
        "# Governance Readiness Certificate",
        "",
        f"**Generated:** {cert['generated']}",
        "",
        "## Coverage",
        f"- Total pills: {cert['total_pills']}",
        f"- Active: {cert['active_pills']}",
        f"- Reserved: {coverage.get('reserved', 0)}",
        f"- Coverage ratio: {cert['coverage_ratio']}",
        "",
        "## Invariants",
        f"- Mandatory pre-execution tests: {', '.join(invariants.get('mandatory_pre_execution_tests', []))}",
        f"- Layer monotonicity: {invariants.get('layer_monotonicity_invariant', '')}",
        "",
        "## Gate stamp",
        f"- Passed: {stamp_passed}",
        f"- Last passed: {stamp_mtime or 'N/A'}",
        "",
    ]
    with open(CERTIFICATE_MD_PATH, "w", encoding="utf-8") as f:
        f.write("\n".join(md_lines))

    return cert


if __name__ == "__main__":
    cert = generate_governance_readiness_certificate()
    print(json.dumps(cert, indent=2))
