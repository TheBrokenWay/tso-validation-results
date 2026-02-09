"""
Governance coverage metric — auto-generated from poison-pill registry.

Purely informational, but powerful for audits. Run after registry changes
or on successful poison-pill gate to keep coverage current.
"""
import json
from pathlib import Path

_GOVERNANCE_DIR = Path(__file__).resolve().parent
REGISTRY_PATH = _GOVERNANCE_DIR / "poison_pills" / "registry.json"
COVERAGE_PATH = _GOVERNANCE_DIR / "governance_coverage.json"


def generate_governance_coverage() -> dict:
    """
    Read registry.json, compute total/active/reserved/coverage, write governance_coverage.json.
    Returns the coverage dict.
    """
    if not REGISTRY_PATH.exists():
        out = {
            "total_pills": 0,
            "active": 0,
            "reserved": 0,
            "coverage": 0.0,
        }
        with open(COVERAGE_PATH, "w", encoding="utf-8") as f:
            json.dump(out, f, indent=2)
        return out

    with open(REGISTRY_PATH, encoding="utf-8") as f:
        registry = json.load(f)

    pills = registry.get("pills", [])
    total = len(pills)
    reserved = sum(1 for p in pills if p.get("status") == "reserved")
    active = total - reserved
    coverage = (active / total) if total else 0.0

    out = {
        "total_pills": total,
        "active": active,
        "reserved": reserved,
        "coverage": round(coverage, 2),
    }

    with open(COVERAGE_PATH, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)

    # Auto-generate Governance Readiness Certificate (coverage + invariants + stamp)
    try:
        from governance.certificate import generate_governance_readiness_certificate
        generate_governance_readiness_certificate()
    except Exception:
        pass

    return out


if __name__ == "__main__":
    d = generate_governance_coverage()
    print(json.dumps(d, indent=2))
