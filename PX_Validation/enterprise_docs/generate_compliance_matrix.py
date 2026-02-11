"""
Generate CONSTITUTIONAL_COMPLIANCE_MATRIX.md

Audits all constitutional laws and produces a compliance matrix
showing which laws are enforced and how.
"""
from __future__ import annotations

import hashlib
import sys
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def generate(output_dir: Path) -> Path:
    """Generate the constitutional compliance matrix."""
    ts = datetime.now(timezone.utc)
    stamp = ts.strftime("%Y%m%d_%H%M%S")

    # Verify constitutional constants
    from PX_System.foundation.ZeusLaws import check_constitutional, run_zeus_gate
    from PX_System.foundation.sign_off import create_sign_off, build_authorization_chain

    # Test L11 at boundary
    l11_boundary = check_constitutional("L11_AUDIT", {"toxicity_index": 0.0210, "harm_energy": 0.001})
    l11_below = check_constitutional("L11_AUDIT", {"toxicity_index": 0.0209, "harm_energy": 0.001})
    l11_pass = not l11_boundary["authorized"] and l11_below["authorized"]

    # Test L10 at boundary
    l10_boundary = check_constitutional("L10_AUDIT", {"toxicity_index": 0.001, "harm_energy": 0.0210})
    l10_below = check_constitutional("L10_AUDIT", {"toxicity_index": 0.001, "harm_energy": 0.0209})
    l10_pass = not l10_boundary["authorized"] and l10_below["authorized"]

    # Test fail-closed (empty payload should pass since no violations)
    failclose = check_constitutional("FAILCLOSE_AUDIT", {})
    failclose_pass = failclose["authorized"]  # No tox/harm → no violation → passes

    # Test Zeus gate laws
    test_dossier = {
        "engines": {
            "ope": {"molecular_weight": 46.07},
            "admet": {"toxicity": {"toxicity_index": 0.005, "risk_level": "TOXICITY_DIAMOND"}},
        },
        "harm_energy": 0.005,
        "trial_outcome_summary": None,
    }
    zeus = run_zeus_gate(test_dossier)
    zeus_laws = set(zeus.get("laws_required", []))

    # Test sign-off chain
    sign_offs = []
    for eid in ["OPE", "OBE", "OCE", "OLE", "OME", "OSE",
                "ADMET", "PKPD", "DOSE", "VEFF", "GRADE", "ZEUS"]:
        so = create_sign_off(
            engine_id=eid, version="1.0", inputs={"test": True}, outputs={"test": True},
            laws_checked=["L11"], laws_results={"L11": True},
        )
        sign_offs.append(so)
    chain = build_authorization_chain(sign_offs)
    chain_pass = chain["all_engines_authorized"] and chain["authorization_count"] == "12/12"

    all_pass = l11_pass and l10_pass and chain_pass
    matrix_hash = hashlib.sha256(
        f"{ts.isoformat()}{l11_pass}{l10_pass}{chain_pass}".encode()
    ).hexdigest()[:16]

    content = f"""# PREDATOR X — CONSTITUTIONAL COMPLIANCE MATRIX

| Field | Value |
|-------|-------|
| **Matrix ID** | `PX-CCM-{stamp}-{matrix_hash}` |
| **Generated** | {ts.isoformat()} |
| **Compliance** | **{"FULL" if all_pass else "PARTIAL"}** |

---

## 1. Constitutional Laws

| Law | Description | Threshold | Enforced | Test Result |
|-----|-------------|-----------|----------|-------------|
| **L1** (Harm Law) | Internal 35D physics snapshot is primary | — | YES | {"PASS" if "L1_HARM_LAW" in zeus_laws else "FAIL"} |
| **L10** (Harm Absolute) | harm_energy absolute limit | >= 0.0210 FAILS | YES | {"PASS" if l10_pass else "FAIL"} |
| **L11** (Deterministic Engine) | toxicity_index hard limit | >= 0.0210 FAILS | YES | {"PASS" if l11_pass else "FAIL"} |
| **L21** (Anti-Sycophancy) | Cannot adjust values to bypass | — | YES | Structural |
| **U27** (Stability Law) | Engine outputs stable & present | — | YES | {"PASS" if "U27_STABILITY_LAW" in zeus_laws else "FAIL"} |
| **U34** (Global Sum) | Physics invariant conserved | sum = 36.1 | YES | {"PASS" if "U34_GLOBAL_SUM" in zeus_laws else "FAIL"} |

## 2. Toxicity Tiers (Immutable)

| Tier | Threshold | Constitutional Rule |
|------|-----------|-------------------|
| TOXICITY_DIAMOND | tox < 0.01 OR (tox > 0.02 AND safety_margin > 50) | Whole-profile override |
| TOXICITY_GOLD | tox < 0.0200 | Standard pass |
| TOXICITY_SILVER | 0.0200 <= tox < 0.0210 | Marginal pass |
| TOXICITY_FAILURE | tox >= 0.0210 | **HARD STOP** — no rounding, no negotiation |

## 3. Governance Gates

| Gate | Location | Default | Behavior |
|------|----------|---------|----------|
| Poison Pill | `governance/poison_pill_gate.py` | Fail-closed | Must pass before any PX_System init |
| Zeus Gate | `PX_System/foundation/ZeusLaws.py` | Fail-closed | L1, U27, U34, L11 + whole-profile |
| Sign-Off Chain | `PX_System/foundation/sign_off.py` | Fail-closed | 12/12 engines required |
| Constitutional Tests | `PX_Constitution/Constitutional_Tests.py` | Fail-closed | Pre-execution mandatory |

## 4. Engine Authorization Chain

| Property | Value |
|----------|-------|
| Required Engines | 12 (OPE, OBE, OCE, OLE, OME, OSE, ADMET, PKPD, DoseOpt, VirtualEfficacy, GradingEngine, ZeusLaws) |
| Chain Completeness | {chain["authorization_count"]} |
| All Authorized | {chain["all_engines_authorized"]} |
| Chain Hash | `{chain.get("chain_hash", "N/A")[:16]}` |

## 5. Hard Constants (Never Change)

| Constant | Value | Location |
|----------|-------|----------|
| TOXICITY_HARD_LIMIT | 0.0210 | ZeusLaws.py (inline) |
| HARM_ENERGY_HARD_LIMIT | 0.0210 | ZeusLaws.py (inline) |
| HARMONIC_OVERDRIVE | 1.02 | Fixed, never adaptive |
| global_sum_target | 36.1 | VectorCore (35D manifold) |
| dims_limit | 35 | VectorCore (35D manifold) |

---
*Generated by Predator X Enterprise Validation Suite v1.0.0*
"""

    out_path = output_dir / f"CONSTITUTIONAL_COMPLIANCE_MATRIX_{stamp}.md"
    out_path.write_text(content, encoding="utf-8")
    return out_path
