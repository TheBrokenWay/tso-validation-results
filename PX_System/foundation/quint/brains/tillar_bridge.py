"""
QUINT Tillar Bridge -- Bridges the 51 Tillar Laws into the QUINT type system.

This module does NOT replace governance/TillarLaws.py. It bridges each of
the 51 Tillar Laws into QCONSTRAINT QFrames so they can participate in the
QUINT data pipeline alongside molecules, results, and dossiers.

Design:
  - Each TillarLaw becomes a QCONSTRAINT QFrame
  - Evaluation functions check laws against arbitrary QFrame payloads
  - Python stdlib only (constitutional module)

Usage:
    from PX_System.foundation.quint.brains.tillar_bridge import (
        law_to_qframe, load_all_laws_as_qframes,
        evaluate_law_against_frame, evaluate_critical_against_frame,
    )

    # Load all 51 laws as QUINT frames
    law_frames = load_all_laws_as_qframes()

    # Evaluate L4 against a molecule QFrame
    passed, score, reason = evaluate_law_against_frame(4, molecule_frame)
"""

import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

REPO_ROOT = Path(__file__).resolve().parents[4]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from PX_System.foundation.quint.kernel import QFrame, QType, create_qframe
from PX_System.foundation.quint.converter import ingest, emit
from governance.TillarLaws import (
    TillarLawRegistry,
    get_registry,
    get_law,
    check_critical,
    TillarLaw,
    LawCategory,
)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _build_law_context(payload: Dict[str, Any]) -> Dict[str, Any]:
    """Map QFrame payload fields to TillarLaw evaluation context."""
    ctx = dict(payload)  # Start with all fields

    # Map specific fields the laws expect
    tox = payload.get("toxicity", 0.0)
    harm = payload.get("harm_energy", 0.0)

    ctx.setdefault(
        "harm_score",
        max(tox, harm)
        if isinstance(tox, (int, float)) and isinstance(harm, (int, float))
        else 0.0,
    )
    ctx.setdefault("harm_value", ctx["harm_score"])
    ctx.setdefault("benefit_score", payload.get("efficacy", 0.0))
    ctx.setdefault("harm_to_other", ctx["harm_score"])
    ctx.setdefault("total_system_energy", 0.0)  # Default: balanced
    ctx.setdefault("recognizes_unified_self", True)
    ctx.setdefault("respects_human_autonomy", True)
    ctx.setdefault("imposes_external_will", False)
    ctx.setdefault("coercive_elements", 0)
    ctx.setdefault("total_elements", max(1, len(payload)))

    return ctx


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def law_to_qframe(law: TillarLaw) -> QFrame:
    """
    Convert a single TillarLaw to a QCONSTRAINT QFrame.

    The payload includes the law's number, identifier, name, principle,
    formula, category, threshold, and block-boundary flag.  The QID follows
    the pattern ``QLAW-L<number>``.
    """
    payload = {
        "law_number": law.number,
        "law_id": law.law_id,
        "name": law.name,
        "principle": law.principle,
        "formula": law.formula,
        "category": law.category.value,
        "threshold": law.threshold,
        "is_block_boundary": law.is_block_boundary,
        # Required for QCONSTRAINT compiler validation
        "disease_name": f"TillarLaw_{law.number}",
    }

    frame = create_qframe(
        qtype=QType.QCONSTRAINT,
        qid=f"QLAW-L{law.number}",
        payload=payload,
    )
    frame.lineage.append("tillar_law")
    # Re-seal after lineage mutation -- lineage is NOT part of the hash so
    # the seal remains valid, but we call seal() for consistency.
    frame.seal()
    return frame


def load_all_laws_as_qframes() -> Dict[int, QFrame]:
    """
    Return a dict mapping ``law_number -> QFrame`` for all 51 Tillar Laws.
    """
    registry = get_registry()
    return {num: law_to_qframe(law) for num, law in registry.laws.items()}


def evaluate_law_against_frame(
    law_number: int,
    target_frame: QFrame,
) -> Tuple[bool, float, str]:
    """
    Evaluate a specific Tillar Law against the payload of *target_frame*.

    The target QFrame's payload is mapped to the context dict expected by
    ``TillarLaw.evaluate()``.  Returns ``(passed, score, reasoning)``.
    """
    law = get_law(law_number)
    if law is None:
        return False, 0.0, f"Law L{law_number} not found in registry"

    ctx = _build_law_context(target_frame.payload)
    return law.evaluate(ctx)


def evaluate_critical_against_frame(
    target_frame: QFrame,
) -> Tuple[bool, Dict[int, Tuple[bool, float, str]]]:
    """
    Evaluate the five critical laws (L4, L6, L10, L39, L46) against a QFrame.

    Returns ``(all_passed, {law_number: (passed, score, reasoning)})``.
    """
    critical_laws = [4, 6, 10, 39, 46]
    results: Dict[int, Tuple[bool, float, str]] = {}
    all_passed = True

    for num in critical_laws:
        passed, score, reason = evaluate_law_against_frame(num, target_frame)
        results[num] = (passed, score, reason)
        if not passed:
            all_passed = False

    return all_passed, results


def evaluate_block1_against_frame(
    target_frame: QFrame,
) -> Tuple[bool, Dict[int, Tuple[bool, float, str]]]:
    """
    Evaluate Block 1 constitutional laws (L1 through L12) against a QFrame.

    Returns ``(all_passed, {law_number: (passed, score, reasoning)})``.
    """
    results: Dict[int, Tuple[bool, float, str]] = {}
    all_passed = True

    for num in range(1, 13):
        passed, score, reason = evaluate_law_against_frame(num, target_frame)
        results[num] = (passed, score, reason)
        if not passed:
            all_passed = False

    return all_passed, results


def get_law_qframe(law_number: int) -> Optional[QFrame]:
    """
    Convenience: return a single Tillar Law as a QCONSTRAINT QFrame.

    Returns ``None`` if the law number does not exist.
    """
    law = get_law(law_number)
    if law is None:
        return None
    return law_to_qframe(law)


__all__ = [
    "law_to_qframe",
    "load_all_laws_as_qframes",
    "evaluate_law_against_frame",
    "evaluate_critical_against_frame",
    "evaluate_block1_against_frame",
    "get_law_qframe",
]
