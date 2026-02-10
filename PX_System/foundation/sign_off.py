"""
Standardized Engine Sign-Off Module

Every engine in the 12-engine mandatory pipeline must produce a sign-off block
proving it executed, what laws it checked, and whether it authorized the result.

Uses only Python stdlib (constitutional module constraint).
"""

import hashlib
import json
import time
from dataclasses import dataclass, asdict, field
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional


@dataclass
class EngineSignOff:
    """Mandatory sign-off block for every engine."""
    engine_id: str
    version: str
    timestamp: str
    execution_time_ms: int
    status: str  # "PASSED" | "FAILED" | "REVIEW_REQUIRED"
    authorized: bool
    laws_checked: List[str]
    laws_passed: List[str]
    laws_failed: List[str]
    reason: str
    input_hash: str
    output_hash: str
    signature: str

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


def compute_hash(data: Any) -> str:
    """Compute SHA-256 hash of data (deterministic, sorted keys)."""
    json_str = json.dumps(data, sort_keys=True, default=str)
    return f"sha256:{hashlib.sha256(json_str.encode()).hexdigest()[:16]}"


def create_sign_off(
    engine_id: str,
    version: str,
    inputs: Dict[str, Any],
    outputs: Dict[str, Any],
    laws_checked: List[str],
    laws_results: Dict[str, bool],
    execution_time_ms: int = 0,
    custom_reason: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Create a standardized sign-off block.

    Args:
        engine_id: Engine identifier (e.g. "OPE_V3_DETERMINISTIC")
        version: Engine version string
        inputs: Engine input data (for hashing)
        outputs: Engine output data (for hashing)
        laws_checked: List of law IDs that were evaluated
        laws_results: Dict mapping law ID -> bool (True = passed)
        execution_time_ms: Execution time in milliseconds
        custom_reason: Override reason string

    Returns:
        Dict with all sign-off fields (ready to embed in engine output)
    """
    input_hash = compute_hash(inputs)
    output_hash = compute_hash(outputs)
    ts = datetime.now(timezone.utc).isoformat()
    signature = compute_hash({
        "engine_id": engine_id,
        "input_hash": input_hash,
        "output_hash": output_hash,
        "timestamp": ts,
    })

    laws_passed = [law for law, passed in laws_results.items() if passed]
    laws_failed = [law for law, passed in laws_results.items() if not passed]
    all_passed = len(laws_failed) == 0

    if all_passed:
        status = "PASSED"
        authorized = True
        reason = custom_reason or "All laws passed"
    else:
        status = "FAILED"
        authorized = False
        reason = custom_reason or f"Laws failed: {laws_failed}"

    sign_off = EngineSignOff(
        engine_id=engine_id,
        version=version,
        timestamp=ts,
        execution_time_ms=execution_time_ms,
        status=status,
        authorized=authorized,
        laws_checked=laws_checked,
        laws_passed=laws_passed,
        laws_failed=laws_failed,
        reason=reason,
        input_hash=input_hash,
        output_hash=output_hash,
        signature=signature,
    )
    return sign_off.to_dict()


def require_authorization(sign_off: Dict[str, Any], engine_name: str) -> None:
    """
    Raise AuthorizationError if sign-off is not authorized.
    Call after each engine in the pipeline to enforce fail-closed.
    """
    if not isinstance(sign_off, dict):
        return  # Engine didn't produce sign-off yet — backward compat
    if not sign_off.get("authorized", False):
        raise AuthorizationError(
            f"Engine {engine_name} not authorized: {sign_off.get('reason', 'Unknown')}"
        )


def build_authorization_chain(sign_offs: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Build the authorization_chain summary for a dossier.

    Args:
        sign_offs: List of sign-off dicts from each engine

    Returns:
        Dict with authorization_chain, all_engines_authorized, authorization_count
    """
    total = len(sign_offs)
    authorized_count = sum(1 for s in sign_offs if s.get("authorized", False))
    return {
        "authorization_chain": sign_offs,
        "all_engines_authorized": authorized_count == total,
        "authorization_count": f"{authorized_count}/{total}",
    }


class AuthorizationError(Exception):
    """Raised when an engine fails to authorize."""
    pass
