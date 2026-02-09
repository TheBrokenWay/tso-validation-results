"""
Gate stamp signing — HMAC(signature, secret) for hostile-modification resistance.

Use when you expect hostile modification of the stamp file; for internal drift
only, hash(GOVERNANCE_SIGNATURE) is sufficient. Secret from env PX_GATE_STAMP_SECRET
or dev default (not for production).
"""
import os
import hmac
import hashlib
from typing import Optional

# Dev default: do not use in production; set PX_GATE_STAMP_SECRET for real signing
_DEV_SECRET = "PX_GATE_STAMP_DEV_DO_NOT_USE_IN_PRODUCTION"


def get_gate_stamp_secret() -> Optional[bytes]:
    """
    Return the gate stamp HMAC secret, or None if not configured.
    Prefer env PX_GATE_STAMP_SECRET (hex or raw); fallback to dev default when
    PX_GATE_STAMP_DEV_ALLOW=1 for local runs.
    """
    raw = os.environ.get("PX_GATE_STAMP_SECRET")
    if raw:
        try:
            return bytes.fromhex(raw)
        except ValueError:
            return raw.encode("utf-8") if isinstance(raw, str) else raw
    if os.environ.get("PX_GATE_STAMP_DEV_ALLOW") == "1":
        return _DEV_SECRET.encode("utf-8")
    return None


def sign_gate_stamp(payload: bytes) -> str:
    """
    HMAC-SHA256(payload, secret). Returns hex digest.
    Use with get_gate_stamp_secret(); if secret is None, caller should skip signing.
    """
    secret = get_gate_stamp_secret()
    if not secret:
        return ""
    return hmac.new(secret, payload, hashlib.sha256).hexdigest()


def verify_gate_stamp(payload: bytes, expected_hmac_hex: str) -> bool:
    """Verify stamp HMAC. Returns True if secret is set and HMAC matches."""
    if not expected_hmac_hex:
        return False
    return hmac.compare_digest(sign_gate_stamp(payload), expected_hmac_hex)


__all__ = ["get_gate_stamp_secret", "sign_gate_stamp", "verify_gate_stamp"]
