"""
Deprecated redirect: use PX_Validation.final_validation_FINAL (hardened suite).
This module preserves the public API so existing callers keep working.
"""
from PX_Validation.final_validation_FINAL import (
    run_hardened_validation,
    run_final_validation,
    seal_for_distribution,
)

__all__ = ["run_final_validation", "run_hardened_validation", "seal_for_distribution"]
