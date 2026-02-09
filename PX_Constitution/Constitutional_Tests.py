"""
Constitutional Tests — Mandatory pre-execution gates.

Formally elevates governance tests from "governance" to constitutional invariant
enforcement. Pipeline runners MUST run these before PX_System init, PX_Warehouse
writes, or PX_Engine invocation.
"""

# Mandatory pre-execution tests: run these at pipeline entry; block if any fail
# Paths relative to repo root (filed under tests/ per Root Filing Plan)
mandatory_pre_execution_tests = [
    "tests/run_poison_pill_test.py",
    "tests/run_temporal_paradox_test.py",
]

# Layer-monotonicity invariant (evolutionary; generalize beyond Nipah when ready):
# A failure may not be reclassified to a later layer once rejected.
# E.g. if PX_Validation rejects, the failure must not later be reported as PX_Engine
# or Byzantium_Council. Preserves architectural purity: ontology stays at PX_Validation.
LAYER_MONOTONICITY_INVARIANT = (
    "A failure may not be reclassified to a later layer once rejected."
)

__all__ = ["mandatory_pre_execution_tests", "LAYER_MONOTONICITY_INVARIANT"]
