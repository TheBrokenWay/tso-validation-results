"""
Layer-monotonicity checks — generalized across all adapters.

Constitutional invariant: A failure may not be reclassified to a later layer
once rejected. Any adapter that exposes last_failure_layer can use these
helpers to assert monotonicity in tests and runners.
"""
from typing import Any, Optional

# Canonical order: earlier index = must reject first; no "reject later" drift
LAYER_ORDER = (
    "PX_Validation",   # 0 — ontology, schema, CFR, temporal
    "PX_Legal",        # 1 — PII, legal
    "PX_Engine",       # 2 — VectorCore, physics
    "Byzantium_Council",  # 3 — governance quorum
)


def layer_index(layer: Optional[str]) -> int:
    """Return index of layer in LAYER_ORDER; unknown layers return -1."""
    if not layer:
        return -1
    try:
        return LAYER_ORDER.index(layer)
    except ValueError:
        return -1


def is_later_layer(actual: Optional[str], expected: Optional[str]) -> bool:
    """
    True if actual is a later layer than expected (higher index).
    Used to detect reclassification drift.
    """
    ia, ie = layer_index(actual), layer_index(expected)
    if ia < 0 or ie < 0:
        return False
    return ia > ie


def assert_layer_monotonicity(
    adapter: Any,
    expected_layer: str,
    *,
    attribute: str = "last_failure_layer",
) -> None:
    """
    Assert that the adapter rejected at the expected layer (no later-layer drift).
    Use in poison-pill tests: failure must occur at PX_Validation (or documented layer),
    not at Engine or Council.
    """
    actual = getattr(adapter, attribute, None)
    assert actual == expected_layer, (
        f"Layer pinning violated: expected {expected_layer!r}, got {actual!r}. "
        "A failure may not be reclassified to a later layer once rejected."
    )
    if is_later_layer(actual, expected_layer):
        raise AssertionError(
            f"Layer monotonicity violated: failure at {actual!r} is later than "
            f"expected {expected_layer!r}. Ontology enforcement must not drift."
        )


def check_layer_monotonicity(adapter: Any, expected_layer: str) -> tuple[bool, Optional[str]]:
    """
    Check (don't assert) layer monotonicity. Returns (passed, error_message).
    """
    actual = getattr(adapter, "last_failure_layer", None)
    if actual != expected_layer:
        return False, f"expected layer {expected_layer!r}, got {actual!r}"
    if is_later_layer(actual, expected_layer):
        return False, f"failure at {actual!r} is later than expected {expected_layer!r}"
    return True, None


def verify_rejection_layer(
    adapter: Any,
    *,
    attribute: str = "last_failure_layer",
) -> None:
    """
    Runtime check: when an adapter rejects (returns None), it must have set
    last_failure_layer to a valid layer in LAYER_ORDER. Call this before
    returning None from mine_candidate (or equivalent) to enforce
    layer-monotonicity in the adapter, not just in tests.
    """
    layer = getattr(adapter, attribute, None)
    if not layer:
        raise AssertionError(
            "Layer pinning required: last_failure_layer must be set when rejecting. "
            "A failure may not be reclassified to a later layer once rejected."
        )
    if layer_index(layer) < 0:
        raise AssertionError(
            f"Layer {layer!r} is not in LAYER_ORDER. "
            "Rejection layer must be one of: " + ", ".join(LAYER_ORDER)
        )


__all__ = [
    "LAYER_ORDER",
    "layer_index",
    "is_later_layer",
    "assert_layer_monotonicity",
    "check_layer_monotonicity",
    "verify_rejection_layer",
]
