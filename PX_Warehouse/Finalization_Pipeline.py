"""Re-export from canonical implementation (placement_gate.Staging)."""
from PX_Warehouse.placement_gate.Staging.Finalization_Pipeline import (
    run_finalization,
    write_finalized_dossier,
    write_discovery_accepted_dossier,
    passes_discovery_bar,
    finalize_and_place,
    validate_finalized_dossier,
    _get_spec,
    _compute_soc_benchmarking_score,
    _compute_novelty_fingerprint_score,
    _compute_synthetic_accessibility_score,
    _validate_finalized_dossier,
)

__all__ = [
    "run_finalization",
    "write_finalized_dossier",
    "write_discovery_accepted_dossier",
    "passes_discovery_bar",
    "finalize_and_place",
    "validate_finalized_dossier",
]
