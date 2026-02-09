# placement_gate: Staging scripts and gate logic
#
# Re-exports from Staging.placement_gate are lazy-loaded to avoid circular
# imports (placement_gate.py -> PX_Warehouse.warehouse_layout -> this package).


def __getattr__(name):
    """Lazy re-export from Staging.placement_gate."""
    _pg_names = {
        "audit_placement",
        "run_placement_audit",
        "is_path_sanctioned",
        "place_queue_file",
        "place_calibration_dossier",
        "place_prv_dossier",
        "place_learning_material_file",
        "get_queue_path_for_reading",
    }
    if name in _pg_names:
        from .Staging import placement_gate as _pg
        return getattr(_pg, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = [
    "audit_placement",
    "run_placement_audit",
    "is_path_sanctioned",
    "place_queue_file",
    "place_calibration_dossier",
    "place_prv_dossier",
    "place_learning_material_file",
    "get_queue_path_for_reading",
]
