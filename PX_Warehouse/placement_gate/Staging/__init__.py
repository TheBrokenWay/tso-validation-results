"""
PX_Warehouse — canonical storage for queues, calibration data, and dossiers.

All file placement MUST go through placement_gate before being sorted.
Use PX_Warehouse.placement_gate.place_* for writes; use warehouse_layout for path resolution.

Note: placement_gate.py imports from PX_Warehouse.warehouse_layout, which re-exports
from this package's warehouse_layout. To avoid a circular import, we defer the import
of placement_gate until first access (lazy loading).
"""
from . import warehouse_layout


def __getattr__(name):
    """Lazy-load placement_gate exports to break circular import chain."""
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
        from . import placement_gate as _pg
        return getattr(_pg, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = [
    "warehouse_layout",
    "audit_placement",
    "run_placement_audit",
    "is_path_sanctioned",
    "place_queue_file",
    "place_calibration_dossier",
    "place_prv_dossier",
    "place_learning_material_file",
    "get_queue_path_for_reading",
]
