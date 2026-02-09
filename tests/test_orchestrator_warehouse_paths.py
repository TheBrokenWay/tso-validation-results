"""
Verify every orchestrator and pipeline writes to the correct warehouse folders.
Run from repo root: python tests/test_orchestrator_warehouse_paths.py
"""
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

# Canonical warehouse layout: Feeder, Calibration_Molecules, Prv_Dossiers, Novel_Dossiers, Learning_Material
WAREHOUSE = ROOT / "PX_Warehouse"
CANONICAL = {
    "Feeder": WAREHOUSE / "Feeder",
    "Calibration_Molecules": WAREHOUSE / "Calibration_Molecules",
    "Calibration_LiveRuns": WAREHOUSE / "Calibration_Molecules" / "LiveRuns",
    "Calibration_BatchRuns": WAREHOUSE / "Calibration_Molecules" / "BatchRuns",
    "Prv_Dossiers": WAREHOUSE / "Prv_Dossiers",
    "Novel_Dossiers": WAREHOUSE / "Novel_Dossiers",
    "Learning_Material": WAREHOUSE / "Learning_Material",
    "Operations_Inputs": WAREHOUSE / "Operations" / "Inputs",
    "WorldLines": WAREHOUSE / "WorldLines",
}


def test_px_live_orchestrator_v2_paths():
    """PX_Live_Orchestrator_v2 must write to PX_Warehouse/Calibration_Molecules/LiveRuns/run_<timestamp>/"""
    from PX_Executive.orchestrators.PX_Live_Orchestrator_v2 import _REPO_ROOT
    from PX_Warehouse.warehouse_layout import get_calibration_molecules_dir
    cal_root = get_calibration_molecules_dir(_REPO_ROOT)
    assert _REPO_ROOT == ROOT, f"_REPO_ROOT should be repo root, got {_REPO_ROOT}"
    assert cal_root == CANONICAL["Calibration_Molecules"], "Calibration_Molecules path mismatch"
    assert (cal_root / "LiveRuns").resolve().is_relative_to(ROOT.resolve()), "LiveRuns must be under repo"


def test_prv_24h_orchestrator_paths():
    """PRV_24H_Orchestrator must write dossiers to Prv_Dossiers/<tier> or Novel_Dossiers/<tier>"""
    from PX_Executive.PRV_24H_Orchestrator import REPO_ROOT
    from PX_Warehouse.warehouse_layout import get_prv_dossier_dir
    assert REPO_ROOT == ROOT, f"REPO_ROOT should be repo root, got {REPO_ROOT}"
    out_rep = get_prv_dossier_dir(False, "Gold", REPO_ROOT)
    out_nov = get_prv_dossier_dir(True, "Gold", REPO_ROOT)
    assert out_rep == WAREHOUSE / "Prv_Dossiers" / "Gold"
    assert out_nov == WAREHOUSE / "Novel_Dossiers" / "Gold"


def test_sovereign_commercial_pipeline_paths():
    """Sovereign_Commercial_Pipeline writes to Prv_Dossiers/<tier> or Novel_Dossiers/<tier>"""
    from PX_Executive.Sovereign_Commercial_Pipeline import _REPO_ROOT
    from PX_Warehouse.warehouse_layout import get_prv_dossier_dir
    assert _REPO_ROOT == ROOT
    gold_rep = get_prv_dossier_dir(False, "Gold", _REPO_ROOT)
    gold_nov = get_prv_dossier_dir(True, "Gold", _REPO_ROOT)
    assert gold_rep == WAREHOUSE / "Prv_Dossiers" / "Gold"
    assert gold_nov == WAREHOUSE / "Novel_Dossiers" / "Gold"


def test_evidence_package_wrap_trial_simulation_default():
    """wrap_trial_simulation default output_dir is Calibration_Molecules (QA vault)"""
    from PX_System.foundation.Evidence_Package import wrap_trial_simulation
    import inspect
    sig = inspect.signature(wrap_trial_simulation)
    default = sig.parameters["output_dir"].default
    assert default == "PX_Warehouse/Calibration_Molecules", f"Expected default PX_Warehouse/Calibration_Molecules, got {default}"


def test_batch_runner_paths():
    """pipeline_batch_runner writes to PX_Warehouse/Calibration_Molecules/BatchRuns"""
    from PX_Executive.batch import pipeline_batch_runner
    from PX_Warehouse.warehouse_layout import get_calibration_molecules_dir
    project_root = Path(pipeline_batch_runner.__file__).resolve().parents[2]
    output_dir = get_calibration_molecules_dir(project_root) / "BatchRuns"
    assert project_root == ROOT
    assert output_dir == WAREHOUSE / "Calibration_Molecules" / "BatchRuns"


def test_universal_pipeline_runner_resolves_warehouse():
    """UniversalPipelineRunner must resolve relative warehouse_root to repo root"""
    from PX_Executive.UniversalPipelineRunner import UniversalPipelineRunner, _REPO_ROOT
    r = UniversalPipelineRunner(warehouse_root="PX_Warehouse")
    assert r.warehouse_root.is_absolute(), "warehouse_root must be absolute"
    assert r.warehouse_root.resolve() == (ROOT / "PX_Warehouse").resolve()
    assert r.sorted_root.resolve().is_relative_to(ROOT.resolve())


def test_placement_gate_audit():
    """Placement gate audit must report 0 misplaced files (all content in sanctioned zones)."""
    from PX_Warehouse.placement_gate import run_placement_audit
    n = run_placement_audit(ROOT, verbose=False)
    assert n == 0, f"Placement audit found {n} misplaced files; all warehouse content must go through placement_gate zones"


def test_canonical_folders_exist_or_creatable():
    """All canonical paths must be under repo and creatable"""
    for name, path in CANONICAL.items():
        path = Path(path)
        assert path.resolve().is_relative_to(ROOT.resolve()), f"{name} must be under repo: {path}"
        path.mkdir(parents=True, exist_ok=True)
        assert path.exists(), f"{name} must be creatable: {path}"


def run_all():
    errors = []
    tests = [
        ("PX_Live_Orchestrator_v2 paths", test_px_live_orchestrator_v2_paths),
        ("PRV_24H_Orchestrator paths", test_prv_24h_orchestrator_paths),
        ("Sovereign_Commercial_Pipeline paths", test_sovereign_commercial_pipeline_paths),
        ("Evidence_Package default output_dir", test_evidence_package_wrap_trial_simulation_default),
        ("pipeline_batch_runner paths", test_batch_runner_paths),
        ("UniversalPipelineRunner warehouse resolution", test_universal_pipeline_runner_resolves_warehouse),
        ("Canonical folders exist/creatable", test_canonical_folders_exist_or_creatable),
        ("Placement gate audit (no misplaced files)", test_placement_gate_audit),
    ]
    for name, fn in tests:
        try:
            fn()
            print(f"  OK  {name}")
        except Exception as e:
            print(f"  FAIL  {name}: {e}")
            errors.append(f"{name}: {e}")
    if errors:
        print("\nFAILED:", errors)
        return 1
    print("\nAll orchestrator/warehouse path checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(run_all())
