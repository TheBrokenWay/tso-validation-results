"""
Edge case tests: validation must FAIL on mixed strain, schema violation, CFR anomaly, temporal drift.
Run from Nipah_Analysis: pytest tests/ -v
"""
import sys
from pathlib import Path
import pytest

# Assume Nipah_Analysis is cwd or parent
SCRIPT_DIR = Path(__file__).resolve().parent.parent
if SCRIPT_DIR.name != "Nipah_Analysis":
    SCRIPT_DIR = SCRIPT_DIR / "Nipah_Analysis"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from validation.validate_raw_data import (
    _load_json_for_validation,
    validate_file_schema,
    validate_forbidden_family,
    validate_cfr_envelope,
    validate_temporal_consistency,
    RawValidationError,
)
from validation.validate_manifest import load_and_validate_manifest, ManifestValidationError

CONFIG_DIR = SCRIPT_DIR / "config"
EDGE_DIR = SCRIPT_DIR / "data" / "raw" / "edge_cases"


def test_manifest_strain_split():
    manifest = load_and_validate_manifest(CONFIG_DIR / "nipah_manifest.json")
    assert "strains" in manifest["pathogen"]
    assert "NiV_Malaysia" in manifest["pathogen"]["strains"]
    assert "NiV_Bangladesh" in manifest["pathogen"]["strains"]


def test_edge_mixed_strain_fails():
    path = EDGE_DIR / "edge_mixed_strain.json"
    if not path.exists():
        pytest.skip("edge_mixed_strain.json not found")
    data = _load_json_for_validation(path)
    with pytest.raises(RawValidationError):
        validate_file_schema(data, path)


def test_edge_schema_violation_fails():
    path = EDGE_DIR / "edge_schema_violation.json"
    if not path.exists():
        pytest.skip("edge_schema_violation.json not found")
    data = _load_json_for_validation(path)
    with pytest.raises(RawValidationError):
        validate_file_schema(data, path)


def test_edge_cfr_anomaly_fails():
    path = EDGE_DIR / "edge_cfr_anomaly.json"
    if not path.exists():
        pytest.skip("edge_cfr_anomaly.json not found")
    data = _load_json_for_validation(path)
    validate_file_schema(data, path)
    with pytest.raises(RawValidationError):
        validate_cfr_envelope(data, path)


def test_edge_temporal_drift_fails():
    path = EDGE_DIR / "edge_temporal_drift.json"
    if not path.exists():
        pytest.skip("edge_temporal_drift.json not found")
    data = _load_json_for_validation(path)
    validate_file_schema(data, path)
    validate_cfr_envelope(data, path)
    with pytest.raises(RawValidationError):
        validate_temporal_consistency(data, path)


def test_valid_malaysia_passes():
    path = SCRIPT_DIR / "data" / "raw" / "valid_malaysia.json"
    if not path.exists():
        pytest.skip("valid_malaysia.json not found")
    data = _load_json_for_validation(path)
    validate_forbidden_family(data, path)
    validate_file_schema(data, path)
    validate_cfr_envelope(data, path)
    validate_temporal_consistency(data, path)


def test_edge_forbidden_family_fails():
    """Forbidden-family tripwire: Orthomyxoviridae (influenza) must FAIL."""
    path = EDGE_DIR / "edge_forbidden_family.json"
    if not path.exists():
        pytest.skip("edge_forbidden_family.json not found")
    data = _load_json_for_validation(path)
    with pytest.raises(RawValidationError) as exc:
        validate_forbidden_family(data, path)
    assert "FORBIDDEN_FAMILY" in str(exc.value) or "Orthomyxoviridae" in str(exc.value)


def test_strain_coverage_assertion_multi_strain_only():
    """When allow_multi_strain=True, full coverage required: missing Bangladesh data must FAIL."""
    from validation.validate_raw_data import validate_all_raw
    manifest = load_and_validate_manifest(CONFIG_DIR / "nipah_manifest.json")
    import tempfile
    import shutil
    with tempfile.TemporaryDirectory() as td:
        raw_dir = Path(td)
        valid_m = SCRIPT_DIR / "data" / "raw" / "valid_malaysia.json"
        if not valid_m.exists():
            pytest.skip("valid_malaysia.json not found")
        shutil.copy(valid_m, raw_dir / "valid_malaysia.json")
        with pytest.raises(RawValidationError) as exc:
            validate_all_raw(
                raw_dir,
                manifest,
                fail_on_no_files=True,
                allow_multi_strain=True,
                require_strain_coverage=True,
            )
        assert "STRAIN_COVERAGE" in str(exc.value) or "no valid raw files" in str(exc.value).lower()


def test_single_strain_mode_no_coverage_requirement():
    """Single-strain mode: only Malaysia data → PASS (no requirement for Bangladesh)."""
    from validation.validate_raw_data import validate_all_raw
    manifest = load_and_validate_manifest(CONFIG_DIR / "nipah_manifest.json")
    import tempfile
    import shutil
    with tempfile.TemporaryDirectory() as td:
        raw_dir = Path(td)
        valid_m = SCRIPT_DIR / "data" / "raw" / "valid_malaysia.json"
        if not valid_m.exists():
            pytest.skip("valid_malaysia.json not found")
        shutil.copy(valid_m, raw_dir / "valid_malaysia.json")
        accepted, _, drift = validate_all_raw(
            raw_dir,
            manifest,
            fail_on_no_files=True,
            allow_multi_strain=False,
            require_strain_coverage=False,
        )
        assert len(accepted) == 1
        assert accepted[0]["strain"] == "NiV_Malaysia"
