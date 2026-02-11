"""
Test finalization spec and pipeline: each step must function correctly.
Run: python tests/test_finalization_spec.py
"""
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def test_spec_load():
    """Step 1: FINALIZATION_SPEC loads and has required keys."""
    from PX_Warehouse.Finalization_Spec import FINALIZATION_SPEC, FINALIZATION_VERSION
    assert "finalized_requires" in FINALIZATION_SPEC
    assert "zeus_gate" in FINALIZATION_SPEC
    assert "warehouse_structure" in FINALIZATION_SPEC
    assert "finalization_checklist" in FINALIZATION_SPEC
    assert "finalization_version" in FINALIZATION_SPEC
    assert FINALIZATION_SPEC["finalization_version"] == FINALIZATION_VERSION
    assert "finalization_version" in FINALIZATION_SPEC["finalized_requires"]
    assert "spec_compliance_report" in FINALIZATION_SPEC["finalized_requires"]
    assert "tier" in FINALIZATION_SPEC["finalized_requires"]
    assert FINALIZATION_SPEC["zeus_gate"]["toxicity_threshold"] == 0.0210
    print("  [PASS] FINALIZATION_SPEC loads and structure OK")


def test_zeus_gate():
    """Step 2: run_zeus_gate returns verdict with laws_required and per-law result."""
    from PX_System.foundation.ZeusLaws import run_zeus_gate
    # Pass dossier
    dossier_ok = {"engines": {"ope": {"mw": 300}, "admet": {"toxicity": {"toxicity_index": 0.01}}}, "harm_energy": 0.01, "causal_trace_log": [{}]}
    v = run_zeus_gate(dossier_ok)
    assert "authorized" in v
    assert "laws_required" in v
    assert "laws_results" in v
    assert v["authorized"] is True
    # Fail dossier (high tox)
    dossier_fail = {"engines": {"ope": {}, "admet": {"toxicity": {"toxicity_index": 0.03}}}, "harm_energy": 0.03}
    v2 = run_zeus_gate(dossier_fail)
    assert v2["authorized"] is False
    print("  [PASS] run_zeus_gate OK")


def test_finalization_pipeline_import():
    """Step 3: Finalization_Pipeline imports and run_finalization signature."""
    from PX_Warehouse.Finalization_Pipeline import run_finalization, validate_finalized_dossier, _get_spec
    spec = _get_spec()
    assert spec is not None
    print("  [PASS] Finalization_Pipeline imports OK")


def test_soc_benchmarking_score():
    """Step 3b: SoC benchmarking score computed and has expected shape."""
    from PX_Warehouse.Finalization_Pipeline import _compute_soc_benchmarking_score
    dossier = {"engines": {"ope": {"binding_affinity_nM": 8.0, "ec50": 0.01}}}
    out = _compute_soc_benchmarking_score(dossier, ["Nipah virus infection"])
    assert isinstance(out, dict)
    assert "score" in out
    assert "candidate_binding_nM" in out
    assert "soc_reference_nM" in out
    assert "better_than_soc" in out
    assert 0 <= out["score"] <= 1
    print("  [PASS] soc_benchmarking_score OK")


def test_novelty_fingerprint_score():
    """Step 3c: Novelty fingerprint score computed and has expected shape."""
    from PX_Warehouse.Finalization_Pipeline import _compute_novelty_fingerprint_score
    dossier = {"candidate": {"smiles": "CCO", "name": "ethanol"}}
    out = _compute_novelty_fingerprint_score(dossier, REPO_ROOT)
    assert isinstance(out, dict)
    assert "score" in out
    assert "reference_set_used" in out
    assert 0 <= out["score"] <= 1
    print("  [PASS] novelty_fingerprint_score OK")


def test_synthetic_accessibility_score():
    """Step 3d: Synthetic accessibility score computed and has expected shape."""
    from PX_Warehouse.Finalization_Pipeline import _compute_synthetic_accessibility_score
    dossier = {"candidate": {"smiles": "CCO"}}
    out = _compute_synthetic_accessibility_score(dossier)
    assert isinstance(out, dict)
    assert "score" in out
    assert "sa_raw" in out
    assert "source" in out
    assert 0 <= out["score"] <= 1
    assert 1 <= out["sa_raw"] <= 10
    print("  [PASS] synthetic_accessibility_score OK")


def test_finalization_version_in_spec():
    """Step 3e: Spec exposes finalization_version 1.0.0 for governance."""
    from PX_Warehouse.Finalization_Spec import FINALIZATION_SPEC, FINALIZATION_VERSION
    assert FINALIZATION_VERSION == "1.0.0"
    assert FINALIZATION_SPEC.get("finalization_version") == "1.0.0"
    print("  [PASS] finalization_version tag in spec OK")


def test_spec_compliance_report_shape():
    """Step 3f: Validator returns spec_compliance_report with required keys."""
    from PX_Warehouse.Finalization_Pipeline import _validate_finalized_dossier, _get_spec
    spec = _get_spec()
    minimal = {"finalization": {"tier": "Bronze"}, "engines": {}, "constitutional_seal": "x", "alcoa_metadata": {}, "causal_trace_log": [], "fda_compliance": "x"}
    ok, missing, report = _validate_finalized_dossier(minimal, spec)
    assert isinstance(report, dict)
    assert "missing_fields" in report
    assert "failed_checks" in report
    assert "passed_checks" in report
    assert "timestamp" in report
    assert not ok
    assert len(missing) > 0
    print("  [PASS] spec_compliance_report shape OK")


def test_finalization_on_real_dossier():
    """Step 4: Run full finalization on a real dossier (low tox) and validate."""
    from PX_Warehouse.Finalization_Pipeline import run_finalization, validate_finalized_dossier
    # Use a dossier that passes Zeus (low toxicity)
    wh = REPO_ROOT / "PX_Warehouse"
    found = None
    for folder in ["Novel_Dossiers", "Prv_Dossiers"]:
        for tier in ["Diamond", "Gold", "Silver", "Bronze"]:
            d = wh / folder / tier
            if not d.exists():
                continue
            for f in d.glob("*.json"):
                if "PATH_TEST" in f.name:
                    continue
                data = json.loads(f.read_text())
                eng = data.get("engines") or {}
                admet = eng.get("admet") or {}
                tox = admet.get("toxicity")
                if isinstance(tox, dict) and tox.get("toxicity_index") is not None:
                    if float(tox["toxicity_index"]) < 0.021:
                        found = (f.stem, data, folder == "Novel_Dossiers")
                        break
            if found:
                break
        if found:
            break
    if not found:
        print("  [SKIP] No low-tox dossier found to test full finalization")
        return
    item_id, dossier, is_novel = found
    finalized, tier, err = run_finalization(dossier, item_id, is_novel, REPO_ROOT)
    assert err is None, err
    assert finalized is not None
    assert "finalization" in finalized
    fin = finalized["finalization"]
    assert "promotion_stage" in fin
    assert "ancestral_trace" in fin
    assert "physics_map_id" in fin
    assert "trial_simulation_id" in fin
    assert "zeus_verdict" in fin
    assert "tier" in fin
    assert "rep_seal" in fin
    assert "soc_benchmarking_score" in fin
    assert "novelty_fingerprint_score" in fin
    assert "synthetic_accessibility_score" in fin
    assert "finalization_version" in fin
    assert fin["finalization_version"] == "1.0.0"
    assert "spec_compliance_report" in fin
    report = fin["spec_compliance_report"]
    assert "missing_fields" in report
    assert "failed_checks" in report
    assert "passed_checks" in report
    assert "timestamp" in report
    assert isinstance(report["passed_checks"], list)
    ok, missing = validate_finalized_dossier(finalized)
    assert ok, f"Missing: {missing}"
    print("  [PASS] Full finalization + validation OK")


def main():
    print("Finalization spec tests (step-by-step)")
    print("-" * 50)
    test_spec_load()
    test_zeus_gate()
    test_finalization_pipeline_import()
    test_soc_benchmarking_score()
    test_novelty_fingerprint_score()
    test_synthetic_accessibility_score()
    test_finalization_version_in_spec()
    test_spec_compliance_report_shape()
    test_finalization_on_real_dossier()
    print("-" * 50)
    print("All steps passed.")


if __name__ == "__main__":
    main()
