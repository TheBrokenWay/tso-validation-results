"""
test_quint_engine_adapter.py
Tests for the QUINT Engine Adapter layer.

Validates that all 20 engine wrappers (12 PRV + 8 auxiliary):
  - Accept both QFrame and plain dict inputs (backward compatibility)
  - Return QFrame with QType.QRESULT
  - Preserve engine_id in payload
  - Track lineage with engine name
"""

import os
import sys
import unittest
import unittest.mock
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from PX_System.foundation.quint.kernel import QFrame, QType, create_qframe
from PX_System.foundation.quint.converter import ingest
from PX_System.foundation.quint.engine_adapter import (
    # PRV pipeline (1-12)
    q_run_ope,
    q_run_obe,
    q_run_oce,
    q_run_ole,
    q_run_ome,
    q_run_ose,
    q_run_admet,
    q_run_pkpd,
    q_run_dose_optimizer,
    q_run_virtual_efficacy,
    q_run_grading,
    q_run_zeus,
    # Auxiliary (13-20)
    q_run_vector_core,
    q_run_simulation,
    q_run_trial,
    q_run_check_constitutional,
    q_run_evidence_dossier,
    q_run_wrap_trial,
    q_run_metabolism_pulse,
    q_run_block_pulse,
    _unwrap,
)


# ---------------------------------------------------------------------------
# Test Data
# ---------------------------------------------------------------------------

TEST_SMILES = "CCO"  # Ethanol — simple, valid SMILES

# Law U34 compliant p_vector: [complexity, energy, dims, valid]
# 4 elements expected by BlockUniverse.project_proposal
TEST_P_VECTOR = [0.1, 0.0, 35.0, 1.0]
TEST_CSA_SCORES = [1.0, 1.0, 1.0, 1.0, 1.0]


def _minimal_dossier():
    """Build a minimal dossier dict with fields needed by GradingEngine and ZeusLaws."""
    return {
        "compound_id": "TEST-CCO-001",
        "smiles": TEST_SMILES,
        "harm_energy": 0.005,
        "causal_trace_log": ["genesis", "ope", "admet"],
        "engines": {
            "ope": {
                "molecular_weight": 46.07,
                "logp": -0.31,
                "hbd": 1,
                "hba": 1,
                "tpsa": 20.23,
                "ec50": 10.0,
                "emax": 0.8,
                "clearance_estimate_L_per_h": 2.0,
                "vd_estimate_L": 40.0,
                "binding_affinity_nM": 8.0,
            },
            "admet": {
                "absorption": {"oral_bioavailability_pct": 85.0},
                "toxicity": {
                    "toxicity_index": 0.005,
                    "risk_level": "TOXICITY_DIAMOND",
                    "safety_margin": 60.0,
                },
                "safety_margins": {"safety_margin": 60.0},
            },
            "obe": {"harm_energy": 0.005},
            "oce": {"coherence": 0.95},
            "ole": {"prv_eligible": True},
            "ome": {"clearance_L_per_h": 2.0, "half_life_h": 24.0},
            "ose": {"toxicity_index": 0.005},
            "virtual_efficacy": {
                "pta": 0.85,
                "responder_rate": 0.75,
                "variability_cv": 20.0,
            },
            "dose_optimizer": {"best_regimen": {"dose_mg": 100, "score": 0.9}},
        },
        "trial_result": {
            "arms": [
                {
                    "arm_id": "active",
                    "patients_enrolled": 21,
                    "exposure_summary": {
                        "auc_mg_h_per_L": {"mean": 300.0, "std": 50.0},
                    },
                    "pd_summary": {
                        "max_effect": {"mean": 0.7, "std": 0.1},
                    },
                }
            ],
        },
    }


# ---------------------------------------------------------------------------
# Tests: _unwrap helper
# ---------------------------------------------------------------------------

class TestUnwrap(unittest.TestCase):
    """Tests for the _unwrap helper function."""

    def test_unwrap_qframe(self):
        """Unwrapping a QFrame returns its dict payload."""
        frame = ingest({"smiles": TEST_SMILES, "engine_id": "test", "status": "OK"})
        result = _unwrap(frame)
        self.assertIsInstance(result, dict)
        self.assertEqual(result.get("smiles"), TEST_SMILES)

    def test_unwrap_dict(self):
        """Unwrapping a dict returns the same dict."""
        d = {"smiles": TEST_SMILES}
        result = _unwrap(d)
        self.assertIs(result, d)

    def test_unwrap_invalid(self):
        """Unwrapping a non-dict/QFrame raises TypeError."""
        with self.assertRaises(TypeError):
            _unwrap(42)
        with self.assertRaises(TypeError):
            _unwrap("not a dict")
        with self.assertRaises(TypeError):
            _unwrap([1, 2, 3])


# ---------------------------------------------------------------------------
# Tests: OPE
# ---------------------------------------------------------------------------

class TestQRunOPE(unittest.TestCase):
    """Tests for the q_run_ope adapter."""

    def test_ope_with_dict(self):
        """q_run_ope with a plain dict returns a QFrame with QRESULT type."""
        frame = q_run_ope({"smiles": TEST_SMILES})
        self.assertIsInstance(frame, QFrame)
        self.assertEqual(frame.qtype, QType.QRESULT)

    def test_ope_with_qframe(self):
        """q_run_ope with a QFrame input returns QFrame with QRESULT type."""
        input_frame = ingest({"smiles": TEST_SMILES})
        frame = q_run_ope(input_frame)
        self.assertIsInstance(frame, QFrame)
        self.assertEqual(frame.qtype, QType.QRESULT)

    def test_ope_result_has_engine_id(self):
        """OPE result payload contains engine_id."""
        frame = q_run_ope({"smiles": TEST_SMILES})
        self.assertEqual(frame.payload.get("engine_id"), "OPE_V3_DETERMINISTIC")

    def test_ope_result_has_molecular_weight(self):
        """OPE result payload contains molecular_weight from RDKit."""
        frame = q_run_ope({"smiles": TEST_SMILES})
        self.assertIn("molecular_weight", frame.payload)
        self.assertIsInstance(frame.payload["molecular_weight"], float)


# ---------------------------------------------------------------------------
# Tests: OBE
# ---------------------------------------------------------------------------

class TestQRunOBE(unittest.TestCase):
    """Tests for the q_run_obe adapter."""

    def test_obe_returns_qframe(self):
        """q_run_obe with dict returns a QFrame."""
        frame = q_run_obe({"smiles": TEST_SMILES})
        self.assertIsInstance(frame, QFrame)
        self.assertEqual(frame.qtype, QType.QRESULT)
        self.assertEqual(frame.payload.get("engine_id"), "OBE_V3_DETERMINISTIC")


# ---------------------------------------------------------------------------
# Tests: OCE
# ---------------------------------------------------------------------------

class TestQRunOCE(unittest.TestCase):
    """Tests for the q_run_oce adapter."""

    def test_oce_returns_qframe(self):
        """q_run_oce with p_vector and csa_scores returns a QFrame."""
        frame = q_run_oce({
            "p_vector": TEST_P_VECTOR,
            "csa_scores": TEST_CSA_SCORES,
        })
        self.assertIsInstance(frame, QFrame)
        self.assertEqual(frame.qtype, QType.QRESULT)
        self.assertEqual(frame.payload.get("engine_id"), "OCE_V3_DETERMINISTIC")


# ---------------------------------------------------------------------------
# Tests: OLE
# ---------------------------------------------------------------------------

class TestQRunOLE(unittest.TestCase):
    """Tests for the q_run_ole adapter."""

    def test_ole_returns_qframe(self):
        """q_run_ole with compound_id and indication returns a QFrame."""
        frame = q_run_ole({
            "compound_id": "TEST-001",
            "indication": "nipah",
            "disease_context": ["nipah virus infection"],
        })
        self.assertIsInstance(frame, QFrame)
        self.assertEqual(frame.qtype, QType.QRESULT)
        self.assertEqual(frame.payload.get("engine_id"), "OLE_V3_DETERMINISTIC")


# ---------------------------------------------------------------------------
# Tests: OME
# ---------------------------------------------------------------------------

class TestQRunOME(unittest.TestCase):
    """Tests for the q_run_ome adapter."""

    def test_ome_returns_qframe(self):
        """q_run_ome with dict returns a QFrame."""
        frame = q_run_ome({"smiles": TEST_SMILES})
        self.assertIsInstance(frame, QFrame)
        self.assertEqual(frame.qtype, QType.QRESULT)
        self.assertEqual(frame.payload.get("engine_id"), "OME_V3_DETERMINISTIC")


# ---------------------------------------------------------------------------
# Tests: OSE
# ---------------------------------------------------------------------------

class TestQRunOSE(unittest.TestCase):
    """Tests for the q_run_ose adapter."""

    def test_ose_returns_qframe(self):
        """q_run_ose with dict returns a QFrame."""
        frame = q_run_ose({"smiles": TEST_SMILES})
        self.assertIsInstance(frame, QFrame)
        self.assertEqual(frame.qtype, QType.QRESULT)
        self.assertEqual(frame.payload.get("engine_id"), "OSE_V3_DETERMINISTIC")


# ---------------------------------------------------------------------------
# Tests: ADMET
# ---------------------------------------------------------------------------

class TestQRunADMET(unittest.TestCase):
    """Tests for the q_run_admet adapter."""

    def test_admet_returns_qframe(self):
        """q_run_admet with smiles and ope_result returns a QFrame."""
        ope_result = {
            "molecular_weight": 46.07,
            "logp": -0.31,
            "hbd": 1,
            "hba": 1,
            "tpsa": 20.23,
            "clearance_estimate_L_per_h": 2.0,
        }
        frame = q_run_admet({"smiles": TEST_SMILES}, ope_result=ope_result)
        self.assertIsInstance(frame, QFrame)
        self.assertEqual(frame.qtype, QType.QRESULT)
        self.assertEqual(frame.payload.get("engine_id"), "ADMET_V3_DETERMINISTIC")

    def test_admet_accepts_qframe_inputs(self):
        """q_run_admet accepts QFrame for all inputs."""
        smiles_frame = ingest({"smiles": TEST_SMILES})
        ope_frame = ingest({
            "molecular_weight": 46.07,
            "logp": -0.31,
            "hbd": 1,
            "hba": 1,
            "tpsa": 20.23,
            "clearance_estimate_L_per_h": 2.0,
            "engine_id": "OPE",
            "status": "PASSED",
        })
        frame = q_run_admet(smiles_frame, ope_result=ope_frame)
        self.assertIsInstance(frame, QFrame)
        self.assertEqual(frame.qtype, QType.QRESULT)


# ---------------------------------------------------------------------------
# Tests: GradingEngine
# ---------------------------------------------------------------------------

class TestQRunGrading(unittest.TestCase):
    """Tests for the q_run_grading adapter."""

    def test_grading_returns_qframe(self):
        """q_run_grading with a minimal dossier returns a QFrame."""
        dossier = _minimal_dossier()
        frame = q_run_grading(dossier)
        self.assertIsInstance(frame, QFrame)
        self.assertEqual(frame.qtype, QType.QRESULT)
        self.assertEqual(frame.payload.get("engine_id"), "GRADING_ENGINE_V2")

    def test_grading_accepts_qframe(self):
        """q_run_grading accepts a QFrame wrapping a dossier."""
        dossier = _minimal_dossier()
        dossier_frame = ingest(dossier)
        frame = q_run_grading(dossier_frame)
        self.assertIsInstance(frame, QFrame)
        self.assertEqual(frame.qtype, QType.QRESULT)


# ---------------------------------------------------------------------------
# Tests: ZeusLaws
# ---------------------------------------------------------------------------

class TestQRunZeus(unittest.TestCase):
    """Tests for the q_run_zeus adapter."""

    def test_zeus_returns_qframe(self):
        """q_run_zeus with a minimal dossier returns a QFrame."""
        dossier = _minimal_dossier()
        frame = q_run_zeus(dossier)
        self.assertIsInstance(frame, QFrame)
        self.assertEqual(frame.qtype, QType.QRESULT)
        self.assertEqual(frame.payload.get("engine_id"), "ZEUS_GATE_V1")

    def test_zeus_authorized_field(self):
        """ZeusLaws result should contain authorized field."""
        dossier = _minimal_dossier()
        frame = q_run_zeus(dossier)
        self.assertIn("authorized", frame.payload)

    def test_zeus_accepts_qframe(self):
        """q_run_zeus accepts a QFrame wrapping a dossier."""
        dossier = _minimal_dossier()
        dossier_frame = ingest(dossier)
        frame = q_run_zeus(dossier_frame)
        self.assertIsInstance(frame, QFrame)


# ===========================================================================
# AUXILIARY ENGINE TESTS (13-20)
# ===========================================================================

# ---------------------------------------------------------------------------
# Test Data — Auxiliary engines
# ---------------------------------------------------------------------------

TEST_PROTOCOL = {
    "trial_id": "TRIAL-QUINT-001",
    "duration_days": 7.0,
    "design": "PARALLEL",
    "arms": [
        {
            "arm_id": "ARM-A",
            "label": "Active",
            "dose_mg": 100.0,
            "dosing_interval_h": 24.0,
            "n_patients": 5,
        }
    ],
}

TEST_ADMET_FULL = {
    "absorption": {"predicted_bioavailability": 1.0, "oral_bioavailability_pct": 85.0},
    "distribution": {"predicted_vd_L_per_kg": 0.7},
    "metabolism": {"predicted_clearance_L_per_h_per_kg": 0.05},
    "toxicity": {"toxicity_index": 0.005, "risk_level": "TOXICITY_DIAMOND", "safety_margin": 60.0},
    "safety_margins": {"safety_margin": 60.0},
}


def _minimal_candidate():
    """Build a minimal candidate dict for Evidence_Package."""
    return {
        "compound_id": "TEST-CCO-002",
        "name": "Ethanol-QUINT",
        "smiles": TEST_SMILES,
    }


def _minimal_engine_results():
    """Build a minimal engine_results dict for Evidence_Package."""
    return {
        "ope": {
            "molecular_weight": 46.07,
            "logp": -0.31,
            "hbd": 1,
            "hba": 1,
            "tpsa": 20.23,
            "ec50": 10.0,
            "emax": 0.8,
            "clearance_estimate_L_per_h": 2.0,
            "vd_estimate_L": 40.0,
            "binding_affinity_nM": 8.0,
        },
        "admet": TEST_ADMET_FULL,
        "obe": {"harm_energy": 0.005},
        "oce": {"coherence": 0.95},
        "ole": {"prv_eligible": True},
        "ome": {"clearance_L_per_h": 2.0, "half_life_h": 24.0},
        "ose": {"toxicity_index": 0.005},
    }


# ---------------------------------------------------------------------------
# Tests: VectorCore (13)
# ---------------------------------------------------------------------------

class TestQRunVectorCore(unittest.TestCase):
    """Tests for the q_run_vector_core adapter."""

    def test_vector_core_with_dict(self):
        """q_run_vector_core with a dict returns a QFrame."""
        frame = q_run_vector_core({"p_vector": TEST_P_VECTOR})
        self.assertIsInstance(frame, QFrame)
        self.assertEqual(frame.qtype, QType.QRESULT)
        self.assertEqual(frame.payload.get("engine_id"), "VECTOR_CORE_V1")

    def test_vector_core_with_list(self):
        """q_run_vector_core accepts a plain list as p_vector."""
        frame = q_run_vector_core(TEST_P_VECTOR)
        self.assertIsInstance(frame, QFrame)
        self.assertEqual(frame.payload.get("engine_id"), "VECTOR_CORE_V1")

    def test_vector_core_with_qframe(self):
        """q_run_vector_core accepts a QFrame input."""
        input_frame = ingest({"p_vector": TEST_P_VECTOR})
        frame = q_run_vector_core(input_frame)
        self.assertIsInstance(frame, QFrame)
        self.assertEqual(frame.qtype, QType.QRESULT)

    def test_vector_core_has_authorized(self):
        """VectorCore result should contain authorized field."""
        frame = q_run_vector_core({"p_vector": TEST_P_VECTOR})
        self.assertIn("authorized", frame.payload)


# ---------------------------------------------------------------------------
# Tests: SimulationEngine (14)
# ---------------------------------------------------------------------------

def _create_temp_worldline(tmpdir, task_num="QSIM01"):
    """Create a minimal .worldline file for SimulationEngine tests."""
    import json
    wl_file = os.path.join(tmpdir, f"{task_num}.worldline")
    wl_data = {
        "worldline_id": f"WL-{task_num}",
        "physical_realization": {
            "toxicity_index": 0.005,
            "binding_affinity_kj": -35.0,
        },
    }
    with open(wl_file, "w") as f:
        json.dump(wl_data, f)
    return f"WL-{task_num}"


class TestQRunSimulation(unittest.TestCase):
    """Tests for the q_run_simulation adapter."""

    def setUp(self):
        import tempfile
        self._tmpdir = tempfile.mkdtemp()
        # Patch DEFAULT_WORLDLINES_PATH at the module level that materialize_candidate imports from
        self._patcher = unittest.mock.patch(
            "PX_Warehouse.WorldLine_Database.DEFAULT_WORLDLINES_PATH",
            self._tmpdir,
        )
        self._patcher.start()

    def tearDown(self):
        self._patcher.stop()
        import shutil
        shutil.rmtree(self._tmpdir, ignore_errors=True)

    def test_simulation_with_dict(self):
        """q_run_simulation with dict returns a QFrame."""
        wl_id = _create_temp_worldline(self._tmpdir, "QSIM01")
        frame = q_run_simulation({"task_id": wl_id, "coherence": 0.95})
        self.assertIsInstance(frame, QFrame)
        self.assertEqual(frame.qtype, QType.QRESULT)
        self.assertEqual(frame.payload.get("engine_id"), "SIMULATION_ENGINE_V1")

    def test_simulation_with_qframe(self):
        """q_run_simulation accepts QFrame input."""
        wl_id = _create_temp_worldline(self._tmpdir, "QSIM02")
        input_frame = ingest({"task_id": wl_id, "coherence": 0.90})
        frame = q_run_simulation(input_frame)
        self.assertIsInstance(frame, QFrame)
        self.assertEqual(frame.qtype, QType.QRESULT)

    def test_simulation_materialization_status(self):
        """SimulationEngine result should contain status field."""
        wl_id = _create_temp_worldline(self._tmpdir, "QSIM03")
        frame = q_run_simulation({"task_id": wl_id, "coherence": 0.95})
        self.assertIn("status", frame.payload)

    def test_simulation_void_low_coherence(self):
        """SimulationEngine returns VOID for coherence < 0.80."""
        frame = q_run_simulation({"task_id": "WL-VOID", "coherence": 0.50})
        self.assertIsInstance(frame, QFrame)
        self.assertEqual(frame.payload.get("status"), "VOID")


# ---------------------------------------------------------------------------
# Tests: TrialEngine (15)
# ---------------------------------------------------------------------------

class TestQRunTrial(unittest.TestCase):
    """Tests for the q_run_trial adapter."""

    def test_trial_with_dict(self):
        """q_run_trial with dict returns a QFrame."""
        frame = q_run_trial(TEST_PROTOCOL, admet=TEST_ADMET_FULL)
        self.assertIsInstance(frame, QFrame)
        self.assertEqual(frame.qtype, QType.QRESULT)
        self.assertEqual(frame.payload.get("engine_id"), "TRIAL_ENGINE_V1")

    def test_trial_with_qframe(self):
        """q_run_trial accepts QFrame inputs."""
        proto_frame = ingest(TEST_PROTOCOL)
        admet_frame = ingest({**TEST_ADMET_FULL, "engine_id": "ADMET", "status": "OK"})
        frame = q_run_trial(proto_frame, admet=admet_frame)
        self.assertIsInstance(frame, QFrame)
        self.assertEqual(frame.qtype, QType.QRESULT)

    def test_trial_has_arms(self):
        """Trial result should contain arms data."""
        frame = q_run_trial(TEST_PROTOCOL, admet=TEST_ADMET_FULL)
        self.assertIn("arms", frame.payload)


# ---------------------------------------------------------------------------
# Tests: CheckConstitutional (16)
# ---------------------------------------------------------------------------

class TestQRunCheckConstitutional(unittest.TestCase):
    """Tests for the q_run_check_constitutional adapter."""

    def test_constitutional_with_safe_data(self):
        """Constitutional check with safe values returns QFrame."""
        data = {"toxicity_index": 0.005, "harm_energy": 0.005}
        frame = q_run_check_constitutional("test_operation", data)
        self.assertIsInstance(frame, QFrame)
        self.assertEqual(frame.qtype, QType.QRESULT)
        self.assertEqual(frame.payload.get("engine_id"), "CONSTITUTIONAL_CHECK_V1")

    def test_constitutional_with_qframe(self):
        """Constitutional check accepts QFrame data."""
        data_frame = ingest({"toxicity_index": 0.005, "harm_energy": 0.005})
        frame = q_run_check_constitutional("test_qframe", data_frame)
        self.assertIsInstance(frame, QFrame)

    def test_constitutional_authorized_field(self):
        """Constitutional result should contain authorized field."""
        data = {"toxicity_index": 0.005, "harm_energy": 0.005}
        frame = q_run_check_constitutional("test_auth", data)
        self.assertIn("authorized", frame.payload)


# ---------------------------------------------------------------------------
# Tests: EvidencePackage (17)
# ---------------------------------------------------------------------------

class TestQRunEvidenceDossier(unittest.TestCase):
    """Tests for the q_run_evidence_dossier adapter."""

    def test_evidence_with_dict(self):
        """Evidence dossier generation with dict returns QFrame."""
        candidate = _minimal_candidate()
        eng = _minimal_engine_results()
        zeus = {"authorized": True, "laws_checked": ["L1", "U27"]}
        frame = q_run_evidence_dossier(candidate, engine_results=eng, zeus_verdict=zeus)
        self.assertIsInstance(frame, QFrame)
        self.assertEqual(frame.qtype, QType.QRESULT)
        self.assertEqual(frame.payload.get("engine_id"), "EVIDENCE_PACKAGE_V1")

    def test_evidence_with_qframe(self):
        """Evidence dossier accepts QFrame inputs."""
        cand_frame = ingest(_minimal_candidate())
        eng_frame = ingest({**_minimal_engine_results(), "status": "OK"})
        frame = q_run_evidence_dossier(cand_frame, engine_results=eng_frame)
        self.assertIsInstance(frame, QFrame)


# ---------------------------------------------------------------------------
# Tests: WrapTrialSimulation (18)
# ---------------------------------------------------------------------------

class TestQRunWrapTrial(unittest.TestCase):
    """Tests for the q_run_wrap_trial adapter."""

    def test_wrap_trial_with_dict(self):
        """Wrap trial simulation with dicts returns QFrame."""
        trial_result = {
            "trial_id": "TRIAL-QUINT-001",
            "arms": [{"arm_id": "ARM-A", "patients_enrolled": 5,
                       "exposure_summary": {"auc_mg_h_per_L": {"mean": 300.0}}}],
        }
        ope = {"molecular_weight": 46.07, "logp": -0.31, "ec50": 10.0, "emax": 0.8}
        import tempfile, os
        out_dir = tempfile.mkdtemp()
        frame = q_run_wrap_trial(
            TEST_PROTOCOL, trial_result=trial_result,
            ope_result=ope, admet_result=TEST_ADMET_FULL,
            output_dir=out_dir,
        )
        self.assertIsInstance(frame, QFrame)
        self.assertEqual(frame.qtype, QType.QRESULT)
        self.assertEqual(frame.payload.get("engine_id"), "WRAP_TRIAL_V1")


# ---------------------------------------------------------------------------
# Tests: Metabolism (19)
# ---------------------------------------------------------------------------

class TestQRunMetabolismPulse(unittest.TestCase):
    """Tests for the q_run_metabolism_pulse adapter."""

    def test_metabolism_with_dict(self):
        """Metabolism pulse with dict returns a QFrame."""
        frame = q_run_metabolism_pulse({
            "task_id": "METAB-QUINT-001",
            "p_vec": TEST_P_VECTOR,
            "csa_s": TEST_CSA_SCORES,
        })
        self.assertIsInstance(frame, QFrame)
        self.assertEqual(frame.qtype, QType.QRESULT)
        self.assertEqual(frame.payload.get("engine_id"), "METABOLISM_V1")

    def test_metabolism_with_string(self):
        """Metabolism pulse accepts a plain string as task_id."""
        frame = q_run_metabolism_pulse("METAB-STRING-001")
        self.assertIsInstance(frame, QFrame)
        self.assertEqual(frame.payload.get("engine_id"), "METABOLISM_V1")

    def test_metabolism_with_qframe(self):
        """Metabolism pulse accepts QFrame input."""
        input_frame = ingest({
            "task_id": "METAB-QF-001",
            "p_vec": TEST_P_VECTOR,
            "csa_s": TEST_CSA_SCORES,
        })
        frame = q_run_metabolism_pulse(input_frame)
        self.assertIsInstance(frame, QFrame)

    def test_metabolism_has_cycle_age(self):
        """Metabolism result should contain cycle_age from pulse tuple."""
        frame = q_run_metabolism_pulse({
            "task_id": "METAB-AGE-001",
            "p_vec": TEST_P_VECTOR,
            "csa_s": TEST_CSA_SCORES,
        })
        self.assertIn("cycle_age", frame.payload)
        self.assertIsInstance(frame.payload["cycle_age"], int)


# ---------------------------------------------------------------------------
# Tests: BlockOrchestrator (20)
# ---------------------------------------------------------------------------

class TestQRunBlockPulse(unittest.TestCase):
    """Tests for the q_run_block_pulse adapter."""

    def test_block_pulse_with_dict(self):
        """Block pulse with dict returns a QFrame."""
        frame = q_run_block_pulse({
            "task_id": "BLOCK-QUINT-001",
            "ctx": {"source_id": "PX-KERNEL-01", "fingerprint": "PX-ROOT-Sovereign"},
            "p_vec": TEST_P_VECTOR,
            "csa_s": TEST_CSA_SCORES,
        })
        self.assertIsInstance(frame, QFrame)
        self.assertEqual(frame.qtype, QType.QRESULT)
        self.assertEqual(frame.payload.get("engine_id"), "BLOCK_ORCHESTRATOR_V1")

    def test_block_pulse_with_qframe(self):
        """Block pulse accepts QFrame input."""
        input_frame = ingest({
            "task_id": "BLOCK-QF-001",
            "ctx": {"source_id": "PX-KERNEL-01"},
            "p_vec": TEST_P_VECTOR,
            "csa_s": TEST_CSA_SCORES,
        })
        frame = q_run_block_pulse(input_frame)
        self.assertIsInstance(frame, QFrame)


# ---------------------------------------------------------------------------
# Tests: Cross-cutting concerns
# ---------------------------------------------------------------------------

class TestBackwardCompatibility(unittest.TestCase):
    """All adapters accept plain dicts without requiring QFrame wrapping."""

    def test_prv_adapters_backward_compatible(self):
        """All 12 PRV adapters work with plain dict input."""
        # OPE
        f = q_run_ope({"smiles": TEST_SMILES})
        self.assertIsInstance(f, QFrame)

        # OBE
        f = q_run_obe({"smiles": TEST_SMILES})
        self.assertIsInstance(f, QFrame)

        # OCE
        f = q_run_oce({"p_vector": TEST_P_VECTOR, "csa_scores": TEST_CSA_SCORES})
        self.assertIsInstance(f, QFrame)

        # OLE
        f = q_run_ole({"compound_id": "T-1", "indication": "nipah"})
        self.assertIsInstance(f, QFrame)

        # OME
        f = q_run_ome({"smiles": TEST_SMILES})
        self.assertIsInstance(f, QFrame)

        # OSE
        f = q_run_ose({"smiles": TEST_SMILES})
        self.assertIsInstance(f, QFrame)

        # ADMET (needs ope_result)
        ope = {"molecular_weight": 46.07, "logp": -0.31, "hbd": 1, "hba": 1, "tpsa": 20.23}
        f = q_run_admet({"smiles": TEST_SMILES}, ope_result=ope)
        self.assertIsInstance(f, QFrame)

        # Grading
        f = q_run_grading(_minimal_dossier())
        self.assertIsInstance(f, QFrame)

        # Zeus
        f = q_run_zeus(_minimal_dossier())
        self.assertIsInstance(f, QFrame)

    def test_auxiliary_adapters_backward_compatible(self):
        """All 8 auxiliary adapters work with plain dict input."""
        # VectorCore
        f = q_run_vector_core({"p_vector": TEST_P_VECTOR})
        self.assertIsInstance(f, QFrame)

        # SimulationEngine (low coherence to avoid needing worldline files)
        f = q_run_simulation({"task_id": "WL-BC-001", "coherence": 0.50})
        self.assertIsInstance(f, QFrame)

        # TrialEngine
        f = q_run_trial(TEST_PROTOCOL, admet=TEST_ADMET_FULL)
        self.assertIsInstance(f, QFrame)

        # CheckConstitutional
        f = q_run_check_constitutional("backcompat_test", {"toxicity_index": 0.005, "harm_energy": 0.005})
        self.assertIsInstance(f, QFrame)

        # EvidencePackage
        f = q_run_evidence_dossier(_minimal_candidate(), engine_results=_minimal_engine_results())
        self.assertIsInstance(f, QFrame)

        # Metabolism
        f = q_run_metabolism_pulse({"task_id": "BC-METAB-001", "p_vec": TEST_P_VECTOR, "csa_s": TEST_CSA_SCORES})
        self.assertIsInstance(f, QFrame)

        # BlockOrchestrator
        f = q_run_block_pulse({"task_id": "BC-BLOCK-001", "ctx": {}, "p_vec": TEST_P_VECTOR, "csa_s": TEST_CSA_SCORES})
        self.assertIsInstance(f, QFrame)


class TestLineageTracking(unittest.TestCase):
    """Lineage in returned QFrames should contain the engine name."""

    def test_prv_lineage_tracking(self):
        """Each PRV adapter appends its engine name to frame lineage."""
        engines_and_inputs = [
            ("OPE", lambda: q_run_ope({"smiles": TEST_SMILES})),
            ("OBE", lambda: q_run_obe({"smiles": TEST_SMILES})),
            ("OCE", lambda: q_run_oce({"p_vector": TEST_P_VECTOR, "csa_scores": TEST_CSA_SCORES})),
            ("OLE", lambda: q_run_ole({"compound_id": "T-1", "indication": "nipah"})),
            ("OME", lambda: q_run_ome({"smiles": TEST_SMILES})),
            ("OSE", lambda: q_run_ose({"smiles": TEST_SMILES})),
            ("GradingEngine", lambda: q_run_grading(_minimal_dossier())),
            ("ZeusLaws", lambda: q_run_zeus(_minimal_dossier())),
        ]

        for engine_name, run_fn in engines_and_inputs:
            with self.subTest(engine=engine_name):
                frame = run_fn()
                lineage_str = " ".join(frame.lineage)
                self.assertIn(
                    f"engine:{engine_name}",
                    lineage_str,
                    f"Lineage should contain 'engine:{engine_name}', got: {frame.lineage}",
                )

    def test_auxiliary_lineage_tracking(self):
        """Each auxiliary adapter appends its engine name to frame lineage."""
        engines_and_inputs = [
            ("VectorCore", lambda: q_run_vector_core({"p_vector": TEST_P_VECTOR})),
            ("SimulationEngine", lambda: q_run_simulation({"task_id": "WL-LIN-001", "coherence": 0.50})),
            ("TrialEngine", lambda: q_run_trial(TEST_PROTOCOL, admet=TEST_ADMET_FULL)),
            ("CheckConstitutional", lambda: q_run_check_constitutional("lineage_test", {"toxicity_index": 0.005, "harm_energy": 0.005})),
            ("EvidencePackage", lambda: q_run_evidence_dossier(_minimal_candidate(), engine_results=_minimal_engine_results())),
            ("Metabolism", lambda: q_run_metabolism_pulse({"task_id": "LIN-METAB-001", "p_vec": TEST_P_VECTOR, "csa_s": TEST_CSA_SCORES})),
            ("BlockOrchestrator", lambda: q_run_block_pulse({"task_id": "LIN-BLOCK-001", "ctx": {}, "p_vec": TEST_P_VECTOR, "csa_s": TEST_CSA_SCORES})),
        ]

        for engine_name, run_fn in engines_and_inputs:
            with self.subTest(engine=engine_name):
                frame = run_fn()
                lineage_str = " ".join(frame.lineage)
                self.assertIn(
                    f"engine:{engine_name}",
                    lineage_str,
                    f"Lineage should contain 'engine:{engine_name}', got: {frame.lineage}",
                )


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    unittest.main()
