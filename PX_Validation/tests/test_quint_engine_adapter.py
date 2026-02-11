"""
test_quint_engine_adapter.py
Tests for the QUINT Engine Adapter layer.

Validates that all 12 PRV pipeline engine wrappers:
  - Accept both QFrame and plain dict inputs (backward compatibility)
  - Return QFrame with QType.QRESULT
  - Preserve engine_id in payload
  - Track lineage with engine name
"""

import sys
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from PX_System.foundation.quint.kernel import QFrame, QType, create_qframe
from PX_System.foundation.quint.converter import ingest
from PX_System.foundation.quint.engine_adapter import (
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


# ---------------------------------------------------------------------------
# Tests: Cross-cutting concerns
# ---------------------------------------------------------------------------

class TestBackwardCompatibility(unittest.TestCase):
    """All adapters accept plain dicts without requiring QFrame wrapping."""

    def test_all_adapters_backward_compatible(self):
        """Every adapter works with plain dict input."""
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


class TestLineageTracking(unittest.TestCase):
    """Lineage in returned QFrames should contain the engine name."""

    def test_lineage_tracking(self):
        """Each adapter appends its engine name to frame lineage."""
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


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    unittest.main()
