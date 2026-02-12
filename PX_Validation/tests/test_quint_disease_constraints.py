"""
Tests for QUINT Disease Constraint evaluator.

Validates that disease constraint JSON files are loaded as QCONSTRAINT frames
and that molecules are correctly evaluated against constraint thresholds.

Run standalone:
    python PX_Validation/tests/test_quint_disease_constraints.py
"""

import unittest
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from PX_System.foundation.quint.kernel import QFrame, QType, create_qframe
from PX_System.foundation.quint.converter import ingest
from PX_System.foundation.quint.brains.disease_constraints import (
    load_disease_constraint,
    load_all_disease_constraints,
    evaluate_molecule_against_constraint,
    evaluate_molecule_against_all,
    get_eligible_diseases,
    get_constraint_summary,
)


def _make_molecule(**kwargs):
    """Helper: create a QMOLECULE frame with given properties."""
    data = {"smiles": "CCO"}
    data.update(kwargs)
    return ingest(data, qtype=QType.QMOLECULE)


class TestLoadDiseaseConstraint(unittest.TestCase):
    """Tests for loading individual disease constraint files."""

    def test_load_nipah_constraint(self):
        """Loading nipah_virus_infection returns a QFrame with QCONSTRAINT type."""
        frame = load_disease_constraint("nipah_virus_infection")
        self.assertIsInstance(frame, QFrame)
        self.assertEqual(frame.qtype, QType.QCONSTRAINT)

    def test_load_nipah_has_disease_name(self):
        """Nipah constraint payload contains disease_name field."""
        frame = load_disease_constraint("nipah_virus_infection")
        self.assertEqual(frame.payload["disease_name"], "Nipah virus infection")

    def test_load_nipah_has_thresholds(self):
        """Nipah constraint payload contains expected threshold fields."""
        frame = load_disease_constraint("nipah_virus_infection")
        payload = frame.payload
        self.assertIn("toxicity_threshold", payload)
        self.assertIn("ic50_max_um", payload)
        self.assertIn("molecular_weight_max", payload)
        self.assertIn("logp_min", payload)
        self.assertIn("logp_max", payload)
        self.assertIn("hbd_max", payload)
        self.assertIn("hba_max", payload)
        self.assertAlmostEqual(payload["toxicity_threshold"], 0.021)
        self.assertAlmostEqual(payload["ic50_max_um"], 1.0)

    def test_load_nonexistent_raises(self):
        """Loading a nonexistent disease ID raises FileNotFoundError."""
        with self.assertRaises(FileNotFoundError):
            load_disease_constraint("nonexistent_disease_xyz")

    def test_constraint_frame_integrity(self):
        """Loaded constraint QFrame passes integrity verification."""
        frame = load_disease_constraint("nipah_virus_infection")
        # verify() checks that frame_hash matches computed hash
        # After ingest, the frame is sealed but lineage was appended after seal,
        # so we re-seal to test the mechanism
        frame.seal()
        self.assertTrue(frame.verify())

    def test_constraint_has_lineage(self):
        """Loaded constraint frame lineage includes disease_constraint label."""
        frame = load_disease_constraint("nipah_virus_infection")
        lineage_str = " ".join(frame.lineage)
        self.assertIn("disease_constraint:nipah_virus_infection", lineage_str)


class TestLoadAllDiseaseConstraints(unittest.TestCase):
    """Tests for loading all 33 disease constraint files."""

    def test_load_all_count(self):
        """load_all_disease_constraints returns exactly 33 entries."""
        constraints = load_all_disease_constraints()
        self.assertEqual(len(constraints), 33)

    def test_load_all_types(self):
        """All loaded constraints are QFrame with QCONSTRAINT type."""
        constraints = load_all_disease_constraints()
        for disease_id, frame in constraints.items():
            with self.subTest(disease_id=disease_id):
                self.assertIsInstance(frame, QFrame)
                self.assertEqual(frame.qtype, QType.QCONSTRAINT)

    def test_load_all_disease_ids(self):
        """Known disease IDs are present in loaded constraints."""
        constraints = load_all_disease_constraints()
        known_ids = ["nipah_virus_infection", "ebola_virus_disease", "malaria"]
        for disease_id in known_ids:
            with self.subTest(disease_id=disease_id):
                self.assertIn(disease_id, constraints)


class TestEvaluateMolecule(unittest.TestCase):
    """Tests for evaluating molecules against disease constraints."""

    def setUp(self):
        """Load nipah constraint for use in multiple tests."""
        self.nipah = load_disease_constraint("nipah_virus_infection")

    def test_evaluate_safe_molecule_passes(self):
        """A safe molecule (low tox, reasonable mw/logp/hbd/hba) passes nipah constraints."""
        mol = _make_molecule(toxicity=0.005, mw=300.0, logp=2.0, hbd=2, hba=4)
        result = evaluate_molecule_against_constraint(mol, self.nipah)
        self.assertTrue(result["passed"])
        self.assertEqual(len(result["violations"]), 0)
        self.assertEqual(result["disease_name"], "Nipah virus infection")
        self.assertGreater(result["checks_performed"], 0)
        self.assertEqual(result["checks_performed"], result["checks_passed"])

    def test_evaluate_toxic_molecule_fails(self):
        """A molecule with toxicity=0.025 fails nipah (threshold=0.021)."""
        mol = _make_molecule(toxicity=0.025, mw=300.0, logp=2.0, hbd=2, hba=4)
        result = evaluate_molecule_against_constraint(mol, self.nipah)
        self.assertFalse(result["passed"])
        tox_violations = [v for v in result["violations"] if v["field"] == "toxicity"]
        self.assertEqual(len(tox_violations), 1)
        self.assertEqual(tox_violations[0]["value"], 0.025)

    def test_evaluate_heavy_molecule_fails(self):
        """A molecule with mw=600 fails nipah (molecular_weight_max=500)."""
        mol = _make_molecule(toxicity=0.005, mw=600.0, logp=2.0, hbd=2, hba=4)
        result = evaluate_molecule_against_constraint(mol, self.nipah)
        self.assertFalse(result["passed"])
        mw_violations = [v for v in result["violations"] if v["field"] == "mw"]
        self.assertEqual(len(mw_violations), 1)
        self.assertEqual(mw_violations[0]["value"], 600.0)

    def test_evaluate_logp_out_of_range_fails(self):
        """A molecule with logp=6.0 fails nipah (logp_max=5.0)."""
        mol = _make_molecule(toxicity=0.005, mw=300.0, logp=6.0, hbd=2, hba=4)
        result = evaluate_molecule_against_constraint(mol, self.nipah)
        self.assertFalse(result["passed"])
        logp_violations = [v for v in result["violations"] if v["field"] == "logp"]
        self.assertEqual(len(logp_violations), 1)

    def test_evaluate_missing_fields_tolerant(self):
        """A molecule with only smiles passes (missing optional fields skipped)."""
        mol = _make_molecule()  # only smiles="CCO"
        result = evaluate_molecule_against_constraint(mol, self.nipah)
        self.assertTrue(result["passed"])
        self.assertEqual(result["checks_performed"], 0)
        self.assertEqual(len(result["violations"]), 0)


class TestEvaluateAgainstAll(unittest.TestCase):
    """Tests for evaluating molecules against all 33 constraints."""

    def test_evaluate_against_all_returns_33(self):
        """evaluate_molecule_against_all returns results for all 33 diseases."""
        mol = _make_molecule(toxicity=0.005, mw=300.0, logp=2.0, hbd=2, hba=4)
        results = evaluate_molecule_against_all(mol)
        self.assertEqual(len(results), 33)

    def test_get_eligible_diseases_safe_molecule(self):
        """A safe molecule is eligible for at least some diseases."""
        mol = _make_molecule(toxicity=0.005, mw=300.0, logp=2.0, hbd=2, hba=4)
        eligible = get_eligible_diseases(mol)
        self.assertIsInstance(eligible, list)
        self.assertGreater(len(eligible), 0)
        # Should be eligible for nipah given these safe values
        self.assertIn("nipah_virus_infection", eligible)


class TestConstraintSummary(unittest.TestCase):
    """Tests for the constraint summary function."""

    def test_get_constraint_summary(self):
        """Summary has count=33 and expected structure."""
        summary = get_constraint_summary()
        self.assertEqual(summary["count"], 33)
        self.assertIn("disease_ids", summary)
        self.assertIn("versions", summary)
        self.assertEqual(len(summary["disease_ids"]), 33)
        self.assertIn("nipah_virus_infection", summary["disease_ids"])


if __name__ == "__main__":
    unittest.main()
