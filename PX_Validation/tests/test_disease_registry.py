"""
Tests for PX_Domain.PRV_Diseases.disease_registry — thin adapter over QUINT disease_constraints.

Validates:
  - Manifest loading (31 diseases from manifest.json)
  - Constraint loading (22 JSON files via QUINT)
  - Candidate evaluation (pass/fail/unknown disease)
  - OLE integration (dynamic PRV disease set)
  - GradingEngine constraint cap
"""

import os
import sys
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)


class TestManifestLoading(unittest.TestCase):
    """Tests for loading the PRV disease manifest."""

    def test_get_prv_diseases_returns_list(self):
        from PX_Domain.PRV_Diseases.disease_registry import get_prv_diseases
        diseases = get_prv_diseases()
        self.assertIsInstance(diseases, list)
        self.assertGreater(len(diseases), 0)

    def test_manifest_has_31_diseases(self):
        from PX_Domain.PRV_Diseases.disease_registry import get_prv_diseases
        diseases = get_prv_diseases()
        self.assertEqual(len(diseases), 31)

    def test_each_disease_has_name(self):
        from PX_Domain.PRV_Diseases.disease_registry import get_prv_diseases
        for d in get_prv_diseases():
            self.assertIn("name", d)
            self.assertTrue(len(d["name"]) > 0)

    def test_get_prv_disease_names(self):
        from PX_Domain.PRV_Diseases.disease_registry import get_prv_disease_names
        names = get_prv_disease_names()
        self.assertEqual(len(names), 31)
        self.assertIn("Nipah virus infection", names)


class TestConstraintLoading(unittest.TestCase):
    """Tests for loading disease constraint files."""

    def test_get_all_constraints(self):
        from PX_Domain.PRV_Diseases.disease_registry import get_all_constraints
        constraints = get_all_constraints()
        self.assertEqual(len(constraints), 22)

    def test_get_constraint_summary(self):
        from PX_Domain.PRV_Diseases.disease_registry import get_constraint_summary
        summary = get_constraint_summary()
        self.assertEqual(summary["count"], 22)
        self.assertIn("nipah_virus_infection", summary["disease_ids"])

    def test_get_single_constraint(self):
        from PX_Domain.PRV_Diseases.disease_registry import get_disease_constraint
        frame = get_disease_constraint("nipah_virus_infection")
        self.assertIsNotNone(frame)
        self.assertEqual(frame.payload["disease_name"], "Nipah virus infection")

    def test_unknown_disease_returns_none(self):
        from PX_Domain.PRV_Diseases.disease_registry import get_disease_constraint
        result = get_disease_constraint("nonexistent_disease_xyz")
        self.assertIsNone(result)


class TestCandidateEvaluation(unittest.TestCase):
    """Tests for evaluating molecules against disease constraints."""

    def test_passing_molecule(self):
        from PX_Domain.PRV_Diseases.disease_registry import evaluate_candidate
        mol = {"mw": 350, "logp": 2.5, "hbd": 2, "hba": 5, "toxicity": 0.01}
        result = evaluate_candidate(mol, "nipah_virus_infection")
        self.assertIsNotNone(result)
        self.assertTrue(result["passed"])
        self.assertEqual(len(result["violations"]), 0)
        self.assertEqual(result["checks_performed"], 5)

    def test_failing_molecule_mw(self):
        from PX_Domain.PRV_Diseases.disease_registry import evaluate_candidate
        mol = {"mw": 800, "logp": 2.5, "hbd": 2, "hba": 5, "toxicity": 0.01}
        result = evaluate_candidate(mol, "nipah_virus_infection")
        self.assertFalse(result["passed"])
        fields = [v["field"] for v in result["violations"]]
        self.assertIn("mw", fields)

    def test_failing_molecule_toxicity(self):
        from PX_Domain.PRV_Diseases.disease_registry import evaluate_candidate
        mol = {"mw": 350, "logp": 2.5, "hbd": 2, "hba": 5, "toxicity": 0.03}
        result = evaluate_candidate(mol, "nipah_virus_infection")
        self.assertFalse(result["passed"])
        fields = [v["field"] for v in result["violations"]]
        self.assertIn("toxicity", fields)

    def test_unknown_disease_returns_none(self):
        from PX_Domain.PRV_Diseases.disease_registry import evaluate_candidate
        mol = {"mw": 350, "logp": 2.5}
        result = evaluate_candidate(mol, "nonexistent_disease_xyz")
        self.assertIsNone(result)

    def test_evaluate_candidate_all(self):
        from PX_Domain.PRV_Diseases.disease_registry import evaluate_candidate_all
        mol = {"mw": 350, "logp": 2.5, "hbd": 2, "hba": 5, "toxicity": 0.01}
        results = evaluate_candidate_all(mol)
        self.assertEqual(len(results), 22)
        # A well-behaved molecule should pass most constraints
        passed = sum(1 for r in results.values() if r["passed"])
        self.assertGreater(passed, 0)


class TestOLEIntegration(unittest.TestCase):
    """Tests that OLE uses the disease registry."""

    def test_ole_loads_dynamic_diseases(self):
        from PX_Engine.operations.OLE import PRV_ELIGIBLE_DISEASES
        # Should include manifest diseases + fallback aliases
        self.assertGreater(len(PRV_ELIGIBLE_DISEASES), 28)
        self.assertIn("nipah virus infection", PRV_ELIGIBLE_DISEASES)
        self.assertIn("nipah", PRV_ELIGIBLE_DISEASES)
        # New diseases from manifest that weren't in old hardcoded set
        self.assertIn("spinal muscular atrophy", PRV_ELIGIBLE_DISEASES)
        self.assertIn("anthrax", PRV_ELIGIBLE_DISEASES)

    def test_ole_execute_prv_eligible(self):
        from PX_Engine.operations.OLE import execute
        result = execute({"compound_id": "TEST", "indication": "Nipah virus infection"})
        self.assertTrue(result["prv_eligible"])
        self.assertEqual(result["regulatory_pathway"], "PRV_ACCELERATED")


if __name__ == "__main__":
    unittest.main()
