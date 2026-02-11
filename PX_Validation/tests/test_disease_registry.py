"""
Tests for PX_Domain.PRV_Diseases.disease_registry — thin adapter over QUINT disease_constraints.

Validates:
  - Manifest loading (31 diseases from manifest.json)
  - Constraint loading (33 JSON files via QUINT)
  - Candidate evaluation (pass/fail/unknown disease)
  - Enhanced schema (molecular_constraints, target_profile, clinical, regulatory)
  - OLE integration (dynamic PRV disease set)
  - All manifest diseases have constraint files (coverage)
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
        self.assertEqual(len(constraints), 33)

    def test_get_constraint_summary(self):
        from PX_Domain.PRV_Diseases.disease_registry import get_constraint_summary
        summary = get_constraint_summary()
        self.assertEqual(summary["count"], 33)
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
        self.assertEqual(len(results), 33)
        # A well-behaved molecule should pass most constraints
        passed = sum(1 for r in results.values() if r["passed"])
        self.assertGreater(passed, 0)


class TestEnhancedSchema(unittest.TestCase):
    """Tests for the enriched constraint schema (v2.0.0)."""

    def test_all_constraints_have_category(self):
        from PX_Domain.PRV_Diseases.disease_registry import get_all_constraints
        for disease_id, frame in get_all_constraints().items():
            self.assertIn("category", frame.payload, f"{disease_id} missing category")

    def test_all_constraints_have_molecular_constraints(self):
        from PX_Domain.PRV_Diseases.disease_registry import get_all_constraints
        for disease_id, frame in get_all_constraints().items():
            mc = frame.payload.get("molecular_constraints")
            self.assertIsNotNone(mc, f"{disease_id} missing molecular_constraints")
            self.assertIn("mw_max", mc)
            self.assertIn("toxicity_max", mc)

    def test_all_constraints_have_target_profile(self):
        from PX_Domain.PRV_Diseases.disease_registry import get_all_constraints
        for disease_id, frame in get_all_constraints().items():
            tp = frame.payload.get("target_profile")
            self.assertIsNotNone(tp, f"{disease_id} missing target_profile")
            self.assertIn("pathogen_type", tp)
            self.assertIn("mechanism_classes", tp)
            self.assertIn("tissue_targets", tp)

    def test_all_constraints_have_clinical_requirements(self):
        from PX_Domain.PRV_Diseases.disease_registry import get_all_constraints
        for disease_id, frame in get_all_constraints().items():
            cr = frame.payload.get("clinical_requirements")
            self.assertIsNotNone(cr, f"{disease_id} missing clinical_requirements")
            self.assertIn("route", cr)
            self.assertIn("treatment_setting", cr)

    def test_all_constraints_have_regulatory(self):
        from PX_Domain.PRV_Diseases.disease_registry import get_all_constraints
        for disease_id, frame in get_all_constraints().items():
            reg = frame.payload.get("regulatory")
            self.assertIsNotNone(reg, f"{disease_id} missing regulatory")
            self.assertEqual(reg["fda_pathway"], "PRV")

    def test_pediatric_have_lower_toxicity(self):
        from PX_Domain.PRV_Diseases.disease_registry import get_all_constraints
        for disease_id, frame in get_all_constraints().items():
            if frame.payload.get("category") == "rare_pediatric":
                tox = frame.payload["molecular_constraints"]["toxicity_max"]
                self.assertLessEqual(tox, 0.018,
                    f"{disease_id}: pediatric toxicity_max {tox} > 0.018")

    def test_all_toxicity_below_constitutional_limit(self):
        from PX_Domain.PRV_Diseases.disease_registry import get_all_constraints
        for disease_id, frame in get_all_constraints().items():
            tox = frame.payload["molecular_constraints"]["toxicity_max"]
            self.assertLessEqual(tox, 0.0210,
                f"{disease_id}: toxicity_max {tox} exceeds constitutional limit 0.0210")

    def test_constraint_version_is_2(self):
        from PX_Domain.PRV_Diseases.disease_registry import get_all_constraints
        for disease_id, frame in get_all_constraints().items():
            self.assertEqual(frame.payload.get("constraint_version"), "2.0.0",
                f"{disease_id} still on old constraint_version")


class TestManifestCoverage(unittest.TestCase):
    """Tests that all manifest diseases have constraint files."""

    def test_all_manifest_diseases_have_constraints(self):
        from PX_Domain.PRV_Diseases.disease_registry import get_prv_diseases, get_all_constraints
        diseases = get_prv_diseases()
        constraints = get_all_constraints()
        for d in diseases:
            disease_id = d["name"].lower().replace(" ", "_")
            self.assertIn(disease_id, constraints,
                f"Missing constraint for {d['name']}")


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
