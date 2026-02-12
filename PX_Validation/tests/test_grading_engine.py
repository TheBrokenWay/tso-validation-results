"""
test_grading_engine.py
Tests for WeightedGradingEngine v2.0 — 100-point MPO grading system

Test categories:
  - Hard limit tests (toxicity, hERG, zeus, constitutional)
  - Weighted scoring tests (DIAMOND, GOLD, SILVER, BRONZE)
  - MPO recovery test (IC50=5nM + tox=0.008 + MW=520 → DIAMOND)
  - Anti-inflation test
  - Score breakdown / reasoning tests
  - Missing data tests
  - Constraint violation test
  - Backward compatibility tests
  - Governance logging tests
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

import unittest
from PX_Engine.operations.GradingEngine import GradingEngine, WeightedGradingEngine


def _make_dossier(*, ec50=50, emax=85, tox=0.008, safety_margin=60,
                  herg_risk=0.1, bioavailability_pct=75, half_life_h=8,
                  clearance_L_h=0.5, vd_L_kg=1.5, mw=400, logp=3.0,
                  hbd=2, hba=5, fto=True, qed=None, synthesis_steps=None,
                  constraint_violations=None, zeus_verdict=None,
                  constitutional_violations=None):
    """Build a dossier dict with controllable fields for testing."""
    d = {
        "engines": {
            "ope": {
                "ec50": ec50,
                "emax": emax,
                "molecular_weight": mw,
                "logp": logp,
                "hbd": hbd,
                "hba": hba,
            },
            "admet": {
                "toxicity": {"toxicity_index": tox, "herg_risk": herg_risk,
                             "safety_margin": safety_margin},
                "absorption": {"oral_bioavailability_percent": bioavailability_pct},
                "metabolism": {"half_life_h": half_life_h,
                               "clearance_L_per_h": clearance_L_h},
                "distribution": {"predicted_vd_L_per_kg": vd_L_kg},
            },
            "ole": {"freedom_to_operate": fto},
        },
    }
    if qed is not None:
        d["qed"] = qed
    if synthesis_steps is not None:
        d["synthesis_steps"] = synthesis_steps
    if constraint_violations is not None:
        d["constraint_violations"] = constraint_violations
    if zeus_verdict is not None:
        d["zeus_verdict"] = zeus_verdict
    if constitutional_violations is not None:
        d["constitutional_violations"] = constitutional_violations
    return d


class TestHardLimits(unittest.TestCase):
    """Hard limit tests — immediate REJECTED, cannot be overridden by scoring."""

    def setUp(self):
        self.engine = GradingEngine(verbose=False)

    def test_toxicity_hard_limit(self):
        """Toxicity >= 0.0210 → REJECTED (Law L11)"""
        d = _make_dossier(tox=0.0210)
        result = self.engine.grade_dossier(d)
        self.assertEqual(result["grade"], "REJECTED")
        self.assertTrue(result["hard_limit_checks"]["toxicity_failed"])
        self.assertIn("L11", result["hard_limit_checks"]["rejection_reasons"][0])

    def test_herg_hard_limit(self):
        """hERG risk > 0.5 → REJECTED"""
        d = _make_dossier(herg_risk=0.6)
        result = self.engine.grade_dossier(d)
        self.assertEqual(result["grade"], "REJECTED")
        self.assertTrue(result["hard_limit_checks"]["herg_failed"])

    def test_zeus_gate_failed(self):
        """Zeus gate not authorized → REJECTED"""
        d = _make_dossier(zeus_verdict={"authorized": False, "rationale": "test fail"})
        result = self.engine.grade_dossier(d)
        self.assertEqual(result["grade"], "REJECTED")
        self.assertTrue(result["hard_limit_checks"]["zeus_failed"])

    def test_constitutional_violations(self):
        """Constitutional violations > 0 → REJECTED"""
        d = _make_dossier(constitutional_violations=["test_violation"])
        result = self.engine.grade_dossier(d)
        self.assertEqual(result["grade"], "REJECTED")
        self.assertTrue(result["hard_limit_checks"]["constitutional_failed"])


class TestWeightedScoring(unittest.TestCase):
    """Weighted scoring tiers — DIAMOND >= 85%, GOLD 70-84%, SILVER 50-69%, BRONZE < 50%."""

    def setUp(self):
        self.engine = GradingEngine(verbose=False)

    def test_diamond_tier(self):
        """Perfect candidate → DIAMOND_TIER (>= 85%)"""
        d = _make_dossier(ec50=5, emax=95, tox=0.003, safety_margin=120,
                          herg_risk=0.05, bioavailability_pct=90, half_life_h=8,
                          clearance_L_h=0.2, vd_L_kg=1.0, mw=350, logp=2.5,
                          hbd=1, hba=4, fto=True, qed=0.8, synthesis_steps=3)
        result = self.engine.grade_dossier(d)
        self.assertEqual(result["grade"], "DIAMOND_TIER")
        self.assertGreaterEqual(result["percentage"], 85.0)

    def test_gold_tier(self):
        """Good candidate → GOLD_TIER (70-84%)"""
        d = _make_dossier(ec50=80, emax=75, tox=0.012, safety_margin=55,
                          herg_risk=0.2, bioavailability_pct=65, half_life_h=6,
                          clearance_L_h=0.8, vd_L_kg=2.0, mw=420, logp=3.5,
                          hbd=3, hba=6, fto=True)
        result = self.engine.grade_dossier(d)
        self.assertEqual(result["grade"], "GOLD_TIER")
        self.assertGreaterEqual(result["percentage"], 70.0)
        self.assertLess(result["percentage"], 85.0)

    def test_silver_tier(self):
        """Moderate candidate → SILVER_TIER (50-69%)"""
        d = _make_dossier(ec50=400, emax=45, tox=0.018, safety_margin=15,
                          herg_risk=0.4, bioavailability_pct=35, half_life_h=1.5,
                          clearance_L_h=2.0, vd_L_kg=8.0, mw=480, logp=4.5,
                          hbd=4, hba=9, fto=True)
        result = self.engine.grade_dossier(d)
        self.assertEqual(result["grade"], "SILVER_TIER")
        self.assertGreaterEqual(result["percentage"], 50.0)
        self.assertLess(result["percentage"], 70.0)

    def test_bronze_tier(self):
        """Poor candidate → BRONZE_TIER (< 50%)"""
        d = _make_dossier(ec50=5000, emax=20, tox=0.0205, safety_margin=5,
                          herg_risk=0.45, bioavailability_pct=10, half_life_h=0.5,
                          clearance_L_h=10.0, vd_L_kg=15.0, mw=600, logp=6.5,
                          hbd=7, hba=12, fto=False)
        result = self.engine.grade_dossier(d)
        self.assertEqual(result["grade"], "BRONZE_TIER")
        self.assertLess(result["percentage"], 50.0)


class TestMPORecovery(unittest.TestCase):
    """MPO recovery guarantee: exceptional efficacy + safety compensates minor property violations."""

    def setUp(self):
        self.engine = GradingEngine(verbose=False)

    def test_ic50_5nM_tox_low_mw_520_is_diamond(self):
        """IC50=5nM + tox=0.008 + MW=520 → DIAMOND (was SILVER/NEEDS_REVIEW in old system)"""
        d = _make_dossier(ec50=5, emax=92, tox=0.008, safety_margin=80,
                          herg_risk=0.08, bioavailability_pct=80, half_life_h=8,
                          clearance_L_h=0.3, vd_L_kg=1.2, mw=520, logp=3.0,
                          hbd=2, hba=6, fto=True)
        result = self.engine.grade_dossier(d)
        self.assertEqual(result["grade"], "DIAMOND_TIER",
                         f"MPO recovery failed: got {result['grade']} ({result['percentage']}%)")
        self.assertGreaterEqual(result["percentage"], 85.0)


class TestAntiInflation(unittest.TestCase):
    """Anti-inflation: one strong metric cannot save overall poor candidate."""

    def setUp(self):
        self.engine = GradingEngine(verbose=False)

    def test_strong_efficacy_poor_everything_else(self):
        """IC50=1nM but terrible safety/PK/druglikeness → not DIAMOND or GOLD"""
        d = _make_dossier(ec50=1, emax=99, tox=0.0205, safety_margin=3,
                          herg_risk=0.45, bioavailability_pct=8, half_life_h=0.3,
                          clearance_L_h=15.0, vd_L_kg=20.0, mw=650, logp=7.0,
                          hbd=8, hba=14, fto=False)
        result = self.engine.grade_dossier(d)
        self.assertNotIn(result["grade"], ("DIAMOND_TIER", "GOLD_TIER"),
                         "One exceptional category should not inflate overall grade")


class TestScoreBreakdown(unittest.TestCase):
    """Score breakdown: sums, percentages, reasoning audit trail."""

    def setUp(self):
        self.engine = GradingEngine(verbose=False)

    def test_category_scores_sum_to_total(self):
        """Sum of category scores must equal total_score"""
        d = _make_dossier()
        result = self.engine.grade_dossier(d)
        cat_sum = sum(cs["score"] for cs in result["category_scores"].values())
        self.assertAlmostEqual(cat_sum, result["total_score"], places=1)

    def test_percentage_matches_total(self):
        """Percentage must equal (total_score / max_score) * 100"""
        d = _make_dossier()
        result = self.engine.grade_dossier(d)
        expected_pct = round((result["total_score"] / result["max_score"]) * 100, 1)
        self.assertAlmostEqual(result["percentage"], expected_pct, places=1)

    def test_reasoning_has_strengths_weaknesses(self):
        """Reasoning must contain strengths and weaknesses lists"""
        d = _make_dossier()
        result = self.engine.grade_dossier(d)
        self.assertIn("strengths", result["reasoning"])
        self.assertIn("weaknesses", result["reasoning"])
        self.assertIn("summary", result["reasoning"])
        self.assertIsInstance(result["reasoning"]["strengths"], list)
        self.assertIsInstance(result["reasoning"]["weaknesses"], list)


class TestMissingData(unittest.TestCase):
    """Missing data handling — partial credit, never zero."""

    def setUp(self):
        self.engine = GradingEngine(verbose=False)

    def test_missing_efficacy_gets_partial_credit(self):
        """No EC50/Emax data → efficacy should score 50% (partial credit)"""
        d = {"engines": {"ope": {}, "admet": {"toxicity": 0.01}}}
        result = self.engine.grade_dossier(d)
        eff = result["category_scores"]["efficacy"]
        self.assertEqual(eff["score"], 15.0)  # 10 (ic50 missing) + 5 (emax missing)
        self.assertEqual(eff["percentage"], 50.0)

    def test_missing_pk_gets_partial_credit(self):
        """No PK data → pk should score 50% (partial credit)"""
        d = {"engines": {"ope": {}, "admet": {"toxicity": 0.01}}}
        result = self.engine.grade_dossier(d)
        pk = result["category_scores"]["pk"]
        self.assertEqual(pk["score"], 10.0)  # 4 + 3 + 1.5 + 1.5
        self.assertEqual(pk["percentage"], 50.0)

    def test_missing_fto_gets_partial_credit(self):
        """No FTO status → fto should score 50% (partial credit)"""
        d = {"engines": {"ope": {}, "admet": {"toxicity": 0.01}}}
        result = self.engine.grade_dossier(d)
        fto = result["category_scores"]["fto"]
        self.assertEqual(fto["score"], 2.5)
        self.assertEqual(fto["percentage"], 50.0)


class TestConstraintViolations(unittest.TestCase):
    """Constraint violations reduce druglikeness score (no tier cap)."""

    def setUp(self):
        self.engine = GradingEngine(verbose=False)

    def test_constraint_violations_reduce_score(self):
        """Each constraint violation → -1 point in druglikeness (max -3)"""
        d_clean = _make_dossier()
        d_dirty = _make_dossier(constraint_violations=["v1", "v2"])

        r_clean = self.engine.grade_dossier(d_clean)
        r_dirty = self.engine.grade_dossier(d_dirty)

        dl_clean = r_clean["category_scores"]["druglikeness"]["score"]
        dl_dirty = r_dirty["category_scores"]["druglikeness"]["score"]
        self.assertEqual(dl_clean - dl_dirty, 2.0)


class TestBackwardCompatibility(unittest.TestCase):
    """Backward compatibility: legacy dossier format and alias."""

    def setUp(self):
        self.engine = GradingEngine(verbose=False)

    def test_legacy_dossier_flat_toxicity(self):
        """Legacy dossier with top-level 'admet': {'toxicity': float} should work"""
        d = {
            "admet": {"toxicity": 0.015},
            "ope": {"ec50": 100, "emax": 70},
        }
        result = self.engine.grade_dossier(d)
        self.assertIn(result["grade"], ("DIAMOND_TIER", "GOLD_TIER", "SILVER_TIER", "BRONZE_TIER"))
        self.assertIn("metrics", result)
        self.assertIn("grading_engine_version", result)

    def test_alias_is_weighted(self):
        """GradingEngine must be an alias for WeightedGradingEngine"""
        self.assertIs(GradingEngine, WeightedGradingEngine)
        self.assertIsInstance(self.engine, WeightedGradingEngine)


class TestGovernanceLogging(unittest.TestCase):
    """Governance logging: required fields, sign-off format."""

    def setUp(self):
        self.engine = GradingEngine(verbose=False)

    def test_required_output_fields(self):
        """Result must contain all required fields from spec section 8"""
        d = _make_dossier()
        result = self.engine.grade_dossier(d)
        required = ["grading_version", "grade", "percentage", "total_score",
                     "max_score", "hard_limit_checks", "category_scores",
                     "reasoning", "metrics", "sign_off", "timestamp",
                     "grading_engine_version"]
        for field in required:
            self.assertIn(field, result, f"Missing required field: {field}")

    def test_sign_off_format(self):
        """Sign-off must have engine_id and version for authorization chain"""
        d = _make_dossier()
        result = self.engine.grade_dossier(d)
        so = result["sign_off"]
        self.assertEqual(so["engine_id"], "GRADING_ENGINE_V2")
        self.assertEqual(so["version"], "v2.0-weighted")
        self.assertIn("timestamp", so)
        self.assertIn("authorized", so)

    def test_statistics_tracking(self):
        """Engine should track grading statistics"""
        d1 = _make_dossier(ec50=5, tox=0.003, emax=95, safety_margin=120,
                           herg_risk=0.05, bioavailability_pct=90, half_life_h=8,
                           qed=0.8, synthesis_steps=3)
        d2 = _make_dossier(ec50=5000, tox=0.0205, emax=20, safety_margin=5,
                           herg_risk=0.45, bioavailability_pct=10, half_life_h=0.5,
                           fto=False)

        self.engine.grade_dossier(d1)
        self.engine.grade_dossier(d2)

        stats = self.engine.get_grade_statistics()
        self.assertEqual(stats["total_graded"], 2)
        self.assertIn("grade_distribution", stats)
        self.assertEqual(stats["engine_version"], "v2.0-weighted")


def run_tests():
    """Run grading engine tests"""
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()

    test_classes = [
        TestHardLimits,
        TestWeightedScoring,
        TestMPORecovery,
        TestAntiInflation,
        TestScoreBreakdown,
        TestMissingData,
        TestConstraintViolations,
        TestBackwardCompatibility,
        TestGovernanceLogging,
    ]

    for cls in test_classes:
        suite.addTests(loader.loadTestsFromTestCase(cls))

    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)

    print("\n" + "="*70)
    print(f"GRADING ENGINE TESTS SUMMARY")
    print("="*70)
    print(f"Tests run:     {result.testsRun}")
    print(f"Successes:     {result.testsRun - len(result.failures) - len(result.errors)}")
    print(f"Failures:      {len(result.failures)}")
    print(f"Errors:        {len(result.errors)}")
    print("="*70)

    if result.wasSuccessful():
        print("ALL GRADING ENGINE TESTS PASSED")
        return 0
    else:
        print("SOME TESTS FAILED")
        return 1


if __name__ == "__main__":
    exit_code = run_tests()
    sys.exit(exit_code)
