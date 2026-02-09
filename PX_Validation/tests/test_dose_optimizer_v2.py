"""
test_dose_optimizer_v2.py
Unit Tests for Dose Optimization v2 - Phase 4

Tests advanced dose optimization with multi-dimensional search.
Constitutional: All tests verify virtual trial-based optimization.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

import unittest
from PX_Engine.operations.DoseOptimizer_v2 import (
    scoring_function,
    is_monotonic_metric,
    evaluate_regimen,
)


class TestScoringFunction(unittest.TestCase):
    """Test regimen scoring logic"""
    
    def test_perfect_score_within_range(self):
        """Score should be 0.0 when mean is within target range"""
        pk_summary = {
            "auc_mg_h_per_L": {
                "mean": 300.0,
                "std": 30.0,
                "cv_percent": 10.0,
                "min": 250.0,
                "max": 350.0,
            }
        }
        
        target_pk_range = {
            "auc_mg_h_per_L": (250.0, 350.0)
        }
        
        score = scoring_function(
            pk_summary=pk_summary,
            pd_summary=None,
            target_pk_range=target_pk_range,
            target_pd_range=None,
        )
        
        self.assertAlmostEqual(score, 0.0, places=2, msg="Perfect score when within range")
    
    def test_penalty_below_range(self):
        """Score should increase when below target range"""
        pk_summary = {
            "auc_mg_h_per_L": {
                "mean": 150.0,  # Below 200
                "std": 15.0,
                "cv_percent": 10.0,
                "min": 130.0,
                "max": 170.0,
            }
        }
        
        target_pk_range = {
            "auc_mg_h_per_L": (200.0, 400.0)
        }
        
        score = scoring_function(
            pk_summary=pk_summary,
            pd_summary=None,
            target_pk_range=target_pk_range,
            target_pd_range=None,
        )
        
        self.assertGreater(score, 0.0, "Should penalize below-range values")
    
    def test_penalty_above_range(self):
        """Score should increase when above target range"""
        pk_summary = {
            "auc_mg_h_per_L": {
                "mean": 500.0,  # Above 400
                "std": 50.0,
                "cv_percent": 10.0,
                "min": 450.0,
                "max": 550.0,
            }
        }
        
        target_pk_range = {
            "auc_mg_h_per_L": (200.0, 400.0)
        }
        
        score = scoring_function(
            pk_summary=pk_summary,
            pd_summary=None,
            target_pk_range=target_pk_range,
            target_pd_range=None,
        )
        
        self.assertGreater(score, 0.0, "Should penalize above-range values")
    
    def test_variability_penalty(self):
        """High CV% should increase score"""
        pk_high_cv = {
            "auc_mg_h_per_L": {
                "mean": 300.0,
                "std": 120.0,
                "cv_percent": 40.0,  # High variability
                "min": 180.0,
                "max": 420.0,
            }
        }
        
        pk_low_cv = {
            "auc_mg_h_per_L": {
                "mean": 300.0,
                "std": 30.0,
                "cv_percent": 10.0,  # Low variability
                "min": 270.0,
                "max": 330.0,
            }
        }
        
        target_pk_range = {
            "auc_mg_h_per_L": (250.0, 350.0)
        }
        
        score_high_cv = scoring_function(
            pk_summary=pk_high_cv,
            pd_summary=None,
            target_pk_range=target_pk_range,
            target_pd_range=None,
        )
        
        score_low_cv = scoring_function(
            pk_summary=pk_low_cv,
            pd_summary=None,
            target_pk_range=target_pk_range,
            target_pd_range=None,
        )
        
        self.assertGreater(score_high_cv, score_low_cv, "High CV% should have worse score")
    
    def test_pd_scoring(self):
        """PD metrics should be scored correctly"""
        pk_summary = {
            "auc_mg_h_per_L": {
                "mean": 300.0,
                "std": 30.0,
                "cv_percent": 10.0,
                "min": 270.0,
                "max": 330.0,
            }
        }
        
        pd_summary = {
            "max_effect": {
                "mean": 0.75,
                "std": 0.05,
                "cv_percent": 6.7,
                "min": 0.70,
                "max": 0.80,
            }
        }
        
        target_pk_range = {
            "auc_mg_h_per_L": (250.0, 350.0)
        }
        
        target_pd_range = {
            "max_effect": (0.70, 0.80)
        }
        
        score = scoring_function(
            pk_summary=pk_summary,
            pd_summary=pd_summary,
            target_pk_range=target_pk_range,
            target_pd_range=target_pd_range,
        )
        
        self.assertAlmostEqual(score, 0.0, places=2, msg="Both PK and PD within range")
    
    def test_multi_metric_scoring(self):
        """Multiple PK metrics should all contribute to score"""
        pk_summary = {
            "auc_mg_h_per_L": {
                "mean": 300.0,
                "std": 30.0,
                "cv_percent": 10.0,
                "min": 270.0,
                "max": 330.0,
            },
            "cmax_mg_per_L": {
                "mean": 15.0,
                "std": 1.5,
                "cv_percent": 10.0,
                "min": 13.5,
                "max": 16.5,
            }
        }
        
        target_pk_range = {
            "auc_mg_h_per_L": (250.0, 350.0),
            "cmax_mg_per_L": (10.0, 20.0),
        }
        
        score = scoring_function(
            pk_summary=pk_summary,
            pd_summary=None,
            target_pk_range=target_pk_range,
            target_pd_range=None,
        )
        
        self.assertAlmostEqual(score, 0.0, places=2, msg="All metrics within range")


class TestMonotonicityCheck(unittest.TestCase):
    """Test monotonic metric identification"""
    
    def test_pk_monotonic_metrics(self):
        """Known PK metrics should be identified as monotonic"""
        self.assertTrue(is_monotonic_metric("auc_mg_h_per_L"))
        self.assertTrue(is_monotonic_metric("cmax_mg_per_L"))
        self.assertTrue(is_monotonic_metric("cmin_steady_state_mg_per_L"))
    
    def test_pd_monotonic_metrics(self):
        """Known PD metrics should be identified as monotonic"""
        self.assertTrue(is_monotonic_metric("max_effect"))
        self.assertTrue(is_monotonic_metric("auec_h"))
        self.assertTrue(is_monotonic_metric("mean_effect"))
    
    def test_non_monotonic_metric(self):
        """Unknown metrics should not be assumed monotonic"""
        self.assertFalse(is_monotonic_metric("unknown_metric"))
        self.assertFalse(is_monotonic_metric("time_above_threshold_h"))  # Can be non-monotonic


class TestEvaluateRegimen(unittest.TestCase):
    """Test regimen evaluation via mini-trial"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.admet = {
            "distribution": {"predicted_vd_L_per_kg": 0.7},
            "metabolism": {"predicted_clearance_L_per_h_per_kg": 0.05},
            "absorption": {"predicted_bioavailability": 1.0}
        }
        
        self.protocol_template = {
            "trial_id": "TEST",
            "duration_days": 3.0,
        }
    
    def test_evaluate_regimen_basic(self):
        """Basic regimen evaluation should return expected structure"""
        result = evaluate_regimen(
            dose_mg=100.0,
            interval_h=24.0,
            admet=self.admet,
            protocol_template=self.protocol_template,
            n_patients=5,
            duration_days=3.0,
        )
        
        # Check structure
        self.assertIn("dose_mg", result)
        self.assertIn("interval_h", result)
        self.assertIn("pk_summary", result)
        self.assertIn("pd_summary", result)
        self.assertIn("trial_result", result)
        
        # Check values
        self.assertEqual(result["dose_mg"], 100.0)
        self.assertEqual(result["interval_h"], 24.0)
        
        # Check PK summary has required fields
        for metric in ["auc_mg_h_per_L", "cmax_mg_per_L", "cmin_steady_state_mg_per_L"]:
            self.assertIn(metric, result["pk_summary"])
            self.assertIn("mean", result["pk_summary"][metric])
            self.assertIn("std", result["pk_summary"][metric])
            self.assertIn("cv_percent", result["pk_summary"][metric])
    
    def test_evaluate_regimen_with_pd(self):
        """Regimen evaluation with PD should include pd_summary"""
        pd_params = {
            "emax": 0.9,
            "ec50": 5.0,
            "hill": 1.5,
        }
        
        result = evaluate_regimen(
            dose_mg=150.0,
            interval_h=12.0,
            admet=self.admet,
            protocol_template=self.protocol_template,
            pd_params=pd_params,
            n_patients=5,
            duration_days=3.0,
        )
        
        # Should have PD summary
        self.assertIsNotNone(result["pd_summary"])
        self.assertIn("max_effect", result["pd_summary"])
        self.assertIn("mean", result["pd_summary"]["max_effect"])
        self.assertIn("cv_percent", result["pd_summary"]["max_effect"])
    
    def test_evaluate_regimen_with_iiv(self):
        """Regimen evaluation with IIV should show variability"""
        variability = {
            "clearance_variation": 0.3,
            "vd_variation": 0.25,
            "n_tiers": 7,
        }
        
        result = evaluate_regimen(
            dose_mg=100.0,
            interval_h=24.0,
            admet=self.admet,
            protocol_template=self.protocol_template,
            variability=variability,
            n_patients=10,
            duration_days=3.0,
        )
        
        # Should have non-zero std
        auc_std = result["pk_summary"]["auc_mg_h_per_L"]["std"]
        self.assertGreater(auc_std, 5.0, "IIV should create variability")
    
    def test_cv_percent_calculation(self):
        """CV% should be calculated correctly"""
        result = evaluate_regimen(
            dose_mg=100.0,
            interval_h=24.0,
            admet=self.admet,
            protocol_template=self.protocol_template,
            n_patients=5,
            duration_days=3.0,
        )
        
        # Check CV% calculation
        for metric, values in result["pk_summary"].items():
            mean = values["mean"]
            std = values["std"]
            cv = values["cv_percent"]
            
            if mean > 0:
                expected_cv = (std / mean) * 100
                self.assertAlmostEqual(cv, expected_cv, places=1,
                                       msg=f"CV% calculation for {metric}")


class TestScoringSensitivity(unittest.TestCase):
    """Test scoring function sensitivity"""
    
    def test_closer_to_target_better_score(self):
        """Regimen closer to target should have better score"""
        target_pk_range = {
            "auc_mg_h_per_L": (300.0, 400.0)
        }
        
        # Regimen 1: Close to target
        pk1 = {
            "auc_mg_h_per_L": {
                "mean": 350.0,
                "std": 35.0,
                "cv_percent": 10.0,
                "min": 315.0,
                "max": 385.0,
            }
        }
        
        # Regimen 2: Far from target
        pk2 = {
            "auc_mg_h_per_L": {
                "mean": 150.0,
                "std": 15.0,
                "cv_percent": 10.0,
                "min": 135.0,
                "max": 165.0,
            }
        }
        
        score1 = scoring_function(pk1, None, target_pk_range, None)
        score2 = scoring_function(pk2, None, target_pk_range, None)
        
        self.assertLess(score1, score2, "Closer to target should have better (lower) score")


def run_tests():
    """Run all dose optimizer v2 unit tests"""
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    # Add all test classes
    suite.addTests(loader.loadTestsFromTestCase(TestScoringFunction))
    suite.addTests(loader.loadTestsFromTestCase(TestMonotonicityCheck))
    suite.addTests(loader.loadTestsFromTestCase(TestEvaluateRegimen))
    suite.addTests(loader.loadTestsFromTestCase(TestScoringSensitivity))
    
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    print("\n" + "="*70)
    print(f"DOSE OPTIMIZER V2 UNIT TESTS SUMMARY")
    print("="*70)
    print(f"Tests run:     {result.testsRun}")
    print(f"Successes:     {result.testsRun - len(result.failures) - len(result.errors)}")
    print(f"Failures:      {len(result.failures)}")
    print(f"Errors:        {len(result.errors)}")
    print("="*70)
    
    if result.wasSuccessful():
        print("✅ ALL DOSE OPTIMIZER V2 UNIT TESTS PASSED")
        return 0
    else:
        print("❌ SOME TESTS FAILED")
        return 1


if __name__ == "__main__":
    exit_code = run_tests()
    sys.exit(exit_code)
