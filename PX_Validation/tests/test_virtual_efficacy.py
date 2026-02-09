"""
test_virtual_efficacy.py
Unit Tests for Virtual Efficacy Analytics - Phase 5

Tests PTA, responder rates, exposure-response curves, and risk assessment.
Constitutional: All tests verify virtual analytics correctness.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

import unittest
from PX_Engine.operations.VirtualEfficacyAnalytics import (
    compute_pta,
    exposure_response_curve,
    virtual_responder_rate,
    effect_variability_risk,
    time_in_therapeutic_window,
)


class TestPTA(unittest.TestCase):
    """Test Probability of Target Attainment"""
    
    def test_pta_all_above_threshold(self):
        """PTA should be ~1.0 when mean >> threshold"""
        trial_result = {
            "arms": [{
                "arm_id": "A1",
                "n_patients": 10,
                "exposure_summary": {
                    "auc_mg_h_per_L": {
                        "mean": 400.0,
                        "std": 40.0,
                        "min": 320.0,
                        "max": 480.0,
                    }
                }
            }]
        }
        
        pta_result = compute_pta(
            trial_result=trial_result,
            metric="auc_mg_h_per_L",
            target_threshold=200.0,  # Well below mean
        )
        
        self.assertGreaterEqual(pta_result["pta"], 0.9, "PTA should be high when mean >> threshold")
        self.assertEqual(pta_result["n_total"], 10)
        self.assertIn("constitutional", pta_result)
    
    def test_pta_all_below_threshold(self):
        """PTA should be ~0.0 when mean << threshold"""
        trial_result = {
            "arms": [{
                "arm_id": "A1",
                "n_patients": 10,
                "exposure_summary": {
                    "auc_mg_h_per_L": {
                        "mean": 150.0,
                        "std": 15.0,
                        "min": 120.0,
                        "max": 180.0,
                    }
                }
            }]
        }
        
        pta_result = compute_pta(
            trial_result=trial_result,
            metric="auc_mg_h_per_L",
            target_threshold=300.0,  # Well above mean
        )
        
        self.assertLess(pta_result["pta"], 0.1, "PTA should be low when mean << threshold")
    
    def test_pta_mid_threshold(self):
        """PTA should be ~0.5 when threshold ≈ mean"""
        trial_result = {
            "arms": [{
                "arm_id": "A1",
                "n_patients": 10,
                "exposure_summary": {
                    "auc_mg_h_per_L": {
                        "mean": 250.0,
                        "std": 25.0,
                        "min": 200.0,
                        "max": 300.0,
                    }
                }
            }]
        }
        
        pta_result = compute_pta(
            trial_result=trial_result,
            metric="auc_mg_h_per_L",
            target_threshold=250.0,  # Equal to mean
        )
        
        self.assertGreater(pta_result["pta"], 0.3)
        self.assertLess(pta_result["pta"], 0.7)


class TestExposureResponseCurve(unittest.TestCase):
    """Test exposure-response curve generation"""
    
    def test_monotonic_curve(self):
        """Monotonic relationship should be detected"""
        trial_results = [
            {
                "arms": [{
                    "arm_id": "A1",
                    "dose_mg": 50.0,
                    "exposure_summary": {"auc_mg_h_per_L": {"mean": 100.0}},
                    "pd_summary": {"max_effect": {"mean": 0.3}},
                }]
            },
            {
                "arms": [{
                    "arm_id": "A2",
                    "dose_mg": 100.0,
                    "exposure_summary": {"auc_mg_h_per_L": {"mean": 200.0}},
                    "pd_summary": {"max_effect": {"mean": 0.5}},
                }]
            },
            {
                "arms": [{
                    "arm_id": "A3",
                    "dose_mg": 200.0,
                    "exposure_summary": {"auc_mg_h_per_L": {"mean": 400.0}},
                    "pd_summary": {"max_effect": {"mean": 0.7}},
                }]
            },
        ]
        
        curve = exposure_response_curve(
            trial_results=trial_results,
            exposure_metric="auc_mg_h_per_L",
            response_metric="max_effect",
        )
        
        self.assertTrue(curve["monotonic"], "Should detect monotonic relationship")
        self.assertEqual(curve["n_points"], 3)
        self.assertIsNotNone(curve["correlation"])
        self.assertGreater(curve["correlation"], 0.9, "Should have strong positive correlation")
    
    def test_non_monotonic_curve(self):
        """Non-monotonic relationship should be detected"""
        trial_results = [
            {
                "arms": [{
                    "arm_id": "A1",
                    "exposure_summary": {"auc_mg_h_per_L": {"mean": 100.0}},
                    "pd_summary": {"max_effect": {"mean": 0.3}},
                }]
            },
            {
                "arms": [{
                    "arm_id": "A2",
                    "exposure_summary": {"auc_mg_h_per_L": {"mean": 200.0}},
                    "pd_summary": {"max_effect": {"mean": 0.7}},
                }]
            },
            {
                "arms": [{
                    "arm_id": "A3",
                    "exposure_summary": {"auc_mg_h_per_L": {"mean": 400.0}},
                    "pd_summary": {"max_effect": {"mean": 0.6}},  # Decreases
                }]
            },
        ]
        
        curve = exposure_response_curve(
            trial_results=trial_results,
            exposure_metric="auc_mg_h_per_L",
            response_metric="max_effect",
        )
        
        self.assertFalse(curve["monotonic"], "Should detect non-monotonic relationship")


class TestResponderRate(unittest.TestCase):
    """Test virtual responder rate computation"""
    
    def test_responder_rate_calculation(self):
        """Responder rate should be computed correctly"""
        trial_result = {
            "arms": [{
                "arm_id": "A1",
                "n_patients": 10,
                "pd_summary": {
                    "max_effect": {
                        "mean": 0.75,
                        "std": 0.05,
                        "min": 0.65,
                        "max": 0.85,
                    }
                }
            }]
        }
        
        result = virtual_responder_rate(
            trial_result=trial_result,
            response_metric="max_effect",
            responder_threshold=0.70,
        )
        
        self.assertIn("responder_rate", result)
        self.assertIn("n_responders", result)
        self.assertIn("n_total", result)
        self.assertIn("constitutional", result)
        
        # Should have high responder rate (mean=0.75, threshold=0.70)
        self.assertGreater(result["responder_rate"], 0.5)


class TestEffectVariabilityRisk(unittest.TestCase):
    """Test effect variability risk assessment"""
    
    def test_low_variability_low_risk(self):
        """Low CV% should result in LOW risk"""
        trial_result = {
            "arms": [{
                "arm_id": "A1",
                "pd_summary": {
                    "max_effect": {
                        "mean": 0.70,
                        "std": 0.07,  # CV% = 10%
                        "min": 0.63,
                        "max": 0.77,
                    }
                }
            }]
        }
        
        risk = effect_variability_risk(
            trial_result=trial_result,
            effect_metric="max_effect",
            safety_threshold=0.90,
            efficacy_threshold=0.50,
        )
        
        arm_risk = risk["risks_by_arm"][0]
        self.assertEqual(arm_risk["variability_risk"], "LOW")
        self.assertEqual(arm_risk["over_response_risk"], "NONE")
        self.assertEqual(arm_risk["under_response_risk"], "NONE")
    
    def test_high_variability_high_risk(self):
        """High CV% should result in HIGH risk"""
        trial_result = {
            "arms": [{
                "arm_id": "A1",
                "pd_summary": {
                    "max_effect": {
                        "mean": 0.70,
                        "std": 0.30,  # CV% = 43%
                        "min": 0.40,
                        "max": 1.00,
                    }
                }
            }]
        }
        
        risk = effect_variability_risk(
            trial_result=trial_result,
            effect_metric="max_effect",
            safety_threshold=0.85,
            efficacy_threshold=0.55,
        )
        
        arm_risk = risk["risks_by_arm"][0]
        self.assertEqual(arm_risk["variability_risk"], "HIGH")


class TestTimeInWindow(unittest.TestCase):
    """Test time-in-therapeutic-window computation"""
    
    def test_pk_window_calculation(self):
        """Should correctly calculate time in PK window"""
        pk_profile = {
            "time_h": [0, 1, 2, 3, 4, 5],
            "concentration_mg_per_L": [0, 5, 8, 6, 4, 2],
        }
        
        result = time_in_therapeutic_window(
            pk_profile=pk_profile,
            pd_profile=None,
            pk_window=(4.0, 9.0),  # Target range
        )
        
        self.assertIn("pk_time_in_window_h", result)
        self.assertIn("pk_window_fraction", result)
        
        # Time above 4.0: hours 1-4 = 3 hours (out of 5 total)
        time_in = result["pk_time_in_window_h"]
        self.assertGreater(time_in, 2.0)
        self.assertLess(time_in, 5.0)
    
    def test_pd_window_calculation(self):
        """Should correctly calculate time in PD window"""
        pk_profile = {
            "time_h": [0, 1, 2, 3, 4, 5],
            "concentration_mg_per_L": [0, 5, 8, 6, 4, 2],
        }
        
        pd_profile = {
            "time_h": [0, 1, 2, 3, 4, 5],
            "effect": [0, 0.4, 0.7, 0.6, 0.5, 0.3],
        }
        
        result = time_in_therapeutic_window(
            pk_profile=pk_profile,
            pd_profile=pd_profile,
            pd_window=(0.5, 0.8),  # Effect window
        )
        
        self.assertIn("pd_time_in_window_h", result)
        self.assertIn("pd_window_fraction", result)
        
        # Effect in [0.5, 0.8]: hours 1-4
        time_in = result["pd_time_in_window_h"]
        self.assertGreater(time_in, 2.0)


def run_tests():
    """Run all virtual efficacy unit tests"""
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    # Add all test classes
    suite.addTests(loader.loadTestsFromTestCase(TestPTA))
    suite.addTests(loader.loadTestsFromTestCase(TestExposureResponseCurve))
    suite.addTests(loader.loadTestsFromTestCase(TestResponderRate))
    suite.addTests(loader.loadTestsFromTestCase(TestEffectVariabilityRisk))
    suite.addTests(loader.loadTestsFromTestCase(TestTimeInWindow))
    
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    print("\n" + "="*70)
    print(f"VIRTUAL EFFICACY ANALYTICS UNIT TESTS SUMMARY")
    print("="*70)
    print(f"Tests run:     {result.testsRun}")
    print(f"Successes:     {result.testsRun - len(result.failures) - len(result.errors)}")
    print(f"Failures:      {len(result.failures)}")
    print(f"Errors:        {len(result.errors)}")
    print("="*70)
    
    if result.wasSuccessful():
        print("✅ ALL VIRTUAL EFFICACY UNIT TESTS PASSED")
        return 0
    else:
        print("❌ SOME TESTS FAILED")
        return 1


if __name__ == "__main__":
    exit_code = run_tests()
    sys.exit(exit_code)
