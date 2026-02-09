"""
test_pkpd.py
Unit Tests for PK/PD Modeling Module (Phase 1 - v2.0)

Tests Emax models, PD metrics, and PK→PD linking.
Constitutional: All tests verify deterministic behavior.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

import unittest
from PX_Engine.operations.PKPD import emax_model, compute_pd_metrics, link_pk_to_pd


class TestEmaxModel(unittest.TestCase):
    """Test the Emax pharmacodynamic model"""
    
    def test_emax_model_monotonicity(self):
        """Effect increases with concentration when C < EC50"""
        emax = 0.9
        ec50 = 5.0
        hill = 1.0
        
        # Test concentrations below EC50
        c1 = 1.0
        c2 = 2.0
        c3 = 3.0
        
        e1 = emax_model(c1, emax, ec50, hill)
        e2 = emax_model(c2, emax, ec50, hill)
        e3 = emax_model(c3, emax, ec50, hill)
        
        # Effect should increase monotonically
        self.assertLess(e1, e2, "Effect should increase from C=1 to C=2")
        self.assertLess(e2, e3, "Effect should increase from C=2 to C=3")
    
    def test_emax_model_ec50_behavior(self):
        """At C = EC50, effect should be 50% of Emax (when hill=1)"""
        emax = 0.9
        ec50 = 5.0
        hill = 1.0
        
        effect = emax_model(ec50, emax, ec50, hill)
        expected = 0.5 * emax  # 0.45
        
        self.assertAlmostEqual(effect, expected, places=6,
                               msg=f"At EC50, effect should be 50% of Emax (expected {expected}, got {effect})")
    
    def test_emax_model_hill_coefficient(self):
        """Hill > 1 produces steeper curve, Hill < 1 produces shallower curve"""
        emax = 0.9
        ec50 = 5.0
        c = 5.0  # At EC50
        
        # Test different Hill coefficients
        effect_hill_05 = emax_model(c, emax, ec50, hill=0.5)
        effect_hill_10 = emax_model(c, emax, ec50, hill=1.0)
        effect_hill_20 = emax_model(c, emax, ec50, hill=2.0)
        
        # At EC50, all should be approximately 0.5*Emax
        # But the slope differs around EC50
        self.assertAlmostEqual(effect_hill_05, 0.5 * emax, places=2)
        self.assertAlmostEqual(effect_hill_10, 0.5 * emax, places=2)
        self.assertAlmostEqual(effect_hill_20, 0.5 * emax, places=2)
        
        # Test at 2*EC50 - steeper curve should have higher effect
        c_high = 2 * ec50
        effect_hill_05_high = emax_model(c_high, emax, ec50, hill=0.5)
        effect_hill_20_high = emax_model(c_high, emax, ec50, hill=2.0)
        
        self.assertLess(effect_hill_05_high, effect_hill_20_high,
                        "Higher Hill coefficient should produce higher effect at 2*EC50")
    
    def test_emax_model_baseline(self):
        """Baseline effect should be returned at zero concentration"""
        emax = 0.9
        ec50 = 5.0
        baseline = 0.1
        
        effect = emax_model(0.0, emax, ec50, hill=1.0, baseline=baseline)
        
        self.assertEqual(effect, baseline,
                         f"At C=0, effect should equal baseline (expected {baseline}, got {effect})")
    
    def test_zero_concentration_no_baseline(self):
        """Zero concentration without baseline should give zero effect"""
        emax = 0.9
        ec50 = 5.0
        
        effect = emax_model(0.0, emax, ec50, hill=1.0)
        
        self.assertEqual(effect, 0.0,
                         f"At C=0 with no baseline, effect should be 0 (got {effect})")
    
    def test_high_concentration_plateau(self):
        """Very high concentration should approach Emax"""
        emax = 0.9
        ec50 = 5.0
        hill = 1.0
        
        # Test at 100x EC50
        c_very_high = 100 * ec50
        effect = emax_model(c_very_high, emax, ec50, hill)
        
        # Should be very close to Emax (within 1%)
        self.assertGreater(effect, 0.99 * emax,
                           f"At C=100*EC50, effect should be >99% of Emax (got {effect/emax*100:.2f}%)")


class TestPDMetrics(unittest.TestCase):
    """Test PD metrics computation"""
    
    def test_pd_metrics_auec_calculation(self):
        """AUEC should be calculated correctly using trapezoidal rule"""
        # Simple rectangular profile: effect=1.0 for 10 hours
        time_h = [0, 5, 10]
        effect = [1.0, 1.0, 1.0]
        
        metrics = compute_pd_metrics(time_h, effect, effect_threshold=0.5)
        
        # AUEC = 1.0 * 10 = 10.0
        expected_auec = 10.0
        self.assertAlmostEqual(metrics["auec_h"], expected_auec, places=2,
                               msg=f"AUEC should be {expected_auec} (got {metrics['auec_h']})")
    
    def test_pd_metrics_time_above_threshold(self):
        """Time above threshold should be calculated correctly"""
        time_h = [0, 1, 2, 3, 4, 5]
        effect = [0.2, 0.4, 0.8, 0.9, 0.7, 0.3]  # Above 0.5 from t=2 to t=5
        threshold = 0.5
        
        metrics = compute_pd_metrics(time_h, effect, effect_threshold=threshold)
        
        # Time above threshold: from t=2 to t=5 = ~3-4 hours
        # (Includes any interval where at least one endpoint is above threshold)
        self.assertGreater(metrics["time_above_threshold_h"], 2.5,
                           "Time above threshold should be >2.5 hours")
        self.assertLessEqual(metrics["time_above_threshold_h"], 4.0,
                             "Time above threshold should be ≤4 hours")
    
    def test_pd_metrics_max_effect_detection(self):
        """Max effect should be detected correctly"""
        time_h = [0, 1, 2, 3, 4]
        effect = [0.1, 0.3, 0.9, 0.5, 0.2]  # Max at t=2
        
        metrics = compute_pd_metrics(time_h, effect, effect_threshold=0.5)
        
        self.assertEqual(metrics["max_effect"], 0.9,
                         f"Max effect should be 0.9 (got {metrics['max_effect']})")
        self.assertEqual(metrics["time_to_max_effect_h"], 2.0,
                         f"Time to max effect should be 2.0h (got {metrics['time_to_max_effect_h']})")
    
    def test_pd_metrics_empty_profile(self):
        """Empty profile should return zero metrics without crashing"""
        time_h = []
        effect = []
        
        metrics = compute_pd_metrics(time_h, effect)
        
        self.assertEqual(metrics["max_effect"], 0.0)
        self.assertEqual(metrics["auec_h"], 0.0)
        self.assertEqual(metrics["time_above_threshold_h"], 0.0)


class TestLinkPKtoPD(unittest.TestCase):
    """Test PK→PD linking function"""
    
    def test_link_pk_to_pd_basic(self):
        """PK profile should convert to PD profile correctly"""
        # Simple PK profile
        pk_profile = {
            "time_grid_h": [0, 1, 2, 3, 4],
            "concentration_mg_per_L": [0, 5.0, 10.0, 5.0, 2.5],
            "summary": {}
        }
        
        pd_params = {
            "emax": 0.9,
            "ec50": 5.0,
            "hill": 1.0
        }
        
        result = link_pk_to_pd(pk_profile, pd_params)
        
        # Verify structure
        self.assertIn("time_h", result)
        self.assertIn("effect", result)
        self.assertIn("pd_summary", result)
        self.assertIn("parameters", result)
        self.assertIn("constitutional", result)
        
        # Verify effect profile length matches PK profile
        self.assertEqual(len(result["effect"]), len(pk_profile["concentration_mg_per_L"]))
        
        # Verify effect at EC50 is approximately 50% of Emax
        # At t=1, C=5.0 (which is EC50)
        effect_at_ec50 = result["effect"][1]
        self.assertAlmostEqual(effect_at_ec50, 0.45, places=2,
                               msg=f"Effect at EC50 should be ~0.45 (got {effect_at_ec50})")
    
    def test_link_pk_to_pd_missing_params(self):
        """Should raise ValueError if emax or ec50 missing"""
        pk_profile = {
            "time_grid_h": [0, 1, 2],
            "concentration_mg_per_L": [0, 5.0, 10.0],
        }
        
        # Missing emax
        pd_params_no_emax = {"ec50": 5.0}
        with self.assertRaises(ValueError):
            link_pk_to_pd(pk_profile, pd_params_no_emax)
        
        # Missing ec50
        pd_params_no_ec50 = {"emax": 0.9}
        with self.assertRaises(ValueError):
            link_pk_to_pd(pk_profile, pd_params_no_ec50)
    
    def test_link_pk_to_pd_pd_summary_metrics(self):
        """PD summary should contain all required metrics"""
        pk_profile = {
            "time_grid_h": [0, 1, 2, 3, 4],
            "concentration_mg_per_L": [0, 5.0, 10.0, 5.0, 2.5],
            "summary": {}
        }
        
        pd_params = {
            "emax": 0.9,
            "ec50": 5.0,
            "hill": 1.0,
            "effect_threshold": 0.5
        }
        
        result = link_pk_to_pd(pk_profile, pd_params)
        pd_summary = result["pd_summary"]
        
        # Verify all metrics present
        required_metrics = [
            "max_effect",
            "time_to_max_effect_h",
            "auec_h",
            "time_above_threshold_h",
            "mean_effect",
            "effect_at_steady_state"
        ]
        
        for metric in required_metrics:
            self.assertIn(metric, pd_summary,
                          f"PD summary should contain '{metric}'")
    
    def test_link_pk_to_pd_constitutional_compliance(self):
        """Result should include constitutional metadata"""
        pk_profile = {
            "time_grid_h": [0, 1, 2],
            "concentration_mg_per_L": [0, 5.0, 10.0],
            "summary": {}
        }
        
        pd_params = {"emax": 0.9, "ec50": 5.0}
        
        result = link_pk_to_pd(pk_profile, pd_params)
        
        # Verify constitutional block
        self.assertIn("constitutional", result)
        const = result["constitutional"]
        
        self.assertEqual(const["status"], "SIMULATED")
        self.assertIn("engine", const)
        self.assertIn("notes", const)
        self.assertIn("clinical validation required", const["notes"].lower())


def run_tests():
    """Run all PKPD unit tests"""
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    # Add all test classes
    suite.addTests(loader.loadTestsFromTestCase(TestEmaxModel))
    suite.addTests(loader.loadTestsFromTestCase(TestPDMetrics))
    suite.addTests(loader.loadTestsFromTestCase(TestLinkPKtoPD))
    
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    print("\n" + "="*70)
    print(f"PKPD UNIT TESTS SUMMARY")
    print("="*70)
    print(f"Tests run:     {result.testsRun}")
    print(f"Successes:     {result.testsRun - len(result.failures) - len(result.errors)}")
    print(f"Failures:      {len(result.failures)}")
    print(f"Errors:        {len(result.errors)}")
    print("="*70)
    
    if result.wasSuccessful():
        print("✅ ALL PKPD UNIT TESTS PASSED")
        return 0
    else:
        print("❌ SOME TESTS FAILED")
        return 1


if __name__ == "__main__":
    exit_code = run_tests()
    sys.exit(exit_code)
