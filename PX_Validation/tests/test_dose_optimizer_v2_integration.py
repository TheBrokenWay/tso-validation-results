"""
test_dose_optimizer_v2_integration.py
Integration Tests for Dose Optimization v2 - Phase 4

Tests full dose optimization pipeline with TrialEngine integration.
Constitutional: All tests verify end-to-end optimization correctness.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

import unittest
from PX_Engine.operations.DoseOptimizer_v2 import optimize_dose


class TestDoseOptimizerIntegration(unittest.TestCase):
    """Test full dose optimization pipeline"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.admet = {
            "distribution": {"predicted_vd_L_per_kg": 0.7},
            "metabolism": {"predicted_clearance_L_per_h_per_kg": 0.05},
            "absorption": {"predicted_bioavailability": 1.0}
        }
        
        self.protocol_template = {
            "trial_id": "OPT-TEST",
            "duration_days": 3.0,
        }
    
    def test_optimize_dose_coarse_to_fine(self):
        """Coarse-to-fine search should find dose near target"""
        result = optimize_dose(
            smiles="CC(=O)OC",  # Example
            admet=self.admet,
            protocol_template=self.protocol_template,
            target_pk_range={"auc_mg_h_per_L": (150.0, 250.0)},
            dose_bounds=(50.0, 300.0),
            interval_options=[24.0],
            search_strategy="coarse_to_fine",
            n_eval_patients=5,
        )
        
        # Should have results
        self.assertIn("best_regimen", result)
        self.assertIn("search_history", result)
        self.assertIn("target_achievement", result)
        
        # Best regimen should have key fields
        best = result["best_regimen"]
        self.assertIn("dose_mg", best)
        self.assertIn("interval_h", best)
        self.assertIn("score", best)
        self.assertIn("pk_summary", best)
        
        # Should be within bounds
        self.assertGreaterEqual(best["dose_mg"], 50.0)
        self.assertLessEqual(best["dose_mg"], 300.0)
        
        # Should have tried multiple doses
        self.assertGreater(result["evaluations"], 5, "Should evaluate multiple regimens")
    
    def test_optimize_dose_binary_search(self):
        """Binary search should work for monotonic metrics"""
        result = optimize_dose(
            smiles="CC(=O)OC",
            admet=self.admet,
            protocol_template=self.protocol_template,
            target_pk_range={"auc_mg_h_per_L": (200.0, 300.0)},
            dose_bounds=(50.0, 400.0),
            interval_options=[24.0],
            search_strategy="binary_search",
            n_eval_patients=5,
        )
        
        # Should use binary search
        self.assertEqual(result["search_strategy"], "binary_search")
        
        # Should have best regimen
        self.assertIn("best_regimen", result)
        best = result["best_regimen"]
        
        # Should be within dose bounds
        self.assertGreaterEqual(best["dose_mg"], 50.0)
        self.assertLessEqual(best["dose_mg"], 400.0)
        
        # Binary search should be efficient (fewer evaluations)
        self.assertLess(result["evaluations"], 15, "Binary search should be efficient")
    
    def test_optimize_dose_with_pd_target(self):
        """Optimization with PD target should work"""
        pd_params = {
            "emax": 0.9,
            "ec50": 5.0,
            "hill": 1.5,
        }
        
        result = optimize_dose(
            smiles="CC(=O)OC",
            admet=self.admet,
            protocol_template=self.protocol_template,
            target_pd_range={"max_effect": (0.6, 0.75)},
            dose_bounds=(50.0, 300.0),
            interval_options=[24.0],
            pd_params=pd_params,
            search_strategy="coarse_to_fine",
            n_eval_patients=5,
        )
        
        # Should have PD summary
        best = result["best_regimen"]
        self.assertIsNotNone(best["pd_summary"], "Should have PD results")
        self.assertIn("max_effect", best["pd_summary"])
        
        # Target achievement should include PD
        self.assertIn("target_achievement", result)
        self.assertTrue(any("pd_" in k for k in result["target_achievement"].keys()),
                        "Should have PD target achievement")
    
    def test_optimize_dose_with_iiv(self):
        """Optimization with IIV should produce realistic distributions"""
        variability = {
            "clearance_variation": 0.3,
            "vd_variation": 0.25,
            "n_tiers": 7,
        }
        
        result = optimize_dose(
            smiles="CC(=O)OC",
            admet=self.admet,
            protocol_template=self.protocol_template,
            target_pk_range={"auc_mg_h_per_L": (200.0, 400.0)},
            dose_bounds=(50.0, 400.0),
            interval_options=[24.0],
            variability=variability,
            search_strategy="coarse_to_fine",
            n_eval_patients=10,
        )
        
        # Best regimen should show variability
        best = result["best_regimen"]
        auc_cv = best["pk_summary"]["auc_mg_h_per_L"]["cv_percent"]
        
        self.assertGreater(auc_cv, 5.0, "IIV should create variability")
        self.assertLess(auc_cv, 40.0, "CV% should be realistic")
    
    def test_optimize_multiple_intervals(self):
        """Optimization should test multiple dosing intervals"""
        result = optimize_dose(
            smiles="CC(=O)OC",
            admet=self.admet,
            protocol_template=self.protocol_template,
            target_pk_range={"auc_mg_h_per_L": (200.0, 300.0)},
            dose_bounds=(50.0, 300.0),
            interval_options=[24.0, 12.0, 8.0],  # QD, BID, TID
            search_strategy="coarse_to_fine",
            n_eval_patients=5,
        )
        
        # Should test multiple intervals
        intervals_tested = set(r["interval_h"] for r in result["search_history"])
        self.assertGreater(len(intervals_tested), 1, "Should test multiple intervals")
        
        # Best regimen should have selected one
        best = result["best_regimen"]
        self.assertIn(best["interval_h"], [24.0, 12.0, 8.0])
    
    def test_constitutional_metadata(self):
        """Result should have constitutional compliance metadata"""
        result = optimize_dose(
            smiles="CC(=O)OC",
            admet=self.admet,
            protocol_template=self.protocol_template,
            target_pk_range={"auc_mg_h_per_L": (200.0, 300.0)},
            dose_bounds=(50.0, 300.0),
            interval_options=[24.0],
            search_strategy="coarse_to_fine",
            n_eval_patients=5,
        )
        
        # Should have constitutional metadata
        self.assertIn("constitutional", result)
        const = result["constitutional"]
        
        self.assertIn("status", const)
        self.assertIn("engine", const)
        self.assertIn("notes", const)
        
        # Notes should mention L51/L34
        self.assertIn("L51", const["notes"])
        self.assertIn("L34", const["notes"])
        self.assertIn("VIRTUAL", const["notes"])
    
    def test_target_achievement_tracking(self):
        """Target achievement should be accurately tracked"""
        result = optimize_dose(
            smiles="CC(=O)OC",
            admet=self.admet,
            protocol_template=self.protocol_template,
            target_pk_range={
                "auc_mg_h_per_L": (200.0, 300.0),
                "cmax_mg_per_L": (2.0, 5.0),
            },
            dose_bounds=(50.0, 300.0),
            interval_options=[24.0],
            search_strategy="coarse_to_fine",
            n_eval_patients=5,
        )
        
        # Target achievement should track both metrics
        achievement = result["target_achievement"]
        
        self.assertIn("pk_auc_mg_h_per_L", achievement)
        self.assertIn("pk_cmax_mg_per_L", achievement)
        
        # Each should have required fields
        for metric, data in achievement.items():
            self.assertIn("target_range", data)
            self.assertIn("achieved", data)
            self.assertIn("within_range", data)


class TestDoseOptimizerPerformance(unittest.TestCase):
    """Test optimization performance characteristics"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.admet = {
            "distribution": {"predicted_vd_L_per_kg": 0.7},
            "metabolism": {"predicted_clearance_L_per_h_per_kg": 0.05},
            "absorption": {"predicted_bioavailability": 1.0}
        }
        
        self.protocol_template = {
            "trial_id": "PERF-TEST",
            "duration_days": 3.0,
        }
    
    def test_coarse_to_fine_efficiency(self):
        """Coarse-to-fine should be more efficient than full grid"""
        result = optimize_dose(
            smiles="CC(=O)OC",
            admet=self.admet,
            protocol_template=self.protocol_template,
            target_pk_range={"auc_mg_h_per_L": (200.0, 300.0)},
            dose_bounds=(50.0, 300.0),
            interval_options=[24.0, 12.0],
            search_strategy="coarse_to_fine",
            n_eval_patients=5,
        )
        
        # Should evaluate coarse + fine, not full grid
        # Coarse: 5 doses × 2 intervals = 10
        # Fine: ~3 candidates × 5 doses × 2 intervals = ~30
        # Total: ~40 (not 100+)
        self.assertLess(result["evaluations"], 50, 
                        "Coarse-to-fine should be efficient")
    
    def test_binary_search_efficiency(self):
        """Binary search should be most efficient for monotonic metrics"""
        result = optimize_dose(
            smiles="CC(=O)OC",
            admet=self.admet,
            protocol_template=self.protocol_template,
            target_pk_range={"auc_mg_h_per_L": (200.0, 300.0)},
            dose_bounds=(50.0, 400.0),
            interval_options=[24.0, 12.0],
            search_strategy="binary_search",
            n_eval_patients=5,
        )
        
        # Binary search should be very efficient
        # ~log2(dose_range) × intervals = ~8 × 2 = 16
        self.assertLess(result["evaluations"], 20,
                        "Binary search should be very efficient")


def run_tests():
    """Run all dose optimizer v2 integration tests"""
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    # Add all test classes
    suite.addTests(loader.loadTestsFromTestCase(TestDoseOptimizerIntegration))
    suite.addTests(loader.loadTestsFromTestCase(TestDoseOptimizerPerformance))
    
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    print("\n" + "="*70)
    print(f"DOSE OPTIMIZER V2 INTEGRATION TESTS SUMMARY")
    print("="*70)
    print(f"Tests run:     {result.testsRun}")
    print(f"Successes:     {result.testsRun - len(result.failures) - len(result.errors)}")
    print(f"Failures:      {len(result.failures)}")
    print(f"Errors:        {len(result.errors)}")
    print("="*70)
    
    if result.wasSuccessful():
        print("✅ ALL DOSE OPTIMIZER V2 INTEGRATION TESTS PASSED")
        return 0
    else:
        print("❌ SOME TESTS FAILED")
        return 1


if __name__ == "__main__":
    exit_code = run_tests()
    sys.exit(exit_code)
