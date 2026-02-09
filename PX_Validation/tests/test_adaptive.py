"""
test_adaptive.py
Unit Tests for Adaptive Trial Logic - Phase 3 v2.0

Tests epoch-based adaptation with REDUCE_DOSE, INCREASE_DOSE, STOP_ARM actions.
Constitutional: All tests verify deterministic, logged behavior.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

import unittest
from PX_Engine.operations.TrialEngine import TrialEngine


class TestAdaptiveRuleEvaluation(unittest.TestCase):
    """Test adaptive rule evaluation logic"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.trial_engine = TrialEngine()
    
    def test_upper_bound_triggers_reduce_dose(self):
        """Upper bound violation should trigger REDUCE_DOSE"""
        adaptive_rules = {
            "metric": "auc_mg_h_per_L",
            "upper_bound": 300.0,
            "action": "REDUCE_DOSE",
            "dose_adjustment_factor": 0.75,
        }
        
        # Simulate high AUC values
        auc_values = [350.0, 360.0, 340.0, 355.0]
        
        decision = self.trial_engine._evaluate_adaptive_rule(
            adaptive_rules=adaptive_rules,
            cmax_values=[],
            auc_values=auc_values,
            cmin_values=[],
            max_effect_values=[],
            auec_values=[],
            time_above_threshold_values=[],
            mean_effect_values=[],
            current_dose_mg=100.0,
            patient_count=4,
        )
        
        self.assertTrue(decision["triggered"], "Rule should trigger for high AUC")
        self.assertEqual(decision["action"], "REDUCE_DOSE")
        self.assertEqual(decision["new_dose_mg"], 75.0, "Dose should be reduced to 75mg (100*0.75)")
        self.assertIn("upper_bound", decision["reason"].lower())
    
    def test_lower_bound_triggers_increase_dose(self):
        """Lower bound violation should trigger INCREASE_DOSE"""
        adaptive_rules = {
            "metric": "auc_mg_h_per_L",
            "lower_bound": 100.0,
            "action": "INCREASE_DOSE",
            "increase_factor": 1.33,
        }
        
        # Simulate low AUC values
        auc_values = [70.0, 75.0, 68.0, 72.0]
        
        decision = self.trial_engine._evaluate_adaptive_rule(
            adaptive_rules=adaptive_rules,
            cmax_values=[],
            auc_values=auc_values,
            cmin_values=[],
            max_effect_values=[],
            auec_values=[],
            time_above_threshold_values=[],
            mean_effect_values=[],
            current_dose_mg=100.0,
            patient_count=4,
        )
        
        self.assertTrue(decision["triggered"], "Rule should trigger for low AUC")
        self.assertEqual(decision["action"], "INCREASE_DOSE")
        self.assertEqual(decision["new_dose_mg"], 133.0, "Dose should increase to 133mg (100*1.33)")
        self.assertIn("lower_bound", decision["reason"].lower())
    
    def test_within_bounds_no_trigger(self):
        """Metric within bounds should not trigger"""
        adaptive_rules = {
            "metric": "auc_mg_h_per_L",
            "lower_bound": 100.0,
            "upper_bound": 300.0,
            "action": "REDUCE_DOSE",
        }
        
        # Simulate values within bounds
        auc_values = [200.0, 210.0, 195.0, 205.0]
        
        decision = self.trial_engine._evaluate_adaptive_rule(
            adaptive_rules=adaptive_rules,
            cmax_values=[],
            auc_values=auc_values,
            cmin_values=[],
            max_effect_values=[],
            auec_values=[],
            time_above_threshold_values=[],
            mean_effect_values=[],
            current_dose_mg=100.0,
            patient_count=4,
        )
        
        self.assertFalse(decision["triggered"], "Rule should not trigger for values within bounds")
        self.assertIsNone(decision["action"])
        self.assertEqual(decision["new_dose_mg"], 100.0, "Dose should remain unchanged")
        self.assertIn("within bounds", decision["reason"].lower())
    
    def test_upper_bound_triggers_stop_arm(self):
        """Upper bound with STOP_ARM action"""
        adaptive_rules = {
            "metric": "cmax_mg_per_L",
            "upper_bound": 10.0,
            "action": "STOP_ARM",
        }
        
        # Simulate high Cmax
        cmax_values = [12.0, 13.0, 11.5, 12.5]
        
        decision = self.trial_engine._evaluate_adaptive_rule(
            adaptive_rules=adaptive_rules,
            cmax_values=cmax_values,
            auc_values=[],
            cmin_values=[],
            max_effect_values=[],
            auec_values=[],
            time_above_threshold_values=[],
            mean_effect_values=[],
            current_dose_mg=100.0,
            patient_count=4,
        )
        
        self.assertTrue(decision["triggered"])
        self.assertEqual(decision["action"], "STOP_ARM")
        self.assertEqual(decision["new_dose_mg"], 100.0, "Dose unchanged when stopping")
    
    def test_pd_metric_evaluation(self):
        """Test adaptive rule on PD metrics (max_effect)"""
        adaptive_rules = {
            "metric": "max_effect",
            "upper_bound": 0.8,
            "action": "REDUCE_DOSE",
            "dose_adjustment_factor": 0.8,
        }
        
        # Simulate high effect
        max_effect_values = [0.85, 0.87, 0.83, 0.86]
        
        decision = self.trial_engine._evaluate_adaptive_rule(
            adaptive_rules=adaptive_rules,
            cmax_values=[],
            auc_values=[],
            cmin_values=[],
            max_effect_values=max_effect_values,
            auec_values=[],
            time_above_threshold_values=[],
            mean_effect_values=[],
            current_dose_mg=100.0,
            patient_count=4,
        )
        
        self.assertTrue(decision["triggered"])
        self.assertEqual(decision["action"], "REDUCE_DOSE")
        self.assertEqual(decision["new_dose_mg"], 80.0)
    
    def test_unknown_metric_no_trigger(self):
        """Unknown metric should not trigger"""
        adaptive_rules = {
            "metric": "invalid_metric",
            "upper_bound": 100.0,
            "action": "REDUCE_DOSE",
        }
        
        decision = self.trial_engine._evaluate_adaptive_rule(
            adaptive_rules=adaptive_rules,
            cmax_values=[],
            auc_values=[200.0],
            cmin_values=[],
            max_effect_values=[],
            auec_values=[],
            time_above_threshold_values=[],
            mean_effect_values=[],
            current_dose_mg=100.0,
            patient_count=1,
        )
        
        self.assertFalse(decision["triggered"])
        self.assertIn("Unknown metric", decision["reason"])
    
    def test_no_data_no_trigger(self):
        """No data should not trigger"""
        adaptive_rules = {
            "metric": "auc_mg_h_per_L",
            "upper_bound": 300.0,
            "action": "REDUCE_DOSE",
        }
        
        decision = self.trial_engine._evaluate_adaptive_rule(
            adaptive_rules=adaptive_rules,
            cmax_values=[],
            auc_values=[],  # Empty
            cmin_values=[],
            max_effect_values=[],
            auec_values=[],
            time_above_threshold_values=[],
            mean_effect_values=[],
            current_dose_mg=100.0,
            patient_count=0,
        )
        
        self.assertFalse(decision["triggered"])
        self.assertIn("No data", decision["reason"])


def run_tests():
    """Run all adaptive unit tests"""
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    # Add all test classes
    suite.addTests(loader.loadTestsFromTestCase(TestAdaptiveRuleEvaluation))
    
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    print("\n" + "="*70)
    print(f"ADAPTIVE UNIT TESTS SUMMARY")
    print("="*70)
    print(f"Tests run:     {result.testsRun}")
    print(f"Successes:     {result.testsRun - len(result.failures) - len(result.errors)}")
    print(f"Failures:      {len(result.failures)}")
    print(f"Errors:        {len(result.errors)}")
    print("="*70)
    
    if result.wasSuccessful():
        print("✅ ALL ADAPTIVE UNIT TESTS PASSED")
        return 0
    else:
        print("❌ SOME TESTS FAILED")
        return 1


if __name__ == "__main__":
    exit_code = run_tests()
    sys.exit(exit_code)
