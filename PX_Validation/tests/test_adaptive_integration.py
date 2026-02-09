"""
test_adaptive_integration.py
Integration Tests for Adaptive Trial Logic - Phase 3 v2.0

Tests full adaptive trial simulations with dose adjustments and arm stopping.
Constitutional: Verifies logged, deterministic adaptive behavior.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

import unittest
from PX_Engine.operations.TrialEngine import TrialEngine


class TestAdaptiveTrialIntegration(unittest.TestCase):
    """Test complete adaptive trial simulations"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.trial_engine = TrialEngine(time_step_h=1.0)
        
        # Standard ADMET
        self.admet = {
            "distribution": {"predicted_vd_L_per_kg": 0.7},
            "metabolism": {"predicted_clearance_L_per_h_per_kg": 0.05},
            "absorption": {"predicted_bioavailability": 1.0}
        }
    
    def test_no_adaptive_rules_backward_compat(self):
        """Trial without adaptive_rules should work as before"""
        protocol = {
            "trial_id": "TEST-NO-ADAPT",
            "duration_days": 7.0,
            "arms": [{
                "arm_id": "A1",
                "label": "100mg QD",
                "dose_mg": 100.0,
                "dosing_interval_h": 24.0,
                "n_patients": 10,
            }]
        }
        
        result = self.trial_engine.run_trial(
            protocol=protocol,
            admet=self.admet,
        )
        
        # Should complete normally
        self.assertEqual(result["trial_id"], "TEST-NO-ADAPT")
        self.assertEqual(len(result["arms"]), 1)
        
        arm = result["arms"][0]
        self.assertEqual(arm["n_patients"], 10)
        self.assertEqual(arm["patients_enrolled"], 10)
        self.assertNotIn("adaptation_log", arm)
    
    def test_adaptive_reduce_dose_triggers(self):
        """High AUC should trigger REDUCE_DOSE"""
        protocol = {
            "trial_id": "TEST-REDUCE",
            "duration_days": 7.0,
            "arms": [{
                "arm_id": "A1",
                "label": "200mg QD (high dose)",
                "dose_mg": 200.0,  # High dose → high AUC
                "dosing_interval_h": 24.0,
                "n_patients": 20,
            }],
            "adaptive_rules": {
                "metric": "auc_mg_h_per_L",
                "upper_bound": 300.0,  # Trigger at 300
                "action": "REDUCE_DOSE",
                "dose_adjustment_factor": 0.75,
                "interim_after_n": 10,  # Check after 10 patients
            }
        }
        
        result = self.trial_engine.run_trial(
            protocol=protocol,
            admet=self.admet,
        )
        
        arm = result["arms"][0]
        
        # Should have adaptation log
        self.assertIn("adaptation_log", arm)
        self.assertGreater(len(arm["adaptation_log"]), 0, "Should have at least one interim")
        
        # Check if adaptation triggered
        triggered_count = arm["adaptations_triggered"]
        self.assertGreater(triggered_count, 0, "Adaptation should trigger for high AUC")
        
        # Dose should be reduced
        self.assertLess(arm["final_dose_mg"], arm["initial_dose_mg"], 
                        "Dose should be reduced after trigger")
        self.assertEqual(arm["final_dose_mg"], 150.0, "Dose should be 150mg (200*0.75)")
    
    def test_adaptive_within_bounds_no_trigger(self):
        """AUC within bounds should not trigger"""
        protocol = {
            "trial_id": "TEST-NO-TRIGGER",
            "duration_days": 7.0,
            "arms": [{
                "arm_id": "A1",
                "label": "100mg QD (normal dose)",
                "dose_mg": 100.0,  # Normal dose
                "dosing_interval_h": 24.0,
                "n_patients": 20,
            }],
            "adaptive_rules": {
                "metric": "auc_mg_h_per_L",
                "lower_bound": 50.0,
                "upper_bound": 500.0,  # Wide bounds
                "action": "REDUCE_DOSE",
                "interim_after_n": 10,
            }
        }
        
        result = self.trial_engine.run_trial(
            protocol=protocol,
            admet=self.admet,
        )
        
        arm = result["arms"][0]
        
        # Should have interim, but not triggered
        self.assertIn("adaptation_log", arm)
        self.assertEqual(arm["adaptations_triggered"], 0, "Should not trigger within bounds")
        
        # Dose unchanged
        self.assertEqual(arm["initial_dose_mg"], arm["final_dose_mg"], 
                         "Dose should remain unchanged")
    
    def test_adaptive_stop_arm(self):
        """STOP_ARM should halt patient enrollment"""
        protocol = {
            "trial_id": "TEST-STOP",
            "duration_days": 7.0,
            "arms": [{
                "arm_id": "A1",
                "label": "300mg QD (very high)",
                "dose_mg": 300.0,  # Very high dose
                "dosing_interval_h": 24.0,
                "n_patients": 30,  # Plan 30 patients
            }],
            "adaptive_rules": {
                "metric": "cmax_mg_per_L",
                "upper_bound": 5.0,  # Low threshold
                "action": "STOP_ARM",
                "interim_after_n": 10,
            }
        }
        
        result = self.trial_engine.run_trial(
            protocol=protocol,
            admet=self.admet,
        )
        
        arm = result["arms"][0]
        
        # Arm should stop
        self.assertTrue(arm["arm_stopped"], "Arm should be stopped")
        self.assertLess(arm["patients_enrolled"], arm["n_patients"], 
                        "Should enroll fewer patients than planned")
        
        # Should have triggered
        self.assertGreater(arm["adaptations_triggered"], 0)
    
    def test_adaptive_with_iiv(self):
        """Adaptive logic should work with IIV"""
        protocol = {
            "trial_id": "TEST-IIV-ADAPT",
            "duration_days": 7.0,
            "arms": [{
                "arm_id": "A1",
                "label": "200mg QD",
                "dose_mg": 200.0,
                "dosing_interval_h": 24.0,
                "n_patients": 21,  # 7-tier coverage
            }],
            "adaptive_rules": {
                "metric": "auc_mg_h_per_L",
                "upper_bound": 350.0,
                "action": "REDUCE_DOSE",
                "dose_adjustment_factor": 0.8,
                "interim_after_n": 14,  # After 2 tier cycles
            }
        }
        
        result = self.trial_engine.run_trial(
            protocol=protocol,
            admet=self.admet,
            variability={
                "clearance_variation": 0.3,
                "vd_variation": 0.25,
                "n_tiers": 7
            }
        )
        
        arm = result["arms"][0]
        
        # Should have adaptation log
        self.assertIn("adaptation_log", arm)
        
        # AUC should show variability
        auc_std = arm["exposure_summary"]["auc_mg_h_per_L"]["std"]
        self.assertGreater(auc_std, 20.0, "Should have IIV variability")
    
    def test_adaptive_with_pkpd(self):
        """Adaptive logic on PD metrics"""
        protocol = {
            "trial_id": "TEST-PD-ADAPT",
            "duration_days": 7.0,
            "arms": [{
                "arm_id": "A1",
                "label": "150mg QD",
                "dose_mg": 150.0,
                "dosing_interval_h": 24.0,
                "n_patients": 20,
            }],
            "adaptive_rules": {
                "metric": "max_effect",  # PD metric
                "upper_bound": 0.75,
                "action": "REDUCE_DOSE",
                "dose_adjustment_factor": 0.85,
                "interim_after_n": 10,
            }
        }
        
        pd_params = {
            "emax": 0.9,
            "ec50": 3.0,
            "hill": 1.5,
            "effect_threshold": 0.6
        }
        
        result = self.trial_engine.run_trial(
            protocol=protocol,
            admet=self.admet,
            pd_params=pd_params,
        )
        
        arm = result["arms"][0]
        
        # Should have PD summary
        self.assertIn("pd_summary", arm)
        
        # Should have adaptation log
        self.assertIn("adaptation_log", arm)
        
        # Check if PD metrics drove adaptation
        if arm["adaptations_triggered"] > 0:
            # If triggered, final dose should be lower
            self.assertLess(arm["final_dose_mg"], arm["initial_dose_mg"])
    
    def test_adaptation_logging_complete(self):
        """Adaptation log should be comprehensive"""
        protocol = {
            "trial_id": "TEST-LOG",
            "duration_days": 7.0,
            "arms": [{
                "arm_id": "A1",
                "label": "180mg QD",
                "dose_mg": 180.0,
                "dosing_interval_h": 24.0,
                "n_patients": 30,
            }],
            "adaptive_rules": {
                "metric": "auc_mg_h_per_L",
                "upper_bound": 320.0,
                "action": "REDUCE_DOSE",
                "dose_adjustment_factor": 0.8,
                "interim_after_n": 10,
            }
        }
        
        result = self.trial_engine.run_trial(
            protocol=protocol,
            admet=self.admet,
        )
        
        arm = result["arms"][0]
        adaptation_log = arm["adaptation_log"]
        
        # Should have entries (30 patients / 10 interim = 2 interims minimum)
        self.assertGreaterEqual(len(adaptation_log), 2, "Should have at least 2 interim analyses")
        
        # Each entry should have required fields
        for entry in adaptation_log:
            self.assertIn("patient_count", entry)
            self.assertIn("metric", entry)
            self.assertIn("mean_value", entry)
            self.assertIn("triggered", entry)
            self.assertIn("reason", entry)
            self.assertIn("action", entry)
            self.assertIn("current_dose_mg", entry)
            self.assertIn("new_dose_mg", entry)


def run_tests():
    """Run all adaptive integration tests"""
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    # Add all test classes
    suite.addTests(loader.loadTestsFromTestCase(TestAdaptiveTrialIntegration))
    
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    print("\n" + "="*70)
    print(f"ADAPTIVE INTEGRATION TESTS SUMMARY")
    print("="*70)
    print(f"Tests run:     {result.testsRun}")
    print(f"Successes:     {result.testsRun - len(result.failures) - len(result.errors)}")
    print(f"Failures:      {len(result.failures)}")
    print(f"Errors:        {len(result.errors)}")
    print("="*70)
    
    if result.wasSuccessful():
        print("✅ ALL ADAPTIVE INTEGRATION TESTS PASSED")
        return 0
    else:
        print("❌ SOME TESTS FAILED")
        return 1


if __name__ == "__main__":
    exit_code = run_tests()
    sys.exit(exit_code)
