"""
test_iiv_integration.py
Integration Tests for IIV (Inter-Individual Variability) - Phase 2 v2.0

Tests population variability in full trial pipeline.
Constitutional: Verifies realistic distributions in trial results.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

import unittest
from PX_Engine.operations.TrialEngine import TrialEngine


class TestTrialEngineIIV(unittest.TestCase):
    """Test TrialEngine with IIV enabled"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.trial_engine = TrialEngine(time_step_h=1.0)
        
        # Standard ADMET
        self.admet = {
            "distribution": {"predicted_vd_L_per_kg": 0.7},
            "metabolism": {"predicted_clearance_L_per_h_per_kg": 0.05},
            "absorption": {"predicted_bioavailability": 1.0}
        }
        
        # Simple trial protocol
        self.protocol = {
            "trial_id": "TEST-IIV-001",
            "duration_days": 7.0,
            "arms": [
                {
                    "arm_id": "A1",
                    "label": "100mg QD",
                    "dose_mg": 100.0,
                    "dosing_interval_h": 24.0,
                    "n_patients": 21,  # Multiple of 7 for good tier coverage
                }
            ]
        }
    
    def test_no_variability_reproduces_phase1(self):
        """With variability=None, results should match v2.0-PHASE1 behavior"""
        # Create protocol without variability parameter
        result = self.trial_engine.run_trial(
            protocol=self.protocol,
            admet=self.admet,
            # No variability parameter
        )
        
        # Verify trial result structure
        self.assertIn("trial_id", result)
        self.assertIn("arms", result)
        
        # Verify arm has exposure_summary
        arm = result["arms"][0]
        self.assertIn("exposure_summary", arm)
        
        # With no variability, std should be very small (weight variation only)
        auc_std = arm["exposure_summary"]["auc_mg_h_per_L"]["std"]
        self.assertGreater(auc_std, 0.0, "std should be > 0 due to weight variation")
        self.assertLess(auc_std, 50.0, "std should be small without PK variability")
    
    def test_with_variability_increases_std(self):
        """With IIV enabled, std should increase significantly"""
        # Trial without IIV
        result_no_iiv = self.trial_engine.run_trial(
            protocol=self.protocol,
            admet=self.admet,
        )
        
        # Trial with IIV
        result_with_iiv = self.trial_engine.run_trial(
            protocol=self.protocol,
            admet=self.admet,
            variability={
                "clearance_variation": 0.3,
                "vd_variation": 0.25,
                "n_tiers": 7
            }
        )
        
        # Extract std values
        std_no_iiv = result_no_iiv["arms"][0]["exposure_summary"]["auc_mg_h_per_L"]["std"]
        std_with_iiv = result_with_iiv["arms"][0]["exposure_summary"]["auc_mg_h_per_L"]["std"]
        
        # With IIV, std should be significantly higher
        self.assertGreater(std_with_iiv, std_no_iiv,
                           f"IIV should increase std: {std_with_iiv} vs {std_no_iiv}")
        
        # With 30% clearance variation, expect substantial std
        self.assertGreater(std_with_iiv, 40.0,
                           f"With IIV, AUC std should be substantial (got {std_with_iiv})")
    
    def test_exposure_distribution_realistic(self):
        """With IIV, exposure distributions should show realistic spread"""
        result = self.trial_engine.run_trial(
            protocol=self.protocol,
            admet=self.admet,
            variability={
                "clearance_variation": 0.3,
                "vd_variation": 0.25,
                "n_tiers": 7
            }
        )
        
        arm = result["arms"][0]
        auc_summary = arm["exposure_summary"]["auc_mg_h_per_L"]
        cmax_summary = arm["exposure_summary"]["cmax_mg_per_L"]
        
        # Verify distribution spread
        auc_range = auc_summary["max"] - auc_summary["min"]
        auc_mean = auc_summary["mean"]
        
        # Range should be substantial (>20% of mean with 30% clearance variation)
        self.assertGreater(auc_range / auc_mean, 0.2,
                           f"AUC range should be >20% of mean (got {auc_range/auc_mean*100:.1f}%)")
        
        # Cmax should also vary
        cmax_range = cmax_summary["max"] - cmax_summary["min"]
        cmax_mean = cmax_summary["mean"]
        
        self.assertGreater(cmax_range / cmax_mean, 0.15,
                           f"Cmax range should be >15% of mean (got {cmax_range/cmax_mean*100:.1f}%)")
    
    def test_pd_variability_propagates(self):
        """IIV should propagate to PD metrics"""
        pd_params = {
            "emax": 0.9,
            "ec50": 5.0,
            "hill": 1.5,
            "effect_threshold": 0.6
        }
        
        # Trial without IIV
        result_no_iiv = self.trial_engine.run_trial(
            protocol=self.protocol,
            admet=self.admet,
            pd_params=pd_params,
        )
        
        # Trial with IIV
        result_with_iiv = self.trial_engine.run_trial(
            protocol=self.protocol,
            admet=self.admet,
            pd_params=pd_params,
            variability={
                "clearance_variation": 0.3,
                "vd_variation": 0.25,
                "n_tiers": 7
            }
        )
        
        # Extract PD std
        pd_no_iiv = result_no_iiv["arms"][0]["pd_summary"]["max_effect"]["std"]
        pd_with_iiv = result_with_iiv["arms"][0]["pd_summary"]["max_effect"]["std"]
        
        # PD std should increase with IIV
        self.assertGreater(pd_with_iiv, pd_no_iiv,
                           f"IIV should increase PD std: {pd_with_iiv} vs {pd_no_iiv}")
    
    def test_std_present_in_all_summaries(self):
        """All distribution summaries should include std"""
        result = self.trial_engine.run_trial(
            protocol=self.protocol,
            admet=self.admet,
            variability={
                "clearance_variation": 0.3,
                "vd_variation": 0.25,
            }
        )
        
        arm = result["arms"][0]
        exposure_summary = arm["exposure_summary"]
        
        # Check all exposure metrics have std
        for metric in ["cmax_mg_per_L", "auc_mg_h_per_L", "cmin_steady_state_mg_per_L"]:
            self.assertIn("std", exposure_summary[metric],
                          f"{metric} should have 'std' field")
            std_val = exposure_summary[metric]["std"]
            self.assertIsInstance(std_val, float,
                                  f"{metric} std should be float")
            self.assertGreaterEqual(std_val, 0.0,
                                    f"{metric} std should be ≥ 0")
    
    def test_multi_arm_iiv(self):
        """IIV should work correctly with multiple arms"""
        protocol_multi = {
            "trial_id": "TEST-IIV-002",
            "duration_days": 7.0,
            "arms": [
                {
                    "arm_id": "A1",
                    "label": "50mg QD",
                    "dose_mg": 50.0,
                    "dosing_interval_h": 24.0,
                    "n_patients": 14,
                },
                {
                    "arm_id": "A2",
                    "label": "200mg QD",
                    "dose_mg": 200.0,
                    "dosing_interval_h": 24.0,
                    "n_patients": 14,
                },
            ]
        }
        
        result = self.trial_engine.run_trial(
            protocol=protocol_multi,
            admet=self.admet,
            variability={
                "clearance_variation": 0.3,
                "vd_variation": 0.25,
                "n_tiers": 7
            }
        )
        
        # Both arms should have variability
        for arm in result["arms"]:
            auc_std = arm["exposure_summary"]["auc_mg_h_per_L"]["std"]
            self.assertGreater(auc_std, 20.0,
                               f"Arm {arm['arm_id']} should have substantial AUC std")
        
        # Higher dose arm should have higher mean but similar relative variability
        arm1_auc = result["arms"][0]["exposure_summary"]["auc_mg_h_per_L"]
        arm2_auc = result["arms"][1]["exposure_summary"]["auc_mg_h_per_L"]
        
        # Higher dose → higher AUC
        self.assertGreater(arm2_auc["mean"], arm1_auc["mean"],
                           "Higher dose should give higher mean AUC")
        
        # Relative std (CV%) should be similar
        cv1 = arm1_auc["std"] / arm1_auc["mean"]
        cv2 = arm2_auc["std"] / arm2_auc["mean"]
        
        # CV% should be within 50% of each other (not exact due to nonlinearity)
        self.assertLess(abs(cv1 - cv2) / cv1, 0.5,
                        f"Relative variability (CV%) should be similar: {cv1:.3f} vs {cv2:.3f}")


def run_tests():
    """Run all IIV integration tests"""
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    # Add all test classes
    suite.addTests(loader.loadTestsFromTestCase(TestTrialEngineIIV))
    
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    print("\n" + "="*70)
    print(f"IIV INTEGRATION TESTS SUMMARY")
    print("="*70)
    print(f"Tests run:     {result.testsRun}")
    print(f"Successes:     {result.testsRun - len(result.failures) - len(result.errors)}")
    print(f"Failures:      {len(result.failures)}")
    print(f"Errors:        {len(result.errors)}")
    print("="*70)
    
    if result.wasSuccessful():
        print("✅ ALL IIV INTEGRATION TESTS PASSED")
        return 0
    else:
        print("❌ SOME TESTS FAILED")
        return 1


if __name__ == "__main__":
    exit_code = run_tests()
    sys.exit(exit_code)
