"""
test_virtual_efficacy_integration.py
Integration Tests for Virtual Efficacy Analytics - Phase 5

Tests full efficacy analytics pipeline with TrialEngine integration.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

import unittest
from PX_Engine.operations.TrialEngine import TrialEngine
from PX_Engine.operations.VirtualEfficacyAnalytics import (
    compute_pta,
    virtual_responder_rate,
    analyze_virtual_efficacy,
)


class TestVirtualEfficacyIntegration(unittest.TestCase):
    """Test efficacy analytics with real trial data"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.trial_engine = TrialEngine()
        self.admet = {
            "distribution": {"predicted_vd_L_per_kg": 0.7},
            "metabolism": {"predicted_clearance_L_per_h_per_kg": 0.05},
            "absorption": {"predicted_bioavailability": 1.0}
        }
        self.pd_params = {
            "emax": 0.9,
            "ec50": 5.0,
            "hill": 1.5,
        }
    
    def test_pta_with_real_trial(self):
        """PTA computation with real TrialEngine output"""
        protocol = {
            "trial_id": "PTA-TEST",
            "duration_days": 7.0,
            "arms": [{
                "arm_id": "A1",
                "label": "150mg QD",
                "dose_mg": 150.0,
                "dosing_interval_h": 24.0,
                "n_patients": 21,
            }]
        }
        
        trial_result = self.trial_engine.run_trial(
            protocol=protocol,
            admet=self.admet,
            variability={
                "clearance_variation": 0.3,
                "vd_variation": 0.25,
            }
        )
        
        pta = compute_pta(
            trial_result=trial_result,
            metric="auc_mg_h_per_L",
            target_threshold=250.0,
        )
        
        self.assertIn("pta", pta)
        self.assertGreaterEqual(pta["pta"], 0.0)
        self.assertLessEqual(pta["pta"], 1.0)
        self.assertEqual(pta["n_total"], 21)
    
    def test_responder_rate_with_pkpd(self):
        """Responder rate with PK/PD trial"""
        protocol = {
            "trial_id": "RESP-TEST",
            "duration_days": 7.0,
            "arms": [{
                "arm_id": "A1",
                "label": "200mg QD",
                "dose_mg": 200.0,
                "dosing_interval_h": 24.0,
                "n_patients": 21,
            }]
        }
        
        trial_result = self.trial_engine.run_trial(
            protocol=protocol,
            admet=self.admet,
            pd_params=self.pd_params,
            variability={
                "clearance_variation": 0.3,
                "vd_variation": 0.25,
            }
        )
        
        responders = virtual_responder_rate(
            trial_result=trial_result,
            response_metric="max_effect",
            responder_threshold=0.65,
        )
        
        self.assertIn("responder_rate", responders)
        self.assertGreaterEqual(responders["responder_rate"], 0.0)
        self.assertLessEqual(responders["responder_rate"], 1.0)
    
    def test_analyze_virtual_efficacy_comprehensive(self):
        """Comprehensive efficacy analysis"""
        protocol = {
            "trial_id": "EFFICACY-TEST",
            "duration_days": 7.0,
            "arms": [{
                "arm_id": "A1",
                "label": "180mg QD",
                "dose_mg": 180.0,
                "dosing_interval_h": 24.0,
                "n_patients": 21,
            }]
        }
        
        trial_result = self.trial_engine.run_trial(
            protocol=protocol,
            admet=self.admet,
            pd_params=self.pd_params,
            variability={
                "clearance_variation": 0.3,
                "vd_variation": 0.25,
            }
        )
        
        analytics = analyze_virtual_efficacy(
            trial_result=trial_result,
            pk_target={"auc_mg_h_per_L": 300.0},
            pd_target={"max_effect": 0.70},
            safety_threshold=0.85,
        )
        
        # Should have all sections
        self.assertIn("pk_pta", analytics)
        self.assertIn("pd_responders", analytics)
        self.assertIn("effect_risk", analytics)
        self.assertIn("constitutional", analytics)
    
    def test_multi_arm_analytics(self):
        """Analytics should work with multiple arms"""
        protocol = {
            "trial_id": "MULTI-ARM",
            "duration_days": 7.0,
            "arms": [
                {
                    "arm_id": "A1",
                    "label": "100mg QD",
                    "dose_mg": 100.0,
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
        
        trial_result = self.trial_engine.run_trial(
            protocol=protocol,
            admet=self.admet,
            pd_params=self.pd_params,
        )
        
        # PTA for all arms
        pta = compute_pta(
            trial_result=trial_result,
            metric="auc_mg_h_per_L",
            target_threshold=250.0,
        )
        
        self.assertEqual(pta["n_total"], 28, "Should include both arms")


def run_tests():
    """Run all virtual efficacy integration tests"""
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    suite.addTests(loader.loadTestsFromTestCase(TestVirtualEfficacyIntegration))
    
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    print("\n" + "="*70)
    print(f"VIRTUAL EFFICACY INTEGRATION TESTS SUMMARY")
    print("="*70)
    print(f"Tests run:     {result.testsRun}")
    print(f"Successes:     {result.testsRun - len(result.failures) - len(result.errors)}")
    print(f"Failures:      {len(result.failures)}")
    print(f"Errors:        {len(result.errors)}")
    print("="*70)
    
    if result.wasSuccessful():
        print("✅ ALL VIRTUAL EFFICACY INTEGRATION TESTS PASSED")
        return 0
    else:
        print("❌ SOME TESTS FAILED")
        return 1


if __name__ == "__main__":
    exit_code = run_tests()
    sys.exit(exit_code)
