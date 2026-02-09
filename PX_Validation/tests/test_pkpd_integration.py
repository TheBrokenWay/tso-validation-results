"""
test_pkpd_integration.py
Integration Tests for PK/PD Pipeline (Phase 1 - v2.0)

Tests end-to-end PK → PK/PD → Trial integration.
Constitutional: Verifies PD summaries in trial results.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

import unittest
from PX_Laboratory.Simulation_Engine import SimulationEngine
from PX_Engine.operations.PKPD import link_pk_to_pd
from PX_Engine.operations.TrialEngine import TrialEngine


class TestPKPDPipeline(unittest.TestCase):
    """Test PK → PK/PD pipeline integration"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.sim_engine = SimulationEngine(time_step_h=0.5)
        
        # Standard ADMET parameters
        self.admet = {
            "distribution": {"predicted_vd_L_per_kg": 0.7},
            "metabolism": {"predicted_clearance_L_per_h_per_kg": 0.05},
            "absorption": {"predicted_bioavailability": 1.0}
        }
        
        # Standard patient
        self.patient = {"weight_kg": 70.0}
        
        # Standard PD parameters
        self.pd_params = {
            "emax": 0.9,
            "ec50": 5.0,
            "hill": 1.5,
            "effect_threshold": 0.6
        }
    
    def test_pk_to_pkpd_pipeline_basic(self):
        """PK simulation → PK/PD linking should work end-to-end"""
        # 1. Run PK simulation
        pk_result = self.sim_engine.simulate_one_compartment(
            dose_mg=100.0,
            duration_h=24.0,
            dosing_interval_h=12.0,
            patient=self.patient,
            admet=self.admet,
        )
        
        # Verify PK result has required fields
        self.assertIn("time_grid_h", pk_result)
        self.assertIn("concentration_mg_per_L", pk_result)
        
        # 2. Run PK/PD linking
        pd_result = link_pk_to_pd(pk_result, self.pd_params)
        
        # Verify PD result structure
        self.assertIn("effect", pd_result)
        self.assertIn("pd_summary", pd_result)
        
        # Verify effect profile matches PK profile length
        self.assertEqual(len(pd_result["effect"]), len(pk_result["time_grid_h"]))
        
        # Verify PD summary has metrics
        pd_summary = pd_result["pd_summary"]
        self.assertIn("max_effect", pd_summary)
        self.assertIn("auec_h", pd_summary)
        self.assertIn("time_above_threshold_h", pd_summary)
        
        # Verify max effect is reasonable (should be < emax)
        self.assertLess(pd_summary["max_effect"], self.pd_params["emax"])
        self.assertGreater(pd_summary["max_effect"], 0.0)
    
    def test_pk_to_pkpd_dose_response(self):
        """Higher dose should produce higher effect"""
        doses = [50.0, 100.0, 200.0]
        max_effects = []
        
        for dose in doses:
            # Run PK
            pk_result = self.sim_engine.simulate_one_compartment(
                dose_mg=dose,
                duration_h=24.0,
                dosing_interval_h=24.0,
                patient=self.patient,
                admet=self.admet,
            )
            
            # Run PK/PD
            pd_result = link_pk_to_pd(pk_result, self.pd_params)
            max_effects.append(pd_result["pd_summary"]["max_effect"])
        
        # Verify dose-response relationship: higher dose → higher effect
        self.assertLess(max_effects[0], max_effects[1],
                        "100mg should produce higher effect than 50mg")
        self.assertLess(max_effects[1], max_effects[2],
                        "200mg should produce higher effect than 100mg")


class TestTrialEnginePD(unittest.TestCase):
    """Test TrialEngine with PD enabled"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.trial_engine = TrialEngine(time_step_h=1.0)
        
        # Standard ADMET
        self.admet = {
            "distribution": {"predicted_vd_L_per_kg": 0.7},
            "metabolism": {"predicted_clearance_L_per_h_per_kg": 0.05},
            "absorption": {"predicted_bioavailability": 1.0}
        }
        
        # Standard PD parameters
        self.pd_params = {
            "emax": 0.9,
            "ec50": 5.0,
            "hill": 1.5,
            "effect_threshold": 0.6
        }
        
        # Simple trial protocol
        self.protocol = {
            "trial_id": "TEST-PKPD-001",
            "duration_days": 7.0,
            "arms": [
                {
                    "arm_id": "A1",
                    "label": "100mg QD",
                    "dose_mg": 100.0,
                    "dosing_interval_h": 24.0,
                    "n_patients": 10,
                }
            ]
        }
    
    def test_trial_with_pd_enabled(self):
        """TrialEngine with pd_params should include PD summary"""
        # Run trial with PD
        result = self.trial_engine.run_trial(
            protocol=self.protocol,
            admet=self.admet,
            pd_params=self.pd_params,
        )
        
        # Verify trial result structure
        self.assertIn("trial_id", result)
        self.assertIn("arms", result)
        self.assertIn("pd_params", result)
        
        # Verify PD params stored for provenance
        self.assertEqual(result["pd_params"], self.pd_params)
        
        # Verify arm has PD summary
        arm = result["arms"][0]
        self.assertIn("pd_summary", arm)
        
        # Verify PD summary structure
        pd_summary = arm["pd_summary"]
        required_metrics = ["max_effect", "auec_h", "time_above_threshold_h", "mean_effect"]
        for metric in required_metrics:
            self.assertIn(metric, pd_summary)
            # Each metric should have distribution stats
            self.assertIn("mean", pd_summary[metric])
            self.assertIn("median", pd_summary[metric])
            self.assertIn("min", pd_summary[metric])
            self.assertIn("max", pd_summary[metric])
    
    def test_trial_without_pd_backward_compat(self):
        """TrialEngine without pd_params should still work (backward compatibility)"""
        # Run trial without PD
        result = self.trial_engine.run_trial(
            protocol=self.protocol,
            admet=self.admet,
            # No pd_params
        )
        
        # Verify trial result structure
        self.assertIn("trial_id", result)
        self.assertIn("arms", result)
        
        # Verify pd_params is None
        self.assertIsNone(result.get("pd_params"))
        
        # Verify arm has exposure_summary but NOT pd_summary
        arm = result["arms"][0]
        self.assertIn("exposure_summary", arm)
        self.assertNotIn("pd_summary", arm)
    
    def test_pd_metrics_change_with_dose(self):
        """PD metrics should change with dose in multi-arm trial"""
        # Multi-arm protocol
        protocol_multi_arm = {
            "trial_id": "TEST-PKPD-002",
            "duration_days": 7.0,
            "arms": [
                {
                    "arm_id": "A1",
                    "label": "50mg QD",
                    "dose_mg": 50.0,
                    "dosing_interval_h": 24.0,
                    "n_patients": 10,
                },
                {
                    "arm_id": "A2",
                    "label": "200mg QD",
                    "dose_mg": 200.0,
                    "dosing_interval_h": 24.0,
                    "n_patients": 10,
                },
            ]
        }
        
        # Run trial with PD
        result = self.trial_engine.run_trial(
            protocol=protocol_multi_arm,
            admet=self.admet,
            pd_params=self.pd_params,
        )
        
        # Extract PD metrics from both arms
        arm1_pd = result["arms"][0]["pd_summary"]
        arm2_pd = result["arms"][1]["pd_summary"]
        
        # Higher dose (arm2) should produce higher effect
        self.assertLess(arm1_pd["max_effect"]["mean"], arm2_pd["max_effect"]["mean"],
                        "Higher dose should produce higher mean max_effect")
        self.assertLess(arm1_pd["auec_h"]["mean"], arm2_pd["auec_h"]["mean"],
                        "Higher dose should produce higher mean AUEC")
    
    def test_trial_constitutional_compliance(self):
        """Trial result with PD should have constitutional metadata"""
        # Run trial with PD
        result = self.trial_engine.run_trial(
            protocol=self.protocol,
            admet=self.admet,
            pd_params=self.pd_params,
        )
        
        # Verify constitutional block
        self.assertIn("constitutional", result)
        const = result["constitutional"]
        
        self.assertEqual(const["status"], "SIMULATED")
        self.assertIn("engine", const)
        self.assertIn("notes", const)
        
        # Should mention PD modeling
        self.assertIn("PK/PD", const["notes"])
        self.assertIn("validation", const["notes"].lower())


def run_tests():
    """Run all PKPD integration tests"""
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    # Add all test classes
    suite.addTests(loader.loadTestsFromTestCase(TestPKPDPipeline))
    suite.addTests(loader.loadTestsFromTestCase(TestTrialEnginePD))
    
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    print("\n" + "="*70)
    print(f"PKPD INTEGRATION TESTS SUMMARY")
    print("="*70)
    print(f"Tests run:     {result.testsRun}")
    print(f"Successes:     {result.testsRun - len(result.failures) - len(result.errors)}")
    print(f"Failures:      {len(result.failures)}")
    print(f"Errors:        {len(result.errors)}")
    print("="*70)
    
    if result.wasSuccessful():
        print("✅ ALL PKPD INTEGRATION TESTS PASSED")
        return 0
    else:
        print("❌ SOME TESTS FAILED")
        return 1


if __name__ == "__main__":
    exit_code = run_tests()
    sys.exit(exit_code)
