"""
test_iiv.py
Unit Tests for Inter-Individual Variability (IIV) - Phase 2 v2.0

Tests realistic population variability with deterministic tiers.
Constitutional: All tests verify deterministic, reproducible behavior.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

import unittest
from PX_Engine.operations.TrialEngine import generate_virtual_population
from PX_Laboratory.Simulation_Engine import SimulationEngine


class TestVirtualPopulationIIV(unittest.TestCase):
    """Test virtual population generation with IIV"""
    
    def test_no_variability_baseline(self):
        """With variability=None, all patients should have default factors (1.0)"""
        pop = generate_virtual_population(n_patients=10, variability=None)
        
        for patient in pop:
            # No IIV factors should be present
            self.assertNotIn("clearance_factor", patient,
                             "clearance_factor should not be present when variability=None")
            self.assertNotIn("vd_factor", patient,
                             "vd_factor should not be present when variability=None")
            self.assertNotIn("ka_factor", patient,
                             "ka_factor should not be present when variability=None")
    
    def test_clearance_factors_applied(self):
        """clearance_factor should be applied with correct tier pattern"""
        variability = {"clearance_variation": 0.3, "n_tiers": 7}
        pop = generate_virtual_population(n_patients=14, variability=variability)
        
        # First 7 patients should have factors cycling through tiers
        # n_tiers=7: -3, -2, -1, 0, 1, 2, 3
        # clearance_variation=0.3: ±30%
        # Expected: 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3
        expected_factors = [0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3]
        
        for i in range(7):
            self.assertIn("clearance_factor", pop[i])
            actual = pop[i]["clearance_factor"]
            expected = expected_factors[i]
            self.assertAlmostEqual(actual, expected, places=6,
                                   msg=f"Patient {i}: expected {expected}, got {actual}")
        
        # Patients 7-13 should repeat the pattern
        for i in range(7, 14):
            self.assertAlmostEqual(pop[i]["clearance_factor"], expected_factors[i - 7], places=6)
    
    def test_vd_factors_applied(self):
        """vd_factor should be applied with correct tier pattern"""
        variability = {"vd_variation": 0.25, "n_tiers": 5}
        pop = generate_virtual_population(n_patients=10, variability=variability)
        
        # n_tiers=5: -2, -1, 0, 1, 2
        # vd_variation=0.25: ±25%
        # Expected: 0.75, 0.875, 1.0, 1.125, 1.25
        expected_factors = [0.75, 0.875, 1.0, 1.125, 1.25]
        
        for i in range(5):
            self.assertIn("vd_factor", pop[i])
            actual = pop[i]["vd_factor"]
            expected = expected_factors[i]
            self.assertAlmostEqual(actual, expected, places=6,
                                   msg=f"Patient {i}: expected {expected}, got {actual}")
    
    def test_ka_factors_applied(self):
        """ka_factor should be applied when ka_variation specified"""
        variability = {"ka_variation": 0.2, "n_tiers": 7}
        pop = generate_virtual_population(n_patients=7, variability=variability)
        
        # n_tiers=7: -3, -2, -1, 0, 1, 2, 3
        # ka_variation=0.2: ±20%
        # Expected: 0.8, 0.867, 0.933, 1.0, 1.067, 1.133, 1.2
        
        for patient in pop:
            self.assertIn("ka_factor", patient)
            factor = patient["ka_factor"]
            self.assertGreaterEqual(factor, 0.5, "ka_factor should be ≥ 0.5 (physiological clamp)")
            self.assertLessEqual(factor, 2.0, "ka_factor should be ≤ 2.0 (physiological clamp)")
    
    def test_deterministic_pattern(self):
        """Same inputs should produce identical patients every time"""
        variability = {
            "clearance_variation": 0.3,
            "vd_variation": 0.25,
            "ka_variation": 0.2,
            "n_tiers": 7
        }
        
        # Generate twice
        pop1 = generate_virtual_population(n_patients=21, variability=variability)
        pop2 = generate_virtual_population(n_patients=21, variability=variability)
        
        # Should be identical
        for i in range(21):
            self.assertEqual(pop1[i]["patient_id"], pop2[i]["patient_id"])
            self.assertAlmostEqual(pop1[i]["weight_kg"], pop2[i]["weight_kg"], places=6)
            self.assertAlmostEqual(pop1[i]["clearance_factor"], pop2[i]["clearance_factor"], places=6)
            self.assertAlmostEqual(pop1[i]["vd_factor"], pop2[i]["vd_factor"], places=6)
            self.assertAlmostEqual(pop1[i]["ka_factor"], pop2[i]["ka_factor"], places=6)
    
    def test_physiological_clamping(self):
        """Factors should be clamped to [0.5, 2.0] range"""
        # Try extreme variation
        variability = {
            "clearance_variation": 10.0,  # Unrealistically high
            "vd_variation": 10.0,
            "ka_variation": 10.0,
            "n_tiers": 7
        }
        
        pop = generate_virtual_population(n_patients=21, variability=variability)
        
        for patient in pop:
            for factor_key in ["clearance_factor", "vd_factor", "ka_factor"]:
                factor = patient[factor_key]
                self.assertGreaterEqual(factor, 0.5,
                                        f"{factor_key} should be clamped to ≥ 0.5 (got {factor})")
                self.assertLessEqual(factor, 2.0,
                                     f"{factor_key} should be clamped to ≤ 2.0 (got {factor})")
    
    def test_factor_equals_one_at_median(self):
        """Middle tier should have factor = 1.0 (no modification)"""
        variability = {
            "clearance_variation": 0.3,
            "vd_variation": 0.25,
            "n_tiers": 7
        }
        
        pop = generate_virtual_population(n_patients=7, variability=variability)
        
        # Middle patient (index 3 for n_tiers=7) should have factor = 1.0
        median_patient = pop[3]
        self.assertAlmostEqual(median_patient["clearance_factor"], 1.0, places=6,
                               msg="Median tier should have clearance_factor = 1.0")
        self.assertAlmostEqual(median_patient["vd_factor"], 1.0, places=6,
                               msg="Median tier should have vd_factor = 1.0")


class TestSimulationEngineIIV(unittest.TestCase):
    """Test SimulationEngine applies IIV factors correctly"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.sim_engine = SimulationEngine(time_step_h=0.5)
        self.admet = {
            "distribution": {"predicted_vd_L_per_kg": 1.0},
            "metabolism": {"predicted_clearance_L_per_h_per_kg": 0.1},
            "absorption": {"predicted_bioavailability": 1.0}
        }
    
    def test_factor_one_reproduces_baseline(self):
        """With factor=1.0, results should match no-factor simulation"""
        patient_no_factor = {"weight_kg": 70.0}
        patient_with_factor_one = {
            "weight_kg": 70.0,
            "clearance_factor": 1.0,
            "vd_factor": 1.0,
            "ka_factor": 1.0
        }
        
        result_no_factor = self.sim_engine.simulate_one_compartment(
            dose_mg=100.0,
            duration_h=24.0,
            dosing_interval_h=24.0,
            patient=patient_no_factor,
            admet=self.admet,
        )
        
        result_with_factor = self.sim_engine.simulate_one_compartment(
            dose_mg=100.0,
            duration_h=24.0,
            dosing_interval_h=24.0,
            patient=patient_with_factor_one,
            admet=self.admet,
        )
        
        # AUC and Cmax should be identical (within floating point precision)
        auc_no = result_no_factor["summary"]["auc_mg_h_per_L"]
        auc_with = result_with_factor["summary"]["auc_mg_h_per_L"]
        self.assertAlmostEqual(auc_no, auc_with, places=4,
                               msg=f"AUC should be identical: {auc_no} vs {auc_with}")
        
        cmax_no = result_no_factor["summary"]["cmax_mg_per_L"]
        cmax_with = result_with_factor["summary"]["cmax_mg_per_L"]
        self.assertAlmostEqual(cmax_no, cmax_with, places=4,
                               msg=f"Cmax should be identical: {cmax_no} vs {cmax_with}")
    
    def test_clearance_factor_affects_auc(self):
        """Higher clearance_factor should decrease AUC"""
        patient_low_cl = {"weight_kg": 70.0, "clearance_factor": 0.7}
        patient_high_cl = {"weight_kg": 70.0, "clearance_factor": 1.3}
        
        result_low = self.sim_engine.simulate_one_compartment(
            dose_mg=100.0,
            duration_h=24.0,
            dosing_interval_h=24.0,
            patient=patient_low_cl,
            admet=self.admet,
        )
        
        result_high = self.sim_engine.simulate_one_compartment(
            dose_mg=100.0,
            duration_h=24.0,
            dosing_interval_h=24.0,
            patient=patient_high_cl,
            admet=self.admet,
        )
        
        auc_low = result_low["summary"]["auc_mg_h_per_L"]
        auc_high = result_high["summary"]["auc_mg_h_per_L"]
        
        # Lower clearance → higher AUC
        self.assertGreater(auc_low, auc_high,
                           f"Lower clearance should give higher AUC: {auc_low} vs {auc_high}")
    
    def test_vd_factor_affects_cmax(self):
        """Higher vd_factor should decrease Cmax"""
        patient_low_vd = {"weight_kg": 70.0, "vd_factor": 0.7}
        patient_high_vd = {"weight_kg": 70.0, "vd_factor": 1.3}
        
        result_low = self.sim_engine.simulate_one_compartment(
            dose_mg=100.0,
            duration_h=24.0,
            dosing_interval_h=24.0,
            patient=patient_low_vd,
            admet=self.admet,
        )
        
        result_high = self.sim_engine.simulate_one_compartment(
            dose_mg=100.0,
            duration_h=24.0,
            dosing_interval_h=24.0,
            patient=patient_high_vd,
            admet=self.admet,
        )
        
        cmax_low = result_low["summary"]["cmax_mg_per_L"]
        cmax_high = result_high["summary"]["cmax_mg_per_L"]
        
        # Lower Vd → higher Cmax
        self.assertGreater(cmax_low, cmax_high,
                           f"Lower Vd should give higher Cmax: {cmax_low} vs {cmax_high}")
    
    def test_ka_factor_affects_tmax(self):
        """Higher ka_factor should decrease Tmax (faster absorption)"""
        patient_low_ka = {"weight_kg": 70.0, "ka_factor": 0.7}
        patient_high_ka = {"weight_kg": 70.0, "ka_factor": 1.3}
        
        result_low = self.sim_engine.simulate_one_compartment(
            dose_mg=100.0,
            duration_h=24.0,
            dosing_interval_h=24.0,
            patient=patient_low_ka,
            admet=self.admet,
        )
        
        result_high = self.sim_engine.simulate_one_compartment(
            dose_mg=100.0,
            duration_h=24.0,
            dosing_interval_h=24.0,
            patient=patient_high_ka,
            admet=self.admet,
        )
        
        tmax_low = result_low["summary"]["tmax_h"]
        tmax_high = result_high["summary"]["tmax_h"]
        
        # Higher ka → faster absorption → earlier Tmax
        self.assertGreater(tmax_low, tmax_high,
                           f"Higher ka should give earlier Tmax: {tmax_low} vs {tmax_high}")


def run_tests():
    """Run all IIV unit tests"""
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    # Add all test classes
    suite.addTests(loader.loadTestsFromTestCase(TestVirtualPopulationIIV))
    suite.addTests(loader.loadTestsFromTestCase(TestSimulationEngineIIV))
    
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    print("\n" + "="*70)
    print(f"IIV UNIT TESTS SUMMARY")
    print("="*70)
    print(f"Tests run:     {result.testsRun}")
    print(f"Successes:     {result.testsRun - len(result.failures) - len(result.errors)}")
    print(f"Failures:      {len(result.failures)}")
    print(f"Errors:        {len(result.errors)}")
    print("="*70)
    
    if result.wasSuccessful():
        print("✅ ALL IIV UNIT TESTS PASSED")
        return 0
    else:
        print("❌ SOME TESTS FAILED")
        return 1


if __name__ == "__main__":
    exit_code = run_tests()
    sys.exit(exit_code)
