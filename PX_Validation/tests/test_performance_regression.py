"""
test_performance_regression.py
Performance Regression Tests for v2.0.0-CORE

Ensures Phase 7 optimizations don't degrade performance.
Tests against baseline thresholds from PERF_BASELINE_v2.0.0-CORE.md
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

import unittest
import time
from PX_Laboratory import SimulationEngine
from PX_Engine.operations import TrialEngine
from PX_Engine.operations.VirtualEfficacyAnalytics import compute_pta
from PX_System.foundation.Evidence_Package import wrap_trial_simulation


class TestPerformanceRegression(unittest.TestCase):
    """Performance regression tests against baseline"""
    
    def setUp(self):
        """Set up common test fixtures"""
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
    
    def test_pk_engine_performance(self):
        """PK engine should complete <2ms per patient (baseline: 1ms)"""
        engine = SimulationEngine(time_step_h=1.0)
        patient = {"weight_kg": 70.0}
        
        # Warm-up
        for _ in range(10):
            engine.simulate_one_compartment(
                dose_mg=100.0,
                duration_h=24.0,
                dosing_interval_h=24.0,
                patient=patient,
                admet=self.admet,
            )
        
        # Measurement
        iterations = 100
        start = time.perf_counter()
        for _ in range(iterations):
            result = engine.simulate_one_compartment(
                dose_mg=100.0,
                duration_h=24.0,
                dosing_interval_h=24.0,
                patient=patient,
                admet=self.admet,
            )
        end = time.perf_counter()
        
        avg_time_ms = ((end - start) / iterations) * 1000
        
        # Should be <2ms (baseline was ~1ms, allow 2× margin)
        self.assertLess(avg_time_ms, 2.0,
                        f"PK engine too slow: {avg_time_ms:.2f}ms (threshold: 2.0ms)")
        
        print(f"✅ PK Engine: {avg_time_ms:.2f}ms per patient (baseline: 1ms, threshold: 2ms)")
    
    def test_trial_engine_performance(self):
        """TrialEngine should complete <200ms for 40 patients (baseline: 100ms)"""
        engine = TrialEngine(time_step_h=2.0)  # Coarser for speed
        
        protocol = {
            "trial_id": "PERF-TEST",
            "duration_days": 3.0,  # Shorter for perf test
            "arms": [
                {
                    "arm_id": "A1",
                    "label": "Arm 1",
                    "dose_mg": 100.0,
                    "dosing_interval_h": 24.0,
                    "n_patients": 20,
                },
                {
                    "arm_id": "A2",
                    "label": "Arm 2",
                    "dose_mg": 150.0,
                    "dosing_interval_h": 24.0,
                    "n_patients": 20,
                },
            ]
        }
        
        # Warm-up
        for _ in range(3):
            engine.run_trial(protocol, self.admet, self.pd_params)
        
        # Measurement
        start = time.perf_counter()
        result = engine.run_trial(protocol, self.admet, self.pd_params)
        end = time.perf_counter()
        
        duration_ms = (end - start) * 1000
        
        # Should be <200ms (baseline was ~100ms, allow 2× margin)
        self.assertLess(duration_ms, 200.0,
                        f"TrialEngine too slow: {duration_ms:.1f}ms (threshold: 200ms)")
        
        print(f"✅ TrialEngine: {duration_ms:.1f}ms for 40 patients (baseline: 100ms, threshold: 200ms)")
    
    def test_virtual_efficacy_performance(self):
        """VirtualEfficacy should complete <20ms (baseline: 10ms)"""
        engine = TrialEngine(time_step_h=2.0)
        
        protocol = {
            "trial_id": "PERF-EFFICACY",
            "duration_days": 3.0,
            "arms": [
                {
                    "arm_id": "A1",
                    "label": "Test",
                    "dose_mg": 100.0,
                    "dosing_interval_h": 24.0,
                    "n_patients": 20,
                }
            ]
        }
        
        trial_result = engine.run_trial(protocol, self.admet, self.pd_params)
        
        # Warm-up
        for _ in range(10):
            compute_pta(trial_result, "auc_mg_h_per_L", 250.0)
        
        # Measurement
        start = time.perf_counter()
        pta_result = compute_pta(trial_result, "auc_mg_h_per_L", 250.0)
        end = time.perf_counter()
        
        duration_ms = (end - start) * 1000
        
        # Should be <20ms (baseline was ~10ms, allow 2× margin)
        self.assertLess(duration_ms, 20.0,
                        f"VirtualEfficacy too slow: {duration_ms:.1f}ms (threshold: 20ms)")
        
        print(f"✅ VirtualEfficacy: {duration_ms:.2f}ms (baseline: 10ms, threshold: 20ms)")
    
    def test_evidence_package_performance(self):
        """Evidence_Package v3 should complete <40ms (baseline: 20ms)"""
        engine = TrialEngine(time_step_h=2.0)
        
        protocol = {
            "trial_id": "PERF-DOSSIER",
            "duration_days": 3.0,
            "arms": [
                {
                    "arm_id": "A1",
                    "label": "Test",
                    "dose_mg": 100.0,
                    "dosing_interval_h": 24.0,
                    "n_patients": 14,
                }
            ]
        }
        
        trial_result = engine.run_trial(protocol, self.admet, self.pd_params)
        ope = {"status": "COMPUTED"}
        
        # Warm-up
        for i in range(5):
            wrap_trial_simulation(
                protocol=protocol,
                trial_result=trial_result,
                ope=ope,
                admet=self.admet,
                output_dir=f"PX_Warehouse/TrialSimulations/TestRuns/perftest_{i}"
            )
        
        # Measurement
        start = time.perf_counter()
        path = wrap_trial_simulation(
            protocol=protocol,
            trial_result=trial_result,
            ope=ope,
            admet=self.admet,
            output_dir="PX_Warehouse/TrialSimulations/TestRuns/perftest_final"
        )
        end = time.perf_counter()
        
        duration_ms = (end - start) * 1000
        
        # Should be <40ms (baseline was ~20ms, allow 2× margin)
        self.assertLess(duration_ms, 40.0,
                        f"Evidence_Package too slow: {duration_ms:.1f}ms (threshold: 40ms)")
        
        print(f"✅ Evidence_Package: {duration_ms:.1f}ms (baseline: 20ms, threshold: 40ms)")


def run_tests():
    """Run performance regression tests"""
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    suite.addTests(loader.loadTestsFromTestCase(TestPerformanceRegression))
    
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    print("\n" + "="*70)
    print(f"PERFORMANCE REGRESSION TESTS SUMMARY")
    print("="*70)
    print(f"Tests run:     {result.testsRun}")
    print(f"Successes:     {result.testsRun - len(result.failures) - len(result.errors)}")
    print(f"Failures:      {len(result.failures)}")
    print(f"Errors:        {len(result.errors)}")
    print("="*70)
    
    if result.wasSuccessful():
        print("✅ ALL PERFORMANCE REGRESSION TESTS PASSED")
        print("⚡ No performance degradation detected")
        return 0
    else:
        print("❌ PERFORMANCE REGRESSION DETECTED")
        return 1


if __name__ == "__main__":
    exit_code = run_tests()
    sys.exit(exit_code)
