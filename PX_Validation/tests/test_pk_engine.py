import unittest
from PX_Laboratory import SimulationEngine


class TestPKEngine(unittest.TestCase):

    def test_basic_pk_profile(self):
        engine = SimulationEngine(time_step_h=1.0)

        patient = {"weight_kg": 70.0}
        admet = {
            "absorption": {"predicted_bioavailability": 1.0},
            "distribution": {"predicted_vd_L_per_kg": 0.7},
            "metabolism": {"predicted_clearance_L_per_h_per_kg": 0.05},
        }

        result = engine.simulate_one_compartment(
            dose_mg=100.0,
            duration_h=24.0,
            dosing_interval_h=24.0,
            patient=patient,
            admet=admet,
        )

        self.assertEqual(result["model"], "ONE_COMPARTMENT_FIRST_ORDER")
        self.assertIn("time_grid_h", result)
        self.assertIn("concentration_mg_per_L", result)
        self.assertGreater(len(result["time_grid_h"]), 1)
        self.assertEqual(len(result["time_grid_h"]), len(result["concentration_mg_per_L"]))

        cmax = result["summary"]["cmax_mg_per_L"]
        auc = result["summary"]["auc_mg_h_per_L"]
        self.assertGreater(cmax, 0.0)
        self.assertGreater(auc, 0.0)

    def test_invalid_inputs(self):
        engine = SimulationEngine(time_step_h=0.5)
        patient = {"weight_kg": 70.0}
        admet = {}

        with self.assertRaises(ValueError):
            engine.simulate_one_compartment(
                dose_mg=0.0,
                duration_h=24.0,
                dosing_interval_h=24.0,
                patient=patient,
                admet=admet,
            )

        with self.assertRaises(ValueError):
            engine.simulate_one_compartment(
                dose_mg=100.0,
                duration_h=0.0,
                dosing_interval_h=24.0,
                patient=patient,
                admet=admet,
            )

        with self.assertRaises(ValueError):
            engine.simulate_one_compartment(
                dose_mg=100.0,
                duration_h=24.0,
                dosing_interval_h=0.0,
                patient=patient,
                admet=admet,
            )

    def test_multiple_doses(self):
        """Test multiple dosing with BID regimen"""
        engine = SimulationEngine(time_step_h=1.0)
        
        patient = {"weight_kg": 70.0}
        admet = {
            "absorption": {"predicted_bioavailability": 1.0},
            "distribution": {"predicted_vd_L_per_kg": 0.7},
            "metabolism": {"predicted_clearance_L_per_h_per_kg": 0.05},
        }
        
        result = engine.simulate_one_compartment(
            dose_mg=50.0,
            duration_h=48.0,
            dosing_interval_h=12.0,  # BID (twice daily)
            patient=patient,
            admet=admet,
        )
        
        # Should have multiple peaks from repeated dosing
        self.assertEqual(result["model"], "ONE_COMPARTMENT_FIRST_ORDER")
        self.assertGreater(result["summary"]["cmax_mg_per_L"], 0.0)
        self.assertGreater(result["summary"]["auc_mg_h_per_L"], 0.0)

    def test_default_parameters(self):
        """Test that engine uses safe defaults when ADMET data is missing"""
        engine = SimulationEngine(time_step_h=1.0)
        
        patient = {"weight_kg": 70.0}
        admet = {}  # Empty ADMET data
        
        result = engine.simulate_one_compartment(
            dose_mg=100.0,
            duration_h=24.0,
            dosing_interval_h=24.0,
            patient=patient,
            admet=admet,
        )
        
        # Should still work with defaults
        self.assertEqual(result["parameters"]["vd_L_per_kg"], 0.7)
        self.assertEqual(result["parameters"]["clearance_L_per_h_per_kg"], 0.05)
        self.assertEqual(result["parameters"]["bioavailability"], 1.0)
        self.assertGreater(result["summary"]["cmax_mg_per_L"], 0.0)

    def test_legacy_materialize_candidate(self):
        """Test materialize_candidate when WorldLine file exists; skip when missing."""
        import os
        import json
        from PX_Warehouse.WorldLine_Database import DEFAULT_WORLDLINES_PATH
        engine = SimulationEngine()
        task_id = "TEST-123456"
        wl_path = os.path.join(DEFAULT_WORLDLINES_PATH, f"{task_id}.worldline")
        os.makedirs(DEFAULT_WORLDLINES_PATH, exist_ok=True)
        try:
            with open(wl_path, "w", encoding="utf-8") as f:
                json.dump({
                    "physical_realization": {
                        "toxicity_index": 0.018,
                        "binding_affinity_kj": -92.5,
                    }
                }, f, indent=2)
            result = engine.materialize_candidate("WL-TEST-123456", 0.85)
            self.assertEqual(result["status"], "READY_FOR_SYNTHESIS")
            self.assertIn("candidate_id", result)
            self.assertIn("binding_affinity_kj", result)
            self.assertIn("toxicity_index", result)
        finally:
            if os.path.exists(wl_path):
                try:
                    os.remove(wl_path)
                except OSError:
                    pass


if __name__ == "__main__":
    unittest.main()
