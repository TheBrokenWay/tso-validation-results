import unittest

from PX_Engine.operations import TrialEngine
from PX_Engine.operations import generate_virtual_population


class TestTrialEngine(unittest.TestCase):

    def test_population_generator(self):
        pop = generate_virtual_population(n_patients=5)
        self.assertEqual(len(pop), 5)
        self.assertIn("patient_id", pop[0])
        self.assertIn("weight_kg", pop[0])

    def test_population_deterministic(self):
        """Test that population generation is deterministic"""
        pop1 = generate_virtual_population(n_patients=10)
        pop2 = generate_virtual_population(n_patients=10)
        
        for i in range(10):
            self.assertEqual(pop1[i]["weight_kg"], pop2[i]["weight_kg"])
            self.assertEqual(pop1[i]["patient_id"], pop2[i]["patient_id"])

    def test_population_weight_bounds(self):
        """Test that weights are within physiological bounds"""
        pop = generate_virtual_population(n_patients=20, base_weight_kg=70.0, weight_sd_kg=10.0)
        
        for patient in pop:
            weight = patient["weight_kg"]
            self.assertGreaterEqual(weight, 40.0)
            self.assertLessEqual(weight, 120.0)

    def test_simple_trial(self):
        engine = TrialEngine(time_step_h=2.0)

        protocol = {
            "trial_id": "TRIAL-PK-001",
            "duration_days": 2.0,
            "arms": [
                {
                    "arm_id": "A1",
                    "label": "QD 100 mg",
                    "dose_mg": 100.0,
                    "dosing_interval_h": 24.0,
                    "n_patients": 5,
                },
                {
                    "arm_id": "A2",
                    "label": "BID 50 mg",
                    "dose_mg": 50.0,
                    "dosing_interval_h": 12.0,
                    "n_patients": 5,
                },
            ],
        }

        admet = {
            "absorption": {"predicted_bioavailability": 1.0},
            "distribution": {"predicted_vd_L_per_kg": 0.7},
            "metabolism": {"predicted_clearance_L_per_h_per_kg": 0.05},
        }

        result = engine.run_trial(protocol, admet)

        self.assertEqual(result["trial_id"], "TRIAL-PK-001")
        self.assertEqual(len(result["arms"]), 2)

        for arm in result["arms"]:
            exp = arm["exposure_summary"]
            for key in ["cmax_mg_per_L", "auc_mg_h_per_L", "cmin_steady_state_mg_per_L"]:
                self.assertIn("mean", exp[key])
                self.assertIn("min", exp[key])
                self.assertIn("max", exp[key])
                self.assertGreaterEqual(exp[key]["max"], exp[key]["min"])

    def test_trial_constitutional_compliance(self):
        """Test that trial results include constitutional tracking"""
        engine = TrialEngine(time_step_h=1.0)
        
        protocol = {
            "trial_id": "TRIAL-CONST-001",
            "duration_days": 1.0,
            "arms": [
                {
                    "arm_id": "A1",
                    "label": "Test Arm",
                    "dose_mg": 100.0,
                    "dosing_interval_h": 24.0,
                    "n_patients": 3,
                },
            ],
        }
        
        admet = {}
        
        result = engine.run_trial(protocol, admet)
        
        self.assertIn("constitutional", result)
        self.assertEqual(result["constitutional"]["status"], "SIMULATED")
        self.assertEqual(result["constitutional"]["engine"], "TRIAL_ENGINE_V1")

    def test_invalid_protocol(self):
        """Test that invalid protocols raise errors"""
        engine = TrialEngine()
        
        # No arms
        with self.assertRaises(ValueError):
            engine.run_trial({"trial_id": "TEST", "arms": []}, {})
        
        # Invalid dose
        with self.assertRaises(ValueError):
            engine.run_trial({
                "trial_id": "TEST",
                "arms": [{"arm_id": "A1", "dose_mg": 0, "n_patients": 5}]
            }, {})

    def test_trial_with_empty_admet(self):
        """Test that trial works with empty ADMET (uses defaults)"""
        engine = TrialEngine(time_step_h=2.0)
        
        protocol = {
            "trial_id": "TRIAL-DEFAULT-001",
            "duration_days": 1.0,
            "arms": [
                {
                    "arm_id": "A1",
                    "label": "Default ADMET",
                    "dose_mg": 100.0,
                    "dosing_interval_h": 24.0,
                    "n_patients": 3,
                },
            ],
        }
        
        result = engine.run_trial(protocol, {})
        
        self.assertEqual(result["trial_id"], "TRIAL-DEFAULT-001")
        self.assertEqual(len(result["arms"]), 1)
        self.assertGreater(result["arms"][0]["exposure_summary"]["cmax_mg_per_L"]["mean"], 0)

    def test_exposure_statistics(self):
        """Test that exposure statistics are reasonable"""
        engine = TrialEngine(time_step_h=1.0)
        
        protocol = {
            "trial_id": "TRIAL-STATS-001",
            "duration_days": 2.0,
            "arms": [
                {
                    "arm_id": "A1",
                    "label": "Stats Test",
                    "dose_mg": 100.0,
                    "dosing_interval_h": 24.0,
                    "n_patients": 10,
                },
            ],
        }
        
        admet = {
            "absorption": {"predicted_bioavailability": 1.0},
            "distribution": {"predicted_vd_L_per_kg": 0.7},
            "metabolism": {"predicted_clearance_L_per_h_per_kg": 0.05},
        }
        
        result = engine.run_trial(protocol, admet)
        arm = result["arms"][0]
        
        # Check that mean is between min and max
        for metric in ["cmax_mg_per_L", "auc_mg_h_per_L", "cmin_steady_state_mg_per_L"]:
            stats = arm["exposure_summary"][metric]
            self.assertLessEqual(stats["min"], stats["mean"])
            self.assertLessEqual(stats["mean"], stats["max"])
            self.assertLessEqual(stats["min"], stats["median"])
            self.assertLessEqual(stats["median"], stats["max"])


if __name__ == "__main__":
    unittest.main()
