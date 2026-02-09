import unittest
from pathlib import Path

from PX_Engine.operations import TrialEngine, run_ope, run_admet
from PX_System.foundation.Evidence_Package import wrap_trial_simulation


class TestSystemTrialEvidencePackage(unittest.TestCase):

    def test_full_trial_dossier_pipeline(self):
        """Test complete trial dossier generation pipeline"""
        smiles = "CC(=O)Oc1ccccc1C(=O)O"

        ope = run_ope(smiles)
        admet = run_admet(smiles, ope)

        protocol = {
            "trial_id": "TRIAL-SYSTEM-001",
            "duration_days": 3.0,
            "arms": [
                {
                    "arm_id": "A1",
                    "label": "QD 100 mg",
                    "dose_mg": 100.0,
                    "dosing_interval_h": 24.0,
                    "n_patients": 10,
                }
            ],
        }

        engine = TrialEngine(time_step_h=1.0)
        trial_result = engine.run_trial(protocol, admet)

        path = wrap_trial_simulation(protocol, trial_result, ope, admet)
        self.assertTrue(Path(path).exists())

        # Basic structural checks
        self.assertIn("TRIAL_SIMULATION_DOSSIER", Path(path).name)


if __name__ == "__main__":
    unittest.main()
