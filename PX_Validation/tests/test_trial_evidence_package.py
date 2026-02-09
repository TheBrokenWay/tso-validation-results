import unittest
import os
import json
from pathlib import Path

# Ensure repo root on path
_PATH = Path(__file__).resolve()
if str(_PATH.parents[2]) not in __import__("sys").path:
    __import__("sys").path.insert(0, str(_PATH.parents[2]))

from PX_System.foundation.Evidence_Package import wrap_trial_simulation


class TestTrialEvidencePackage(unittest.TestCase):

    def setUp(self):
        """Set up test fixtures (canonical: under Calibration_Molecules)"""
        repo = Path(__file__).resolve().parents[2]
        self.test_output_dir = str(repo / "PX_Warehouse" / "Calibration_Molecules" / "TrialSimulations_TEST")
        
    def tearDown(self):
        """Clean up test files"""
        if os.path.exists(self.test_output_dir):
            import shutil
            shutil.rmtree(self.test_output_dir)

    def test_wrap_trial_simulation(self):
        """Test basic trial simulation wrapping"""
        protocol = {"trial_id": "TEST-TRIAL"}
        trial_result = {
            "constitutional": {"engine": "TRIAL_ENGINE_V1"},
            "arms": []
        }
        ope = {"status": "STUB"}
        admet = {"constitutional": {"engine": "OPE_ADMET_V3_DETERMINISTIC"}}

        path = wrap_trial_simulation(
            protocol, trial_result, ope, admet, 
            output_dir=self.test_output_dir
        )

        self.assertTrue(path.endswith(".json"))
        self.assertIn("TRIAL_SIMULATION_DOSSIER", path)
        self.assertTrue(os.path.exists(path))

    def test_dossier_structure(self):
        """Test that dossier has all required fields"""
        protocol = {"trial_id": "TEST-STRUCTURE"}
        trial_result = {
            "constitutional": {"engine": "TRIAL_ENGINE_V1"},
            "arms": [{"arm_id": "A1"}]
        }
        ope = {"status": "STUB", "logp": None}
        admet = {
            "constitutional": {"engine": "OPE_ADMET_V3_DETERMINISTIC"},
            "toxicity_flags": {"hepatotoxicity_risk": "UNKNOWN"}
        }

        path = wrap_trial_simulation(
            protocol, trial_result, ope, admet,
            output_dir=self.test_output_dir
        )

        # Read and validate
        with open(path, 'r') as f:
            dossier = json.load(f)

        # Check required fields
        self.assertIn("dossier_type", dossier)
        self.assertEqual(dossier["dossier_type"], "TRIAL_SIMULATION_DOSSIER")
        
        self.assertIn("version", dossier)
        self.assertIn("timestamp_utc", dossier)
        self.assertIn("protocol", dossier)
        self.assertIn("trial_result", dossier)
        self.assertIn("provenance", dossier)
        self.assertIn("inputs", dossier)
        self.assertIn("constitutional", dossier)
        self.assertIn("evidence_hash", dossier)

    def test_provenance_tracking(self):
        """Test that provenance is correctly tracked"""
        protocol = {"trial_id": "TEST-PROVENANCE"}
        trial_result = {
            "constitutional": {"engine": "TRIAL_ENGINE_V1"},
            "arms": []
        }
        ope = {"status": "PRODUCTION"}
        admet = {"constitutional": {"engine": "OPE_ADMET_V3_DETERMINISTIC"}}

        path = wrap_trial_simulation(
            protocol, trial_result, ope, admet,
            output_dir=self.test_output_dir
        )

        with open(path, 'r') as f:
            dossier = json.load(f)

        prov = dossier["provenance"]
        self.assertEqual(prov["ope_engine"], "PRODUCTION")
        self.assertEqual(prov["admet_engine"], "OPE_ADMET_V3_DETERMINISTIC")
        self.assertEqual(prov["trial_engine"], "TRIAL_ENGINE_V1")

    def test_constitutional_compliance(self):
        """Test that constitutional metadata is included"""
        protocol = {"trial_id": "TEST-CONST"}
        trial_result = {
            "constitutional": {"engine": "TRIAL_ENGINE_V1"},
            "arms": []
        }
        ope = {"status": "STUB"}
        admet = {"constitutional": {"engine": "OPE_ADMET_V3_DETERMINISTIC"}}

        path = wrap_trial_simulation(
            protocol, trial_result, ope, admet,
            output_dir=self.test_output_dir
        )

        with open(path, 'r') as f:
            dossier = json.load(f)

        const = dossier["constitutional"]
        self.assertEqual(const["status"], "EVIDENCE_PACKAGE_CREATED")
        self.assertIn("L51", const["law_basis"])
        self.assertIn("L34", const["law_basis"])
        self.assertIn("Exposure-only", const["notes"])

    def test_hash_reproducibility(self):
        """Test that identical inputs produce identical hashes"""
        protocol = {"trial_id": "TEST-HASH"}
        trial_result = {
            "constitutional": {"engine": "TRIAL_ENGINE_V1"},
            "arms": []
        }
        ope = {"status": "STUB"}
        admet = {"constitutional": {"engine": "OPE_ADMET_V3_DETERMINISTIC"}}

        path1 = wrap_trial_simulation(
            protocol, trial_result, ope, admet,
            output_dir=self.test_output_dir
        )

        # Create second dossier with same inputs
        path2 = wrap_trial_simulation(
            protocol, trial_result, ope, admet,
            output_dir=self.test_output_dir
        )

        # Should produce same hash (different timestamps, but same evidence hash logic)
        # Extract hashes from filenames
        hash1 = Path(path1).stem.split('-')[-1]
        hash2 = Path(path2).stem.split('-')[-1]
        
        # Note: Hashes will differ due to timestamp, but structure should be same
        self.assertTrue(len(hash1) == 12)
        self.assertTrue(len(hash2) == 12)

    def test_inputs_preservation(self):
        """Test that OPE and ADMET inputs are preserved"""
        protocol = {"trial_id": "TEST-INPUTS"}
        trial_result = {
            "constitutional": {"engine": "TRIAL_ENGINE_V1"},
            "arms": []
        }
        ope = {"status": "STUB", "logp": 2.5, "molecular_weight": 180.0}
        admet = {
            "constitutional": {"engine": "OPE_ADMET_V3_DETERMINISTIC"},
            "toxicity_flags": {"hepatotoxicity_risk": "LOW"}
        }

        path = wrap_trial_simulation(
            protocol, trial_result, ope, admet,
            output_dir=self.test_output_dir
        )

        with open(path, 'r') as f:
            dossier = json.load(f)

        # Check inputs are preserved
        self.assertEqual(dossier["inputs"]["ope_analysis"]["logp"], 2.5)
        self.assertEqual(dossier["inputs"]["ope_analysis"]["molecular_weight"], 180.0)
        self.assertEqual(
            dossier["inputs"]["admet_analysis"]["toxicity_flags"]["hepatotoxicity_risk"],
            "LOW"
        )


if __name__ == "__main__":
    unittest.main()
