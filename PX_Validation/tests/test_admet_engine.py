"""
Unit tests for ADMET Engine
Tests deterministic hepatotoxicity risk assessment
"""

import unittest
import sys
import os

# Add foundation to path
ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if ROOT_DIR not in sys.path:
    sys.path.append(ROOT_DIR)

from PX_Engine.operations.ADMET import run_admet

class TestADMET(unittest.TestCase):
    """Test suite for ADMET engine"""

    def test_hepatotoxicity_low(self):
        """Test LOW hepatotoxicity risk for logP < 3.5"""
        ope = {"logp": 2.5}
        admet = run_admet("CCO", ope)
        self.assertEqual(admet["toxicity_flags"]["hepatotoxicity_risk"], "LOW")

    def test_hepatotoxicity_medium(self):
        """Test MEDIUM hepatotoxicity risk for 3.5 < logP < 4.5"""
        ope = {"logp": 3.8}
        admet = run_admet("CCO", ope)
        self.assertEqual(admet["toxicity_flags"]["hepatotoxicity_risk"], "MEDIUM")

    def test_hepatotoxicity_high(self):
        """Test HIGH hepatotoxicity risk for logP > 4.5"""
        ope = {"logp": 5.0}
        admet = run_admet("CCO", ope)
        self.assertEqual(admet["toxicity_flags"]["hepatotoxicity_risk"], "HIGH")

    def test_unknown_logp(self):
        """Test UNKNOWN hepatotoxicity risk when logP is not provided"""
        ope = {}
        admet = run_admet("CCO", ope)
        self.assertEqual(admet["toxicity_flags"]["hepatotoxicity_risk"], "UNKNOWN")

    def test_admet_structure(self):
        """Test that ADMET returns expected structure"""
        ope = {"logp": 3.0}
        admet = run_admet("CCO", ope)
        
        # Check all required sections exist
        self.assertIn("absorption", admet)
        self.assertIn("distribution", admet)
        self.assertIn("metabolism", admet)
        self.assertIn("excretion", admet)
        self.assertIn("toxicity_flags", admet)
        self.assertIn("constitutional", admet)
        
        # Check constitutional compliance (engine may report VERIFIED or PARTIAL)
        self.assertIn(admet["constitutional"]["status"], ("PARTIAL", "VERIFIED"))
        self.assertIn("L11", admet["constitutional"]["law_basis"])
        # 4-tier risk_level (Order 66 taxonomy)
        self.assertIn(admet["toxicity"]["risk_level"], ("TOXICITY_DIAMOND", "TOXICITY_GOLD", "TOXICITY_SILVER", "TOXICITY_FAILURE"))

    def test_none_values_compliance(self):
        """Test that required keys exist and toxicity_flags are set"""
        ope = {"logp": 3.0}
        admet = run_admet("CCO", ope)
        
        self.assertIn("predicted_bioavailability", admet["absorption"])
        self.assertIn("bcs_class", admet["absorption"])
        self.assertIn("predicted_vd_L_per_kg", admet["distribution"])
        self.assertIn("bbb_penetration", admet["distribution"])
        self.assertIn(admet["toxicity_flags"]["cardiotoxicity_risk"], ("HIGH", "LOW", "UNKNOWN"))

if __name__ == "__main__":
    unittest.main()
