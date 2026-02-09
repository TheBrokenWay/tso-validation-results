"""
test_orchestrator_v2.py
Integration Tests for PX_Live_Orchestrator_v2

Tests full v2.0 computational pipeline end-to-end.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

import unittest
import json
from PX_Executive.orchestrators.PX_Live_Orchestrator_v2 import PredatorXOrchestratorV2


class TestOrchestratorV2(unittest.TestCase):
    """Test Orchestrator v2 pipeline"""
    
    def test_full_pipeline_aspirin(self):
        """Full pipeline should complete all stages for Aspirin; Zeus gate may reject warehouse write"""
        orchestrator = PredatorXOrchestratorV2(verbose=False)

        smiles = "CC(=O)Oc1ccccc1C(=O)O"  # Aspirin
        metadata = {"id": "TEST-ASPIRIN", "name": "Aspirin"}

        results = orchestrator.run_pipeline(smiles, metadata)

        # Check all 6 computational stages completed
        stages = results.get("stages", {})
        self.assertIn("ope", stages)
        self.assertIn("admet", stages)
        self.assertIn("trial_result", stages)
        self.assertIn("dose_optimization", stages)
        self.assertIn("virtual_efficacy", stages)

        # Evidence package is always built (analysis result)
        self.assertIn("evidence_package", results)

        # Pipeline either writes dossier (Zeus pass) or rejects (Zeus fail)
        if results.get("zeus_rejected"):
            # Toxic compound correctly rejected — no dossier written
            self.assertIn("zeus_verdict", results)
            self.assertNotIn("dossier_path", results)
        else:
            # Compound passed governance — dossier written
            self.assertIn("dossier_path", results)
            dossier_path = Path(results["dossier_path"])
            self.assertTrue(dossier_path.exists(), f"Dossier not found: {dossier_path}")
            with open(dossier_path, 'r') as f:
                dossier = json.load(f)
            self.assertIn("metadata", dossier)
            self.assertIn("ope", dossier)
            self.assertIn("admet", dossier)
    
    def test_pipeline_results_structure(self):
        """Pipeline results should have expected structure"""
        orchestrator = PredatorXOrchestratorV2(verbose=False)

        smiles = "CCO"  # Ethanol (simple test case)
        metadata = {"id": "TEST-CCO", "name": "Ethanol"}
        results = orchestrator.run_pipeline(smiles, metadata)

        # Top-level fields (v2 API)
        self.assertIn("version", results)
        self.assertIn("smiles", results)
        self.assertIn("run_id", results)
        self.assertIn("timestamp", results)
        self.assertIn("stages", results)
        self.assertIn("evidence_package", results)

        # Stage results under stages
        stages = results["stages"]
        self.assertIn("ope", stages)
        self.assertIn("admet", stages)
        self.assertIn("trial_result", stages)
        self.assertIn("dose_optimization", stages)
        self.assertIn("virtual_efficacy", stages)

        # Either dossier_path (Zeus pass) or zeus_rejected (Zeus fail)
        if results.get("zeus_rejected"):
            self.assertNotIn("dossier_path", results)
        else:
            self.assertIn("dossier_path", results)
            self.assertTrue(Path(results["dossier_path"]).exists())
    
    def test_dose_optimization_integration(self):
        """Dose optimization should integrate correctly"""
        orchestrator = PredatorXOrchestratorV2(verbose=False)
        
        smiles = "CC(=O)Oc1ccccc1C(=O)O"
        metadata = {"id": "TEST-ASPIRIN", "name": "Aspirin"}
        results = orchestrator.run_pipeline(smiles, metadata)
        
        dose_opt = results["stages"]["dose_optimization"]
        
        # Check structure
        self.assertIn("best_regimen", dose_opt)
        self.assertIn("search_history", dose_opt)
        self.assertIn("evaluations", dose_opt)
        
        # Check best regimen
        best = dose_opt["best_regimen"]
        self.assertIn("dose_mg", best)
        self.assertIn("interval_h", best)
        self.assertIn("score", best)
        
        # Check reasonable values
        self.assertGreater(best["dose_mg"], 0)
        self.assertGreater(best["interval_h"], 0)
    
    def test_virtual_efficacy_integration(self):
        """Virtual efficacy analytics should integrate correctly"""
        orchestrator = PredatorXOrchestratorV2(verbose=False)
        
        smiles = "CC(=O)Oc1ccccc1C(=O)O"
        metadata = {"id": "TEST-ASPIRIN", "name": "Aspirin"}
        results = orchestrator.run_pipeline(smiles, metadata)
        
        efficacy = results["stages"]["virtual_efficacy"]
        
        # Check structure
        self.assertIn("trial_id", efficacy)
        self.assertIn("constitutional", efficacy)
        
        # Check PTA results if present
        if "pk_pta" in efficacy:
            pta = efficacy["pk_pta"]
            self.assertIn("auc_mg_h_per_L", pta)
            
            auc_pta = pta["auc_mg_h_per_L"]
            self.assertIn("pta", auc_pta)
            
            # PTA should be between 0 and 1
            self.assertGreaterEqual(auc_pta["pta"], 0.0)
            self.assertLessEqual(auc_pta["pta"], 1.0)
    
    def test_evidence_package_v3_schema(self):
        """Evidence package should conform to v3 schema (regardless of Zeus gate outcome)"""
        orchestrator = PredatorXOrchestratorV2(verbose=False)

        smiles = "CC(=O)Oc1ccccc1C(=O)O"
        metadata = {"id": "TEST-ASPIRIN", "name": "Aspirin"}
        results = orchestrator.run_pipeline(smiles, metadata)

        # Evidence package is always in results (analysis is complete before Zeus gate)
        dossier = results["evidence_package"]

        # Evidence Package v3: metadata, ope, admet, pkpd, pd, dose_optimization, oce_gate, virtual_efficacy
        self.assertIn("metadata", dossier)
        self.assertIn("ope", dossier)
        self.assertIn("admet", dossier)
        self.assertIn("pkpd", dossier)
        self.assertIn("dose_optimization", dossier)
        self.assertIn("virtual_efficacy", dossier)
        self.assertIn("zeus", dossier)

        meta = dossier.get("metadata", {})
        self.assertIn("version", meta)
        self.assertEqual(meta.get("version"), "3.0")


def run_tests():
    """Run orchestrator v2 tests"""
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    suite.addTests(loader.loadTestsFromTestCase(TestOrchestratorV2))
    
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    print("\n" + "="*70)
    print(f"ORCHESTRATOR V2 TESTS SUMMARY")
    print("="*70)
    print(f"Tests run:     {result.testsRun}")
    print(f"Successes:     {result.testsRun - len(result.failures) - len(result.errors)}")
    print(f"Failures:      {len(result.failures)}")
    print(f"Errors:        {len(result.errors)}")
    print("="*70)
    
    if result.wasSuccessful():
        print("✅ ALL ORCHESTRATOR V2 TESTS PASSED")
        return 0
    else:
        print("❌ SOME TESTS FAILED")
        return 1


if __name__ == "__main__":
    exit_code = run_tests()
    sys.exit(exit_code)
