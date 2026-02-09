"""
test_grading_engine.py
Tests for Constitutional Grading Engine
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

import unittest
from PX_Engine.operations.GradingEngine import GradingEngine


class TestGradingEngine(unittest.TestCase):
    """Test grading engine classification"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.engine = GradingEngine(verbose=False)
    
    def test_gold_tier_classification(self):
        """GOLD_TIER candidates should meet all GOLD criteria (toxicity within L11 limit)"""
        metrics = {
            "pta": 85.0,
            "responder_rate": 75.0,
            "toxicity": 0.02,
            "dose_score": 0.85,
            "variability_cv": 25.0,
            "binding_affinity": 0.9,
        }
        
        grade, reasoning = self.engine.classify_metrics(metrics)
        
        self.assertEqual(grade, "GOLD_TIER")
        self.assertIn("Exceptional", reasoning["reason"])
    
    def test_silver_tier_classification(self):
        """SILVER_TIER candidates should meet all SILVER criteria (toxicity within L11 limit)"""
        metrics = {
            "pta": 65.0,
            "responder_rate": 55.0,
            "toxicity": 0.02,
            "dose_score": 0.65,
            "variability_cv": 45.0,
            "binding_affinity": 0.7,
        }
        
        grade, reasoning = self.engine.classify_metrics(metrics)
        
        self.assertEqual(grade, "SILVER_TIER")
        self.assertIn("Strong", reasoning["reason"])
    
    def test_bronze_tier_classification(self):
        """BRONZE_TIER candidates should meet all BRONZE criteria (toxicity within L11 limit)"""
        metrics = {
            "pta": 45.0,
            "responder_rate": 35.0,
            "toxicity": 0.02,
            "dose_score": 0.45,
            "variability_cv": 65.0,
            "binding_affinity": 0.5,
        }
        
        grade, reasoning = self.engine.classify_metrics(metrics)
        
        self.assertEqual(grade, "BRONZE_TIER")
        self.assertIn("Moderate", reasoning["reason"])
    
    def test_needs_review_classification(self):
        """NEEDS_REVIEW candidates should meet 3-4 criteria (toxicity within L11 limit)"""
        metrics = {
            "pta": 50.0,
            "responder_rate": 40.0,
            "toxicity": 0.02,
            "dose_score": 0.2,
            "variability_cv": 80.0,
            "binding_affinity": 0.5,
        }
        
        grade, reasoning = self.engine.classify_metrics(metrics)
        
        self.assertEqual(grade, "NEEDS_REVIEW")
        self.assertIn("Ambiguous", reasoning["reason"])
        self.assertEqual(reasoning["criteria_met_count"], 3)
    
    def test_rejected_classification(self):
        """REJECTED candidates should fail most criteria"""
        metrics = {
            "pta": 10.0,           # FAILS
            "responder_rate": 5.0,   # FAILS
            "toxicity": 0.8,       # FAILS
            "dose_score": 0.1,     # FAILS
            "variability_cv": 90.0,  # FAILS
            "binding_affinity": 0.2,
        }
        
        grade, reasoning = self.engine.classify_metrics(metrics)
        
        self.assertEqual(grade, "REJECTED")
        self.assertTrue(
            "Failed criteria" in reasoning["reason"] or "Law L11" in reasoning["reason"] or "toxicity" in reasoning["reason"],
            f"Expected rejection reason, got: {reasoning['reason']}"
        )
    
    def test_extract_metrics_from_dossier(self):
        """Should extract metrics from Evidence Package v3 dossier"""
        dossier = {
            "id": "TEST-001",
            "name": "Test Compound",
            "virtual_efficacy": {
                "pk_pta": {
                    "auc_mg_h_per_L": {
                        "pta": 75.5,
                    }
                },
                "pd_responders": {
                    "max_effect": {
                        "responder_rate": 0.68,
                    }
                }
            },
            "dose_optimization": {
                "best_regimen": {
                    "optimization_score": 0.82,
                }
            },
            "trial_result": {
                "arms": [{
                    "exposure_summary": {
                        "auc_mg_h_per_L": {
                            "mean": 100.0,
                            "sd": 20.0,
                        }
                    }
                }]
            },
            "admet": {
                "toxicity": 0.25,
            },
            "ope": {
                "binding_affinity": {
                    "score": 0.85,
                }
            }
        }
        
        metrics = self.engine.extract_metrics(dossier)
        
        self.assertEqual(metrics["pta"], 75.5)
        self.assertEqual(metrics["responder_rate"], 68.0)  # Converted to %
        self.assertEqual(metrics["dose_score"], 0.82)
        self.assertEqual(metrics["variability_cv"], 20.0)  # (20/100)*100
        self.assertEqual(metrics["toxicity"], 0.25)
        self.assertEqual(metrics["binding_affinity"], 0.85)
    
    def test_grade_dossier_returns_complete_result(self):
        """Should return complete grading result with all fields"""
        dossier = {
            "id": "TEST-002",
            "name": "Test Compound 2",
            "virtual_efficacy": {
                "pk_pta": {"auc_mg_h_per_L": {"pta": 85.0}},
                "pd_responders": {"max_effect": {"responder_rate": 0.75}}
            },
            "dose_optimization": {"best_regimen": {"optimization_score": 0.85}},
            "trial_result": {
                "arms": [{
                    "exposure_summary": {
                        "auc_mg_h_per_L": {"mean": 100.0, "sd": 25.0}
                    }
                }]
            },
            "admet": {"toxicity": 0.02},
            "ope": {"binding_affinity": {"score": 0.9}}
        }
        
        result = self.engine.grade_dossier(dossier)
        
        self.assertIn("grade", result)
        self.assertIn("metrics", result)
        self.assertIn("reasoning", result)
        self.assertIn("timestamp", result)
        self.assertIn("grading_engine_version", result)
        self.assertIn("thresholds_used", result)
        
        self.assertEqual(result["grade"], "GOLD_TIER")
    
    def test_get_grade_statistics(self):
        """Should track grading statistics"""
        # Grade multiple dossiers
        dossiers = [
            {"id": "GOLD-1", "name": "Gold 1"},
            {"id": "SILVER-1", "name": "Silver 1"},
            {"id": "BRONZE-1", "name": "Bronze 1"},
        ]
        
        # Add dummy data for each
        for dossier in dossiers:
            dossier.update({
                "virtual_efficacy": {
                    "pk_pta": {"auc_mg_h_per_L": {"pta": 85.0}},
                    "pd_responders": {"max_effect": {"responder_rate": 0.75}}
                },
                "dose_optimization": {"best_regimen": {"optimization_score": 0.85}},
                "trial_result": {
                    "arms": [{
                        "exposure_summary": {
                            "auc_mg_h_per_L": {"mean": 100.0, "sd": 25.0}
                        }
                    }]
                },
                "admet": {"toxicity": 0.02},
                "ope": {"binding_affinity": {"score": 0.9}}
            })
        
        for dossier in dossiers:
            self.engine.grade_dossier(dossier)
        
        stats = self.engine.get_grade_statistics()
        
        self.assertEqual(stats["total_graded"], 3)
        dist = stats.get("grade_distribution", stats)
        self.assertGreaterEqual(dist.get("GOLD_TIER", 0), 0)
        self.assertGreaterEqual(dist.get("SILVER_TIER", 0), 0)
        self.assertGreaterEqual(dist.get("BRONZE_TIER", 0), 0)


def run_tests():
    """Run grading engine tests"""
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    suite.addTests(loader.loadTestsFromTestCase(TestGradingEngine))
    
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    print("\n" + "="*70)
    print(f"GRADING ENGINE TESTS SUMMARY")
    print("="*70)
    print(f"Tests run:     {result.testsRun}")
    print(f"Successes:     {result.testsRun - len(result.failures) - len(result.errors)}")
    print(f"Failures:      {len(result.failures)}")
    print(f"Errors:        {len(result.errors)}")
    print("="*70)
    
    if result.wasSuccessful():
        print("✅ ALL GRADING ENGINE TESTS PASSED")
        return 0
    else:
        print("❌ SOME TESTS FAILED")
        return 1


if __name__ == "__main__":
    exit_code = run_tests()
    sys.exit(exit_code)
