"""Tests for TSO_Validator reproducibility checks."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import unittest
from audit.reproducibility import check_determinism, check_records_deterministic


class TestReproducibility(unittest.TestCase):

    def test_clean_record_passes(self):
        rec = {"record_id": "R001", "toxicity_index": 0.01}
        failures = check_determinism(rec)
        self.assertEqual(failures, [])

    def test_random_without_seed_fails(self):
        rec = {"record_id": "R002", "generation_method": "random"}
        failures = check_determinism(rec)
        self.assertTrue(any("rng_seed" in f for f in failures))

    def test_random_with_seed_passes(self):
        rec = {"record_id": "R003", "generation_method": "random", "rng_seed": 42}
        failures = check_determinism(rec)
        self.assertEqual(failures, [])

    def test_nan_detected(self):
        rec = {"record_id": "R004", "toxicity_index": float("nan")}
        failures = check_determinism(rec)
        self.assertTrue(any("NaN" in f for f in failures))

    def test_inf_detected(self):
        rec = {"record_id": "R005", "binding_affinity_kj": float("inf")}
        failures = check_determinism(rec)
        self.assertTrue(any("Infinity" in f for f in failures))

    def test_batch_check(self):
        records = [
            {"record_id": "R006", "toxicity_index": 0.01},
            {"record_id": "R007", "toxicity_index": float("nan")},
        ]
        result = check_records_deterministic(records)
        self.assertFalse(result["passed"])
        self.assertEqual(result["records_failed"], 1)


if __name__ == "__main__":
    unittest.main()
