"""Tests for TSO_Validator ethical guard."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import unittest
from audit.ethical_guard import (
    check_no_harm, check_no_mock, run_ethical_guard,
    TOXICITY_HARD_LIMIT, HARM_ENERGY_HARD_LIMIT, HARMONIC_OVERDRIVE,
)


class TestEthicalGuard(unittest.TestCase):

    def test_clean_record_passes(self):
        rec = {"record_id": "R001", "toxicity_index": 0.005, "harm_energy": 0.003}
        result = run_ethical_guard(rec)
        self.assertTrue(result["passed"])
        self.assertEqual(result["failures"], [])

    def test_toxicity_over_limit_fails(self):
        rec = {"record_id": "R002", "toxicity_index": 0.025}
        failures = check_no_harm(rec)
        self.assertTrue(any("L11" in f for f in failures))

    def test_toxicity_at_boundary_fails(self):
        """0.0210 is >= limit — must fail (no rounding)."""
        rec = {"record_id": "R003", "toxicity_index": 0.0210}
        failures = check_no_harm(rec)
        self.assertTrue(any("L11" in f for f in failures))

    def test_toxicity_below_boundary_passes(self):
        rec = {"record_id": "R004", "toxicity_index": 0.0209}
        failures = check_no_harm(rec)
        self.assertEqual(failures, [])

    def test_harm_energy_over_limit_fails(self):
        rec = {"record_id": "R005", "harm_energy": 0.0210}
        failures = check_no_harm(rec)
        self.assertTrue(any("L10" in f for f in failures))

    def test_mock_record_id_detected(self):
        rec = {"record_id": "MOCK-001", "toxicity_index": 0.001}
        failures = check_no_mock(rec)
        self.assertTrue(any("Mock record_id" in f for f in failures))

    def test_string_value_detected_as_mock(self):
        rec = {"record_id": "R006", "toxicity_index": "low"}
        failures = check_no_mock(rec)
        self.assertTrue(any("Mock data" in f for f in failures))

    def test_placeholder_value_detected(self):
        rec = {"record_id": "R007", "toxicity_index": -1}
        failures = check_no_mock(rec)
        self.assertTrue(any("Placeholder" in f for f in failures))

    def test_harmonic_overdrive_is_fixed(self):
        self.assertEqual(HARMONIC_OVERDRIVE, 1.02)

    def test_thresholds_match_constitution(self):
        self.assertEqual(TOXICITY_HARD_LIMIT, 0.0210)
        self.assertEqual(HARM_ENERGY_HARD_LIMIT, 0.0210)

    def test_missing_fields_pass(self):
        """Records without physics fields pass (fields are optional)."""
        rec = {"record_id": "R008", "timestamp": "2026-01-01T00:00:00Z"}
        result = run_ethical_guard(rec)
        self.assertTrue(result["passed"])


class TestNoSycophancy(unittest.TestCase):
    """Law 21: thresholds are absolute. No adjustments."""

    def test_cannot_round_past_threshold(self):
        """0.02099999 rounds to 0.021 but raw value is < 0.0210 — MUST pass."""
        rec = {"record_id": "R010", "toxicity_index": 0.02099999}
        result = run_ethical_guard(rec)
        self.assertTrue(result["passed"])

    def test_at_threshold_fails_no_exception(self):
        """Exactly 0.0210 MUST fail. No special cases."""
        rec = {"record_id": "R011", "toxicity_index": 0.0210}
        result = run_ethical_guard(rec)
        self.assertFalse(result["passed"])


if __name__ == "__main__":
    unittest.main()
