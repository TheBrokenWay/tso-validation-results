"""Basic module tests for Lymphatic filariasis."""

import os
import sys
import unittest

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)


class TestLymphaticFilariasisModule(unittest.TestCase):

    def test_module_imports(self):
        from PX_Domain.PRV_Diseases.Lymphatic_Filariasis import adapters, analytics, validation
        self.assertIsNotNone(adapters)
        self.assertIsNotNone(analytics)
        self.assertIsNotNone(validation)

    def test_constraint_file_exists(self):
        config_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "config")
        constraint_file = os.path.join(config_dir, "lymphatic_filariasis_constraints.json")
        self.assertTrue(os.path.isfile(constraint_file), f"Missing: {constraint_file}")

    def test_adapter_class_exists(self):
        from PX_Domain.PRV_Diseases.Lymphatic_Filariasis.adapters.lymphatic_filariasis_adapter import LymphaticFilariasisAdapter
        adapter = LymphaticFilariasisAdapter()
        self.assertEqual(adapter.disease_id, "lymphatic_filariasis")


if __name__ == "__main__":
    unittest.main()
