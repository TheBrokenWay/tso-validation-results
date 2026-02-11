"""Basic module tests for Duchenne muscular dystrophy."""

import os
import sys
import unittest

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)


class TestDuchenneMuscularDystrophyModule(unittest.TestCase):

    def test_module_imports(self):
        from PX_Domain.PRV_Diseases.Duchenne_Muscular_Dystrophy import adapters, analytics, validation
        self.assertIsNotNone(adapters)
        self.assertIsNotNone(analytics)
        self.assertIsNotNone(validation)

    def test_constraint_file_exists(self):
        config_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "config")
        constraint_file = os.path.join(config_dir, "duchenne_muscular_dystrophy_constraints.json")
        self.assertTrue(os.path.isfile(constraint_file), f"Missing: {constraint_file}")

    def test_adapter_class_exists(self):
        from PX_Domain.PRV_Diseases.Duchenne_Muscular_Dystrophy.adapters.duchenne_muscular_dystrophy_adapter import DuchenneMuscularDystrophyAdapter
        adapter = DuchenneMuscularDystrophyAdapter()
        self.assertEqual(adapter.disease_id, "duchenne_muscular_dystrophy")


if __name__ == "__main__":
    unittest.main()
