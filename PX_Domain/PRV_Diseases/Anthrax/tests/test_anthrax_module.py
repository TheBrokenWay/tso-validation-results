"""Basic module tests for Anthrax."""

import os
import sys
import unittest

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)


class TestAnthraxModule(unittest.TestCase):

    def test_module_imports(self):
        from PX_Domain.PRV_Diseases.Anthrax import adapters, analytics, validation
        self.assertIsNotNone(adapters)
        self.assertIsNotNone(analytics)
        self.assertIsNotNone(validation)

    def test_constraint_file_exists(self):
        config_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "config")
        constraint_file = os.path.join(config_dir, "anthrax_constraints.json")
        self.assertTrue(os.path.isfile(constraint_file), f"Missing: {constraint_file}")

    def test_adapter_class_exists(self):
        from PX_Domain.PRV_Diseases.Anthrax.adapters.anthrax_adapter import AnthraxAdapter
        adapter = AnthraxAdapter()
        self.assertEqual(adapter.disease_id, "anthrax")


if __name__ == "__main__":
    unittest.main()
