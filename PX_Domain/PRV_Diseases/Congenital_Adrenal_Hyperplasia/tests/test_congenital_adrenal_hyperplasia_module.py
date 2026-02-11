"""Basic module tests for Congenital adrenal hyperplasia."""

import os
import sys
import unittest

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)


class TestCongenitalAdrenalHyperplasiaModule(unittest.TestCase):

    def test_module_imports(self):
        from PX_Domain.PRV_Diseases.Congenital_Adrenal_Hyperplasia import adapters, analytics, validation
        self.assertIsNotNone(adapters)
        self.assertIsNotNone(analytics)
        self.assertIsNotNone(validation)

    def test_constraint_file_exists(self):
        config_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "config")
        constraint_file = os.path.join(config_dir, "congenital_adrenal_hyperplasia_constraints.json")
        self.assertTrue(os.path.isfile(constraint_file), f"Missing: {constraint_file}")

    def test_adapter_class_exists(self):
        from PX_Domain.PRV_Diseases.Congenital_Adrenal_Hyperplasia.adapters.congenital_adrenal_hyperplasia_adapter import CongenitalAdrenalHyperplasiaAdapter
        adapter = CongenitalAdrenalHyperplasiaAdapter()
        self.assertEqual(adapter.disease_id, "congenital_adrenal_hyperplasia")


if __name__ == "__main__":
    unittest.main()
