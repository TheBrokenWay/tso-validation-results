"""
Create folder structures for all PRV diseases (except Nipah).

Reads PX_Domain/PRV_Diseases/manifest.json, skips Nipah (already has full module),
and generates standardised module folders for all 30 remaining diseases.
"""

import json
import os
import shutil
import sys
import textwrap

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PRV_DIR = os.path.join(_REPO_ROOT, "PX_Domain", "PRV_Diseases")
MANIFEST_PATH = os.path.join(PRV_DIR, "manifest.json")


# ── helpers ──────────────────────────────────────────────────────────


def name_to_folder(name: str) -> str:
    """Convert disease name to TitleCase_Underscore folder name.
    e.g. 'Ebola virus disease' -> 'Ebola_Virus_Disease'
    """
    return "_".join(word.capitalize() for word in name.split())


def name_to_id(name: str) -> str:
    """Convert disease name to lowercase_underscore id.
    e.g. 'Ebola virus disease' -> 'ebola_virus_disease'
    """
    return "_".join(word.lower() for word in name.split())


def name_to_class(name: str) -> str:
    """Convert disease name to PascalCase class name.
    e.g. 'Ebola virus disease' -> 'EbolaVirusDisease'
    """
    return "".join(word.capitalize() for word in name.split())


def write_file(path: str, content: str) -> None:
    """Write content to path, creating parent dirs as needed."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(content)


def write_gitkeep(dir_path: str) -> None:
    """Create directory with .gitkeep."""
    os.makedirs(dir_path, exist_ok=True)
    gk = os.path.join(dir_path, ".gitkeep")
    if not os.path.exists(gk):
        with open(gk, "w", encoding="utf-8") as f:
            pass


# ── template generators ──────────────────────────────────────────────


def gen_top_init(disease_name: str) -> str:
    return f'"""{disease_name} analysis module."""\n'


def gen_adapters_init(disease_name: str) -> str:
    return f'"""{disease_name} data adapters."""\n'


def gen_adapter(disease_name: str, disease_id: str, disease_class: str) -> str:
    return textwrap.dedent(f'''\
        """
        {disease_name} data adapter.

        TODO: Implement disease-specific data mining and ingestion.
        Follows the pattern established by Nipah adapter.
        """

        import os
        import sys

        _REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))
        if _REPO_ROOT not in sys.path:
            sys.path.insert(0, _REPO_ROOT)


        class {disease_class}Adapter:
            """Data adapter for {disease_name} research data."""

            def __init__(self):
                self.disease_id = "{disease_id}"
                self.disease_name = "{disease_name}"

            def mine(self, filepath: str) -> dict:
                """Mine raw data file and return structured results.

                TODO: Implement disease-specific data mining.
                """
                raise NotImplementedError(f"{{self.disease_name}} adapter not yet implemented")

            def validate(self, data: dict) -> bool:
                """Validate mined data against disease constraints.

                TODO: Implement disease-specific validation.
                """
                raise NotImplementedError(f"{{self.disease_name}} validation not yet implemented")
    ''')


def gen_analytics_init(disease_name: str) -> str:
    return f'"""{disease_name} analytics."""\n'


def gen_analytics(disease_name: str, disease_id: str, disease_class: str) -> str:
    return textwrap.dedent(f'''\
        """
        {disease_name} analytics module.

        TODO: Implement disease-specific analysis pipelines.
        """


        class {disease_class}Analytics:
            """Analytics engine for {disease_name} compound evaluation."""

            def __init__(self):
                self.disease_id = "{disease_id}"

            def analyze(self, compound_data: dict) -> dict:
                """Run disease-specific analysis on compound data.

                TODO: Implement analysis pipeline.
                """
                return {{
                    "disease_id": self.disease_id,
                    "status": "NOT_IMPLEMENTED",
                    "message": f"{{self.disease_id}} analytics pending implementation",
                }}
    ''')


def gen_tests_init(disease_name: str) -> str:
    return f'"""{disease_name} tests."""\n'


def gen_test_module(disease_name: str, disease_id: str, disease_class: str, folder_name: str) -> str:
    return textwrap.dedent(f'''\
        """Basic module tests for {disease_name}."""

        import os
        import sys
        import unittest

        _REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))
        if _REPO_ROOT not in sys.path:
            sys.path.insert(0, _REPO_ROOT)


        class Test{disease_class}Module(unittest.TestCase):

            def test_module_imports(self):
                from PX_Domain.PRV_Diseases.{folder_name} import adapters, analytics, validation
                self.assertIsNotNone(adapters)
                self.assertIsNotNone(analytics)
                self.assertIsNotNone(validation)

            def test_constraint_file_exists(self):
                config_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "config")
                constraint_file = os.path.join(config_dir, "{disease_id}_constraints.json")
                self.assertTrue(os.path.isfile(constraint_file), f"Missing: {{constraint_file}}")

            def test_adapter_class_exists(self):
                from PX_Domain.PRV_Diseases.{folder_name}.adapters.{disease_id}_adapter import {disease_class}Adapter
                adapter = {disease_class}Adapter()
                self.assertEqual(adapter.disease_id, "{disease_id}")


        if __name__ == "__main__":
            unittest.main()
    ''')


def gen_validation_init(disease_name: str) -> str:
    return f'"""{disease_name} data validation."""\n'


def gen_validation(disease_name: str, disease_id: str, disease_class: str) -> str:
    return textwrap.dedent(f'''\
        """
        {disease_name} data validation.

        TODO: Implement disease-specific data validation rules.
        Follows the pattern established by Nipah validator.
        """

        import json
        import os


        class {disease_class}ValidationError(Exception):
            pass


        def validate_data(data: dict) -> bool:
            """Validate data against {disease_name} constraints.

            TODO: Implement disease-specific validation rules.
            """
            if not data:
                raise {disease_class}ValidationError("Empty data")
            return True


        def load_constraints() -> dict:
            """Load disease constraint file."""
            config_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "config")
            constraint_path = os.path.join(config_dir, "{disease_id}_constraints.json")
            with open(constraint_path, "r", encoding="utf-8") as f:
                return json.load(f)
    ''')


# ── main ─────────────────────────────────────────────────────────────


def main():
    # Load manifest
    with open(MANIFEST_PATH, "r", encoding="utf-8") as f:
        manifest = json.load(f)

    diseases = manifest["prv_diseases"]
    print(f"Manifest contains {len(diseases)} diseases.")

    created = []
    skipped = []
    errors = []

    for entry in diseases:
        disease_name = entry["name"]

        # Skip Nipah
        if "nipah" in disease_name.lower():
            skipped.append(disease_name)
            print(f"  SKIP: {disease_name} (already has full module)")
            continue

        folder_name = name_to_folder(disease_name)
        disease_id = name_to_id(disease_name)
        disease_class = name_to_class(disease_name)
        base = os.path.join(PRV_DIR, folder_name)

        # Check if folder already exists
        if os.path.isdir(base):
            skipped.append(disease_name)
            print(f"  SKIP: {disease_name} (folder already exists: {folder_name}/)")
            continue

        try:
            # Create base directory
            os.makedirs(base, exist_ok=True)

            # __init__.py (top level)
            write_file(os.path.join(base, "__init__.py"), gen_top_init(disease_name))

            # adapters/
            adapters_dir = os.path.join(base, "adapters")
            os.makedirs(adapters_dir, exist_ok=True)
            write_file(os.path.join(adapters_dir, "__init__.py"), gen_adapters_init(disease_name))
            write_file(
                os.path.join(adapters_dir, f"{disease_id}_adapter.py"),
                gen_adapter(disease_name, disease_id, disease_class),
            )

            # analytics/
            analytics_dir = os.path.join(base, "analytics")
            os.makedirs(analytics_dir, exist_ok=True)
            write_file(os.path.join(analytics_dir, "__init__.py"), gen_analytics_init(disease_name))
            write_file(
                os.path.join(analytics_dir, f"{disease_id}_analytics.py"),
                gen_analytics(disease_name, disease_id, disease_class),
            )

            # config/  (copy constraint JSON)
            config_dir = os.path.join(base, "config")
            os.makedirs(config_dir, exist_ok=True)
            src_constraint = os.path.join(PRV_DIR, f"{disease_id}.json")
            dst_constraint = os.path.join(config_dir, f"{disease_id}_constraints.json")
            if os.path.isfile(src_constraint):
                shutil.copy2(src_constraint, dst_constraint)
            else:
                print(f"  WARNING: No constraint file found at {src_constraint}")

            # data/raw/ and data/processed/
            write_gitkeep(os.path.join(base, "data", "raw"))
            write_gitkeep(os.path.join(base, "data", "processed"))

            # results/
            write_gitkeep(os.path.join(base, "results"))

            # tests/
            tests_dir = os.path.join(base, "tests")
            os.makedirs(tests_dir, exist_ok=True)
            write_file(os.path.join(tests_dir, "__init__.py"), gen_tests_init(disease_name))
            write_file(
                os.path.join(tests_dir, f"test_{disease_id}_module.py"),
                gen_test_module(disease_name, disease_id, disease_class, folder_name),
            )

            # validation/
            validation_dir = os.path.join(base, "validation")
            os.makedirs(validation_dir, exist_ok=True)
            write_file(os.path.join(validation_dir, "__init__.py"), gen_validation_init(disease_name))
            write_file(
                os.path.join(validation_dir, f"validate_{disease_id}_data.py"),
                gen_validation(disease_name, disease_id, disease_class),
            )

            created.append(disease_name)
            print(f"  CREATED: {folder_name}/")

        except Exception as exc:
            errors.append((disease_name, str(exc)))
            print(f"  ERROR: {disease_name} -> {exc}", file=sys.stderr)

    # ── Summary ──
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"  Total in manifest : {len(diseases)}")
    print(f"  Created           : {len(created)}")
    print(f"  Skipped           : {len(skipped)}")
    print(f"  Errors            : {len(errors)}")
    print()

    if created:
        print("Created modules:")
        for name in sorted(created):
            folder = name_to_folder(name)
            print(f"  PX_Domain/PRV_Diseases/{folder}/")

    if errors:
        print("\nErrors:")
        for name, err in errors:
            print(f"  {name}: {err}")

    print()
    return 0 if not errors else 1


if __name__ == "__main__":
    sys.exit(main())
