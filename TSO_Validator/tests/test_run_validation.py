"""Tests for TSO_Validator end-to-end run_validation."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import json
import shutil
import tempfile
import unittest
from unittest.mock import patch
from run_validation import run_validation, RAW_DIR, RESULTS_DIR, TSO_ROOT


class TestRunValidation(unittest.TestCase):

    def setUp(self):
        """Back up and clear raw/results dirs for isolated tests."""
        self._backup_raw = []
        if RAW_DIR.exists():
            for f in RAW_DIR.glob("*.json"):
                self._backup_raw.append((f.name, f.read_text()))
                f.unlink()
        self._orig_results = list(RESULTS_DIR.glob("run_*")) if RESULTS_DIR.exists() else []

    def tearDown(self):
        """Restore original raw files."""
        RAW_DIR.mkdir(parents=True, exist_ok=True)
        for name, content in self._backup_raw:
            (RAW_DIR / name).write_text(content)

    def test_empty_raw_returns_validated(self):
        summary = run_validation()
        self.assertEqual(summary["status"], "TSO_VALIDATED")
        self.assertEqual(summary["files_checked"], 0)

    def test_clean_record_validates(self):
        RAW_DIR.mkdir(parents=True, exist_ok=True)
        rec = {"record_id": "CLEAN-001", "timestamp": "2026-01-01T00:00:00Z", "toxicity_index": 0.005}
        (RAW_DIR / "clean.json").write_text(json.dumps(rec))
        summary = run_validation()
        self.assertEqual(summary["status"], "TSO_VALIDATED")
        self.assertEqual(summary["files_checked"], 1)
        self.assertEqual(summary["failures"], [])
        # Cleanup
        (RAW_DIR / "clean.json").unlink()

    def test_toxic_record_fails(self):
        RAW_DIR.mkdir(parents=True, exist_ok=True)
        rec = {"record_id": "TOX-001", "timestamp": "2026-01-01T00:00:00Z", "toxicity_index": 0.05}
        (RAW_DIR / "toxic.json").write_text(json.dumps(rec))
        summary = run_validation()
        self.assertEqual(summary["status"], "TSO_FAILED")
        self.assertTrue(any("L11" in f for f in summary["failures"]))
        # Cleanup
        (RAW_DIR / "toxic.json").unlink()

    def test_run_creates_output_files(self):
        RAW_DIR.mkdir(parents=True, exist_ok=True)
        rec = {"record_id": "OUT-001", "timestamp": "2026-01-01T00:00:00Z"}
        (RAW_DIR / "out.json").write_text(json.dumps(rec))
        summary = run_validation()
        run_dir = RESULTS_DIR / summary["run_id"]
        self.assertTrue((run_dir / "run_summary.json").exists())
        self.assertTrue((run_dir / "provenance.json").exists())
        self.assertTrue((run_dir / "historical_context.json").exists())
        # Cleanup
        (RAW_DIR / "out.json").unlink()
        shutil.rmtree(run_dir, ignore_errors=True)

    def test_malformed_json_fails(self):
        RAW_DIR.mkdir(parents=True, exist_ok=True)
        (RAW_DIR / "bad.json").write_text("{this is not valid json")
        summary = run_validation()
        self.assertEqual(summary["status"], "TSO_FAILED")
        self.assertTrue(any("Parse error" in f for f in summary["failures"]))
        # Cleanup
        (RAW_DIR / "bad.json").unlink()

    def test_no_px_imports(self):
        """TSO_Validator must not import any PX_* module."""
        tso_files = list(TSO_ROOT.rglob("*.py"))
        for f in tso_files:
            content = f.read_text(encoding="utf-8")
            in_docstring = False
            for line in content.splitlines():
                stripped = line.strip()
                if stripped.startswith("#"):
                    continue
                # Skip docstrings
                if '"""' in stripped or "'''" in stripped:
                    if stripped.count('"""') == 1 or stripped.count("'''") == 1:
                        in_docstring = not in_docstring
                    continue
                if in_docstring:
                    continue
                # Check actual import statements only
                if stripped.startswith("from PX_") or stripped.startswith("import PX_"):
                    self.fail(f"TSO_Validator file {f.name} imports PX_* module: {stripped}")


if __name__ == "__main__":
    unittest.main()
