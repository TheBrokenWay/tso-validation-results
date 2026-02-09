"""Tests for TSO_Validator provenance tracking."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import json
import tempfile
import unittest
from audit.provenance import compute_checksum, build_provenance, verify_provenance


class TestProvenance(unittest.TestCase):

    def test_checksum_deterministic(self):
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump({"test": 1}, f)
            p = Path(f.name)
        try:
            h1 = compute_checksum(p)
            h2 = compute_checksum(p)
            self.assertEqual(h1, h2)
            self.assertEqual(len(h1), 64)  # SHA-256 hex length
        finally:
            p.unlink()

    def test_build_provenance_empty(self):
        with tempfile.TemporaryDirectory() as d:
            result = build_provenance(Path(d))
            self.assertEqual(result["file_count"], 0)
            self.assertEqual(result["checksums"], {})

    def test_build_and_verify_provenance(self):
        with tempfile.TemporaryDirectory() as d:
            raw = Path(d)
            (raw / "a.json").write_text('{"x":1}')
            (raw / "b.json").write_text('{"y":2}')
            prov = build_provenance(raw)
            self.assertEqual(prov["file_count"], 2)
            failures = verify_provenance(prov, raw)
            self.assertEqual(failures, [])

    def test_verify_detects_modification(self):
        with tempfile.TemporaryDirectory() as d:
            raw = Path(d)
            (raw / "a.json").write_text('{"x":1}')
            prov = build_provenance(raw)
            (raw / "a.json").write_text('{"x":2}')  # Modify
            failures = verify_provenance(prov, raw)
            self.assertTrue(any("Checksum mismatch" in f for f in failures))

    def test_verify_detects_missing_file(self):
        with tempfile.TemporaryDirectory() as d:
            raw = Path(d)
            (raw / "a.json").write_text('{"x":1}')
            prov = build_provenance(raw)
            (raw / "a.json").unlink()  # Remove
            failures = verify_provenance(prov, raw)
            self.assertTrue(any("Missing file" in f for f in failures))


if __name__ == "__main__":
    unittest.main()
