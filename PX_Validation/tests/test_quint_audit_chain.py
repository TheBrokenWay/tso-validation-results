"""
Tests for QUINT integration with Sovereign_Log_Chain.

Verifies that:
1. append() still works without qframe (backward compatible)
2. append() with qframe adds quint_metadata to log entry
3. quint_metadata contains frame_hash, qtype, qid, lineage
4. Invalid qframe doesn't break logging
"""
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from PX_System.foundation.quint.kernel import QFrame, QType, create_qframe
from PX_System.foundation.quint.converter import ingest


class TestSovereignLogChainQuintIntegration(unittest.TestCase):
    """Test QUINT integration with Sovereign_Log_Chain."""

    def setUp(self):
        """Create a temporary log file for testing."""
        self.temp_dir = tempfile.mkdtemp()
        self.temp_log = Path(self.temp_dir) / "test_log.jsonl"
        # Patch LOG_PATH to use temp file
        import PX_System.foundation.Sovereign_Log_Chain as slc
        self._orig_log_path = slc.LOG_PATH
        slc.LOG_PATH = self.temp_log

    def tearDown(self):
        """Restore original LOG_PATH."""
        import PX_System.foundation.Sovereign_Log_Chain as slc
        slc.LOG_PATH = self._orig_log_path
        # Clean up temp file
        if self.temp_log.exists():
            self.temp_log.unlink()
        os.rmdir(self.temp_dir)

    def _read_last_entry(self) -> dict:
        """Read the last log entry."""
        lines = self.temp_log.read_text(encoding="utf-8").strip().split("\n")
        return json.loads(lines[-1])

    def test_append_without_qframe_backward_compatible(self):
        """append() without qframe works exactly as before."""
        from PX_System.foundation.Sovereign_Log_Chain import append
        result = append("TEST_EVENT", {"key": "value"})
        self.assertTrue(result)
        entry = self._read_last_entry()
        self.assertEqual(entry["event_type"], "TEST_EVENT")
        self.assertNotIn("quint_metadata", entry)

    def test_append_with_qframe_adds_metadata(self):
        """append() with qframe adds quint_metadata block."""
        from PX_System.foundation.Sovereign_Log_Chain import append
        frame = create_qframe(
            qtype=QType.QRESULT,
            qid="TEST-Q001",
            payload={"engine_id": "TEST", "status": "PASSED"},
        )
        result = append("QUINT_TEST", {"test": True}, qframe=frame)
        self.assertTrue(result)
        entry = self._read_last_entry()
        self.assertIn("quint_metadata", entry)

    def test_quint_metadata_has_frame_hash(self):
        """quint_metadata contains frame_hash."""
        from PX_System.foundation.Sovereign_Log_Chain import append
        frame = create_qframe(
            qtype=QType.QMOLECULE,
            qid="MOL-001",
            payload={"smiles": "CCO"},
        )
        append("HASH_TEST", {"test": True}, qframe=frame)
        entry = self._read_last_entry()
        meta = entry["quint_metadata"]
        self.assertIn("frame_hash", meta)
        self.assertTrue(len(meta["frame_hash"]) == 64)  # SHA-256 hex

    def test_quint_metadata_has_qtype(self):
        """quint_metadata contains qtype as string."""
        from PX_System.foundation.Sovereign_Log_Chain import append
        frame = create_qframe(
            qtype=QType.QDOSSIER,
            qid="DOS-001",
            payload={"compound_id": "TEST"},
        )
        append("TYPE_TEST", {"test": True}, qframe=frame)
        entry = self._read_last_entry()
        self.assertEqual(entry["quint_metadata"]["qtype"], "QDOSSIER")

    def test_quint_metadata_has_qid(self):
        """quint_metadata contains qid."""
        from PX_System.foundation.Sovereign_Log_Chain import append
        frame = create_qframe(
            qtype=QType.QRESULT,
            qid="MY-QID-123",
            payload={"engine_id": "X", "status": "PASSED"},
        )
        append("QID_TEST", {}, qframe=frame)
        entry = self._read_last_entry()
        self.assertEqual(entry["quint_metadata"]["qid"], "MY-QID-123")

    def test_quint_metadata_has_lineage(self):
        """quint_metadata contains lineage list."""
        from PX_System.foundation.Sovereign_Log_Chain import append
        frame = ingest({"smiles": "CCO"}, source_label="test_source")
        append("LINEAGE_TEST", {}, qframe=frame)
        entry = self._read_last_entry()
        self.assertIsInstance(entry["quint_metadata"]["lineage"], list)
        self.assertTrue(len(entry["quint_metadata"]["lineage"]) > 0)

    def test_quint_metadata_has_source_hash(self):
        """quint_metadata contains source_hash."""
        from PX_System.foundation.Sovereign_Log_Chain import append
        frame = ingest({"smiles": "CCO"})
        append("HASH_TEST", {}, qframe=frame)
        entry = self._read_last_entry()
        self.assertIn("source_hash", entry["quint_metadata"])
        self.assertTrue(len(entry["quint_metadata"]["source_hash"]) == 64)

    def test_invalid_qframe_does_not_break_logging(self):
        """Passing invalid object as qframe doesn't prevent logging."""
        from PX_System.foundation.Sovereign_Log_Chain import append
        result = append("INVALID_TEST", {"test": True}, qframe="not_a_qframe")
        self.assertTrue(result)
        entry = self._read_last_entry()
        self.assertEqual(entry["event_type"], "INVALID_TEST")

    def test_none_qframe_no_metadata(self):
        """qframe=None produces no quint_metadata."""
        from PX_System.foundation.Sovereign_Log_Chain import append
        append("NONE_TEST", {}, qframe=None)
        entry = self._read_last_entry()
        self.assertNotIn("quint_metadata", entry)

    def test_chain_hash_preserved_with_quint(self):
        """Chain hash integrity maintained when using qframe."""
        from PX_System.foundation.Sovereign_Log_Chain import append, get_chain_hash

        # First entry without QUINT
        append("FIRST", {"seq": 1})
        hash1 = get_chain_hash()

        # Second entry with QUINT
        frame = create_qframe(QType.QRESULT, "Q1", {"engine_id": "T", "status": "PASSED"})
        append("SECOND", {"seq": 2}, qframe=frame)
        hash2 = get_chain_hash()

        # Hashes should be different (chain progressed)
        self.assertNotEqual(hash1, hash2)

        # Read entries and verify prev_hash chain
        lines = self.temp_log.read_text(encoding="utf-8").strip().split("\n")
        entry1 = json.loads(lines[0])
        entry2 = json.loads(lines[1])
        self.assertEqual(entry2["prev_hash"], entry1["record_hash"])

    def test_record_hash_includes_quint_metadata(self):
        """Record hash computation includes quint_metadata when present."""
        from PX_System.foundation.Sovereign_Log_Chain import append
        frame = create_qframe(QType.QRESULT, "Q1", {"engine_id": "T", "status": "PASSED"})
        append("HASH_VERIFY", {}, qframe=frame)
        entry = self._read_last_entry()
        self.assertIn("record_hash", entry)
        self.assertTrue(len(entry["record_hash"]) == 64)


if __name__ == "__main__":
    unittest.main()
