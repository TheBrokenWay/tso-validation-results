"""
QUINT Internal Language — Comprehensive Test Suite

Tests:
  1. Kernel: QFrame creation, seal, verify, round-trip serialization, overhead
  2. Compiler: type detection, alias normalization, field blocking, validation
  3. Converter: ingest (dict/str/bytes), emit, batch, rejection, stats
  4. Runtime: validate, transform, merge, extract, annotate, chain
  5. Registry: schema lookup, payload check, freeze
  6. Resonance: standard vectors, aggregate, overhead target
  7. Integration: full pipeline molecule → QFrame → engine results → dossier

Acceptance criteria:
  - Internal overhead < 500 bytes per QFrame envelope
  - Round-trip resonance >= 99% on standard vectors
  - Per-frame creation < 100 us
  - Zero regressions on existing 175 tests
"""

import json
import os
import sys
import time
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from PX_System.foundation.quint.kernel import (
    QFrame, QType, SCHEMA_VERSION, create_qframe, _hash_payload,
)
from PX_System.foundation.quint.compiler import (
    detect_qtype, normalize, compile_batch, CompilationError, REQUIRED_FIELDS,
)
from PX_System.foundation.quint.converter import (
    ingest, emit, emit_json, ingest_batch, emit_batch,
    round_trip, get_stats, reset_stats, ConversionError,
)
from PX_System.foundation.quint.runtime import (
    validate, transform, merge, extract, annotate, chain,
)
from PX_System.foundation.quint.runtime import RuntimeError as QRuntimeError
from PX_System.foundation.quint.registry import (
    get_schema, registered_types, check_payload,
)
from PX_System.foundation.quint.resonance import (
    measure, measure_aggregate, run_standard_resonance_test,
    STANDARD_VECTORS, ResonanceResult,
)


class TestKernel(unittest.TestCase):
    """QFrame core type tests."""

    def test_create_and_seal(self):
        f = create_qframe(QType.QMOLECULE, "K-001", {"smiles": "CCO", "mw": 46.07})
        self.assertEqual(f.qtype, QType.QMOLECULE)
        self.assertEqual(f.qid, "K-001")
        self.assertNotEqual(f.frame_hash, "")
        self.assertTrue(f.verify())

    def test_tamper_detection(self):
        f = create_qframe(QType.QMOLECULE, "K-002", {"smiles": "CCO"})
        self.assertTrue(f.verify())
        f.payload["smiles"] = "TAMPERED"
        self.assertFalse(f.verify())

    def test_round_trip_serialization(self):
        f = create_qframe(QType.QRESULT, "K-003", {"engine_id": "OPE", "status": "OK"})
        d = f.to_dict()
        f2 = QFrame.from_dict(d)
        self.assertEqual(f2.qtype, f.qtype)
        self.assertEqual(f2.qid, f.qid)
        self.assertEqual(f2.payload, f.payload)
        self.assertEqual(f2.frame_hash, f.frame_hash)
        self.assertTrue(f2.verify())

    def test_overhead_under_500_bytes(self):
        f = create_qframe(QType.QMOLECULE, "K-004", {"smiles": "CCO"})
        self.assertLess(f.size_bytes(), 500)

    def test_all_qtypes(self):
        for qt in QType:
            f = create_qframe(qt, f"K-{qt.value}", {"test": True})
            self.assertTrue(f.verify())

    def test_deterministic_hash(self):
        f1 = create_qframe(QType.QMOLECULE, "K-005", {"smiles": "CCO", "mw": 46.07})
        f2 = create_qframe(QType.QMOLECULE, "K-005", {"smiles": "CCO", "mw": 46.07})
        self.assertEqual(f1.frame_hash, f2.frame_hash)


class TestCompiler(unittest.TestCase):
    """Compiler normalization tests."""

    def test_type_detection_molecule(self):
        self.assertEqual(detect_qtype({"smiles": "CCO"}), QType.QMOLECULE)

    def test_type_detection_result(self):
        self.assertEqual(detect_qtype({"engine_id": "OPE", "status": "OK"}), QType.QRESULT)

    def test_type_detection_constraint(self):
        self.assertEqual(detect_qtype({"disease_name": "Nipah"}), QType.QCONSTRAINT)

    def test_type_detection_raw(self):
        self.assertEqual(detect_qtype({"random_key": 1}), QType.QRAW)

    def test_alias_normalization(self):
        qt, payload = normalize({"SMILES": "CCO", "molecular_weight": 46.07})
        self.assertEqual(payload["smiles"], "CCO")
        self.assertEqual(payload["mw"], 46.07)

    def test_blocked_field_stripping(self):
        qt, payload = normalize({"smiles": "CCO", "password": "x", "api_key": "y"})
        self.assertNotIn("password", payload)
        self.assertNotIn("api_key", payload)
        self.assertIn("smiles", payload)

    def test_missing_required_fields(self):
        with self.assertRaises(CompilationError):
            normalize({}, QType.QMOLECULE)

    def test_batch_compilation(self):
        results = compile_batch([{"smiles": "CCO"}, {"smiles": "CC"}])
        self.assertEqual(len(results), 2)
        self.assertEqual(results[0][0], QType.QMOLECULE)


class TestConverter(unittest.TestCase):
    """Single-gateway converter tests."""

    def setUp(self):
        reset_stats()

    def test_ingest_dict(self):
        f = ingest({"smiles": "CCO", "mw": 46.07})
        self.assertIsInstance(f, QFrame)
        self.assertEqual(f.qtype, QType.QMOLECULE)
        self.assertTrue(f.verify())

    def test_ingest_json_string(self):
        f = ingest('{"smiles": "CC", "mw": 30.07}')
        self.assertEqual(f.payload["smiles"], "CC")

    def test_ingest_bytes(self):
        f = ingest(b'{"engine_id": "OPE", "status": "PASSED"}')
        self.assertEqual(f.qtype, QType.QRESULT)

    def test_emit_payload_only(self):
        f = ingest({"smiles": "CCO", "mw": 46.07})
        out = emit(f)
        self.assertEqual(out, {"smiles": "CCO", "mw": 46.07})

    def test_emit_with_meta(self):
        f = ingest({"smiles": "CCO"})
        out = emit(f, include_meta=True)
        self.assertIn("_qtype", out)
        self.assertIn("_frame_hash", out)

    def test_round_trip_preserves_data(self):
        original = {"smiles": "c1ccccc1", "mw": 78.11, "logp": 1.56}
        result = round_trip(original)
        self.assertEqual(result, original)

    def test_reject_bad_json(self):
        with self.assertRaises(ConversionError):
            ingest("not json")

    def test_reject_bad_type(self):
        with self.assertRaises(ConversionError):
            ingest(12345)

    def test_batch_ingest_emit(self):
        frames = ingest_batch([{"smiles": "CCO"}, {"smiles": "CC"}])
        self.assertEqual(len(frames), 2)
        outs = emit_batch(frames)
        self.assertEqual(outs[0]["smiles"], "CCO")

    def test_stats(self):
        ingest({"smiles": "CCO"})
        f = ingest({"smiles": "CC"})
        emit(f)
        stats = get_stats()
        self.assertEqual(stats["ingested"], 2)
        self.assertEqual(stats["emitted"], 1)


class TestRuntime(unittest.TestCase):
    """QFrame operation tests."""

    def test_validate_good(self):
        f = create_qframe(QType.QMOLECULE, "R-001", {"smiles": "CCO"})
        self.assertTrue(validate(f))

    def test_validate_tampered(self):
        f = create_qframe(QType.QMOLECULE, "R-002", {"smiles": "CCO"})
        f.payload["smiles"] = "BAD"
        with self.assertRaises(QRuntimeError):
            validate(f)

    def test_transform(self):
        f = create_qframe(QType.QMOLECULE, "R-003", {"smiles": "CCO", "mw": 46.07})
        f2 = transform(f, lambda p: {**p, "logp": 1.0}, label="add_logp")
        self.assertEqual(f2.payload["logp"], 1.0)
        self.assertNotIn("logp", f.payload)
        self.assertTrue(f2.verify())

    def test_merge(self):
        fa = create_qframe(QType.QRESULT, "R-A", {"engine_id": "OPE", "status": "OK", "s": 0.9})
        fb = create_qframe(QType.QRESULT, "R-B", {"engine_id": "OBE", "status": "OK", "b": 0.8})
        m = merge([fa, fb], qtype=QType.QDOSSIER, qid="R-MERGE")
        self.assertEqual(m.payload["s"], 0.9)
        self.assertEqual(m.payload["b"], 0.8)
        self.assertTrue(m.verify())

    def test_extract(self):
        f = create_qframe(QType.QMOLECULE, "R-005", {"smiles": "CCO", "mw": 46.07, "logp": 1.0})
        e = extract(f, {"smiles", "mw"})
        self.assertIn("smiles", e.payload)
        self.assertNotIn("logp", e.payload)

    def test_annotate(self):
        f = create_qframe(QType.QMOLECULE, "R-006", {"smiles": "CCO"})
        a = annotate(f, {"tier": "GOLD"})
        self.assertEqual(a.payload["tier"], "GOLD")
        self.assertNotIn("tier", f.payload)

    def test_chain(self):
        f = create_qframe(QType.QMOLECULE, "R-007", {"smiles": "CCO", "mw": 46.07})
        result = chain(f, [
            lambda x: transform(x, lambda p: {**p, "logp": 1.0}),
            lambda x: annotate(x, {"processed": True}),
        ])
        self.assertEqual(result.payload["logp"], 1.0)
        self.assertTrue(result.payload["processed"])


class TestRegistry(unittest.TestCase):
    """Schema registry tests."""

    def test_all_types_registered(self):
        types = registered_types()
        self.assertEqual(len(types), len(QType))

    def test_schema_lookup(self):
        s = get_schema(QType.QMOLECULE)
        self.assertIsNotNone(s)
        self.assertIn("smiles", s.required_fields)

    def test_payload_check_missing(self):
        missing = check_payload(QType.QMOLECULE, {"mw": 46.07})
        self.assertIn("smiles", missing)

    def test_payload_check_ok(self):
        missing = check_payload(QType.QMOLECULE, {"smiles": "CCO"})
        self.assertEqual(len(missing), 0)


class TestResonance(unittest.TestCase):
    """Resonance measurement tests."""

    def test_single_measurement(self):
        r = measure({"smiles": "CCO", "mw": 46.07})
        self.assertEqual(r.resonance_pct, 100.0)
        self.assertTrue(r.passed)

    def test_standard_vectors(self):
        agg = run_standard_resonance_test()
        self.assertTrue(agg["all_passed"],
                        f"Standard resonance failed: min={agg['min_resonance']:.2f}%")
        self.assertGreaterEqual(agg["mean_resonance"], 99.0)

    def test_large_payload(self):
        big = {"smiles": "C" * 200}
        for i in range(100):
            big[f"field_{i}"] = i * 0.1
        r = measure(big)
        self.assertGreaterEqual(r.resonance_pct, 99.0)

    def test_nested_payload(self):
        data = {"smiles": "CCO", "meta": {"source": "chembl", "version": 2}}
        r = measure(data)
        self.assertGreaterEqual(r.resonance_pct, 99.0)


class TestPerformance(unittest.TestCase):
    """Overhead and throughput tests."""

    def test_qframe_overhead(self):
        f = create_qframe(QType.QMOLECULE, "P-001", {"smiles": "CCO"})
        self.assertLess(f.size_bytes(), 500, "QFrame envelope overhead exceeds 500 bytes")

    def test_creation_throughput(self):
        t0 = time.perf_counter()
        for i in range(5000):
            create_qframe(QType.QRESULT, f"P-{i}", {"engine_id": "OPE", "status": "OK"})
        elapsed = time.perf_counter() - t0
        per_us = (elapsed / 5000) * 1_000_000
        self.assertLess(per_us, 100, f"QFrame creation too slow: {per_us:.1f} us/frame")

    def test_round_trip_throughput(self):
        reset_stats()
        t0 = time.perf_counter()
        for i in range(5000):
            f = ingest({"smiles": f"C{i}", "mw": 100 + i})
            emit(f)
        elapsed = time.perf_counter() - t0
        per_us = (elapsed / 5000) * 1_000_000
        self.assertLess(per_us, 200, f"Round-trip too slow: {per_us:.1f} us")


class TestIntegration(unittest.TestCase):
    """Full pipeline integration: molecule → engine results → dossier."""

    def test_molecule_pipeline(self):
        reset_stats()

        # Step 1: External molecule enters through the gateway
        mol_frame = ingest(
            {"smiles": "CC(=O)Oc1ccccc1C(=O)O", "mw": 180.16, "logp": 1.19},
            source_label="external_chembl",
        )
        self.assertEqual(mol_frame.qtype, QType.QMOLECULE)
        self.assertTrue(mol_frame.verify())

        # Step 2: Engine processes the molecule (simulate OPE result)
        ope_result = transform(
            mol_frame,
            lambda p: {**p, "ec50": 0.5, "emax": 0.95, "binding_affinity": 8.2},
            label="OPE_V3",
            new_qtype=QType.QRESULT,
        )

        # Step 3: Another engine adds ADMET data
        admet_result = transform(
            ope_result,
            lambda p: {**p, "toxicity": 0.015, "safety_margin": 65.0, "clearance": 5.2},
            label="ADMET_V3",
        )

        # Step 4: Merge into dossier
        dossier = annotate(
            admet_result,
            {
                "compound_id": "PRV_NOV_aspirin_analog",
                "tier": "GOLD",
                "zeus_verdict": {"authorized": True, "laws_passed": ["L1", "U27", "U34", "L11"]},
            },
            label="finalization",
        )

        # Step 5: Emit through the gateway
        output = emit(dossier, include_meta=True, target_label="warehouse")

        # Verify everything made it through
        self.assertEqual(output["smiles"], "CC(=O)Oc1ccccc1C(=O)O")
        self.assertEqual(output["mw"], 180.16)
        self.assertEqual(output["ec50"], 0.5)
        self.assertEqual(output["toxicity"], 0.015)
        self.assertEqual(output["tier"], "GOLD")
        self.assertTrue(output["zeus_verdict"]["authorized"])
        self.assertIn("_frame_hash", output)

        # Verify lineage tracked all steps
        self.assertIn("OPE_V3", dossier.lineage)
        self.assertIn("ADMET_V3", dossier.lineage)
        self.assertIn("finalization", dossier.lineage)

        # Verify stats
        stats = get_stats()
        self.assertEqual(stats["ingested"], 1)
        self.assertEqual(stats["emitted"], 1)


if __name__ == "__main__":
    unittest.main(verbosity=2)
