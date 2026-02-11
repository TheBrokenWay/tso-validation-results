"""
PREDATOR X — Enterprise Document Generator

Master script that generates all 5 enterprise validation documents.

Usage:
  python -m PX_Validation.enterprise_docs.document_generator
  python PX_Validation/enterprise_docs/document_generator.py
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

OUTPUT_DIR = REPO_ROOT / "PX_Validation" / "enterprise_docs" / "output"


def main() -> int:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("PREDATOR X — ENTERPRISE DOCUMENT GENERATOR")
    print("=" * 70)
    print(f"Output: {OUTPUT_DIR}\n")

    generators = [
        ("SYSTEM_VALIDATION_CERTIFICATE", "generate_certificate"),
        ("PENETRATION_TEST_REPORT", "generate_pentest_report"),
        ("CONSTITUTIONAL_COMPLIANCE_MATRIX", "generate_compliance_matrix"),
        ("AUDIT_TRAIL_INTEGRITY_REPORT", "generate_audit_report"),
        ("PERFORMANCE_BENCHMARK_REPORT", "generate_benchmark_report"),
    ]

    results = []
    total_t0 = time.perf_counter()

    for doc_name, module_name in generators:
        print(f"  Generating {doc_name}...", end=" ", flush=True)
        t0 = time.perf_counter()
        try:
            mod = __import__(
                f"PX_Validation.enterprise_docs.{module_name}",
                fromlist=["generate"],
            )
            out_path = mod.generate(OUTPUT_DIR)
            elapsed = time.perf_counter() - t0
            print(f"OK ({elapsed:.1f}s) -> {out_path.name}")
            results.append((doc_name, True, out_path))
        except Exception as e:
            elapsed = time.perf_counter() - t0
            print(f"FAIL ({elapsed:.1f}s): {e}")
            results.append((doc_name, False, None))

    total_elapsed = time.perf_counter() - total_t0

    print(f"\n{'=' * 70}")
    passed = sum(1 for _, ok, _ in results if ok)
    print(f"RESULT: {passed}/{len(results)} documents generated in {total_elapsed:.1f}s")

    if passed == len(results):
        print(f"\nAll enterprise documents written to:")
        for _, _, path in results:
            if path:
                print(f"  {path}")
    else:
        print("\nFailed documents:")
        for name, ok, _ in results:
            if not ok:
                print(f"  FAIL: {name}")

    print("=" * 70)
    return 0 if passed == len(results) else 1


if __name__ == "__main__":
    sys.exit(main())
