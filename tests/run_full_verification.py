#!/usr/bin/env python3
r"""
Full platform verification: every layer, critical files, E2E.
Run from repo root with conda env: conda run -p E:\foundation\.conda python tests/run_full_verification.py
Outputs: PX_LOGS/VERIFICATION_MANIFEST_<timestamp>.md
"""
from __future__ import annotations

import os
import sys
import subprocess
from pathlib import Path
from datetime import datetime, timezone

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

LOG_DIR = REPO_ROOT / "PX_LOGS"
LOG_DIR.mkdir(parents=True, exist_ok=True)
STAMP = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
MANIFEST_PATH = LOG_DIR / f"VERIFICATION_MANIFEST_{STAMP}.md"
RESULTS: list[tuple[str, bool, str]] = []


def run(name: str, ok: bool, detail: str = ""):
    RESULTS.append((name, ok, detail))
    sym = "PASS" if ok else "FAIL"
    print(f"  [{sym}] {name}" + (f" — {detail}" if detail else ""))


def main() -> int:
    print("=" * 60)
    print("PREDATOR X — FULL PLATFORM VERIFICATION")
    print("=" * 60)

    # --- 1. CRITICAL IMPORTS ---
    print("\n1. Critical layer imports")
    r0 = subprocess.run(
        [sys.executable, str(REPO_ROOT / "tests" / "verify_all_imports.py")],
        cwd=str(REPO_ROOT), capture_output=True, text=True, timeout=30
    )
    run("verify_all_imports.py (OPE, ADMET, GradingEngine, Finalization, Zeus, Legal, Patent)", r0.returncode == 0, r0.stderr[:150] if r0.returncode else "")

    # --- 2. SPEC & PIPELINE TESTS ---
    print("\n2. Finalization spec tests")
    r = subprocess.run(
        [sys.executable, str(REPO_ROOT / "tests" / "test_finalization_spec.py")],
        cwd=str(REPO_ROOT), capture_output=True, text=True, timeout=90
    )
    run("test_finalization_spec.py", r.returncode == 0, r.stderr[:200] if r.returncode else "")

    print("\n3. Orchestrator/warehouse path tests")
    r2 = subprocess.run(
        [sys.executable, str(REPO_ROOT / "tests" / "test_orchestrator_warehouse_paths.py")],
        cwd=str(REPO_ROOT), capture_output=True, text=True, timeout=60
    )
    run("test_orchestrator_warehouse_paths.py", r2.returncode == 0, r2.stderr[:200] if r2.returncode else "")

    # --- 4. VALIDATION SUITE (key tests) ---
    print("\n4. GradingEngine tests")
    r3 = subprocess.run(
        [sys.executable, str(REPO_ROOT / "PX_Validation" / "tests" / "test_grading_engine.py")],
        cwd=str(REPO_ROOT), capture_output=True, text=True, timeout=30
    )
    run("test_grading_engine.py", r3.returncode == 0)

    print("5. ADMET engine tests")
    r4 = subprocess.run(
        [sys.executable, str(REPO_ROOT / "PX_Validation" / "tests" / "test_admet_engine.py")],
        cwd=str(REPO_ROOT), capture_output=True, text=True, timeout=30
    )
    run("test_admet_engine.py", r4.returncode == 0)

    # --- 6. MONITOR (import + one count) ---
    print("\n6. Monitor warehouse (import + count_files)")
    try:
        from PX_Executive.monitor_warehouse import count_files, PATHS
        n = count_files(PATHS["Novel"])
        r_count = count_files(PATHS["Repurposed"])
        run("monitor_warehouse count_files", True, f"N={n} R={r_count}")
    except Exception as e:
        run("monitor_warehouse", False, str(e))

    # --- 7. FTO / Legal ---
    print("\n7. FTO diagnosis (script runs)")
    r5 = subprocess.run(
        [sys.executable, str(REPO_ROOT / "PX_Executive" / "diagnose_fto_failure_rate.py")],
        cwd=str(REPO_ROOT), capture_output=True, text=True, timeout=30
    )
    run("diagnose_fto_failure_rate.py", r5.returncode == 0)

    # --- 8. One-cycle E2E (Feeder → Orchestrator → Finalization) ---
    print("\n8. One-cycle E2E (Feeder, Engine, Finalization)")
    r6 = subprocess.run(
        [sys.executable, str(REPO_ROOT / "PX_Executive" / "run_one_cycle_test.py")],
        cwd=str(REPO_ROOT), capture_output=True, text=True, timeout=300
    )
    run("run_one_cycle_test.py (full cycle)", r6.returncode == 0, r6.stderr[:300] if r6.returncode else "")

    # --- 9. Orchestrator 2-item E2E ---
    print("\n9. Orchestrator 2-item E2E (IN→PP→E2E→LG→OK)")
    env = os.environ.copy()
    env["PRV_MAX_ITEMS"] = "2"
    env["PRV_LIVE_RESEARCH"] = "1"
    env["PRV_API_PACING_SEC"] = "1"
    if str(REPO_ROOT) not in env.get("PYTHONPATH", ""):
        env["PYTHONPATH"] = str(REPO_ROOT) + (os.pathsep + env.get("PYTHONPATH", "") if env.get("PYTHONPATH") else "")
    r7 = subprocess.run(
        [sys.executable, str(REPO_ROOT / "PX_Executive" / "PRV_24H_Orchestrator.py")],
        cwd=str(REPO_ROOT), env=env, capture_output=True, text=True, timeout=120
    )
    out = (r7.stdout or "") + (r7.stderr or "")
    done_ok = "DONE=2" in out or "DONE=1" in out
    run("PRV_24H_Orchestrator 2 items", r7.returncode == 0 and ("OK" in out or done_ok), "exit=%s" % r7.returncode)

    # --- WRITE MANIFEST ---
    passed = sum(1 for _, ok, _ in RESULTS if ok)
    total = len(RESULTS)
    with open(MANIFEST_PATH, "w", encoding="utf-8") as f:
        f.write("# PREDATOR X — Verification Manifest\n\n")
        f.write(f"**Generated:** {datetime.now(timezone.utc).isoformat()}\n\n")
        f.write(f"**Result:** {passed}/{total} checks passed\n\n")
        f.write("## Per-check status\n\n")
        for name, ok, detail in RESULTS:
            status = "PASS" if ok else "FAIL"
            f.write(f"- **{status}** `{name}`")
            if detail:
                f.write(f" — {detail[:200]}")
            f.write("\n")
        f.write("\n## Files / layers verified\n\n")
        f.write("- `tests/verify_all_imports.py` — Engine, Warehouse, Zeus, Legal, Patent_Local_Index\n")
        f.write("- `tests/test_finalization_spec.py` — Finalization_Spec, Finalization_Pipeline, Zeus gate\n")
        f.write("- `tests/test_orchestrator_warehouse_paths.py` — Paths, Evidence_Package, UniversalPipelineRunner\n")
        f.write("- `PX_Validation/tests/test_grading_engine.py` — GradingEngine\n")
        f.write("- `PX_Validation/tests/test_admet_engine.py` — ADMET\n")
        f.write("- `PX_Executive/monitor_warehouse.py` — count_files, PATHS\n")
        f.write("- `PX_Executive/diagnose_fto_failure_rate.py` — FTO diagnosis\n")
        f.write("- `PX_Executive/run_one_cycle_test.py` — Feeder → Orchestrator → Finalization\n")
        f.write("- `PX_Executive/PRV_24H_Orchestrator.py` — 2-item E2E (IN→PP→E2E→LG→OK)\n")

    print("\n" + "=" * 60)
    print(f"VERIFICATION: {passed}/{total} passed. Manifest: {MANIFEST_PATH}")
    print("=" * 60)
    return 0 if passed == total else 1


if __name__ == "__main__":
    sys.exit(main())
