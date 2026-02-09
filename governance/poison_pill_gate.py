"""
Poison-pill execution gate: enforces run order at the runner.

Before any of:
  - PX_System initialization
  - PX_Warehouse writes
  - PX_Engine invocation

the runner MUST run the poison pill test and set poison_pill_test_ran = True.
This makes the poison pill a true entry gate, not just an early test.
"""
import os
import sys
import subprocess
from pathlib import Path

# Repo root (parent of governance/)
_GATE_DIR = Path(__file__).resolve().parent
_ROOT_DIR = _GATE_DIR.parent

# Set by run_poison_pill_gate() when the test passes; asserted before pipeline work
poison_pill_test_ran = False


def run_poison_pill_gate(timeout_seconds: int = 60) -> bool:
    """
    Run all mandatory pre-execution poison pill tests as subprocesses.
    On success for all, set poison_pill_test_ran = True.
    Call this at the very start of any main pipeline runner.
    """
    global poison_pill_test_ran
    try:
        from PX_Constitution import mandatory_pre_execution_tests
    except ImportError:
        mandatory_pre_execution_tests = ["run_poison_pill_test.py"]
    for test_script_name in mandatory_pre_execution_tests:
        test_script = _ROOT_DIR / test_script_name
        if not test_script.exists():
            sys.stderr.write(
                f"CRITICAL: Poison pill test not found at {test_script}. "
                "Constitutional gate cannot run.\n"
            )
            return False
        result = subprocess.run(
            [sys.executable, str(test_script)],
            cwd=str(_ROOT_DIR),
            timeout=timeout_seconds,
            capture_output=False,
        )
        if result.returncode != 0:
            return False
    poison_pill_test_ran = True
    return True


def require_poison_pill_gate_before_pipeline() -> None:
    """
    Assert that the poison pill gate has been run and passed.
    Call this immediately after run_poison_pill_gate() and before:
      - PX_System initialization
      - PX_Warehouse writes
      - PX_Engine invocation
    """
    assert poison_pill_test_ran is True, (
        "Poison pill gate must run before PX_System init, PX_Warehouse writes, or PX_Engine invocation. "
        "Call run_poison_pill_gate() at pipeline entry and assert this before proceeding."
    )
