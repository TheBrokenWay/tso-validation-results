#!/usr/bin/env python3
"""
Verify every engine, orchestrator, governance, and research .py file in the repo.
For each file: resolve module name, import it, record pass/fail.
Excludes: 99_LEGACY_CODE, __pycache__, virtual envs.
Run from repo root: python tests/verify_every_file.py
"""
from __future__ import annotations

import os
import sys
import importlib
import traceback
from pathlib import Path
from datetime import datetime, timezone

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

SKIP_DIRS = {"99_LEGACY_CODE", "archive", "__pycache__", ".conda", "venv", "env", ".git", "node_modules"}
SKIP_FILES = {"verify_every_file.py"}  # self

RESULTS: list[tuple[str, bool, str]] = []
FAILED: list[tuple[str, str]] = []


def path_to_module(rel_path: Path) -> str | None:
    """Convert path like PX_Engine/operations/OPE.py to PX_Engine.operations.OPE."""
    parts = rel_path.parts
    if not parts or rel_path.suffix != ".py":
        return None
    name = rel_path.stem
    if name.startswith("_"):
        return None  # __init__ etc. are imported via package
    return ".".join(parts[:-1] + (name,))


def collect_py_files() -> list[Path]:
    out: list[Path] = []
    for root, dirs, files in os.walk(REPO_ROOT):
        dirs[:] = [d for d in dirs if d not in SKIP_DIRS]
        root_path = Path(root)
        try:
            rel = root_path.relative_to(REPO_ROOT)
        except ValueError:
            continue
        if any(p in rel.parts for p in SKIP_DIRS):
            continue
        for f in files:
            if f.endswith(".py") and f not in SKIP_FILES:
                out.append(rel / f)
    return sorted(out)


def verify_one(rel_path: Path) -> tuple[bool, str]:
    module_name = path_to_module(rel_path)
    if not module_name:
        return True, "(skip init)"
    try:
        importlib.import_module(module_name)
        return True, ""
    except Exception as e:
        return False, traceback.format_exc()


def main() -> int:
    print("=" * 70)
    print("PREDATOR X — VERIFY EVERY FILE (import each module)")
    print("=" * 70)
    py_files = collect_py_files()
    print(f"Found {len(py_files)} .py files (excluding 99_LEGACY_CODE, __pycache__)\n")

    for i, rel_path in enumerate(py_files):
        mod_name = path_to_module(rel_path) or str(rel_path)
        ok, err = verify_one(rel_path)
        path_str = str(rel_path).replace("\\", "/")
        RESULTS.append((path_str, ok, err))
        if ok:
            print(f"  [PASS] {path_str}")
        else:
            print(f"  [FAIL] {path_str}")
            FAILED.append((path_str, err))

    passed = sum(1 for _, ok, _ in RESULTS if ok)
    total = len(RESULTS)
    print("\n" + "=" * 70)
    print(f"RESULT: {passed}/{total} files verified (import OK)")
    if FAILED:
        print(f"\n{len(FAILED)} FAILED (first 5 errors below):")
        for path_str, err in FAILED[:5]:
            print(f"\n--- {path_str} ---")
            print(err[:800])
        out_path = REPO_ROOT / "PX_LOGS" / "verify_every_file_failures.txt"
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with open(out_path, "w", encoding="utf-8") as f:
            for path_str, err in FAILED:
                f.write(f"\n{'='*60}\n{path_str}\n{'='*60}\n{err}\n")
        print(f"\nFull failures written to: {out_path}")
    print("=" * 70)
    return 0 if not FAILED else 1


if __name__ == "__main__":
    sys.exit(main())
