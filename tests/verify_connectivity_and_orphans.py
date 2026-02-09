#!/usr/bin/env python3
"""
PREDATOR X — Verify every engine/orchestrator is connected (no zombies/orphans).
Run from repo root. Output: PX_LOGS/CONNECTIVITY_REPORT_<timestamp>.md
"""
from __future__ import annotations

import ast
import os
import sys
from pathlib import Path
from datetime import datetime, timezone
from collections import defaultdict

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

LOG_DIR = REPO / "PX_LOGS"
LOG_DIR.mkdir(parents=True, exist_ok=True)
STAMP = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
REPORT_PATH = LOG_DIR / f"CONNECTIVITY_REPORT_{STAMP}.md"

# Canonical engines per LEAN_REPO / ENGINE_INTEGRATION_SPEC
ENGINES_PX_ENGINE = [
    "Genesis_Engine", "Trajectory_Predictor", "Metabolism", "Vector_Core",
    "Block_Orchestrator", "Engine_Orchestrator", "Stress_Test",
]
ENGINES_OPERATIONS = [
    "OBE", "OCE", "OLE", "OME", "OPE", "OSE", "ADMET", "TrialEngine", "PKPD",
    "DoseOptimizer_v2", "VirtualEfficacyAnalytics", "GradingEngine",
]
ENGINES_LAB = ["Simulation_Engine", "Manufacturing_Manifest"]

ORCHESTRATORS_DIRS = [REPO / "PX_Executive", REPO / "run_e2e_layers.py"]
TEST_DIRS = [REPO / "PX_Validation" / "tests", REPO / "tests"]
AUDIT_OR_OTHER = [REPO / "PX_Audit", REPO / "Nipah_Analysis", REPO / "PX_System"]


def grep_imports(root: Path, pattern: str, glob_pat: str = "*.py") -> set[str]:
    found = set()
    for py in root.rglob(glob_pat):
        if "__pycache__" in str(py) or "99_LEGACY_CODE" in str(py):
            continue
        try:
            text = py.read_text(encoding="utf-8", errors="replace")
            if pattern in text:
                found.add(str(py.relative_to(REPO)))
        except Exception:
            pass
    return found


def main() -> int:
    lines = [
        "# PREDATOR X — Connectivity & Orphan Report",
        f"Generated: {datetime.now(timezone.utc).isoformat()}",
        "",
        "## 1. Engine connectivity (Section 5)",
        "",
    ]

    all_ok = True
    for name in ENGINES_PX_ENGINE + ENGINES_OPERATIONS:
        mod = name if name in ENGINES_PX_ENGINE else f"PX_Engine.operations.{name}" if name != "ADMET" else "PX_Engine.operations.ADMET"
        needle = name.replace("_", ".") if name in ("Genesis_Engine", "Trajectory_Predictor", "Vector_Core", "Block_Orchestrator", "Engine_Orchestrator", "Stress_Test", "Metabolism") else name
        refs = set()
        for d in ORCHESTRATORS_DIRS + TEST_DIRS + AUDIT_OR_OTHER:
            if d.is_file():
                refs |= grep_imports(d.parent, needle, d.name)
            elif d.exists():
                refs |= grep_imports(d, needle)
        if not refs:
            # Try generic PX_Engine
            refs = grep_imports(REPO, f"PX_Engine.{needle}") or grep_imports(REPO, needle)
        status = "OK" if refs else "ORPHAN"
        if not refs:
            all_ok = False
        lines.append(f"- **{name}**: {status}" + (f" (refs: {len(refs)})" if refs else ""))

    for name in ENGINES_LAB:
        refs = grep_imports(REPO, name)
        status = "OK" if refs else "ORPHAN"
        if not refs:
            all_ok = False
        lines.append(f"- **{name}**: {status}" + (f" (refs: {len(refs)})" if refs else ""))

    lines.extend(["", "## 2. Root scripts (have __main__ or are imported)", ""])
    for root_py in ["run_e2e_layers.py"]:
        p = REPO / root_py
        if not p.exists():
            continue
        refs = grep_imports(REPO, root_py.replace(".py", ""))
        has_main = "__main__" in (p.read_text(encoding="utf-8", errors="replace"))
        lines.append(f"- **{root_py}**: " + ("has __main__" if has_main else "no __main__") + f", refs: {len(refs)}")

    lines.extend(["", "## 3. Result", ""])
    lines.append("**ALL CONNECTED**" if all_ok else "**ORPHANS DETECTED — do not archive further.**")
    report = "\n".join(lines)
    REPORT_PATH.write_text(report, encoding="utf-8")
    print(report)
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
