# PREDATOR X — Lean Repo Instruction Set (IDE Execution Version)

Four instruction groups merged into one master list. The IDE must enforce these rules when performing cleanup, archive, or structural changes.

---

## SECTION 1 — PROTECT THE WAREHOUSE (DO NOT TOUCH)

The IDE must **lock** the following directories and files. They must **never** be moved, deleted, renamed, or modified.

**Protected paths:**

- `PX_Warehouse/`
- `PX_Warehouse/WorldLines/`
- `PX_Warehouse/Prv_Dossiers/`
- `PX_Warehouse/Finalized_Dossiers/`
- `PX_Warehouse/Novel_Dossiers/`
- `PX_Warehouse/Physics_Maps/` (if present; reserved for DVC/tracked physics)
- `PX_LOGS/`
- `PX_LOGS/current/` (if present)
- `PX_Warehouse/warehouse_layout.py`

**Rule:** If any operation attempts to modify, move, delete, or rename these paths, the IDE must **abort immediately**.

---

## SECTION 2 — VALIDATION BEFORE ANY ACTION

Before archiving or moving anything, the IDE must run:

```bash
python PX_Validation/tests/run_all_tests.py
python run_e2e_layers.py
```

**Both must pass 100%.** If either fails → **STOP** and restore the repo state. This ensures the system is fully operational before cleanup.

---

## SECTION 3 — IDENTIFY AND ARCHIVE NON‑ESSENTIAL FILES

The IDE must move the following categories into `archive/<category>/`.

### 3.1 — Legacy Code (Safe to Archive)

Move everything under:

- `99_LEGACY_CODE/`
- `99_WAREHOUSE_ARCHIVE/` (if present)
- `PX_LOGS/old/` (if present)
- `PX_Executive/docs/` (entire folder: FILE_PLACEMENT_POLICY.md, Manufacturing_Manifest.py, USAGE.md, etc.)

### 3.2 — Developer‑Only Files (Safe to Archive)

Move (if present):

- `PX_Validation/manual_tests/`
- `PX_Validation/final_validation.py`
- `PX_Validation/final_validation_FINAL.py`
- `PX_SYSTEM_TEST_REPORT.json` (root; duplicate in PX_Audit/reports/ may exist)

### 3.3 — Duplicate or Deprecated Scripts

Move **if present** and **only after reference check** (see Rule below):

- `extract_reprocess_candidates.py` → under `PX_Warehouse/Operations/scripts/`
- `extract_reprocess_candidates_v2.py` → under `PX_Warehouse/Operations/scripts/`
- `reprocess_pipeline.py` → under `PX_Warehouse/Operations/scripts/`
- `reprocess_pipeline_direct.py` → under `PX_Warehouse/Operations/scripts/`
- `PX_LOGS/extract_fto_blocked_candidates.py`
- `PX_Warehouse/test_35d_ignition.py`
- `TSO_Validator/TSO_PUBLIC_RELEASE/run_poison_pill_test.py` (if present)
- `PX_Executive/docs/Manufacturing_Manifest.py`

**Rule:** Before archiving each file, the IDE must:

1. Search the entire repo for imports or references.
2. Confirm no orchestrator depends on it.
3. Confirm no engine depends on it.
4. Confirm no warehouse script depends on it.

**If referenced → do not archive.**

---

## SECTION 4 — VALIDATION AFTER EACH ARCHIVE ACTION

After each file or folder is moved to `archive/`, the IDE must re‑run:

```bash
python PX_Validation/tests/run_all_tests.py
python run_e2e_layers.py
```

If either fails:

- Restore the file.
- Mark it as **required**.
- Do not attempt to archive it again.

This guarantees zero regressions.

---

## SECTION 5 — ENFORCE ENGINE CONNECTIVITY (NO ORPHAN ENGINES)

The IDE must verify that every engine in `PX_Engine/`, `PX_Engine/operations/`, and `PX_Laboratory/` is imported by at least one orchestrator in `PX_Executive/` (or by a test in `PX_Validation/tests/` or `tests/`).

**Engines that must be connected:**

- **PX_Engine root:** Genesis_Engine, Trajectory_Predictor, Metabolism, Vector_Core, Block_Orchestrator, Engine_Orchestrator, Stress_Test
- **PX_Engine.operations:** OBE, OCE, OLE, OME, OPE, OSE, ADMET, TrialEngine, PKPD, PKPD_Simple, DoseOptimizer, DoseOptimizer_v2, DoseOptimizer_Simple, VirtualEfficacy_Simple, VirtualEfficacyAnalytics, GradingEngine, GradingSchema_Discovery.json
- **PX_Laboratory:** Simulation_Engine, Manufacturing_Manifest

If any engine is not imported by any orchestrator (or test): **IDE must flag it**, **halt cleanup**, and **not archive anything further**. This ensures 100% pipeline connectivity.

---

## SECTION 6 — ENFORCE ORCHESTRATOR COMPLETENESS

Each orchestrator must call the correct engines per stage:

| Stage | Required engines |
|-------|------------------|
| **Feed** | Genesis_Engine, Vector_Core, Metabolism, Trajectory_Predictor |
| **PRV** | OBE, OCE, OLE, OME, OPE, OSE, ADMET, PKPD or PKPD_Simple, DoseOptimizer or DoseOptimizer_v2, VirtualEfficacy_Simple, GradingEngine, GradingSchema_Discovery.json |
| **Trial** | TrialEngine, Simulation_Engine, VirtualEfficacyAnalytics, Manufacturing_Manifest |
| **Finalization** | Evidence_Package, WorldLine_Database, GradingEngine, GradingSchema_Discovery.json |

If any orchestrator is missing required calls → **IDE must halt cleanup.**

---

## SECTION 7 — ENFORCE GOVERNANCE HOOKS

Every orchestrator must import or call:

- **ZeusLaws** (e.g. `run_zeus_gate` or `check_constitutional` from `PX_System.foundation.ZeusLaws`)
- **Hydra** (configuration governance, where used)
- **run_e2e_layers.py** (governance stress test; orchestrators must be covered by it)

If any orchestrator bypasses governance → **IDE must halt cleanup.**

---

## SECTION 8 — ENFORCE DATA FLOW INTEGRITY

The IDE must verify:

- **Feed → PRV → Trial → Finalization** (data passed forward).
- No stage regenerates data that should be passed forward.
- No stage bypasses governance.
- No stage bypasses validation.

If any violation is detected → **IDE must halt cleanup.**

---

## SECTION 9 — FINAL CLEAN REPO STRUCTURE (TARGET STATE)

After cleanup, the repo should contain **only**:

**Directories:**

- `PX_Executive/`
- `PX_Engine/`
- `PX_Engine/operations/`
- `PX_System/`
- `PX_Laboratory/`
- `PX_Warehouse/` **(protected)**
- `PX_LOGS/` **(protected)**
- `governance/`
- `scripts/`
- `PX_Validation/` (required for tests)
- `PX_Warehouse/Operations/` (canonical warehouse scripts; do not archive wholesale without reference check)
- Other modules referenced by ENGINE_INTEGRATION_SPEC.md (e.g. Nipah_Analysis, PX_Constitution, tests, 99_LEGACY_CODE as archive only)

**Root files:**

- `flake.nix`, `.envrc`, `justfile`
- `README.md`, `STATUS.md`, `SYSTEM_FLOW.md`, `ARCHITECTURE_DIAGRAM.md`
- `GOVERNANCE_VERSION_LOCK.md`, `CONTRIBUTING.md`, `EXEC_SUMMARY_ONEPAGER.md`, `EXECUTIVE_SUMMARY.md`
- `ENGINE_INTEGRATION_SPEC.md`, `LEAN_REPO_INSTRUCTION_SET.md`, `DETERMINISTIC_SETUP.md`, `JUSTFILE.md`, `REPO_STRUCTURE.md`

**Everything else** (that is non-essential and not referenced) goes to `archive/`.

---

## SECTION 10 — FINAL VALIDATION (MANDATORY)

When cleanup is complete, the IDE must run:

```bash
python PX_Validation/tests/run_all_tests.py
python run_e2e_layers.py
```

If both pass:

- ✔ The repo is lean  
- ✔ The repo is deterministic  
- ✔ The repo is governed  
- ✔ The repo is partner‑ready  
- ✔ The repo is safe to publish  

If either fails → **IDE must restore the last known good state.**

---

## Repo alignment notes

- **Protected:** `PX_Warehouse/Physics_Maps/` may not exist yet; treat as reserved. `PX_LOGS/current/` may not exist; protect if created.
- **Archive candidates:** `extract_reprocess_candidates*.py`, `reprocess_pipeline*.py` live under `PX_Warehouse/Operations/scripts/`; archive only after reference check. `test_35d_ignition.py` is in `PX_Warehouse/` root.
- **99_WAREHOUSE_ARCHIVE:** Not present at repo root; may be under PX_Warehouse or created in archive; leave as specified.
- **PX_Validation:** Required for Sections 2, 4, 10; do not archive `PX_Validation/tests/` or `run_all_tests.py`.
