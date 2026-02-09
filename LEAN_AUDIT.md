# Predator X — Lean Repo Audit (Legacy & Unused Files)

**Purpose:** Identify legacy scripts, obsolete tests, and unused files so the repo stays lean **without losing functionality or data**. Every removal must be verified so the pipeline does not break.

**Verification after any removal:** Run from repo root:
1. `python PX_Validation/tests/run_all_tests.py`
2. `python run_e2e_layers.py`
3. If either fails, do **not** commit the removal; restore or fix first.

---

## 1. Already archived (do not use for execution)

| Location | Notes |
|----------|--------|
| **99_LEGACY_CODE/** | Orchestrators and pipelines superseded by canonical paths. **Do not import or run from here.** See `99_LEGACY_CODE/README.md` and `.cursor/rules/05_governance_version_lock.mdc`. |
| 99_LEGACY_CODE/orchestrators/Predator_X_v3_Orchestrator.py | BAN: use PRV_24H_Orchestrator.py |
| 99_LEGACY_CODE/orchestrators/PX_Live_Orchestrator.py | BAN: use PX_Live_Orchestrator_v2.py |
| 99_LEGACY_CODE/pipelines/PRV_Master_Pipeline*.py | Reference only; canonical pipeline is run via PRV_24H + PX_Live_Orchestrator_v2 |
| 99_LEGACY_CODE/warehouse_scripts/*.py | Old v3 stratification/ranking; not used by active pipeline |

**Confirmed:** No active code imports from `99_LEGACY_CODE`. `tests/run_novel_pipeline_test.py` defines PRV_TARGETS inline and does **not** import from legacy.

---

## 2. Root and critical paths (keep)

| File / folder | Reason |
|---------------|--------|
| **run_e2e_layers.py** | Required by governance (RULE_3_VALIDATION_FIRST). Do not remove. |
| **tests/** | Contains constitutional pre-execution tests (`run_poison_pill_test.py`, `run_temporal_paradox_test.py`) and finalization/orchestrator tests. Referenced by `PX_Constitution/Constitutional_Tests.py` and `governance/poison_pill_gate.py`. |
| **PX_Executive/run_finalize_dossiers.py** | Canonical backfill/finalization. |
| **PX_Executive/orchestrators/PX_Live_Orchestrator_v2.py** | ORCHESTRATOR_LIVE per governance. |
| **PX_Executive/PRV_24H_Orchestrator.py** | ORCHESTRATOR_24H per governance. |

---

## 3. Fixes applied (no removal; pipeline preserved)

| Item | Change |
|------|--------|
| **PX_Validation/final_validation_FINAL.py** | Replaced broken `from PX_Laboratory.PX_Materialize import materialize_candidate` with `Simulation_Engine.SimulationEngine().materialize_candidate(...)`. PX_Materialize does not exist; materialization is canonical on Simulation_Engine. |
| **PX_Warehouse/test_35d_ignition.py** | Added `random.seed(42)` for determinism (RULE_2_DETERMINISTIC_PHYSICS). |

---

## 4. Validation suite (keep; optional consolidation)

| File | Status |
|------|--------|
| **PX_Validation/final_validation.py** | Stub suite (all `results.append(True)`). Referenced in system_manifest.json and ARCHITECTURE_FILE_LISTING. Safe to keep as lightweight entry point or later replace with a single runner that calls PX_Validation/tests/run_all_tests.py + run_e2e_layers. |
| **PX_Validation/final_validation_FINAL.py** | Hardened integration suite; referenced by warehouse simulation reports. **Fixed** (see §3). Keep. |
| **PX_Validation/tests/run_all_tests.py** | Canonical test runner. Keep. |
| **PX_Validation/manual_tests/** | Manual/adjudication tests. Keep unless you explicitly retire manual flows. |

---

## 5. Scripts to review (consolidate or keep, do not delete without verification)

- **extract_reprocess_candidates.py** vs **extract_reprocess_candidates_v2.py**  
  No other code imports these. If only one is used in operations, consider deprecating the other and documenting the canonical script.

- **reprocess_pipeline.py** vs **reprocess_pipeline_direct.py**  
  `reprocess_pipeline.py` invokes PX_Live_Orchestrator_v2 (canonical). `reprocess_pipeline_direct.py` writes to its own log. Decide which path is canonical; document and optionally retire the other.

- **PX_LOGS/extract_fto_blocked_candidates.py**  
  Standalone utility. Can stay in PX_LOGS or move to `PX_Warehouse/Operations/scripts/` for consistency. Moving is optional; no need to delete.

- **PX_Warehouse/test_35d_ignition.py**  
  Dev/test ignition script. Now deterministic. Can stay; or move to `tests/` if you want all ad‑hoc tests under one tree.

---

## 6. Data and logs (do not remove without policy)

- **PX_LOGS/**  
  Contains run summaries, discovery reports, FTO audit, warehouse_simulation reports. **Do not bulk-delete** without a retention policy. Old Batch_* stdout/stderr were already removed (per git status); current contents are operational artifacts.

- **PX_Warehouse/**  
  WorldLines, Novel_Dossiers, Prv_Dossiers, Finalized_Dossiers, etc. are pipeline data. **Do not remove** except via documented archival/purge scripts (e.g. purge_hallucinations, run_warehouse_execution).

- **99_WAREHOUSE_ARCHIVE**  
  Referenced in `unified_warehouse_consolidation.py` and `purge_hallucinations.py` for quarantine/archive. Keep path and references unless you formally retire that archive location.

---

## 7. Duplicate or optional copies (safe to archive after verification)

- **TSO_Validator/TSO_PUBLIC_RELEASE/run_poison_pill_test.py**  
  Duplicate of `tests/run_poison_pill_test.py`. If TSO release is self-contained, keep for the release bundle; otherwise you can document that the canonical test is `tests/run_poison_pill_test.py` and remove the copy after running both test paths once.

- **PX_Executive/docs/Manufacturing_Manifest.py**  
  Different from `PX_Laboratory/Manufacturing_Manifest.py` (worldline production order vs lab manifest). Not duplicate; keep unless you retire the doc/script.

---

## 8. Safe removal checklist (before deleting anything)

1. **Grep for imports and references**  
   Search repo for the file name and any `import`/`from` that would load it. If anything outside 99_LEGACY_CODE or tests references it, do not delete without updating callers.

2. **Run validation**  
   - `python PX_Validation/tests/run_all_tests.py`  
   - `python run_e2e_layers.py`  

3. **If removing a script that might be run by hand**  
   Document the removal in CHANGELOG or LEAN_AUDIT.md and point to the canonical replacement (if any).

4. **Prefer archive over delete**  
   When unsure, move to `99_LEGACY_CODE/` (or an archive subfolder) instead of deleting, so functionality and data are preserved.

---

## 9. Summary

- **No active code** imports from `99_LEGACY_CODE`; legacy is correctly isolated.
- **Root and tests:** `run_e2e_layers.py` and `tests/` (including poison pill and temporal paradox) are required; keep them.
- **Fixes applied:** `final_validation_FINAL.py` now uses Simulation_Engine for materialization; `test_35d_ignition.py` is deterministic.
- **Lean actions that are safe after verification:** (1) Consolidate extract_reprocess_candidates and reprocess_pipeline variants once a canonical choice is documented; (2) optionally move or retire duplicate TSO run_poison_pill_test copy; (3) leave PX_LOGS and warehouse data intact unless you have a retention/archival policy.

**Last audit:** 2026-02-06
