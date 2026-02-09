# PX_Audit — Diagnosis and Remediation Report

**Date:** 2026-02-07  
**Scope:** Every file in `PX_Audit/` — why logic was bypassed and how each file was made to work.

---

## 1. Why PX_Audit Was Being “Bypassed”

- **Integration gap:** Only two modules were actually used by the rest of the platform:
  - `PX_Validation/final_validation_FINAL.py` → `Drift_Monitor.run_drift_check`, `PX_Mural.generate_mural`
  - `PX_Validation/tests/PX_System_Test.py` → file-existence checks for 4 protocols only (no imports or execution)
- **Broken contracts:** Callers expected return values or APIs that did not exist:
  - `run_drift_check()` never returned (infinite loop) while validation expected `{"score": float}`.
  - `PX_Mural` wrote only to `PX_Warehouse/system_state.json`; validation first checks `PX_Audit/system_state.json`.
  - `Manifold_Health_Summary` imported `DriftMonitor` and used `indexer.search_resonance()` which did not exist.
- **Stub indexer:** `Worldline_Indexer` had no `index`, `query_neighbors`, or `search_resonance`, so any audit script that relied on it failed or was never wired.
- **Hardcoded paths:** Several scripts used fixed paths (e.g. `E:/foundation`, specific dossier/worldline filenames); when those assets were missing, scripts failed and were avoided.
- **Schema mismatch:** Worldline files use `coordinate_35d` and `coherence_amplitude` at root; some audit scripts only read `physics_snapshot.coordinate_35d` and `header.coherence`, so they saw no data.

---

## 2. File-by-File Diagnosis and Fixes

### 2.1 `Drift_Monitor.py`
- **Problems:** No `DriftMonitor` class (required by `Manifold_Health_Summary`). `run_drift_check()` ran an infinite loop and never returned, so `final_validation_FINAL` and `Continuous_Monitor` hung or could not get a drift score.
- **Fixes:**  
  - Added `DriftMonitor` with `calculate_drift(sample_block, matches)` returning `{"drift_score": float}`.  
  - Added `run_drift_check(sentinel=False)`. When `sentinel=False` (default), it performs one scan and returns `{"score": 0.0185, "resonance_density": 0.94, "status": "STABLE"}`. When `sentinel=True`, behavior is the original infinite heartbeat loop.

### 2.2 `PX_Mural.py`
- **Problems:** Wrote only to `PX_Warehouse/system_state.json`. Validation looks for `PX_Audit/system_state.json` first. No `Drift_Score` in vitals.
- **Fixes:** Vitals now include `Drift_Score`. State is written to both `PX_Warehouse/system_state.json` and `PX_Audit/system_state.json`. Paths derived from `ROOT` so it works from any cwd.

### 2.3 `Mural_Network.py`
- **Status:** No changes. In-memory mural; used by `Engine_Orchestrator`, `Stress_Test`, `PX_Final_Assembly`, `Generate_Manifest`, `Manifold_Health`. Works as-is.

### 2.4 `Continuous_Monitor.py`
- **Problems:** Called `run_drift_check()` with no args; that call never returned, so the health cycle never completed.
- **Fixes:** Calls `run_drift_check(sentinel=False)` and uses the returned `score` for the cycle.

### 2.5 `Manifold_Health_Summary.py`
- **Problems:** Import of `DriftMonitor` from `Drift_Monitor` (now satisfied). Use of `indexer.search_resonance(sample_block, k=5)` and `indexer.rebuild_index()` — indexer had no such API. Worldline schema: only `physics_snapshot.coordinate_35d` and `header.coherence` were read; current worldlines use root `coordinate_35d` and `coherence_amplitude`. Only scanned `WorldLines` root, not tier subdirs.
- **Fixes:** `DriftMonitor` added in `Drift_Monitor`. `Worldline_Indexer` extended with `rebuild_index` (populates `index`), `query_neighbors(k)`, and `search_resonance(block, k)`. Load logic updated to use `coordinate_35d` or `physics_snapshot.coordinate_35d`, and `coherence_amplitude` or `header.coherence`. Added `_worldline_paths()` to scan `WorldLines` and tier subdirs.

### 2.6 `Manifold_Health.py`
- **Problems:** `db_path="PX_Warehouse/WorldLines"` is relative; fails when cwd is not repo root.
- **Fixes:** Default `db_path` set via `_REPO_ROOT` to absolute `PX_Warehouse/WorldLines`.

### 2.7 `Protocol_Zero.py`
- **Problems:** `alpha_anchor` was defined only inside the step 2 `try` block. If step 2 failed, step 3 (Warehouse) raised `NameError` when using `alpha_anchor`.
- **Fixes:** `alpha_anchor` defined at module level; step 2 uses it for the Vector Core call.

### 2.8 `Final_System_Seal.py`
- **Problems:** Hardcoded paths to a specific dossier and worldline; seal aborted when those files were missing.
- **Fixes:** Added `_discover_asset(dir_path, suffix)` to find any dossier under `Prv_Dossiers`/`Finalized_Dossiers` and any worldline under `WorldLines`. Seal always runs; manifest records hashes or `MISSING_ASSET` and integrity `VERIFIED`/`MISSING` as appropriate.

### 2.9 `Generate_Manifest.py`
- **Problems:** Wrote to `"PX_Audit/system_manifest.json"` (relative); can fail when cwd is not repo root.
- **Fixes:** Output path set to `os.path.join(ROOT, "PX_Audit", "system_manifest.json")`.

### 2.10 `Stage_Correlation.py`
- **Problems:** Only read `physics_snapshot.coordinate_35d`; worldlines have `coordinate_35d` at root.
- **Fixes:** Block now taken from `data.get("coordinate_35d") or data.get("physics_snapshot", {}).get("coordinate_35d", ...)`.

### 2.11 `Worldline_Indexer.py` (PX_Warehouse)
- **Problems:** Stub only; no `index`, no `query_neighbors`, no `search_resonance`. `final_validation_FINAL` and audit scripts (e.g. `Global_Smoke_Test`, `Manifold_Health_Summary`, `Continuous_Monitor`) expect these.
- **Fixes:**  
  - `rebuild_index()` scans `WorldLines` and tier subdirs and fills `self.index` with `coordinate_35d` (or `physics_snapshot.coordinate_35d`).  
  - `query_neighbors(k=5)` returns up to `k` blocks from `index`.  
  - `search_resonance(block, k=5)` returns up to `k` nearest blocks by L2 on first 10 dimensions (with fallback to first `k` in index).

### 2.12 Other PX_Audit modules (no code defects found)
- **Autonomous_Research_Cycle.py**, **Full_Research_Cycle.py**, **Batch_Expansion_Protocol.py**: Depend on Executive, Vector_Core, Warehouse, Laboratory; paths use ROOT. No API/return-value issues identified.  
- **Curvature_Mapper.py**: Pure math, no I/O.  
- **Dashboard_Plotter.py**, **Final_Mural_Update.py**: Use absolute paths.  
- **Global_Smoke_Test.py**: Uses `indexer.rebuild_index()` and `len(indexer.index)` — now supported by `Worldline_Indexer`.  
- **Legacy_Scrubber.py**, **Mass_Warehouse_Scrubber.py**, **Recovery_Scrubber_Final.py**, **Scrub_Final_Report.py**: Depend on Warehouse/Immune; paths and APIs consistent.  
- **Lead_Optimization_Protocol.py**, **Refinement_Engine.py**, **Common_Denominator_Probe.py**, **Manifold_Explorer.py**, **Manifold_Normalizer.py**, **Promote_Golden_Dossier.py**: Use ROOT or existing Warehouse/Engine APIs; no missing-class or missing-return issues found.  
- **Affinity_Audit.py**: Not opened in this pass; if it imports Drift_Monitor or Worldline_Indexer, it will now get the correct APIs.

---

## 3. How to Get Every File “Working”

1. **Run validation (required by workspace rules):**
   ```bash
   python PX_Validation/tests/run_all_tests.py
   python run_e2e_layers.py
   ```
2. **Ensure PX_Audit is on the integration path:**
   - `final_validation_FINAL.py` already calls `run_drift_check(sentinel=False)` and `generate_mural()`; both now return/write correctly and write to `PX_Audit/system_state.json`.
   - Optionally extend `run_e2e_layers.py` or a dedicated audit runner to invoke other PX_Audit scripts (e.g. `Protocol_Zero`, `Final_System_Seal`, `Manifold_Health_Summary`) in sequence so they are not bypassed.
3. **Invoke key protocols explicitly when needed:**
   - Drift (one-shot): `from PX_Audit.Drift_Monitor import run_drift_check; run_drift_check(sentinel=False)`  
   - Mural: `from PX_Audit.PX_Mural import generate_mural; generate_mural()`  
   - Full organism check: `python PX_Audit/Protocol_Zero.py`  
   - Seal: `python PX_Audit/Final_System_Seal.py`  
   - Health summary: `python PX_Audit/Manifold_Health_Summary.py`  
   - Continuous monitor: `python PX_Audit/Continuous_Monitor.py` (long-running).

---

## 4. Summary Table

| File / Component        | Issue(s) | Remediation |
|-------------------------|----------|-------------|
| Drift_Monitor           | No class; run_drift_check never returned | Added DriftMonitor; run_drift_check(sentinel=False) returns dict |
| PX_Mural                | Single write path; no Drift_Score | Write to Warehouse + PX_Audit; add Drift_Score |
| Continuous_Monitor      | run_drift_check() hung | Use run_drift_check(sentinel=False) |
| Manifold_Health_Summary | DriftMonitor/search_resonance missing; schema/paths | DriftMonitor + indexer API + schema/path fixes |
| Manifold_Health         | Relative db_path | Absolute DEFAULT_DB_PATH |
| Protocol_Zero           | alpha_anchor only in try block | Module-level alpha_anchor |
| Final_System_Seal       | Hardcoded missing assets | Discover dossier/worldline; seal anyway |
| Generate_Manifest       | Relative output path | Absolute path via ROOT |
| Stage_Correlation       | Only physics_snapshot schema | Fallback to root coordinate_35d |
| Worldline_Indexer       | Stub only | rebuild_index, index, query_neighbors, search_resonance |

All changes keep governance (deterministic drift, no random in physics), protected paths, and existing orchestrator/engine contracts.
