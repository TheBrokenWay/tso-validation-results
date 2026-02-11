# ENTERPRISE CLEANUP AUDIT REPORT
Generated: 2026-02-11T15:30:00Z

## EXECUTIVE SUMMARY

| Metric | Value |
|--------|-------|
| Total Python files | ~359 (excluding .venv, .conda, __pycache__) |
| Total lines of Python code | ~49,883 |
| Legacy folders to remove | 4 (01_Executive, 02_Audit, PX_STATE, MANIFESTS) |
| Folders to keep separate (NOT merge) | 2 (Nipah_Analysis, TSO_Validator) |
| Duplicate files found | 78 (PX_Warehouse root shadows of placement_gate/Staging/) |
| Orphan files found | ~17 (PX_Audit utilities, deprecated wrappers) |
| Deprecated engines to archive | 4 (PKPD_Simple, DoseOptimizer_Simple, DoseOptimizer, VirtualEfficacy_Simple) |
| Deprecated orchestrators to archive | 6 (run_*.py wrappers + PRV_24H_Orchestrator) |
| Drive root (E:/) | CLEAN — no duplicates, no orphans |
| Estimated file reduction | ~100 files |

---

## SECTION 1: LEGACY FOLDER CONSOLIDATION

### 01_Executive/ → PX_Executive/

| File | Action | Reason | Conflicts |
|------|--------|--------|-----------|
| `CSA_Pentarchy.py` | **MERGE** to PX_Executive/ | Dynamically loaded by PX_System/benchmark.py and PX_System/console.py | Update importlib load paths to standard imports |
| `__init__.py` | DELETE | Empty package marker | None |

**Risk: HIGH** — Must update dynamic load paths in benchmark.py and console.py before removing.

### 02_Audit/ → PX_Audit/ (or PX_System/)

| File | Action | Reason | Conflicts |
|------|--------|--------|-----------|
| `AAS_Verification.py` | **MERGE** to PX_Audit/ | Dynamically loaded by PX_System/benchmark.py and PX_System/console.py | Update importlib load paths |
| `__init__.py` | DELETE | Empty package marker | None |

**Risk: HIGH** — Same dynamic loading issue as 01_Executive.

### tests/ (root) → PX_Validation/tests/

| File | Action | Reason | Conflicts |
|------|--------|--------|-----------|
| `run_poison_pill_test.py` | **MERGE** (carefully) | **CONSTITUTIONAL** — called by PX_Constitution/Constitutional_Tests.py and governance/poison_pill_gate.py | Must update hardcoded paths in Constitutional_Tests.py line 12 and poison_pill_gate.py line 37 |
| `run_temporal_paradox_test.py` | **MERGE** (carefully) | **CONSTITUTIONAL** — called by PX_Constitution/Constitutional_Tests.py line 13 | Same path update required |
| `run_full_verification.py` | MERGE | Meta-test runner | No conflicts |
| `run_genesis_e2e_test.py` | MERGE | E2E cycle test | No conflicts |
| `run_novel_pipeline_test.py` | MERGE | PRV pipeline test | No conflicts |
| `run_warehouse_path_test.py` | MERGE | Warehouse structure test | No conflicts |
| `test_finalization_spec.py` | MERGE | Finalization spec test | No conflicts |
| `verify_all_imports.py` | MERGE | Import smoke test | No conflicts |
| `test_orchestrator_warehouse_paths.py` | ARCHIVE | Low priority, no active references | None |
| `verify_connectivity_and_orphans.py` | ARCHIVE | Superseded by other tools | None |
| `verify_every_file.py` | ARCHIVE | Redundant with lineage_status.py | None |
| `__init__.py` | DELETE | Package marker | None |

**Risk: CRITICAL** — Poison pill tests are constitutional entry gates. Must update ALL path references before moving.

### scripts/ → Split

| File | Action | Destination | Reason |
|------|--------|-------------|--------|
| `lineage_status.py` | **KEEP** | scripts/ (or governance/) | Active — called by `just lineage` |
| `pre_commit_check.py` | **KEEP** | scripts/ (or governance/) | Active — wired into .git/hooks/pre-commit |
| `predator_x_deterministic_setup.sh` | ARCHIVE | governance/setup_guides/ | One-time setup, not runtime |
| `fix_git_and_nix_wsl.sh` | ARCHIVE | governance/setup_guides/ | Manual recovery, not runtime |

**Risk: MEDIUM** — pre_commit_check.py is wired into git hooks; lineage_status.py into justfile. If moving, update both.

### MANIFESTS/ → governance/

| File | Action | Reason | Conflicts |
|------|--------|--------|-----------|
| `SMART_Antiviral_Fork.manifest.json` | MERGE to governance/ or PX_Domain/ | Standalone BARDA spec, no imports | Hardcoded output path to PX_Warehouse/SMART_Antiviral_Dossiers |

**Risk: LOW** — No code imports this file.

### Nipah_Analysis/ — KEEP SEPARATE

| Decision | Reason |
|----------|--------|
| **DO NOT MERGE into PX_Domain/** | Active E2E layer 2 integration. Imported by run_e2e_layers.py (NipahMinerAdapter). Complementary to PX_Domain (epidemiology vs. pharmaceutical constraints). Has own tests, data, and analytics pipeline. |

**Risk: NONE** — Leave as-is.

### TSO_Validator/ — KEEP SEPARATE

| Decision | Reason |
|----------|--------|
| **DO NOT MERGE into PX_Validation/** | Architecturally isolated by design (zero PX_* imports, stdlib only). Runs as E2E layer 7 via subprocess. Tests already integrated into run_all_tests.py harness. Standalone safety gate that can validate ANY data. |

**Risk: NONE** — Leave as-is.

### PX_STATE/ → DELETE (after archiving value)

| File | Action | Reason | Conflicts |
|------|--------|--------|-----------|
| `current_state.json` | DELETE | Ephemeral state snapshot (batch 170, Jan 23). Zero imports, zero consumers. Superseded by PX_LOGS/ audit trails. | None |

**Risk: LOW** — No code references this directory.

---

## SECTION 2: PX_WAREHOUSE DEDUPLICATION

**Critical finding:** 78 Python files in `PX_Warehouse/` root are identical shadow copies of files in `PX_Warehouse/placement_gate/Staging/` (the canonical location per CLAUDE.md).

| Category | Count | Action |
|----------|-------|--------|
| Identical root shadows of Staging/ files | ~78 | DELETE root copies, keep Staging/ |
| Canonical files (warehouse_layout.py, Finalization_Pipeline.py, WorldLine_Database.py) | 3 | KEEP — these re-export from Staging/ |
| Finalization_Spec.py | 1 | KEEP |

**Key duplicates to remove (sample):**
purge_hallucinations.py, run_batch_pipeline.py, sort_live_results.py, one_hour_pipeline.py, execute_pipeline_batch.py, process_learning_material.py, archive_state.py, generate_analytics.py, rollback_audit.py, SystemTelemetry.py, extract_reprocess_candidates.py, unified_warehouse_consolidation.py, bootstrap_canonical_warehouse.py, flatten_tier_folders.py, Worldline_Indexer.py, run_intake_harvester.py, implement_research_plan_v2.py, migrate_legacy_labels.py, master_fix_calibration_and_prv_recovery.py, warehouse_architect_phase2.py, recover_misplaced_warehouse_data.py, run_warehouse_execution.py, and ~56 more.

**Risk: MEDIUM** — Must verify no imports reference `PX_Warehouse.filename` (vs `PX_Warehouse.placement_gate.Staging.filename`) before deleting.

---

## SECTION 3: DEPRECATED ENGINES & ORCHESTRATORS

### Deprecated Engines (zero imports)

| File | Replaced By | Action |
|------|-------------|--------|
| `PX_Engine/operations/PKPD_Simple.py` | PKPD.py | DELETE |
| `PX_Engine/operations/DoseOptimizer_Simple.py` | DoseOptimizer_v2.py | DELETE |
| `PX_Engine/operations/DoseOptimizer.py` | DoseOptimizer_v2.py | DELETE |
| `PX_Engine/operations/VirtualEfficacy_Simple.py` | VirtualEfficacyAnalytics.py | DELETE |

### Deprecated Orchestrators (backward compat wrappers, no direct callers)

| File | Delegates To | Action |
|------|-------------|--------|
| `PX_Executive/run_genesis_feed.py` | px_feed.py | DELETE |
| `PX_Executive/run_repurposed_feed.py` | px_feed.py | DELETE |
| `PX_Executive/run_prv_novel.py` | PRV_24H_Orchestrator.py | DELETE |
| `PX_Executive/run_prv_repurposed.py` | PRV_24H_Orchestrator.py | DELETE |
| `PX_Executive/run_finalize_dossiers.py` | px_finalize.py | DELETE |
| `PX_Executive/PRV_24H_Orchestrator.py` | px_prv.py (12-engine) | DELETE |

---

## SECTION 4: DRIVE ROOT CLEANUP

| Path | Contents | Action | Reason |
|------|----------|--------|--------|
| E:/PX_Warehouse/ | Does not exist | N/A | No duplicate |
| E:/archive/ | ChimeraX 1.11 (external tool) | IGNORE | Not repo-related |
| E:/audit/ | Does not exist | N/A | Clean |
| E:/tmp/ | Does not exist | N/A | Clean |
| E:/OLYMPUS_GOLD/ | Does not exist | N/A | Clean |
| E:/foundation-full-history.bundle | 226MB git bundle | KEEP | Documented backup |
| E:/*.py, *.json, *.md | None found | N/A | Clean |

**Drive root is CLEAN. No action required.**

---

## SECTION 5: IMPORT UPDATES REQUIRED

Files that import from legacy folders (need path updates after migration):

| File | Current Import | New Import |
|------|----------------|------------|
| `PX_System/benchmark.py` (line 20) | `importlib.util.spec_from_file_location("...", "01_Executive/CSA_Pentarchy.py")` | `from PX_Executive.CSA_Pentarchy import CSAPentarchy` |
| `PX_System/benchmark.py` (line 21) | `importlib.util.spec_from_file_location("...", "02_Audit/AAS_Verification.py")` | `from PX_Audit.AAS_Verification import AASVerification` |
| `PX_System/console.py` (line 11) | `importlib.util.spec_from_file_location("...", "01_Executive/CSA_Pentarchy.py")` | `from PX_Executive.CSA_Pentarchy import CSAPentarchy` |
| `PX_System/console.py` (line 12) | `importlib.util.spec_from_file_location("...", "02_Audit/AAS_Verification.py")` | `from PX_Audit.AAS_Verification import AASVerification` |
| `PX_Constitution/Constitutional_Tests.py` (line 12) | `tests/run_poison_pill_test.py` (hardcoded path) | `PX_Validation/tests/run_poison_pill_test.py` |
| `PX_Constitution/Constitutional_Tests.py` (line 13) | `tests/run_temporal_paradox_test.py` (hardcoded path) | `PX_Validation/tests/run_temporal_paradox_test.py` |
| `governance/poison_pill_gate.py` (line 37) | `tests/run_poison_pill_test.py` (subprocess path) | `PX_Validation/tests/run_poison_pill_test.py` |

---

## SECTION 6: ORPHAN FILES

Files not imported anywhere (candidates for review/deletion):

### PX_Audit Orphans (standalone utilities, zero imports)

| File | Lines | Recommendation |
|------|-------|----------------|
| `PX_Audit/Curvature_Mapper.py` | ~100 | REVIEW — analysis utility |
| `PX_Audit/Scrub_Final_Report.py` | ~80 | REVIEW — cleanup utility |
| `PX_Audit/Common_Denominator_Probe.py` | ~90 | REVIEW — analysis utility |
| `PX_Audit/Affinity_Audit.py` | ~70 | REVIEW — analysis utility |
| `PX_Audit/Mass_Warehouse_Scrubber.py` | ~100 | REVIEW — cleanup utility |
| `PX_Audit/Recovery_Scrubber_Final.py` | ~80 | REVIEW — recovery utility |
| `PX_Audit/Manifold_Normalizer.py` | ~90 | REVIEW — data normalization |
| `PX_Audit/Legacy_Scrubber.py` | ~80 | REVIEW — legacy cleanup |
| `PX_Audit/Stage_Correlation.py` | ~70 | REVIEW — analysis utility |
| `PX_Audit/Refinement_Engine.py` | ~90 | REVIEW — analysis utility |
| `PX_Audit/Final_System_Seal.py` | ~60 | REVIEW — sealing utility |
| `PX_Audit/Lead_Optimization_Protocol.py` | ~100 | REVIEW — optimization protocol |
| `PX_Audit/Promote_Golden_Dossier.py` | ~90 | REVIEW — promotion utility |

### PX_Executive Orphans

| File | Lines | Recommendation |
|------|-------|----------------|
| `PX_Executive/debug_factory.py` | ~50 | REVIEW — dev utility |
| `PX_Executive/Byzantium_Council.py` | ~100 | REVIEW — appears unused |
| `PX_Executive/harvest_leads.py` | ~80 | REVIEW — lead harvesting |
| `PX_Executive/PX_Refinery.py` | ~90 | REVIEW — dossier refinement |
| `PX_Executive/validate_v2_release.py` | ~70 | REVIEW — v2 test utility |

---

## SECTION 7: ENVIRONMENT CACHES

| Path | Action | Reason |
|------|--------|--------|
| `.conda/` | DELETE | Regenerates; use .venv only |
| `.direnv/` | DELETE | Regenerates from .envrc |
| `.pytest_cache/` | DELETE | Regenerates on test run |
| `nul` (root) | DELETE | Windows null device artifact |

Already in .gitignore: `.venv/` (KEEP but not tracked).

---

## SECTION 8: RECOMMENDED EXECUTION ORDER

| Batch | Actions | Test Command | Risk |
|-------|---------|-------------|------|
| 0 | Review this report together | — | — |
| 1 | Delete caches (.pytest_cache, .direnv, .conda, nul) | `just test` | LOW |
| 2 | Delete PX_STATE/ | `just test` | LOW |
| 3 | Move CSA_Pentarchy.py → PX_Executive/, AAS_Verification.py → PX_Audit/, update importlib paths | `just test` + `just govern` | HIGH |
| 4 | Delete empty 01_Executive/, 02_Audit/ | `just test` | LOW |
| 5 | Move tests/ → PX_Validation/tests/, update Constitutional_Tests.py + poison_pill_gate.py paths | `just test` + `just govern` | CRITICAL |
| 6 | Move scripts/ utilities to governance/setup_guides/ (keep lineage + pre-commit in scripts/) | `just test` | LOW |
| 7 | Move MANIFESTS/ → governance/manifest_forks/ | `just test` | LOW |
| 8 | Delete deprecated engines (PKPD_Simple, DoseOptimizer_Simple, DoseOptimizer, VirtualEfficacy_Simple) | `just test` | LOW |
| 9 | Delete deprecated orchestrators (run_*.py wrappers, PRV_24H_Orchestrator) | `just test` + `just govern` | MEDIUM |
| 10 | Delete PX_Warehouse root shadow files (78 duplicates of Staging/) | `just test` | MEDIUM |
| 11 | Review and handle PX_Audit orphan utilities (13 files) | `just test` | LOW |
| 12 | Update governance docs (CLAUDE.md, README.md, .cursor/rules) | `just test` + `just govern` | LOW |
| 13 | Full verification | `just cycle` (Feed → PRV → Finalize) | — |
| 14 | Commit and push | `git status` clean | — |

---

## SECTION 9: GOVERNANCE DOCUMENTATION UPDATES

After cleanup, update these files:
- [ ] **CLAUDE.md** — Update directory structure table (remove 01_Executive, 02_Audit, tests/, MANIFESTS/, PX_STATE; note Nipah_Analysis and TSO_Validator stay)
- [ ] **README.md** — Update for enterprise-ready presentation
- [ ] **.cursor/rules/06_engine_integration.mdc** — Update canonical orchestrator list (remove deprecated run_*.py references)
- [ ] **.cursor/rules/07_lean_repo.mdc** — Update target state to match new structure
- [ ] **justfile** — Update any paths that reference moved files
- [ ] **.git/hooks/pre-commit** — Update if scripts/pre_commit_check.py moves

---

## KEY DECISIONS FOR FOUNDER

1. **Nipah_Analysis/** — Audit recommends KEEP SEPARATE (not merge into PX_Domain). It's actively integrated in E2E governance layer 2. Confirm?

2. **TSO_Validator/** — Audit recommends KEEP SEPARATE (not merge into PX_Validation). It's architecturally isolated by design (zero PX_* imports). Confirm?

3. **PX_Audit orphans (13 files)** — These are standalone analysis/cleanup utilities with zero imports. DELETE all, KEEP all, or REVIEW individually?

4. **PX_Warehouse shadow files (78 files)** — All identical to placement_gate/Staging/ canonical copies. Safe to delete but need import verification first. Proceed?

5. **Root-level .md files** (ARCHITECTURE_DIAGRAM.md, CONTRIBUTING.md, ENGINE_INTEGRATION_SPEC.md, GOVERNANCE_VERSION_LOCK.md, JUSTFILE.md, LEAN_AUDIT.md, LEAN_REPO_INSTRUCTION_SET.md, REPO_STRUCTURE.md, STATUS.md, SYSTEM_FLOW.md) — Move to governance/ or keep at root?
