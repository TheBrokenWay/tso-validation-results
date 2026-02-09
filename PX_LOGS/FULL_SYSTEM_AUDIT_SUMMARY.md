# FULL SYSTEM CONNECTIVITY AUDIT SUMMARY

**Date:** 2026-02-08
**Scope:** 26 directory areas, ~6,100 files (~280 .py)
**Method:** Parallel sub-agent audit per directory — tracing imports, importers, tests, orchestrator connections, warehouse writes, and gating status

---

## EXECUTIVE SUMMARY

| Metric | Count |
|--------|-------|
| Total files audited | ~6,100 |
| Python files | ~280 |
| Connected (live code) | ~198 |
| Disconnected (dead code) | 82 |
| Warehouse write paths | 28 |
| Gated writes (placement_gate/Zeus) | 5 |
| **Ungated writes (bypass)** | **23** |
| Bypass violations | 22 |
| Broken imports | 1 |
| Rule violations (CRITICAL) | 1 |
| Tests in harness | 7 of 22 |
| Misplaced warehouse root files | 5,631 |

**Bottom line:** Only 18% of warehouse writes are gated. 82 files are dead code. 5,631 dossiers sit at warehouse root instead of tier subdirs. The finalization pipeline exists but never executes because of the schema mismatch (toxicity_index field location).

---

## CONNECTIVITY MATRIX — DIRECTORY COMMUNICATION GRAPH

```
                    PX_Executive
                   /     |     \
                  v      v      v
            PX_Engine  PX_System  PX_Laboratory
               |    \     |         |
               v     v    v         v
        PX_Constitution  PX_Warehouse <----+---- PX_Audit
               ^          |    ^           |
               |          v    |           |
            governance   WorldLines     PX_Security
               ^          |
               |          v
          Nipah_Analysis  Finalized_Dossiers (EMPTY)
               |
               v
           PX_Domain (manifest.json)
```

### Who talks to whom:

| Source | Destinations |
|--------|-------------|
| PX_Executive | PX_Engine, PX_Warehouse, PX_System, PX_Laboratory, PX_Constitution |
| PX_Engine | PX_Constitution, PX_Security, PX_Warehouse, PX_Laboratory |
| PX_Audit | PX_Executive, PX_Engine, PX_Laboratory, PX_Warehouse, PX_Constitution, PX_Security |
| PX_System | PX_Warehouse, PX_Audit |
| PX_Security | PX_Constitution, PX_Warehouse, PX_Audit |
| PX_Validation | PX_Engine, PX_Executive, PX_Laboratory, PX_System, PX_Warehouse, PX_Security |
| PX_Laboratory | PX_Warehouse, PX_Engine |
| Nipah_Analysis | PX_Engine, PX_Executive, governance, PX_Warehouse |
| governance | PX_Constitution |
| tests_root | PX_Engine, PX_Executive, PX_Warehouse, PX_System, PX_Constitution, governance, Nipah_Analysis |

### Chain Breaks (no outbound connections):

- **PX_Constitution** — talked to by 4 dirs, talks to none (pure dependency leaf)
- **PX_Discovery** — stub, returns empty list
- **PX_Data** — data-only, no code imports outward
- **PX_Domain** — data-only, manifest.json consumed by 6+ dirs
- **MANIFESTS** — fully orphaned
- **scripts** — setup-only (WSL/Nix), no code dependencies

### Isolated Islands:

- **01_Executive** — thin wrapper around PX_Executive.GAIP_Gateway, only connected to PX_System
- **02_Audit** — thin wrapper, only connected to PX_System
- **MANIFESTS** — zero connections in or out

---

## WAREHOUSE WRITE PATHS — COMPLETE INVENTORY

### GATED (5 files) — pass through placement_gate and/or Zeus gate:

| File | Gate Method |
|------|-------------|
| `PX_Executive/PRV_24H_Orchestrator.py` | placement_gate.place_prv_dossier() |
| `PX_Executive/run_finalize_dossiers.py` | Finalization_Pipeline + Zeus gate |
| `PX_Executive/reprocess_warehouse_lifecycle.py` | GradingEngine + placement_gate |
| `PX_Engine/Block_Orchestrator.py` | WorldLine_Database (materialization gated) |
| `PX_Warehouse/placement_gate/Staging/placement_gate.py` | IS the gate (canonical) |

### UNGATED (23 files) — write directly, bypassing grading + placement:

| # | File | Target |
|---|------|--------|
| 1 | `PX_Executive/PX_Live_Orchestrator_v2.py` | PX_Warehouse root |
| 2 | `PX_Executive/generators/SMART_Antiviral_Fork.py` | PX_Warehouse root |
| 3 | `PX_Executive/Rank_Diamonds.py` | PX_Warehouse root |
| 4 | `PX_Executive/Generate_BARDA_Brief.py` | PX_Warehouse root |
| 5 | `PX_Executive/run_genesis_feed.py` | PX_Warehouse/Feeder |
| 6 | `PX_Executive/Sovereign_Commercial_Pipeline.py` | PX_Warehouse root |
| 7 | `PX_Executive/Gold_Rush_Miner.py` | PX_Warehouse root |
| 8 | `PX_Executive/batch/pipeline_batch_runner.py` | PX_Warehouse root |
| 9 | `PX_Executive/Acquisition/olympus_bridge_v2.py` | PX_Warehouse root |
| 10 | `PX_Executive/debug_factory.py` | PX_Warehouse root |
| 11 | `PX_Engine/Metabolism.py` | PX_Warehouse (no placement_gate) |
| 12 | `PX_Laboratory/Manufacturing_Manifest.py` | PX_Warehouse/Orders/ |
| 13 | `PX_System/foundation/Evidence_Package.py` | Calibration_Molecules (relative path default) |
| 14 | `PX_Security/AAS_CircuitBreaker.py` | aas_circuit_breaker.json |
| 15 | `PX_Security/LogForensics.py` | logs/forensic_manifest.json |
| 16 | `PX_Security/Immune_Test.py` | system_state.json |
| 17 | `PX_Audit/Autonomous_Research_Cycle.py` | PX_Warehouse root |
| 18 | `PX_Audit/Batch_Expansion_Protocol.py` | PX_Warehouse root |
| 19 | `PX_Audit/Full_Research_Cycle.py` | PX_Warehouse root |
| 20 | `PX_Audit/Lead_Optimization_Protocol.py` | PX_Warehouse root |
| 21 | `PX_Audit/Legacy_Scrubber.py` | PX_Warehouse root |
| 22 | `PX_Audit/Manifold_Normalizer.py` | PX_Warehouse root |
| 23 | `PX_Audit/Protocol_Zero.py` | PX_Warehouse root |

**Root cause of 5,631 misplaced files:** Items 1-10 and 17-23 write dossiers directly to `PX_Warehouse/` root without calling `placement_gate.place_prv_dossier()`. The placement gate exists and works — it's just never called by most code paths.

---

## ORPHANED DATA

### Misplaced warehouse root files (5,631 total):

| Pattern | Count | Correct Location |
|---------|-------|-----------------|
| `PRV_REP_*.json` | ~2,500 | `Prv_Dossiers/<DIAMOND\|GOLD\|SILVER>` |
| `PRV_NOV_*.json` | ~3,000 | `Novel_Dossiers/<DIAMOND\|GOLD\|SILVER>` |
| `TRIAL_SIMULATION_DOSSIER_*.json` | 31 | `Calibration_Molecules/` |
| Policy/config artifacts | 6 | Various or remove |

### Orphaned warehouse subdirectories (empty, no code references):

1. `Archive_Novel/` — empty
2. `Archive_Primary/` — empty
3. `Backup_Pre_Refinery/` — empty
4. `CommercialAssets/` — empty
5. `TrialSimulations/` — empty

### Other orphaned data:

- `MANIFESTS/SMART_Antiviral_Fork.manifest.json` — no code references it

---

## DEAD CODE INVENTORY (82 files)

### PX_Audit (15 disconnected):
- Affinity_Audit.py, Common_Denominator_Probe.py, Curvature_Mapper.py
- Final_System_Seal.py, Final_Mural_Update.py, Lead_Optimization_Protocol.py
- Legacy_Scrubber.py, Manifold_Explorer.py, Mass_Warehouse_Scrubber.py
- Promote_Golden_Dossier.py, PX_Mural.py, Recovery_Scrubber_Final.py
- Refinement_Engine.py, Scrub_Final_Report.py, Stage_Correlation.py

### PX_System (14 disconnected):
- foundation/core.py, foundation/Disease_Constraint_Model.py
- foundation/Emergency_Stop.py, foundation/Research_Checkpoint.py
- foundation/Novelty_Budget_Engine.py, foundation/integrations/smiles_security.py
- foundation/integrations/net_policy.py, foundation/integrations/retry.py
- PX_Handshake.py, auto_patch.py, benchmark.py, console.py
- Sign_Off_GAIP.py, PX_Final_Assembly.py

### PX_Executive (10 disconnected):
- Patent_Local_Index.py, PX_Refinery.py
- Acquisition/patent_harvester.py, Acquisition/harvest_failed_trials.py
- patent_data_refresh.py, monitor_warehouse.py, purge_fto_error_logs.py
- harvest_leads.py, validate_v2_release.py, diagnose_fto_failure_rate.py

### PX_Engine (4 disconnected):
- Stress_Test.py
- operations/DoseOptimizer_Simple.py **(BANNED by version lock)**
- operations/PKPD_Simple.py **(BANNED by version lock)**
- operations/VirtualEfficacy_Simple.py **(BANNED by version lock)**

### PX_Constitution (3 disconnected):
- Block_Universe.py, Filing_Rules.py, gate_stamp.py

### PX_Validation (3 disconnected):
- benchmarks/dashboard.py, benchmarks/report_generator.py, benchmarks/run_benchmarks.py

### tests_root (4 disconnected):
- run_novel_pipeline_test.py, run_warehouse_path_test.py
- verify_every_file.py, verify_connectivity_and_orphans.py

### governance (3 disconnected):
- Design_Principles.py, poison_pills/README.md, ARCHITECTURE_AND_RULES_REFERENCE.md

### Other (6 disconnected):
- PX_Data: drugbank/drugbank.csv, drugbank/README.md, patents/README.md
- scripts: 2 shell scripts (WSL/Nix setup)
- MANIFESTS: SMART_Antiviral_Fork.manifest.json

---

## RULE VIOLATIONS

### CRITICAL: Rule 2 — No random physics

| File | Lines | Detail |
|------|-------|--------|
| `PX_Executive/generators/SMART_Antiviral_Fork.py` | 73, 82 | Uses `np.random.uniform()` for `target_potency` and `host_potency` in physics calculations |

### Broken Imports

| File | Import | Status |
|------|--------|--------|
| `PX_Audit/Refinement_Engine.py` | `PX_Warehouse.RAG_Query_Engine` | Module does not exist; uses stub fallback |

### Test Harness Gap

Only **7 of 22** test files are included in `run_all_tests.py`. **15 tests are excluded:**

test_adaptive.py, test_adaptive_integration.py, test_dose_optimizer_v2.py,
test_dose_optimizer_v2_integration.py, test_grading_engine.py, test_iiv.py,
test_iiv_integration.py, test_orchestrator_v2.py, test_performance_regression.py,
test_pkpd.py, test_pkpd_integration.py, test_system_trial_evidence_package.py,
test_virtual_efficacy.py, test_virtual_efficacy_integration.py, test_warehouse_integrity.py

---

## RECOMMENDATIONS

### REMOVE (safe to delete):

| Target | Reason |
|--------|--------|
| `PX_Engine/operations/DoseOptimizer_Simple.py` | Banned by version lock |
| `PX_Engine/operations/PKPD_Simple.py` | Banned by version lock |
| `PX_Engine/operations/VirtualEfficacy_Simple.py` | Banned by version lock |
| `PX_Warehouse/Archive_Novel/` | Empty, orphaned |
| `PX_Warehouse/Archive_Primary/` | Empty, orphaned |
| `PX_Warehouse/Backup_Pre_Refinery/` | Empty, orphaned |
| `PX_Warehouse/CommercialAssets/` | Empty, orphaned |
| `PX_Warehouse/TrialSimulations/` | Empty, orphaned |
| `MANIFESTS/` directory | Fully orphaned |

### RELOCATE (through placement_gate):

| Files | From | To |
|-------|------|----|
| ~2,500 `PRV_REP_*.json` | `PX_Warehouse/` root | `Prv_Dossiers/<tier>/` via placement_gate |
| ~3,000 `PRV_NOV_*.json` | `PX_Warehouse/` root | `Novel_Dossiers/<tier>/` via placement_gate |
| 31 `TRIAL_SIMULATION_DOSSIER_*.json` | `PX_Warehouse/` root | `Calibration_Molecules/` |

### FIX (code changes required):

| Priority | Fix |
|----------|-----|
| P0 | Route all 23 ungated warehouse write paths through `placement_gate.place_prv_dossier()` |
| P0 | Fix schema mismatch: GradingEngine expects top-level `toxicity_index`, dossiers nest it inside `harm_energy`/`engines` |
| P0 | Fix `SMART_Antiviral_Fork.py` Rule 2 violation — replace `np.random` with deterministic physics |
| P1 | Add 15 excluded test files to `run_all_tests.py` harness |
| P1 | Fix `Evidence_Package.wrap_trial_simulation` default `output_dir` (relative path causes root-level dumps) |
| P2 | Remove broken import in `Refinement_Engine.py` (PX_Warehouse.RAG_Query_Engine) |
| P2 | Review 82 dead code files for archival to `99_LEGACY_CODE/` or deletion |

### KEEP (connected, functional):

All files not listed above — approximately 198 connected Python files across 12 active departments plus supporting directories (governance, Nipah_Analysis, tests).

---

## PHASE 2 ACTION SEQUENCE

1. **Fix the schema mismatch** (toxicity_index field location) — this unblocks finalization
2. **Route the 23 ungated write paths** through placement_gate — this prevents future root dumps
3. **Run placement_gate on 5,631 root files** — this clears the warehouse root
4. **Fix Rule 2 violation** in SMART_Antiviral_Fork.py
5. **Add 15 tests to harness** and verify they pass
6. **Delete banned _Simple files** and empty warehouse subdirs
7. **Archive dead code** (82 files) to 99_LEGACY_CODE or delete
8. **Verify with `just test` and `just govern`**
