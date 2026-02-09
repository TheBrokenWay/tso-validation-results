# PX Platform — IDE Deep-Indexing Reference

**Purpose:** Enable the IDE to build dependency graphs, execution DAGs, state machines, governance lattice, and scientific-integrity boundaries. Answers topology, interface contracts, dependency rules, data flow, lifecycle, determinism, failure propagation, governance, and self-awareness.

---

## 1. Repository Topology & Ownership

### 1.1 Top-Level Execution Roots (Entry Points)

| Entry Point | Path | Role |
|-------------|------|------|
| TSO Validator | `TSO_Validator/run_validation.py` | Physics validation pipeline; produces run_YYYYMMDD_HHMMSS |
| Nipah Analysis | `Nipah_Analysis/run_nipah_analysis.py` | Strain-aware validation → normalization → analytics |
| Nipah (PowerShell) | `Nipah_Analysis/run_nipah_analysis.ps1` | Wraps Python; forwards --allow-multi-strain |
| Universal Pipeline | `PX_Executive/UniversalPipelineRunner.py` | Runs molecule through GradingEngine + dossier |
| Batch Pipeline | `PX_Executive/batch/universal_pipeline_batch.py`, `pipeline_batch_runner.py` | Batch discovery and run |
| Live Orchestrator | `PX_Executive/orchestrators/PX_Live_Orchestrator_v2.py` (main) | Orchestration entry |
| PRV Master | `PX_Executive/PRV_Master_Pipeline_V2.py`, `PRV_Master_Pipeline.py` | PRV/master engine |
| Harvest Failed | `harvest_failed_trials.py` (root) | Trial harvest script |
| Predator X v3 | `predator_x_v3_stratification.py`, `predator_x_v3_monetization.py` (root) | Stratification/monetization |

### 1.2 Pure Libraries vs Active Runtime

| Type | Directories |
|------|-------------|
| **Pure libraries** (no direct entry; used by others) | `PX_System/foundation` (api, Disease_Constraint_Model, Evidence_Package, ZeusLaws, etc.), `PX_Constitution` (Block_Universe, Virtual_Machine), `PX_Engine/operations` (OBE, OCE, OLE, OME, OPE, OSE, TrialEngine, PKPD, ADMET, GradingEngine, etc.) |
| **Active runtime** (entry points or orchestration) | `PX_Executive`, `TSO_Validator`, `Nipah_Analysis`, `PX_Engine` (Engine_Orchestrator, Vector_Core, Metabolism, Trajectory_Predictor), `PX_Laboratory` (Simulation_Engine) |

### 1.3 Runtime-Critical vs Development-Only vs Audit-Only

| Classification | Modules |
|----------------|---------|
| **Runtime-critical** | PX_System (api, Evidence_Package, Disease_Constraint_Model), PX_Engine (Vector_Core, Engine_Orchestrator, operations), PX_Executive (UniversalPipelineRunner, Sovereign_Commercial_Pipeline, orchestrators), PX_Constitution (ZeusLaws, Block_Universe), PX_Warehouse (WorldLine_Database) |
| **Development-only** | PX_Executive/demos, PX_Executive/validate_v2_release.py, PX_Validation manual_tests, TSO_Validator/tests, Nipah_Analysis/tests |
| **Audit-only** | PX_Audit (Continuous_Monitor, Drift_Monitor, Final_System_Seal, Protocol_Zero, Mural_Network, etc.), TSO_Validator/audit (traceability, reproducibility, drift, ethical_guard, history, provenance) |

### 1.4 Immutable Governance vs Mutable Research

| Layer | Directories | Notes |
|-------|-------------|--------|
| **Immutable governance** | `.cursor/rules` (00_constitution, 01_research_audit, 02_physics_integrity, 03_anti_sycophancy), `PX_Constitution`, `PX_System/foundation/ZeusLaws.py`, `TSO_Validator/audit/ethical_guard.py`, `MANIFESTS` | Constitutional rules; ethical lock; manifest constraints |
| **Mutable research** | `PX_Warehouse`, `PX_Discovery`, `PX_Executive/generators`, `Nipah_Analysis/data`, `TSO_Validator/data` | Data, dossiers, worldlines, discovery outputs |

---

## 2. Module Interface Contracts

### 2.1 PX_System

| Question | Answer |
|----------|--------|
| **Public interfaces** | `api.py`: run_boot_sequence, verify_integrity, adjudicate_package, enforce_law_l1, append_audit_record, get_system_status. `Disease_Constraint_Model`: DiseaseConstraintModel, get_disease_constraints, create_chagas_dcm, create_nipah_malaysia_dcm, create_nipah_bangladesh_dcm, get_nipah_dcm. `Evidence_Package`: generate_dossier, validate_dossier, wrap_trial_simulation. `ZeusLaws`: check_constitutional. `Sovereign_Log_Chain`: append, get_chain_hash. `Data_Sources`: get_data_sources. `Emergency_Stop`: EmergencyStop. `Novelty_Budget_Engine`: NoveltyBudgetEngine. |
| **DCM representation** | Dict[str, Any] from create_*_dcm() / get_disease_constraints(disease_name). Keys: disease_id, disease_name, pathogen, strain (Nipah), target_proteins, cfr_range, ic50_max_um, selectivity_index_min, etc. |
| **DCM construction/validation/serialization** | Constructed in Disease_Constraint_Model.py (create_chagas_dcm, create_nipah_malaysia_dcm, create_nipah_bangladesh_dcm). Validated via get_disease_constraints and check_constitutional (ZeusLaws). Serialization: DCMs are dicts; persisted in MANIFESTS and Nipah_Analysis config. |
| **Formal I/O contract** | generate_dossier(candidate, engine_results, zeus_verdict?, context?) → dossier dict. adjudicate_package(evidence) → {verdict, input_hash}. check_constitutional(operation, payload) → {authorized, operation, rationale, timestamp}. |

### 2.2 PX_Engine

| Question | Answer |
|----------|--------|
| **Computational primitives** | Vector_Core: execute(p_vector) → state (amplitude, authorized, trace_id). Metabolism, Trajectory_Predictor. operations: OBE, OCE, OLE, OME, OPE, OSE, TrialEngine, PKPD, ADMET, GradingEngine, DoseOptimizer, VirtualEfficacy. |
| **Solver(s)** | Vector physics (Vector_Core); PK/PD models (PKPD, DoseOptimizer); grading (GradingEngine). No external solver package required; NumPy/SciPy. |
| **Constraint violations** | ZeusLaws.check_constitutional: toxicity_index >= 0.0210 or harm_energy >= 0.0210 → authorized False. Evidence_Package computes harm_energy from engine results; dossier carries constitutional evidence. |
| **Determinism** | TSO_Validator sets np.random.seed(RNG_SEED); PX_Engine/Executive do not globally set seed in all paths (see §6). |

### 2.3 PX_Executive

| Question | Answer |
|----------|--------|
| **Workflow mechanism** | Pipeline/DAG-style: discover assets → run pipeline per molecule (GradingEngine, dossier generation) → optional TSO gate. UniversalPipelineRunner.run() drives flow. PX_Live_Orchestrator_v2, PRV_Master_Pipeline_V2 class-based. |
| **Retry/recovery** | PX_System.foundation.integrations.retry: retry_on_transient, retry_transient_http. No global FSM retry in Executive. |
| **Response to validation failures** | Pipeline returns result dict; caller interprets. Nipah: validation failure → exit 1, no partial runs. TSO: TSO_FAILED → run_summary status; Nipah --tso-validate fails run if TSO_FAILED. |
| **Response to engine divergence / audit violations** | Constitutional check (ZeusLaws) before/after; audit via append_audit_record, Sovereign_Log_Chain. No explicit “interrupt execution” from PX_Audit in all code paths. |

### 2.4 PX_Validation

| Question | Answer |
|----------|--------|
| **Authoritative schemas** | PX_Validation/system_manifest.json (structure snapshot). GradingSchema_Discovery.json in PX_Engine/operations. Nipah: config/nipah_schema.json, nipah_ontology.json. TSO: physics/classification (report_schema, RUN_SUMMARY_SCHEMA_VERSION). |
| **Ontology derivation** | Nipah: config/nipah_ontology.json (allowed_strains, forbidden_families, cfr_bounds_by_strain). TSO: config/TSO_CONSTANTS.json, TSO_SECOND_PASS_SPEC.json. PX_System: DCMs encode disease/constraint ontology. |
| **Blocking vs advisory** | Nipah: validation is blocking (exit 1 on any reject). TSO: run produces PASS/FAIL/INDETERMINATE; --tso-validate in Nipah makes TSO blocking. PX_Validation tests are advisory (CI). |
| **Validation state → PX_System** | Nipah_Analysis feeds strain-specific DCMs via get_disease_constraints("NiV_Malaysia") etc. TSO outputs (run_summary, falsification) are not automatically consumed by PX_System; optional gate only. |

### 2.5 PX_Audit

| Question | Answer |
|----------|--------|
| **Events captured** | audit_trail.jsonl (api.append_audit_record). Sovereign_Log_Chain append. health_monitor_log.jsonl (Continuous_Monitor.log_event). Mural_Network: update_node, update_edge, get_mural. WorldLine_Database: record_materialization (worldline artifacts). TSO: run_summary, provenance, historical_context, run.log. |
| **Immutability** | Append-only audit_trail.jsonl and log_event writes. TSO run dirs immutable (history reads only; writes only to current run). Nipah validation/checksums.json updated per run (mutable registry). |
| **Passive vs active** | Predominantly passive (logging, trace). check_constitutional can deny operation (active). TSO ethical_guard: pre/post gate can raise (active). |
| **Can interrupt execution** | Yes where constitutional/ethical check raises or returns authorized=False and caller treats as fatal. PX_Audit itself does not globally halt; callers decide. |

### 2.6 PX_Security

| Question | Answer |
|----------|--------|
| **Security policies** | RedSurface: decoy layer (generate_decoy_response). PX_System.foundation.integrations: check_external_access_allowed, get_rate_limit, online_sync_enabled (net_policy); validate_smiles_string, sanitize_smiles (smiles_security). |
| **Threat model** | External access and rate limits; SMILES injection/validation; decoy responses for adversarial probing. |
| **Sandboxing** | No process-level sandbox documented; security is policy checks and validated inputs. |

### 2.7 PX_Warehouse

| Question | Answer |
|----------|--------|
| **Data models** | WorldLine_Database: worldline artifacts (header, physics_snapshot, physical_realization). Operational_Data/WorldLine_Records; commercial dossiers; research assets. |
| **Lifecycle state** | Persisted per task_id (.worldline files). Batch/sort scripts (batch_run, sort_live_results). |
| **Append-only vs mutable** | New worldlines written; existing files overwritten on update. Not strictly append-only. |

### 2.8 PX_Laboratory

| Question | Answer |
|----------|--------|
| **Normalization location** | Simulation_Engine: PK simulation (one-compartment). Nipah: normalize_data.py (units → copies/mL, PII anonymization). TSO: normalization/tso_bridge.py (delta_energy, delta_entropy, control_parameter, time). |
| **Feature normalization guarantee** | Nipah: validation gate before normalization; schema + ontology. TSO: tso_bridge produces canonical schema; no interpolation. |
| **Reversible transforms** | PII anonymization (Nipah) is one-way hash. Unit conversion (viral_load → copies_per_mL) is reversible if formula stored. |

### 2.9 PX_Discovery

| Question | Answer |
|----------|--------|
| **Hypothesis generation** | candidate_discovery_engine.discover_candidates. AutonomousResearchController (start_research_cycle) stubbed. |
| **Feeds constraints upstream** | Discovery outputs feed into pipeline/dossier; constraints come from PX_System DCMs and MANIFESTS. |
| **Scientific novelty** | Novelty_Budget_Engine in PX_System; discovery generators (e.g. SMART_Antiviral_Fork) use variance/coherence. |

### 2.10 PX_Constitution

| Question | Answer |
|----------|--------|
| **Invariants** | ZeusLaws (PX_System.foundation.ZeusLaws): toxicity_index < 0.0210, harm_energy < 0.0210. L1: internal physics snapshot primary; external cannot override. PX_Constitution: Block_Universe (35D manifold, coherence = 1/(1+distance^2)), Virtual_Machine (get_vm_fingerprint). |
| **Violation detection** | check_constitutional(operation, payload) returns authorized False + rationale. Evidence_Package computes harm_energy from engines. |
| **Passive vs active** | Active: check_constitutional used by callers; can block. enforce_law_l1 used in dossier generation. |
| **Centralized vs distributed** | ZeusLaws in PX_System.foundation; Block_Universe, Virtual_Machine in PX_Constitution. Constitutional checks invoked from Evidence_Package, Sovereign_Commercial_Pipeline, api. Engine_Orchestrator uses get_vm_fingerprint from PX_Constitution. |

### 2.11 Nipah_Analysis

| Question | Answer |
|----------|--------|
| **Input → output contract** | Input: data/raw/*.json (strain, records with outbreak_year, cfr). Output: results/summary.json (mode, gates, normalized_files), results/validation_report.json (provenance, constraint_drift), validation/checksums.json, data/normalized/*.json, results/cross_strain_comparative.json, constraint_divergence.json, evolutionary_drift.json, intervention_sensitivity.json. Exit 0 only if all validation passes. |
| **Biological logic allowed** | Strain-specific CFR bounds, temporal bounds, forbidden_family checks; DCM wiring to PX_System. No protocols, synthesis, or experimental procedures (ethical_guard in TSO; Nipah has no procedural output). |
| **Strictly forbidden** | No biological protocols, synthesis routes, or experimental instructions in outputs (same principle as TSO ethical lock). Forbidden families in raw data (Orthomyxoviridae, etc.) → reject. |
| **Injection into PX_System** | Via get_disease_constraints("NiV_Malaysia") / ("NiV_Bangladesh") and get_nipah_dcm(strain). Nipah does not write into PX_System codebase; PX_System imports and exposes Nipah DCMs. |

### 2.12 TSO_Validator

| Question | Answer |
|----------|--------|
| **Semantic checks** | Physics: scaling (ΔE ∝ ΔS), sigmoid decoherence, divergence (or second-pass: crossover, stretched exponential, tanh, alpha, topology, avalanche). Classification: PASS/FAIL/INDETERMINATE. Ethical guard: no procedural/protocol text. |
| **Failure propagation** | Any test FAIL → falsification.status TSO_FAILED. run_validation exits 0; status in run_summary. Nipah --tso-validate: if TSO_FAILED → Nipah run exit 1. |
| **Trust metrics** | run_summary.json (canonical), provenance.json (git, config hash, dataset hashes), historical_context.json (drift, lineage, falsification_log), confidence per test, falsification matrix. |

### 2.13 PX_STATE

| Question | Answer |
|----------|--------|
| **State machine** | current_state.json: status (e.g. ACTIVE), timestamp, batch_id, total_molecules, total_survivors. Boot sequence checks existence for SYSTEM_READY_V2. |
| **Atomic transitions** | Single JSON file; no formal atomicity (file overwrite). |
| **Replay** | State is summary only; full replay would require audit logs + deterministic engine (see §6). |

### 2.14 PX_LOGS

| Question | Answer |
|----------|--------|
| **Event taxonomy** | Batch_*_stdout.log, Batch_*_stderr.log; health_monitor_log.jsonl; audit_trail.jsonl; constitutional_trace.log; TSO run.log. |
| **Failure taxonomy** | Per-module (e.g. TSO classification reason_code; Nipah RawValidationError code). No single global failure code enum. |
| **Cryptographic verification** | Sovereign_Log_Chain.get_chain_hash(). TSO provenance.json (hashes). No global log signing. |

---

## 3. Dependency Graph Validation

### 3.1 Allowed Dependency Directions

- **PX_Validation** → PX_System (tests, evidence packages).
- **PX_System** → PX_Engine (orchestration, dossier from engine results).
- **PX_Engine** → PX_Audit (Mural_Network), PX_Constitution (Virtual_Machine fingerprint).
- **PX_Executive** → PX_Engine, PX_System, PX_Audit, PX_Warehouse.
- **PX_Warehouse** ↔ PX_Laboratory (data flow; Warehouse stores; Lab can read/write).
- **Nipah_Analysis** → PX_System (DCMs only via get_disease_constraints; no Nipah → PX_Engine).
- **TSO_Validator** → standalone (no PX_* imports); optional second gate for Nipah (Nipah runs TSO subprocess).
- **PX_Constitution** → used by PX_System, PX_Engine, PX_Executive (governs all).

### 3.2 Forbidden Cross-Layer Calls

- PX_Engine must not depend on PX_Executive.
- PX_Constitution must not depend on PX_Warehouse or PX_Discovery.
- TSO_Validator must not import PX_* or Olympus (standalone).
- Nipah_Analysis must not contain biological protocols or synthesis logic.

### 3.3 Cyclic Imports / Layer Violations / Backward Dependencies

- **Potential cycles:** PX_Executive imports PX_Engine, PX_System; PX_System.foundation.Evidence_Package may import Sovereign_Log_Chain; PX_Engine.Engine_Orchestrator imports PX_Audit.Mural_Network, PX_Constitution. Audit/Constitution are intended downstream of Engine/Executive; no deep cycle identified in the explored files.
- **Backward dependency:** PX_Engine → PX_Audit (Mural_Network) is Engine → Audit, which is allowed (audit records engine activity).

---

## 4. Data Flow Integrity

| Question | Answer |
|----------|--------|
| **Single canonical pipeline** | No single global DAG. Pipelines: (1) Nipah: raw → validate → normalize → analytics. (2) TSO: raw → normalize → physics tests → falsification → run_summary. (3) Executive: warehouse/assets → GradingEngine → dossier → audit. |
| **Raw data outside PX_Validation** | Nipah raw validated in Nipah_Analysis; TSO raw in TSO_Validator. PX_Validation (tests) operate on fixtures/snapshots. Executive can read warehouse assets without going through TSO/Nipah. |
| **Normalized data bypass** | Nipah: analytics consume normalized only. TSO: physics tests consume normalized only. Executive pipeline uses engine outputs (not necessarily “normalized” in the same sense). |
| **Biological logic duplication** | DCMs centralized in PX_System. Nipah ontology (CFR bounds, strains) in Nipah config; Nipah DCMs in PX_System. Avoid duplicating constraint bounds in multiple repos. |

---

## 5. Lifecycle State Model

| Concept | Details |
|---------|---------|
| **Execution states** | BOOT (INITIALIZING_CORE / SYSTEM_READY_V2), ACTIVE (current_state.json), run_YYYYMMDD_HHMMSS (TSO), validation PASS/FAIL (Nipah). |
| **Transitions** | Boot → ready; batch run → state update; TSO run → new run dir; Nipah run → summary + registry. |
| **Who can mutate state** | PX_STATE: boot/runner scripts. PX_Warehouse: pipeline writes. TSO: run_* dirs. Nipah: checksums.json, results. PX_Audit: append-only logs. |
| **Rollback** | No formal rollback. Restore from backup or recompute from raw + deterministic run. |

---

## 6. Execution Determinism

| Question | Answer |
|----------|--------|
| **Bit-for-bit reproducibility** | TSO_Validator: RNG_SEED set at start; config + data hashes in provenance. Nipah: deterministic validation and hashes. PX_Engine/Executive: SMART_Antiviral_Fork and some audit code use np.random without global seed in all entry points. |
| **Nondeterministic elements** | np.random in PX_Executive/generators, PX_Audit/Manifold_Normalizer unless seeded by caller. Time stamps (datetime.now()) in logs and state. |
| **Random seeds** | TSO: config.yaml rng_seed, _set_global_rng_seed(seed). Nipah: no RNG in validation/normalization. |
| **Solver nondeterminism** | SciPy/NumPy fits; TSO seeds numpy. PX_Engine uses NumPy; no explicit seed in Engine_Orchestrator. |

---

## 7. Failure Propagation Semantics

| Failure | Recoverable? | Fatal? | Logged? | State-tracked? |
|---------|-------------|--------|---------|----------------|
| **Validation fails (Nipah)** | No | Yes (exit 1) | stderr | No partial run |
| **Normalization fails (Nipah)** | No | Yes (exit 1) | stderr | — |
| **Analytics fail (Nipah)** | Yes | No (warning) | stderr | summary still written |
| **TSO run** | — | No (exit 0) | run.log, run_summary | status in run_summary |
| **Nipah --tso-validate + TSO_FAILED** | No | Yes (exit 1) | stderr | — |
| **Engine divergence** | Caller-dependent | If constitutional fail | Audit/dossier | zeus_verdict |
| **Audit violation** | If check_constitutional returns False | Caller can treat as fatal | rationale in result | — |
| **Governance violated** | No | When enforced | rationale / ethical_guard | — |

---

## 8. Governance Enforcement

| Question | Answer |
|----------|--------|
| **Where enforced** | ZeusLaws.check_constitutional (PX_System); Evidence_Package harm_energy + zeus; enforce_law_l1 in api and Evidence_Package. TSO ethical_guard at startup and before reports. Nipah: validation gate (schema, strain, CFR, forbidden family). |
| **Bypass possible?** | Only if callers skip check_constitutional or ethical_guard. No hard process boundary; governance is invoked at defined call sites. |
| **Passive vs active** | Active where checks return/raise and caller acts. |
| **If violated** | check_constitutional: authorized False, rationale. ethical_guard: EthicalGuardViolation exception. Nipah: RawValidationError, exit 1. |

---

## 9. Scientific Integrity Controls

| Question | Answer |
|----------|--------|
| **Domain truth encoded** | PX_System DCMs (Chagas, Nipah Malaysia/Bangladesh). Nipah config (ontology, schema). TSO config (constants, second-pass spec). MANIFESTS. |
| **Duplication prevented** | Single DCM source (PX_System); Nipah ontology in one config; TSO constants in config. |
| **Drift detected** | TSO: historical_context.json (drift vs prior run). Nipah: constraint_drift in report; evolutionary_drift.json; checksums.json for file drift. PX_Audit Drift_Monitor. |
| **Corrections propagated** | Manual config/DCM updates; no automatic propagation. |

---

## 10. Platform Self-Awareness

| Question | Answer |
|----------|--------|
| **Introspect dependencies?** | No single tool. system_manifest.json (PX_Validation) is a structure snapshot. Dependency graph can be built by static analysis (imports). |
| **Introspect constraint model?** | DCMs are dicts; get_disease_constraints, create_*_dcm expose them. Nipah manifest + ontology in config. |
| **Introspect drift?** | Yes: TSO historical_context; Nipah evolutionary_drift, checksums registry; PX_Audit Drift_Monitor. |
| **Detect architectural violations?** | Not automated. Allowed/forbidden dependency directions (§3) can be checked by static analysis or custom lint. |

---

## Appendix A: Dependency Diagram (Allowed Flow)

```
PX_Constitution (ZeusLaws, Block_Universe, Virtual_Machine)
         ↑
PX_System (api, DCM, Evidence_Package, Sovereign_Log_Chain)
         ↑
PX_Engine (Vector_Core, operations, Engine_Orchestrator) → PX_Audit (Mural)
         ↑
PX_Executive (UniversalPipelineRunner, orchestrators, batch)
         ↑
PX_Validation (tests) ──→ PX_System
PX_Warehouse ←→ PX_Laboratory
Nipah_Analysis ──→ PX_System (DCM only)
TSO_Validator (standalone)
Optional: Nipah --tso-validate → subprocess TSO_Validator
```

## Appendix B: Execution DAG (Simplified)

- **Nipah:** Load manifest → Validate raw → Normalize → [Optional TSO] → Analytics → Summary/Reports.
- **TSO:** Ingest → Normalize → Physics tests → Falsification → Provenance/History → run_summary + reports.
- **Executive:** Discover assets → For each: GradingEngine → Dossier (Evidence_Package) → Zeus check → Audit/Warehouse.

## Appendix C: Governance Lattice

- **Constitutional:** ZeusLaws (L10 harm, L11 toxicity), L1 (internal snapshot primary).
- **Ethical:** TSO ethical_guard (no procedural/protocol output); Nipah no forbidden families/protocols.
- **Validation:** Nipah schema + ontology; TSO physics + classification; PX_Validation tests.

## Appendix D: Failure Propagation Tree

- **Nipah validation fail** → exit 1 (no normalization/analytics).
- **Nipah normalization fail** → exit 1.
- **Nipah analytics fail** → warning; summary written.
- **TSO test FAIL** → falsification TSO_FAILED; run_summary.status.
- **Nipah --tso-validate and TSO_FAILED** → Nipah exit 1.
- **check_constitutional authorized False** → caller can abort or log.
- **ethical_guard raise** → exception; abort if not caught.

---

*Document generated for IDE deep-indexing. Update when topology or contracts change.*
