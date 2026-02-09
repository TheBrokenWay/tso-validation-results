# PX Enterprise Layout — IDE-Ready

**Enterprise Edition — Final, Locked Directory Structure**

This file defines the **authoritative, frozen** directory structure for the PX Enterprise Edition.  
Any directory or file **not listed below** must be placed under `archive/` and is **not part of the enterprise runtime**.

**The IDE must enforce this layout.**

---

## 1. Constitutional & Governance Layer (MANDATORY)

| Path | Purpose |
|------|---------|
| `PX_Constitution/` | Constitutional invariants, VM fingerprint, mandatory pre-execution tests, gate stamp signing |
| `governance/` | Poison-pill registry, gate enforcement, layer monotonicity, coverage, certificate |
| `PX_PLATFORM_FULL_DEEP_INDEX_AND_GOVERNANCE_MAP.md` | Platform index and governance map |
| `PX_ENTERPRISE_EDITION.md` | Enterprise edition specification |
| `PX_ENTERPRISE_LAYOUT.md` | This file — frozen directory layout |
| `PX_WAREHOUSE_ENTERPRISE_MODEL.md` | Warehouse model: ONE CommercialAssets, ONE WorldLines, ONE Archive, ONE terminology |
| `PX_WAREHOUSE_MIGRATION_INSTRUCTION.md` | Authoritative migration directive: mandatory steps to align PX_Warehouse with the enterprise model; IDE must treat as enforceable |
| `PX_WAREHOUSE_IDE_ENFORCEMENT_SPEC.md` | IDE execution spec: folder/file enforcement, lifecycle, tiering, cleanup, script updates, verification, lock |
| `WAREHOUSE_DRIFT_REPORT_TEMPLATE.md` | Template for IDE-generated drift reports when verification fails |
| `WAREHOUSE_DRIFT_REPORT_EXECUTION_LOGIC.md` | Verification algorithm: load canonical model, scan warehouse, compare, generate actions, re-verify, certification |
| `bootstrap_canonical_warehouse.py` | Idempotent script: creates canonical PX_Warehouse folder structure (run from repo root) |
| `PX_WAREHOUSE_REDUNDANCY_SCANNER_SPEC.md` | Safe-mode redundancy scanner: Discovery (read-only) → Simulation (propose actions) → Execution (validated only); classification ACTIVE/REFERENCED/LEGACY/DUPLICATE/UNUSED/UNKNOWN |
| `PX_WAREHOUSE_FULL_REBUILD_DIRECTIVE.md` | Full rebuild command: recursive scan, reclassification, renaming, legacy removal, **system-wide reference update**; Simulation (dry run) → Execution → Post-rebuild verification; supersedes all legacy warehouse structures and references |
| `DOSSIER_COMMERCIALASSETS_FILING_RULES.md` | **ENFORCE_FINAL_PRODUCT_FILING:** Only assets that pass all scientific, legal, commercial, and governance checks may enter Dossier_Final; all commercial-grade outputs must be routed to correct tier (Diamond/Gold/Silver/Bronze/Executive_Summary/Audit_Trail); fail-closed; canonical path CommercialAssets/\<Tier\>/\<AssetID\>/ |
| `PX_Constitution/Filing_Rules.py` | Reference implementation: `may_enter_dossier_final`, `get_tier_for_asset`, `validate_filing_path`; call before any write into CommercialAssets or Dossier_Final |

---

## 2. Core Runtime Organs (MANDATORY)

| Path | Purpose |
|------|---------|
| `PX_System/` | Foundation: ZeusLaws, Evidence_Package, Disease_Constraint_Model, Sovereign_Log_Chain, api |
| `PX_Engine/` | Vector_Core, Metabolism, Trajectory_Predictor, operations (OPE, ADMET, TrialEngine, etc.) |
| `PX_Executive/` | Gold_Rush_Miner, Sovereign_Commercial_Pipeline, orchestrators, GAIP_Gateway, PX_Legal_Check, Acquisition (Harvest/Bridge) |
| `PX_Validation/` | Test suite, benchmarks, warehouse integrity, system inventory |
| `PX_Audit/` | Mural network, pipeline runs, sovereign log chain, reports |
| `PX_Security/` | RedSurface, circuit breaker, PredatorImmune_Block |
| `PX_Warehouse/` | **ONE** warehouse: CommercialAssets, WorldLines, Archive, TrialSimulations, Operations; governed by **PX_WAREHOUSE_ENTERPRISE_MODEL.md** |
| `PX_Laboratory/` | Simulation_Engine, Manufacturing_Manifest, Synthetic_Expansion |
| `PX_Discovery/` | Candidate discovery engine (stub), adapters |
| `PX_STATE/` | Current state (e.g. current_state.json) |
| `PX_LOGS/` | Runtime logs; `PX_LOGS/archive/` for run-specific logs (reprocess, etc.) |

---

## 3. Domain Slice (REFERENCE IMPLEMENTATION)

| Path | Purpose |
|------|---------|
| `Nipah_Analysis/` | Strain-aware Nipah validation, adapters, analytics, poison-pill payloads |

---

## 4. Active Constitutional Second Gate (MANDATORY)

| Path | Purpose |
|------|---------|
| `TSO_Validator/` | TSO physics validation — second gate; independent of PX runtime |

---

## 5. Evidence & Trial Layer (MANDATORY)

These modules are **mandatory** for the enterprise runtime; they live in the paths below:

| Module | Actual Path |
|--------|-------------|
| Evidence_Package | `PX_System/foundation/Evidence_Package.py` |
| TrialEngine | `PX_Engine/operations/TrialEngine.py` |
| VectorCore | `PX_Engine/Vector_Core.py` |

---

## 6. Orchestrator Layer (MANDATORY)

| Path | Purpose |
|------|---------|
| `PX_Executive/orchestrators/PX_Live_Orchestrator_v2.py` | Predator X v2.0-CORE orchestrator (corrected pipeline order) |
| `run_e2e_layers.py` | End-to-end layer verification script (root) |

---

## 7. Enterprise Test Suite (MANDATORY)

| Path | Purpose |
|------|---------|
| `run_poison_pill_test.py` | Governance poison-pill canary (CFR envelope) — root (constitutional) |
| `run_temporal_paradox_test.py` | Governance poison-pill canary (temporal paradox) — root (constitutional) |
| `run_e2e_layers.py` | End-to-end layer verification — root (constitutional) |
| `tests/manual/` | Manual/debug scripts (test_pipeline_single, test_subprocess) |
| `PX_Validation/tests/` | Engine, trial, warehouse, orchestrator tests |

---

## 8. Optional Runtime Dependencies (KEEP ONLY IF IMPORTED)

| Path | Purpose |
|------|---------|
| `MANIFESTS/` | Pipeline/manifest config (e.g. SMART_Antiviral_Fork) |
| `utils/` | Shared utilities — retain only if referenced by enterprise modules |
| `common/` | Common code — retain only if referenced by enterprise modules |

The IDE must verify that these folders are **only retained if referenced** by enterprise modules.

---

## 9. Archive (EVERYTHING NOT LISTED ABOVE)

All legacy, deprecated, experimental, or unused modules must be moved under:

```
archive/
    legacy_code/
    legacy_scripts/
    experiments/
    deprecated/
```

**Examples of items that must be archived** (when not part of Sections 1–8):

- `Olympus_Research/`
- `tso-validation-results/`
- `archive/v3_core/` (predator_x_v3 scripts — preserved for reference; enterprise uses v2 orchestrator)
- `reprocess_*` (scripts not in PX_Warehouse/Operations)
- `run_pipeline*` (duplicate/legacy runners)
- `.pytest_cache/`
- `__pycache__/` (build artifacts; may be excluded from layout enforcement)

**Clean Desk (Filing Protocol):** Harvest/Bridge → `PX_Executive/Acquisition/`; debug scripts → `tests/manual/`; v3 root scripts → `archive/v3_core/`. Root `run_poison_pill_test.py`, `run_temporal_paradox_test.py`, `run_e2e_layers.py` remain at root per Constitution.

`.cursor/` is IDE-specific and may be excluded from archive enforcement.

---

## 10. IDE Enforcement Rules

The IDE **must** enforce:

1. **No new top-level directories** unless added to this file.
2. **No runtime code** outside the directories listed in Sections 1–8.
3. **All unlisted directories** must be moved to `archive/` (or documented as optional).
4. **Architecture drift must be zero** (Type-A..E).
5. **Poison-pill gate must pass** before any pipeline execution.
6. **Governance Readiness Certificate** must be valid.
7. **TSO_Validator** must remain available as a second gate.
8. **PX_Constitution invariants** must be loaded and enforced.

---

## 11. Enterprise Edition Lock Clause

**This layout is frozen.**

Any modification requires:

1. Constitutional amendment  
2. Update to this file  
3. Passing all governance tests  
4. Regeneration of the Governance Readiness Certificate  
5. Zero architecture drift  
