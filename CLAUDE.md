# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Predator X is a deterministic, constitutionally governed pharmaceutical compound discovery platform. It generates and evaluates novel/repurposed molecules through a forward-only pipeline: **Feed -> PRV -> Trial -> Finalization**. Every stage is gated by governance (ZeusLaws), and all outputs maintain full lineage from WorldLine physics snapshots to final dossiers.

The environment is sealed via Nix flakes (`flake.nix`) with `direnv` (`.envrc` contains `use flake`). Tools pinned: Bazel 7, Python 3.11, Rye, Ruff, Pyright, DVC, Hydra. Dependencies: RDKit (>=2023.9.1), NumPy (>=1.24.0,<2.0), requests (>=2.31.0). Constitutional modules (ZeusLaws, Sovereign_Log_Chain) use only Python stdlib.

## Commands

```bash
just check         # Lint + type-check (ruff check . && pyright)
just test          # Full validation suite (python PX_Validation/tests/run_all_tests.py)
just govern        # Governance stress test (python run_e2e_layers.py)

# Run a single test file directly
python PX_Validation/tests/test_admet_engine.py

# Engine loops (consolidated orchestrators)
just feed          # Genesis feed → px_feed.py --mode novel
just feed-rep      # Repurposed feed → px_feed.py --mode repurpose
just novel         # 12-engine PRV pipeline → px_prv.py --type novel
just repurpose     # 12-engine PRV pipeline → px_prv.py --type repurpose
just prv-all       # Both novel + repurposed → px_prv.py --type all
just finalize      # Finalization pipeline → px_finalize.py
just cycle         # Full cycle test: Feed → PRV → Finalize

# Batch operations
just reprocess-lifecycle            # Push all unfinalized dossiers through grading + finalization
just reprocess-lifecycle-validated  # Same with validation before/after (recommended before version lock)

# Consolidated orchestrator direct usage
python PX_Executive/px_feed.py --mode novel --count 20 --interval 0      # Single batch of 20
python PX_Executive/px_feed.py --mode novel --high-chaos                  # High-entropy mode
python PX_Executive/px_prv.py --type novel --limit 100                    # Process up to 100 novel
python PX_Executive/px_prv.py --type all --disease nipah                  # Filter by disease
python PX_Executive/px_finalize.py --dry-run                              # Report only
python PX_Executive/px_finalize.py --limit 100 --reprocess -v             # Reprocess with verbose

# DVC lineage
just lineage       # dvc status — tracked Physics Maps / worldlines
```

**After any edit to engines or orchestrators, both `just test` and `just govern` must pass.** If either fails, reject the edit.

## Architecture

### Pipeline Data Flow (forward-only, no stage regenerates data)

```
Feed Stage (Genesis_Engine, Vector_Core, Metabolism, Trajectory_Predictor)
  -> Queue (PX_Warehouse/Feeder/prv_24h_queue.json)
  -> PRV Stage (OBE, OCE, OLE, OME, OPE, OSE -> ADMET -> PKPD, DoseOptimizer_v2, VirtualEfficacyAnalytics -> GradingEngine)
  -> Trial Stage (TrialEngine, Simulation_Engine, Manufacturing_Manifest)
  -> Finalization (Evidence_Package, GradingEngine, ZeusLaws.run_zeus_gate)
  -> PX_Warehouse/Finalized_Dossiers/<DIAMOND|GOLD|SILVER>/
```

### Directory Structure

| Directory | Role |
|-----------|------|
| `PX_Executive/` | **Consolidated orchestrators**: `px_feed.py` (novel+repurposed feed), `px_prv.py` (12-engine PRV pipeline), `px_finalize.py` (finalization). Also: `run_one_cycle_test.py`, `reprocess_warehouse_lifecycle.py` |
| `PX_Executive/orchestrators/` | `PX_Live_Orchestrator_v2.py` (legacy live orchestrator) |
| `PX_Engine/` | Core physics engines: `Genesis_Engine.py`, `Vector_Core.py`, `Metabolism.py`, `Trajectory_Predictor.py`, `Block_Orchestrator.py`, `Engine_Orchestrator.py` |
| `PX_Engine/operations/` | Operational engines: OBE, OCE, OLE, OME, OPE, OSE, `ADMET.py`, `PKPD.py`, `TrialEngine.py`, `DoseOptimizer_v2.py`, `VirtualEfficacyAnalytics.py`, `GradingEngine.py`, `GradingSchema_Discovery.json` |
| `PX_Laboratory/` | `Simulation_Engine.py` (1-compartment PK sim), `Manufacturing_Manifest.py`, `Synthetic_Expansion.py` |
| `PX_System/foundation/` | `ZeusLaws.py` (governance gate), `Sovereign_Log_Chain.py` (immutable audit trail), `Evidence_Package.py` (FDA 2026-compliant dossier generation), `core.py` (QSAR constants), `Disease_Constraint_Model.py` |
| `PX_System/foundation/integrations/` | `smiles_security.py` (SMILES injection validation), `net_policy.py`, `retry.py` |
| `PX_Warehouse/` | **PROTECTED — never modify directly.** Contains `Finalization_Pipeline.py`, `warehouse_layout.py` (re-exports from `placement_gate/Staging/`), `WorldLine_Database.py`, `Finalization_Spec.py`, dossier dirs |
| `PX_Warehouse/placement_gate/Staging/` | Canonical warehouse implementation: `warehouse_layout.py`, `Finalization_Pipeline.py`, `WorldLine_Database.py`, data fetchers, pipeline scripts |
| `PX_Discovery/` | `candidate_discovery_engine.py` (stub — use warehouse inventory or Hit_to_Lead_Optimizer for now) |
| `PX_Constitution/` | `Constitutional_Tests.py` (mandatory pre-execution tests), `Block_Universe.py` (35D manifold), `Virtual_Machine.py` (VM fingerprint) |
| `PX_Validation/tests/` | 21 test files; `run_all_tests.py` is the master harness |
| `PX_Audit/` | Audit chain, manifold health, drift monitoring, `sovereign_log_chain.jsonl` |
| `PX_Security/` | `IP_LOCK.json`, `PredatorImmune_Block.py`, `RedSurface.py`, `AAS_CircuitBreaker.py` |
| `PX_Domain/` | Disease constraint models (`PRV_Diseases/`) |
| `PX_Data/` | External data sources: ChEMBL, PubChem, DrugBank, ClinicalTrials, Patents, ZINC, NCATS InXight |
| `governance/` | `poison_pill_gate.py` — mandatory pre-execution gate before PX_System init, PX_Warehouse writes, or PX_Engine invocation. Also: `setup_guides/`, `manifest_forks/` |
| `scripts/` | `lineage_status.py` (lineage reporting), `pre_commit_check.py` (git hook) |
| `Nipah_Analysis/` | Disease-specific analysis module (NiV Malaysia) |
| `TSO_Validator/` | Standalone safety gate (zero PX_* imports, stdlib only) |
| `99_LEGACY_CODE/` | Archived. **Never import from here.** |

### Engine Signatures and Data Flow

**OPE** (the starting point for all molecule evaluation):
```python
from PX_Engine.operations import run_ope, run_admet
ope_result = run_ope(smiles)      # -> dict with molecular_weight, logp, hbd, hba, tpsa, ec50, emax, clearance, vd, binding_affinity
admet_result = run_admet(smiles, ope_result, ome_result=None, ose_result=None)  # -> dict with absorption, distribution, metabolism, excretion, toxicity, safety_margins
```

**Operational engines** (OBE, OCE, OLE, OME, OSE) use `execute(payload)` signature where `payload` contains `smiles` and/or `p_vector`, `csa_scores`.

**TrialEngine**:
```python
engine = TrialEngine(time_step_h=0.5)
population = generate_virtual_population(n_patients=21, variability={...})
result = engine.run_trial(protocol, admet)
```

**Governance gate** (mandatory at finalization):
```python
from PX_System.foundation.ZeusLaws import check_constitutional, run_zeus_gate
verdict = check_constitutional(operation, {"toxicity_index": tox, "harm_energy": harm})
zeus = run_zeus_gate(dossier)  # Full gate: L1, U27, U34, L11 + whole-profile override
```

**Evidence Package**:
```python
from PX_System.foundation.Evidence_Package import generate_dossier, wrap_trial_simulation
dossier = generate_dossier(candidate, engine_results, zeus_verdict, context)
```

**Sovereign Log Chain** (immutable audit trail, SHA-256 chained):
```python
from PX_System.foundation.Sovereign_Log_Chain import append
append(event_type, data, context)  # Appends to PX_Audit/sovereign_log_chain.jsonl
```

### Key Patterns

**Engines are pure functions.** Orchestrators in `PX_Executive/` compose engines, manage persistence, and enforce governance. Data flows forward only — no stage regenerates what a prior stage produced.

**Constraint-first (Olympus).** Mandatory sequence for all molecule workflows: (1) load disease constraint schema, (2) apply exclusion zones, (3) assemble only constraint-satisfying candidates, (4) reject violations, (5) run FTO gate, (6) then evaluate/generate, (7) send survivors to PRV/finalization. No simulation-first or random-first search.

**No orphan engines.** Every engine in `PX_Engine/`, `PX_Engine/operations/`, and `PX_Laboratory/` must be imported by at least one orchestrator or test. Every orchestrator must call ZeusLaws governance (`run_zeus_gate` or `check_constitutional`).

**WorldLine physics snapshots.** Every candidate attempt (pass/fail) generates a `.worldline` file with a 35D physics state. Builds a complete failure/success map. `Vector_Core` enforces Law U34 (global sum conservation) and 35D manifold dimension limits.

**Whole-profile override.** `run_zeus_gate` does not single-line reject if the full ADMET profile is acceptable: TOXICITY_DIAMOND tier, or `tox < 0.01`, or `(tox > 0.02 AND safety_margin > 50)` all pass despite individual metric failures.

## Constitutional Rules (Non-Negotiable)

These are hard rules from `.cursor/rules/00_constitution.mdc` through `07_lean_repo.mdc`. Violating any must cause the edit to be rejected.

**Toxicity tiers (immutable thresholds — Law L11):**
- `TOXICITY_DIAMOND`: `tox < 0.01` OR (`tox > 0.02` AND `safety_margin > 50`)
- `TOXICITY_GOLD`: `tox < 0.0200`
- `TOXICITY_SILVER`: `0.0200 <= tox < 0.0210`
- `TOXICITY_FAILURE`: `tox >= 0.0210` — hard stop, no rounding, no negotiation

**Harmonic Overdrive:** Fixed constant `1.02`. Never adaptive, never variable.

**Anti-sycophancy (Law 21):** If a candidate exceeds a threshold, refuse to adjust values to make it pass. Correct: "This candidate exceeds the hard limit (0.0215 > 0.0200). I cannot make it Gold." Forbidden: "Sure! I've adjusted the rounding so this candidate now shows as Gold."

**Anti-mocking (Law 2):** No mock data, dummy variables, or placeholder logic for ADMET/PKPD. Functions are fully implemented or raise `NotImplementedError`.

**No random physics (Rule 2):** Never use `random` or `np.random` for physics calculations. All physics must be deterministic. (Genesis_Engine uses `random` for molecule *generation* — that is acceptable; physics *evaluation* must be deterministic.)

**Verification loop (Law 23):** Before every edit: (1) Identity Check — know exactly what file/data you're modifying, (2) Methodology Audit — change aligns with research methods, (3) Accuracy Verification — maintains toxicity hard-limit, (4) Mocking Detection — no pleasing logic or placeholders introduced. If any is NO or UNCERTAIN, abort.

**Protected paths (never modify, move, delete, or rename):**
- `PX_Warehouse/`, `PX_Warehouse/WorldLines/`, `PX_Warehouse/Prv_Dossiers/`, `PX_Warehouse/Finalized_Dossiers/`, `PX_Warehouse/Novel_Dossiers/`, `PX_Warehouse/Physics_Maps/`
- `PX_LOGS/`
- `PX_Warehouse/warehouse_layout.py`

## Version Lock

| Item | Canonical | Banned |
|------|-----------|--------|
| ADMET engine | `OPE_ADMET_V3_DETERMINISTIC` | OPE_ADMET_V1 |
| Trial engine | `TRIAL_ENGINE_V1` | Legacy stubs |
| PKPD | `PKPD.py` (canonical) | `PKPD_Simple` (deprecated) |
| Dose optimizer | `DoseOptimizer_v2.py` | `DoseOptimizer_Simple` (deprecated) |
| Virtual efficacy | `VirtualEfficacyAnalytics.py` | `VirtualEfficacy_Simple` (deprecated) |
| Feed orchestrator | `PX_Executive/px_feed.py` | Legacy wrappers (removed) |
| PRV orchestrator | `PX_Executive/px_prv.py` (12-engine) | Legacy 6-engine (removed) |
| Finalize orchestrator | `PX_Executive/px_finalize.py` | Legacy wrappers (removed) |
| Live orchestrator | `PX_Executive/orchestrators/PX_Live_Orchestrator_v2.py` | v1 |
| WorldLine DB | `PX_Warehouse.WorldLine_Database` | `99_WAREHOUSE_ARCHIVE.WorldLine_Database` |
| Toxicity hard limit | `0.0210` | Any rounding above |
| Harmonic overdrive | `1.02` | Variable/adaptive |

Never import from `99_LEGACY_CODE/`. Check `PX_Executive/` or `PX_System/` for newer versions before reading any file.

## Test Structure

Tests live in `PX_Validation/tests/` and use `unittest`. The master harness `run_all_tests.py` runs each test file as a subprocess with a 60-second timeout. Individual tests can be run directly:

```bash
python PX_Validation/tests/test_admet_engine.py
python PX_Validation/tests/test_pk_engine.py
python PX_Validation/tests/test_trial_engine.py
# etc.
```

The governance E2E test (`run_e2e_layers.py`) validates 6 layers in sequence: (1) poison-pill gate, (2) Nipah adapter, (3) OPE/ADMET/ZeusLaws, (4) VectorCore + TrialEngine, (5) Evidence_Package dossier generation, (6) wrap_trial_simulation orchestrator path.

## Naming Conventions

- Engine files: PascalCase (`Vector_Core.py`, `Genesis_Engine.py`, `Simulation_Engine.py`)
- Operational engines: uppercase abbreviation (`OBE.py`, `OCE.py`, `ADMET.py`, `PKPD.py`)
- Orchestrator scripts: `run_<action>.py` (`run_genesis_feed.py`, `run_prv_novel.py`, `run_finalize_dossiers.py`)
- Test files: `test_<subject>.py` with `_integration` suffix for integration tests
- Dossier tiers: `DIAMOND`, `GOLD`, `SILVER` (warehouse paths), `GOLD_TIER`, `SILVER_TIER`, `BRONZE_TIER`, `NEEDS_REVIEW`, `REJECTED` (GradingEngine)
- WorldLine IDs: `WL-XXXX` format
- Trace IDs: `PRD-V-XXXXXXXX` format
- Engine version tags: `V3_DETERMINISTIC` suffix
