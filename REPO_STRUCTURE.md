# Predator X — Repository Structure

This document provides a clear overview of the Predator X repository layout and how each directory participates in the deterministic, governed scientific workflow.

---

## Top-Level Structure

| Directory / File | Purpose |
|------------------|---------|
| **PX_Executive/** | Engine loop entrypoints (Feed, Novel, Repurposed, Finalization, Cycle) |
| **PX_System/** | Core system logic, Evidence_Package, ZeusLaws, Intake_Policy |
| **PX_Engine/** | OPE, ADMET, TrialEngine, PK, Vector_Core, Genesis, Metabolism |
| **PX_Validation/** | Full validation suite (25/25 tests) |
| **PX_Warehouse/** | Dossiers, WorldLines, Feeder; Physics Maps versioned via DVC |
| **governance/** | Poison-pill gate, policy, constitutional checks |
| **scripts/** | Setup scripts (e.g. deterministic environment installer) |
| **flake.nix** | Environment ledger (pinned toolchain) |
| **.envrc** | Environment trigger (direnv) |
| **justfile** | Ergonomic command panel |
| **README.md** | Top-level overview |
| **STATUS.md** | Deterministic & governed system status |
| **GOVERNANCE_VERSION_LOCK.md** | Version lock and constitutional enforcement |
| **DETERMINISTIC_SETUP.md** | Nix + direnv setup |
| **ARCHITECTURE_DIAGRAM.md** | Layered architecture diagram |
| **JUSTFILE.md** | Justfile command reference |

---

## Engine Entrypoints (PX_Executive)

| Script | Purpose |
|--------|---------|
| `run_genesis_feed.py` | Genesis feed engine |
| `run_repurposed_feed.py` | Repurposed feed engine |
| `run_prv_novel.py` | Novel mode wrapper for PRV orchestrator |
| `run_prv_repurposed.py` | Repurposed mode wrapper |
| `run_finalize_dossiers.py` | Finalization engine |
| `run_one_cycle_test.py` | Full-cycle test |
| `reprocess_warehouse_lifecycle.py` | One-time batch: push all unfinalized Prv/Novel dossiers through grading + finalization → Finalized_Dossiers/\<tier\>; use `--validate-before` / `--validate-after` before version lock |

---

## Governance Entrypoint

- **run_e2e_layers.py** (repo root)  
  Performs poison-pill, Nipah, OPE/ADMET, PK, Trial, and Evidence governance checks.

---

## Validation Entrypoint

- **PX_Validation/tests/run_all_tests.py**  
  Runs the full 25/25 validation suite.

---

## Deterministic Environment Files

- **.envrc** → auto-loads flake  
- **flake.nix** → pins Bazel, Rye, Ruff, Pyright, OPA, Hydra, DVC, Python 3.11, RDKit deps  
- **justfile** → ergonomic command panel  

---

## Notes

- All engines run under the deterministic environment.  
- No global system tools are used.  
- Governance and validation pipelines are fully aligned with repo structure.  
- TrialEngine, OPE, ADMET, PK, and Evidence_Package live inside **PX_Engine/operations/** and **PX_System/foundation/**; see “Key modules” or GOVERNANCE_VERSION_LOCK.md and PX_Warehouse/warehouse_layout.py for details.
