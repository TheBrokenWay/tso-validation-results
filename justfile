# Predator X — Justfile command layer
# All commands run from repo root. Just installed via WinGet.
#
# Windows: just defaults to sh (often missing). Use cmd.exe so just cycle etc. work.
# Linux/Mac: unchanged (sh/bash from environment).
set windows-shell := ["cmd.exe", "/c"]

default:
    @just --list

# ---- Engine loops (consolidated orchestrators) ----

# Genesis feeder: invent novel SMILES and write to Feeder queue
feed:
    python PX_Executive/px_feed.py --mode novel

# Repurposed feed: populate queue from repurposed sources (run before repurpose if needed)
feed-rep:
    python PX_Executive/px_feed.py --mode repurpose

# Novel engine: full 12-engine PRV pipeline (processes type N from queue)
novel:
    python PX_Executive/px_prv.py --type novel

# Repurposed engine: full 12-engine PRV pipeline (processes type R from queue)
repurpose:
    python PX_Executive/px_prv.py --type repurpose

# All PRV candidates: novel + repurposed
prv-all:
    python PX_Executive/px_prv.py --type all

# Finalization: run finalization pipeline on unfinalized dossiers
finalize:
    python PX_Executive/px_finalize.py

# Full cycle: Feed → PRV (1 novel + 1 repurposed) → Finalize — for testing
cycle:
    python PX_Executive/run_one_cycle_test.py

# One-time batch: push all unfinalized Prv/Novel dossiers through grading + finalization → Finalized_Dossiers/<tier>
reprocess-lifecycle:
    python PX_Executive/reprocess_warehouse_lifecycle.py

# Same with validation before and after (recommended before version lock)
reprocess-lifecycle-validated:
    python PX_Executive/reprocess_warehouse_lifecycle.py --validate-before --validate-after

# ---- Development & governance ----

# Lint and type-check using pinned Ruff + Pyright
check:
    ruff check .
    pyright

# CI-ready enterprise compliance checklist (file-level rules, fails build on violation)
comply:
    python scripts/ci_compliance_check.py

# Governance stress test (poison-pill gate + Nipah + OPE/ADMET + Evidence)
govern:
    python run_e2e_layers.py

# Enterprise build (test suite removed)
test:
    python PX_Validation/tests/run_all_tests.py

# TSO_Validator standalone safety gate
tso:
    python TSO_Validator/run_validation.py

# TSO_Validator tests
tso-test:
    python -m pytest TSO_Validator/tests/ -q

# Full warehouse simulation — structural audit, classification, proposed actions
warehouse-audit:
    python PX_Warehouse/run_warehouse_simulation.py --enforce

# Lineage status: Sovereign Log Chain, WorldLines, DataLineage graph
lineage:
    python scripts/lineage_status.py
