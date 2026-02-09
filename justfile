# Predator X — Justfile command layer
# All commands run from repo root. Use with Nix + direnv (just is in flake.nix).

default:
    @just --list

# ---- Engine loops (canonical entry points) ----

# Genesis feeder: invent novel SMILES and write to Feeder queue
feed:
    python PX_Executive/run_genesis_feed.py

# Repurposed feed: populate queue from repurposed sources (run before repurpose if needed)
feed-rep:
    python PX_Executive/run_repurposed_feed.py

# Novel engine: 24H orchestrator in NOVEL mode (processes type N from queue)
novel:
    python PX_Executive/run_prv_novel.py

# Repurposed engine: 24H orchestrator in REPURPOSED mode (processes type R from queue)
repurpose:
    python PX_Executive/run_prv_repurposed.py

# Finalization: run finalization pipeline on unfinalized dossiers
finalize:
    python PX_Executive/run_finalize_dossiers.py

# Full cycle: Feeder → Novel (1 item) → Repurposed (1 item) → Finalize (1) — for testing
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

# Governance stress test (poison-pill gate + Nipah + OPE/ADMET + Evidence)
govern:
    python run_e2e_layers.py

# Full validation suite (25 tests)
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
