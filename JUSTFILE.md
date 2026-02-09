# Predator X — Justfile Command Layer

This file defines ergonomic commands for operating the four Predator X engine loops and common development tasks. These commands rely on the deterministic environment provided by Nix + direnv + flake.nix.

## Engine Loops

### Feed Engine
```bash
just feed
```
Runs the Feeder loop using pinned Bazel, Python, and Physics Maps.

### Novel Engine
```bash
just novel
```
Executes the Novel Discovery loop with deterministic physics and governance checks.

### Repurposed Engine
```bash
just repurpose
```
Runs the Repurposed loop with full ADMET, PK, and governance validation.

### Finalization Engine
```bash
just finalize
```
Executes the Finalization loop and prepares the Evidence Package.

### Full Cycle
```bash
just cycle
```
Runs Feeder → Novel → Repurposed → Finalization in sequence.

### Batch lifecycle reprocess (one-time backfill)
```bash
just reprocess-lifecycle
```
Pushes all unfinalized dossiers in Prv_Dossiers/Novel_Dossiers through grading and finalization into Finalized_Dossiers/\<tier\>. Use `just reprocess-lifecycle-validated` to run validation before and after the batch (recommended before version lock).

## Development Utilities

### Lint & Type Check
```bash
just check
```
Runs Ruff + Pyright using the pinned toolchain.

### Governance Check
```bash
just govern
```
Runs the governance stress test (poison-pill gate, Nipah adapter, OPE/ADMET, Evidence Package). In this repo the canonical run is `run_e2e_layers.py`; OPA/Hydra enforce rules at design time and via that pipeline.

### Physics Lineage Snapshot
```bash
just lineage
```
Captures a DVC snapshot of Physics Maps and worldlines. Requires a configured `dvc.yaml`.

## Notes

- All commands run inside the deterministic environment. No global system tools are used.
- Optional: `just feed-rep` populates the queue from repurposed sources; `just test` runs the full 25-test validation suite.
