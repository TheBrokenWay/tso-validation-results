# Predator X — Deterministic Pharmaceutical Compound Discovery Platform

Predator X is a fully deterministic, constitutionally governed pharmaceutical compound discovery platform. It generates and evaluates novel/repurposed molecules through a forward-only pipeline: **Feed -> PRV -> Trial -> Finalization**. Every stage is gated by governance (ZeusLaws), and all outputs maintain full lineage from WorldLine physics snapshots to final dossiers.

---

## Environment

| Component | Tool | Version |
|-----------|------|---------|
| Language | Python | 3.13 (Miniconda3, Windows) |
| Build | setuptools + wheel | via pyproject.toml |
| Lint | Ruff | ruff.toml at repo root |
| Type-check | Pyright | pyrightconfig.json at repo root |
| Task runner | Just | justfile at repo root |
| Dependencies | pip / conda | requirements-lock.txt for pinned versions |
| Version control | Git + Git LFS | .gitattributes tracks *.db, *.tar.gz, *.tar.bz2 |

Core dependencies: RDKit (>=2023.9.1), NumPy (>=1.24.0), requests (>=2.31.0), reportlab (>=4.0).

Constitutional modules (ZeusLaws, Sovereign_Log_Chain, QUINT) use only Python stdlib.

---

## Quick Start

```bash
just test          # Run full test suite (447 tests across 40 files)
just govern        # Governance stress test (7 layers)
just check         # Lint + type-check (ruff + pyright)
just feed          # Genesis feed (novel candidates)
just novel         # 12-engine PRV pipeline (novel)
just finalize      # Finalization pipeline
just cycle         # Full cycle: Feed -> PRV -> Finalize
```

---

## Architecture

```
Feed Stage (Genesis_Engine, Vector_Core, Metabolism, Trajectory_Predictor)
  -> Queue (PX_Warehouse/Feeder/prv_24h_queue.json)
  -> PRV Stage (OPE -> Constraints -> OBE -> OCE -> OLE -> OME -> OSE -> ADMET -> PKSim -> PKPD -> DoseOpt -> Trial -> VirtualEfficacy -> GradingEngine -> ZeusLaws)
  -> Finalization (Evidence_Package, GradingEngine, ZeusLaws.run_zeus_gate)
  -> PX_Warehouse/Finalized_Dossiers/<DIAMOND|GOLD|SILVER>/
```

See `CLAUDE.md` for full architecture reference and constitutional rules.

---

## Governance

- **ZeusLaws**: Constitutional governance gate (Laws L1, L11, U27, U34)
- **Toxicity hard limit**: 0.0210 (Law L11) — no rounding, no negotiation
- **Harmonic Overdrive**: Fixed constant 1.02 — never adaptive
- **Fail-closed**: All gates default to `authorized=False`
- **Sovereign Log Chain**: SHA-256 chained immutable audit trail

---

## Validation Status

- **447/447** PX_Validation tests passed (40 test files)
- **7/7** governance E2E layers passed
- **TSO_Validator**: Standalone safety gate (zero PX_* imports, stdlib only)
- **65** preflight integration checks across 3 diseases, all 12 engines

---

## Current Status

**v1.0.0-enterprise** — Predator X is fully operational with 12-engine PRV pipeline, QUINT internal language, pharma-grade dossier system, and 33 disease constraint models.
