# Predator X — System Status

**Version**: 1.0.0-enterprise
**Last validated**: February 2026
**Platform**: Python 3.13.11, Windows, Miniconda3

## 1. Environment

The environment uses standard Python tooling:
- **setuptools + wheel** for packaging (pyproject.toml)
- **Just** for task automation (justfile)
- **Ruff** for linting (ruff.toml)
- **Pyright** for type checking (pyrightconfig.json)
- **Git + Git LFS** for version control
- **requirements-lock.txt** for pinned dependency versions

Core dependencies: RDKit 2025.9.3, NumPy 2.3.5, requests, reportlab, scipy 1.16.3.

## 2. Pipeline Architecture

Forward-only 12-engine PRV pipeline with QUINT internal language routing:

OPE -> ConstraintCheck -> OBE -> OCE -> OLE -> OME -> OSE -> ADMET -> PKSim -> PKPD -> DoseOpt -> Trial -> VirtualEfficacy -> GradingEngine -> ZeusLaws

All engines routed through QUINT — zero direct engine imports.

## 3. Governance Layer

- **ZeusLaws**: Python-native constitutional gate (Laws L1, L11, U27, U34)
- **Sovereign Log Chain**: SHA-256 chained audit trail (stdlib only)
- **Engine Sign-Off**: All 12 engines produce signed authorization chains
- **Fail-closed**: Gates default to `authorized=False`

## 4. Data Lineage

- **WorldLines**: 6,393 physics snapshots in PX_Warehouse/WorldLines/
- **Sovereign Log Chain**: 1,160+ audit entries in PX_Audit/sovereign_log_chain.jsonl
- **Lineage script**: `python scripts/lineage_status.py`

## 5. Validation Status

- **447/447** tests passed across 40 test files
- **7/7** governance E2E layers passed
- **TSO_Validator**: Standalone safety gate (zero PX_* imports)
- **65** preflight integration checks across 3 diseases

## 6. Warehouse Status

- **47 active dossiers**: 3 DIAMOND, 26 GOLD, 18 SILVER
- **871 REJECTED dossiers**: Archived as learning material
- **6 pharma packages**: DIAMOND-tier PDF dossier packages
- **Queue**: Clean (0 stale entries)

## 7. Disease Coverage

- **33 disease constraint files** (v2.0.0 schema)
- **31 diseases** in manifest (+ 2 extras: cholera, yellow_fever)
- **30 disease module folders** with adapters, analytics, validators
- Categories: tropical, rare pediatric, medical countermeasures (MCM)

## 8. Current Status

**Predator X is fully operational.**
No fictional build tools. No mock data. No placeholder logic.
