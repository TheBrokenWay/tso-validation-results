# Predator X — System Flow

This document explains the full operational lifecycle of the Predator X platform, from environment initialization to governance and lineage.

---

## 1. Initialize (Deterministic Environment)

When entering the repo:

- `direnv` loads the Nix flake  
- Bazel, Rye, Ruff, Pyright, OPA, Hydra, DVC are pinned  
- Global system paths are ignored  

This ensures reproducibility.

---

## 2. Configure (Hydra)

Hydra governs:

- Discovery parameters  
- PRV modes  
- Experimental seeds  
- Safety thresholds  

All configuration is versioned and auditable.

---

## 3. Execute (Engine Loops)

Use the Justfile:

- `just feed`  
- `just novel`  
- `just repurpose`  
- `just finalize`  
- `just cycle`  

These call the real engine scripts in `PX_Executive/`.

---

## 4. Govern (OPA + Governance Pipeline)

Governance is enforced through:

- `run_e2e_layers.py`  
- Poison-pill paradox defense  
- Nipah timeline defense  
- OPE/ADMET  
- PK engine  
- Trial engine  
- Evidence package  

This ensures constitutional compliance.

---

## 5. Audit (DVC Lineage)

Physics Maps and worldlines are versioned via DVC.

- `just lineage`  
- Tracks changes  
- Enables time-travel  
- Ensures reproducibility of physics states  

---

## 6. Validate (PX_Validation)

Full validation suite:

- `just test`  
- 25/25 tests  
- Integration tests  
- Engine tests  
- Governance tests  

This confirms system integrity.

---

## Summary

The Predator X system flow is:

**Initialize → Configure → Execute → Govern → Audit → Validate**

Every step is deterministic, governed, and reproducible.
