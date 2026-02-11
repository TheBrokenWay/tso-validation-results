# Contributing to Predator X

Thank you for contributing to the Predator X deterministic scientific platform.  
This document explains how to work inside the governed, reproducible environment.

---

## 1. Environment Setup

Before contributing:

1. Install Nix + direnv  
2. Ensure `.envrc` and `flake.nix` are present  
3. Run:

```bash
direnv allow
```

This loads the pinned toolchain (Bazel, Rye, Ruff, Pyright, OPA, Hydra, DVC).

---

## 2. Running Engine Loops

Use the Justfile commands:

```bash
just feed
just novel
just repurpose
just finalize
just cycle
```

These commands call the real engine scripts in `PX_Executive/`.

---

## 3. Running Validation

Before submitting changes:

```bash
just test
```

This runs the full 25/25 validation suite.

---

## 4. Governance Check

Every contribution must pass the governance pipeline:

```bash
just govern
```

This runs:

- Poison-pill paradox defense  
- Nipah timeline defense  
- OPE/ADMET  
- PK engine  
- Trial engine  
- Evidence package  

---

## 5. Static Analysis

All code must pass:

```bash
just check
```

This runs Ruff + Pyright using the pinned versions.

---

## 6. Physics Lineage

If modifying Physics Maps:

```bash
just lineage
```

(DVC must be configured with a `dvc.yaml`.)

---

## 7. Pull Request Requirements

Every PR must include:

- Passing validation suite  
- Passing governance pipeline  
- Passing static analysis  
- Clear description of changes  
- No bypass of deterministic or governance layers  

---

## 8. Notes

- Do not install global Python packages.  
- Do not run engines outside the deterministic environment.  
- All contributions must maintain reproducibility and constitutional compliance.
