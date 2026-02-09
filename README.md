# Predator X — Deterministic Scientific Platform

Predator X is a fully deterministic, constitutionally governed scientific research platform. Every tool, engine, and physics layer is pinned, reproducible, and validated end‑to‑end.

This README provides the top‑level overview and links to the core documents that define the system.

---

## 🔒 Deterministic Environment
Predator X uses:
- **Nix** for version locking  
- **direnv** for automatic environment activation  
- **flake.nix** as the environment ledger  

This ensures identical results on any machine.

See: `DETERMINISTIC_SETUP.md`

---

## 🧱 Hermetic Build System
All engines run under:
- **Bazel 7** (hermetic execution, caching)

This guarantees reproducible simulations and prevents drift.

---

## 🧬 Python Architecture & Safety
The Python layer is governed by:
- **Rye** (project structure)
- **Ruff + Pyright** (static safety)

This eliminates runtime type errors and ensures molecular engine stability.

---

## 🏛 Governance & Constitution
Predator X enforces:
- **OPA** (Zeus Laws)
- **Hydra** (configuration governance)

Governance is externalized and cannot be bypassed.

See: `GOVERNANCE_VERSION_LOCK.md`

---

## 🗂 Physics Lineage
All Physics Maps and worldlines are versioned using:
- **DVC**

This provides full lineage and time‑travel capability.

---

## 🧪 Validation Status
- **25/25 PX_Validation tests passed**
- **Governance stress tests passed**
  - Thermodynamic paradox → rejected  
  - Temporal paradox → rejected  

See: `STATUS.md`

---

## 🎛 Ergonomic Control Panel
A Justfile provides simple commands for engine loops:

- `just feed`
- `just novel`
- `just repurpose`
- `just finalize`
- `just cycle`

See: `JUSTFILE.md`

---

## 📐 Architecture Diagram
A full textual architecture diagram is available in:

`ARCHITECTURE_DIAGRAM.md`

---

## ✔ Current Status
**Predator X is fully installed, deterministic, governed, and operational.  
No additional configuration is required.**
