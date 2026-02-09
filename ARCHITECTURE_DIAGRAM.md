# Predator X — System Architecture Diagram (Textual)

Below is a clear, layered architecture diagram describing how the Predator X deterministic platform is structured.

```
+-------------------------------------------------------------+
|                     EXECUTIVE LAYER                         |
|                 (Justfile Command Panel)                    |
+-------------------------------------------------------------+
                              |
                              v
+-------------------------------------------------------------+
|                 GOVERNANCE & CONFIG LAYER                   |
|   - OPA (Zeus Laws, constitutional enforcement)             |
|   - Hydra (parameter governance, config ledger)             |
+-------------------------------------------------------------+
                              |
                              v
+-------------------------------------------------------------+
|                PYTHON & DATA MANAGEMENT LAYER               |
|   - Rye (Python project manager)                            |
|   - Ruff + Pyright (static safety)                          |
|   - DVC (Physics Map lineage, worldline versioning)         |
+-------------------------------------------------------------+
                              |
                              v
+-------------------------------------------------------------+
|                BUILD & EXECUTION LAYER                      |
|   - Bazel 7 (hermetic execution, caching)                   |
|   - RDKit/Physics dependencies (zlib, expat)                |
+-------------------------------------------------------------+
                              |
                              v
+-------------------------------------------------------------+
|                DETERMINISTIC ENVIRONMENT LAYER              |
|   - Nix (version lock, reproducibility)                     |
|   - direnv (auto-activation)                                |
|   - flake.nix (environment ledger)                          |
+-------------------------------------------------------------+
                              |
                              v
+-------------------------------------------------------------+
|                     HOST SYSTEM (WSL2)                      |
+-------------------------------------------------------------+
```

This diagram shows the flow from the deterministic substrate up through governance and into the ergonomic control layer.
