# PX Enterprise Edition — Final, Locked Specification

**Constitutionally-Governed Deterministic Scientific Computing Platform**

This document defines the **official, enterprise-ready edition** of the Predator X platform. It describes the platform's purpose, guarantees, runtime surface, and governance requirements. It is paired with **PX_ENTERPRISE_LAYOUT.md**, which defines the frozen directory structure.

**This edition is locked.**  
Any modification requires a constitutional amendment.

---

## 1. Purpose of the Enterprise Edition

The PX Enterprise Edition is a **deterministic, governed** scientific computing platform. It provides:

- Constitutional invariants  
- Dual-DAG governance overlay  
- Mandatory poison-pill gate  
- Optional but active TSO second gate  
- Deterministic validation and normalization  
- Physics-based engine execution  
- System-level disease modeling  
- Evidence and trial simulation  
- Full auditability and provenance  
- Governance readiness certification  

This edition is designed for **regulated, multi-team, enterprise environments**.

---

## 2. Enterprise Guarantees

The platform guarantees:

| Guarantee | Description |
|-----------|-------------|
| **Determinism** | All computations follow fixed seeds, fixed execution order, and fixed engine semantics. |
| **Reproducibility** | Every run produces a reproducibility bundle and a reproducibility ID where applicable. |
| **Governance Enforcement** | Execution is only allowed after: poison-pill gate passes; Governance Readiness Certificate is valid; architecture drift is zero; constitutional invariants load successfully. |
| **Auditability** | All layers emit structured logs and provenance artifacts. |
| **Security** | Runtime modules enforce layer-monotonicity and reject reclassified failures. |
| **Scientific Validity** | The platform satisfies the Scientific Validity Framework defined in the governance map. |

---

## 3. Runtime Surface (What the Enterprise Edition Provides)

### Core Organs

- PX_System  
- PX_Engine  
- PX_Executive  
- PX_Validation  
- PX_Audit  
- PX_Security  
- PX_Warehouse  
- PX_Laboratory  
- PX_Discovery  
- PX_STATE  
- PX_LOGS  

### Domain Slice

- **Nipah_Analysis** (reference implementation)  

### Second Gate

- **TSO_Validator** (active constitutional gate)  

### Evidence & Trial Layer

- Evidence_Package  
- TrialEngine  
- VectorCore  

### Orchestrator

- **PredatorXOrchestratorV2** (PX_Executive/orchestrators/PX_Live_Orchestrator_v2.py)  
- **run_e2e_layers.py**  

### Governance

- PX_Constitution  
- governance/  
- Poison-pill tests  
- Gate stamp signing  
- Governance Readiness Certificate  

---

## 4. Execution Requirements

The IDE and all enterprise environments **must** enforce:

### 4.1 Mandatory Pre-Execution Steps

1. Run poison-pill gate  
2. Verify gate stamp signature (when configured)  
3. Regenerate governance coverage (on gate pass)  
4. Regenerate Governance Readiness Certificate  
5. Confirm architecture drift = 0  
6. Load constitutional invariants  
7. Validate domain input  
8. (Optional) Run TSO second gate  

### 4.2 Runtime Invariants

- Layer-monotonicity  
- No governance bypass  
- No shadow ontology  
- No dual-truth collision  
- No temporal paradox  
- No energy paradox (reserved class)  

### 4.3 Post-Execution Requirements

- Evidence dossier generation  
- Trial simulation  
- Provenance bundle  
- Reproducibility bundle (where applicable)  

---

## 5. Directory Structure

The Enterprise Edition directory structure is defined in:

**PX_ENTERPRISE_LAYOUT.md**

The IDE must enforce that:

- Only directories listed in the layout file exist at the root (or are explicitly optional/archive).  
- All unlisted directories are moved to `archive/`.  
- No new runtime modules are allowed without constitutional amendment.  

---

## 6. Governance Readiness Certificate

The certificate includes:

- Governance coverage  
- Active vs reserved poison pills  
- Constitutional invariants  
- Gate stamp status  
- Timestamp  
- Drift status  
- TSO gate status (when run)  
- Reproducibility metadata (when applicable)  

**Execution is not permitted unless the certificate is valid** (stamp passed, invariants loaded, coverage generated).

---

## 7. Amendment Process

Any change to:

- Directory structure  
- Runtime modules  
- Constitutional invariants  
- Governance tests  
- Engine semantics  
- Validation ontology  
- Orchestrator behavior  

**requires**:

1. Constitutional amendment  
2. Update to PX_ENTERPRISE_LAYOUT.md  
3. Update to this file  
4. Passing all governance tests  
5. Zero drift  
6. Regenerated readiness certificate  

---

## 8. Enterprise Edition Lock

**This edition is locked.**

It represents the clean, governed, enterprise-ready version of the platform. It is suitable for deployment into any enterprise environment.
