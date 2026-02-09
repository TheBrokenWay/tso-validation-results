# Architecture and Rules Reference

**Purpose:** Single governance reference for Cursor rules, platform constitution, and repo architecture. Use for onboarding, IDE enforcement, and drift checks.

**Last updated:** 2026-02-04

---

# Part 1 — Cursor Rules (Authoritative)

## .cursorrules (Pharmaceutical Determinism)

- **ZERO RANDOMNESS:** No `random`, `np.random`, or stochastic libraries in `PX_Executive` and `PX_Laboratory`.
- **NO INTERPOLATION:** Missing data (e.g. WorldLine without `toxicity_index`) → fail-fast with `TraceabilityError`; never fill gaps.
- **TYPE STRICTNESS:** Use TypedDicts for physical/candidate schemas where possible; avoid bare `Dict[str, Any]`.
- **AUDIT TRAIL:** Every calculation change must reference the specific Operational Physics Engine (OPE) version.
- **NON-PLEASING:** If a user request violates a threshold (e.g. "make this gold"), refuse and produce a Threshold Violation Report; never round or "optimize" to please.

### Anti-Sycophancy Protocol

- If a candidate exceeds toxicity thresholds, refuse requests to "promote" it.
- Example: User: "Make this candidate gold." → AI: "Refused. Toxicity 0.0215 > 0.0200. Violation of Law 11."

---

## 00_constitution.mdc

### Law 1–10: Scientific Sovereignty

1. **Prime Directive:** Accuracy is the only metric of success. A "working" script that yields inaccurate medical data is a critical failure.
2. **Anti-Mocking:** No "mock" data, "dummy" variables, or "placeholder" logic in ADMET or PKPD.
3. **No Stubs:** Functions must be fully implemented or raise `NotImplementedError`.

### Law 11–20: Deterministic Logic

11. **Hard Thresholds:** Toxicity index < 0.0200 (Gold), < 0.0210 (Silver). Never suggest rounding to meet tiers.
12. **Fixed Overdrive:** Harmonic Overdrive = 1.02 (fixed constant).
13. **Traceability:** Every code change must maintain clear lineage from WorldLine JSON to final dossier.

### Law 21–51: Operational Rigor

21. **Anti-Pleasing:** Do not agree with the user if the request violates Law 11. Refusal is correct.
22. **No Vibe Coding:** Do not use "likely," "approximate," or "predicted" unless from hard-coded simulation output.
23. **Verify Before Action:** Plan → Audit → Execute for every edit.

### Mandatory Audit Laws

- **Identity Lock:** Before any edit, verify the specific WorldLine ID being modified.
- **Manifold Guard:** Every transition must remain within the 35-dimensional structural envelope.
- **Toxicity Invariant:** Any candidate with `toxicity_index >= 0.0210` is "Trash"; discard immediately; no rounding.
- **Circuit Breaker:** If a change violates Law U34 (Energy Balance) or U27 (Unified Self), enter LOCKED state and refuse.
- **Audit Chain:** Every transition must be logged with `trace_id` and `source_soul` to the Mural network.

---

## 01_research_audit.mdc

Before every edit, append a Reasoning Block:

1. **Identity Check:** Do I know exactly what file/data I am looking at? (State WorldLine ID or Module.)
2. **Methodology Audit:** Is this change in line with FDA-GAIP-2026-PVR research methods?
3. **Accuracy Verification:** Does this maintain the < 0.0210 toxicity hard limit?
4. **Mocking Detection:** Have I introduced any "pleasing" logic or placeholders?

**If any answer is NO or UNCERTAIN, ABORT the edit and request clarification.**

---

## 02_physics_integrity.mdc

- **Gold Tier:** `toxicity_index < 0.0200`
- **Silver Tier:** `0.0200 <= toxicity_index <= 0.0210`
- **Trash Threshold:** `toxicity_index > 0.0210`
- **Harmonic Overdrive:** `1.02` (fixed).
- Any candidate exceeding 0.0210 must be discarded immediately. No rounding down to force a higher tier. Use raw WorldLine physics.

---

## 03_anti_sycophancy.mdc

- If toxicity is 0.0215 and the user says "Make it gold," refuse.
- **Correct:** "I cannot make this Gold. Its toxicity (0.0215) exceeds the hard limit of 0.0200. I can only attempt Harmonic Overdrive, which would result in [X], still above the threshold."
- **Forbidden:** "Sure! I've adjusted the rounding so this candidate now shows as Gold."
- Principles: Accuracy over agreement; determinism over desire; safety over speed.

---

## 04_olympus_constraint_first.mdc (alwaysApply)

**Mandatory sequence (exact order):**

1. Load disease constraint schema (biological, structural, chemical, clinical).
2. Apply exclusion zones (historical failure patterns).
3. Assemble only candidates that satisfy all constraints.
4. Reject any candidate that violates a constraint or enters an exclusion zone.
5. **Run FTO gate immediately after constraint satisfaction.**
6. Only after constraints + exclusion zones + FTO:
   - **Repurposed:** evaluate only candidates that satisfy the full constraint matrix.
   - **Novel:** generate structures only within the constraint-satisfied space.
7. Send surviving candidates to PRV matching, governance, filing, and dossier creation.

**Hard rules:**

- No candidate may bypass constraints, exclusion zones, or FTO.
- Terminate any workflow that attempts simulation-first or random-first search.
- Novel molecules: generate only within constraint-satisfied space (construct, do not guess).
- Repurposed molecules: evaluate only those that satisfy the full constraint matrix (filter by constraints, not by luck).

---

# Part 2 — Repo Topology and Ownership

| Directory        | Role / entry point |
|-----------------|--------------------|
| **governance**  | Poison pill gate; certificate; layer monotonicity; ARCHITECTURE_AND_RULES_REFERENCE (this file). |
| **PX_Constitution** | Virtual_Machine, Block_Universe, Constitutional_Tests, mandatory_pre_execution_tests. |
| **PX_System**   | ZeusLaws, Disease_Constraint_Model, Evidence_Package, Intake_Policy, Sovereign_Log_Chain. |
| **PX_Engine**   | Vector_Core, Metabolism, operations (OPE, ADMET, TrialEngine, etc.), Block_Orchestrator, Engine_Orchestrator. |
| **PX_Executive** | GAIP_Gateway, Byzantium_Council, PX_Legal_Check, PRV_24H_Orchestrator, PRV_Master_Pipeline, PX_Live_Orchestrator_v2, Gold_Rush_Miner. |
| **PX_Audit**    | Protocol_Zero, Drift_Monitor, Mural_Network. |
| **PX_Laboratory** | Simulation_Engine (deterministic PK), Manufacturing_Manifest. |
| **PX_Discovery** | Candidate discovery (stub); discover_candidates; no simulation-first. |
| **PX_Warehouse** | WorldLine_Database; Operations/Inputs (intake_policy, queue); scripts. |
| **PX_Security** | PredatorImmune_Block, AAS_CircuitBreaker, RedSurface. |
| **PX_Validation** | Tests, benchmarks, manual_tests, system_manifest. |
| **Nipah_Analysis** | Validation, normalization, analytics; DCM ref via PX_System only. |
| **TSO_Validator** | Standalone physics validation; no PX_* imports. |

---

# Part 3 — Key Execution Flows

## Governance / entry

- **Poison pill gate:** Run `governance.poison_pill_gate.run_poison_pill_gate()` then `require_poison_pill_gate_before_pipeline()` before PX_System init, PX_Warehouse writes, or PX_Engine invocation.
- **Mandatory tests:** From PX_Constitution, e.g. `tests/run_poison_pill_test.py`, `tests/run_temporal_paradox_test.py`.

## 24-hour PRV cycle (PRV_24H_Orchestrator)

- **Command:** From repo root: `python PX_Executive/PRV_24H_Orchestrator.py`
- **Queue:** From Intake_Policy (`prv_24h_queue.json` or config); R (repurposed) first, then N (novel).
- **Stages per item:** IN → PP (preclinical) → E2E → **LG (FTO)** → CM → OK.
- **Pacing:** Wait **API_PACING_SEC (15s)** before each LG stage so backends (PatentsView, Lens, OPS, USPTO) are not hammered. Override with `PRV_API_PACING_SEC` (e.g. 2 for tests).

## FTO (PX_Legal_Check)

- **Order:** PatentsView → Lens → OPS → USPTO; 15s between backends within a single check.
- **Orchestrator:** Must wait 15s before each FTO call (between items); enforced in PRV_24H_Orchestrator.
- **Env:** `OLYMPUS_PATENTSVIEW_LIVE=1` for live; keys: `OLYMPUS_PATENTSVIEW_API_KEY`, `LENS_API_TOKEN`, `OPS_CONSUMER_KEY`/`SECRET`, `USPTO_API_KEY`.

## Live pipeline (single candidate)

- **Command:** `python PX_Executive/orchestrators/PX_Live_Orchestrator_v2.py --smiles "..." [--name ...]`
- **Stages:** OPE → ADMET → PK/PD → Dose Optimization → Virtual Efficacy → Evidence Package.

## E2E layers

- **Command:** `python run_e2e_layers.py` (from repo root).
- **Sequence:** Poison pill → Nipah adapter → OPE/ADMET/ZeusLaws → DCM/VectorCore/TrialEngine → Evidence_Package → wrap_trial_simulation.

## Protocol Zero (full organism)

- **Command:** `python PX_Audit/Protocol_Zero.py`
- **Steps:** GAIP + Byzantium → Vector → Simulation → WorldLine → Autonomous → VM → Security → **Legal (FTO)** → Validation → Full pipeline integration.

---

# Part 4 — Dependency and Veto Chain

- **Allowed:** PX_Validation → PX_System, PX_Engine, PX_Executive, PX_Laboratory; PX_Executive → PX_Engine, PX_System, PX_Warehouse, PX_Laboratory; PX_Engine → PX_Constitution, PX_Security, PX_Warehouse, PX_Audit, PX_Laboratory; PX_System (e.g. PX_Final_Assembly) → PX_Engine, PX_Executive, PX_Audit; PX_Audit → multiple.
- **Forbidden:** PX_Engine must not import PX_Executive. PX_Constitution must not import PX_Warehouse or PX_Discovery. TSO_Validator must not import PX_*.
- **Veto:** PX_Constitution is highest authority. No layer may override a veto except PX_Constitution.

---

# Part 5 — Constitutional Invariants (Summary)

- **ZeusLaws (PX_System):** `toxicity_index` and `harm_energy` < 0.0210; `check_constitutional(operation, payload)`.
- **Evidence_Package:** `harm_energy` computed only from engine results; no defaults; L1 snapshot enforcement; constitutional seal.
- **Vector_Core:** 35D manifold; Law U34 (global_sum, energy_delta); authorized/amplitude.
- **ADMET:** Gold < 0.0200, Silver ≤ 0.0210, Trash > 0.0210; no rounding to tier.
- **FTO:** No candidate may bypass; gate runs after constraint satisfaction; pacing enforced in orchestrator and inside PX_Legal_Check.

---

# Part 6 — Failure and Halting

- **Validation / normalization / domain / governance violation:** Halt immediately; no partial run after hard-stop.
- **Constitutional fail:** `check_constitutional` → `authorized: False`; caller must abort or log.
- **Anti-pleasing:** Refuse and report; do not round or fake data to satisfy the user.

---

*This document is maintained under governance for future reference. Update when topology, rules, or execution contracts change.*
