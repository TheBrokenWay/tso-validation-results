# PREDATOR X — Engine Integration Specification (IDE Enforcement Edition)

Ensures 100% engine connectivity, no bypasses, no dead modules, no partial throughput.

---

## 1. Canonical Engine Registry (All Engines That Must Be Connected)

The system contains three engine families, all of which must be wired into the pipeline.

### A. PX_Engine Root Engines

| Engine | Purpose |
|--------|---------|
| Genesis_Engine | Novel SMILES generation |
| Trajectory_Predictor | Next-state safety |
| Metabolism | Pulse / heartbeat |
| Vector_Core | Physics gate (p_vector) |
| Block_Orchestrator | Block-level orchestration |
| Engine_Orchestrator | Engine-level orchestration |
| Stress_Test | Stress testing |

### B. PX_Engine.operations Engines

| Engine | Purpose |
|--------|---------|
| OBE | Operational Blocker Engine |
| OCE | Operational Coherence Engine |
| OLE | Operational Logic Engine |
| OME | Operational Momentum Engine |
| OPE | Operational Physics Engine |
| OSE | Operational Status Engine |
| ADMET | ADMET evaluation |
| TrialEngine | 7-tier trial simulation |
| PKPD | PK/PD linking (canonical) |
| PKPD_Simple | PK/PD (deprecated; governance prefers PKPD) |
| DoseOptimizer | Dose optimization |
| DoseOptimizer_v2 | Dose optimization (canonical) |
| DoseOptimizer_Simple | Deprecated |
| VirtualEfficacy_Simple | Virtual efficacy (deprecated) |
| VirtualEfficacyAnalytics | Virtual efficacy (canonical) |
| GradingEngine | Discovery grading |
| GradingSchema_Discovery.json | Schema dependency for GradingEngine |

### C. PX_Laboratory Engines

| Engine | Purpose |
|--------|---------|
| Simulation_Engine | Simulation |
| Manufacturing_Manifest | Production order / manifest |

---

## 2. Canonical Orchestrators (Must Call Engines Above)

| Orchestrator | Location |
|--------------|----------|
| PRV_24H_Orchestrator.py | PX_Executive/ |
| PX_Live_Orchestrator_v2.py | PX_Executive/orchestrators/ |
| run_genesis_feed.py | PX_Executive/ |
| run_repurposed_feed.py | PX_Executive/ |
| run_prv_novel.py | PX_Executive/ |
| run_prv_repurposed.py | PX_Executive/ |
| run_finalize_dossiers.py | PX_Executive/ |
| run_one_cycle_test.py | PX_Executive/ |

These orchestrators must:

- Import the engines they require
- Call the engines in the correct stage order
- Pass validated data between engines
- Enforce governance (ZeusLaws: `run_zeus_gate` or `check_constitutional`)
- Produce outputs consumed by downstream engines

---

## 3. Engine → Orchestrator Connectivity Matrix

### A. Feed Stage

**Orchestrators:** `run_genesis_feed.py`, `run_repurposed_feed.py`

**Must call:**

- Genesis_Engine
- Vector_Core (via PRV path or Refinery; Genesis feed uses OPE + Vector Core for every novel SMILES)
- Metabolism
- Trajectory_Predictor (downstream in 24H after WorldLine persist)

### B. Novel / Repurposed Stage (PRV)

**Orchestrators:** `run_prv_novel.py`, `run_prv_repurposed.py`, `PRV_24H_Orchestrator.py`

**Must call:**

- OBE, OCE, OLE, OME, OPE, OSE
- ADMET
- PKPD or PKPD_Simple (governance: prefer PKPD)
- DoseOptimizer or DoseOptimizer_v2 (governance: DoseOptimizer_v2)
- VirtualEfficacy_Simple or VirtualEfficacyAnalytics (governance: VirtualEfficacyAnalytics)
- GradingEngine (via Finalization_Pipeline)
- GradingSchema_Discovery.json (loaded by GradingEngine)

### C. Trial Stage

**Orchestrators:** `PX_Live_Orchestrator_v2.py`, `run_one_cycle_test.py` (indirect via PRV → Finalization)

**Must call:**

- TrialEngine
- Simulation_Engine (e.g. Gold_Rush_Miner; demos)
- VirtualEfficacyAnalytics
- Manufacturing_Manifest (reachable; validated by PX_System_Test)

### D. Finalization Stage

**Orchestrator:** `run_finalize_dossiers.py` (uses `PX_Warehouse.Finalization_Pipeline`)

**Must call:**

- Evidence_Package (PX_System.foundation.Evidence_Package)
- WorldLine_Database (PX_Warehouse.WorldLine_Database)
- GradingEngine
- GradingSchema_Discovery.json (via GradingEngine)
- ZeusLaws.run_zeus_gate (PX_System.foundation.ZeusLaws)

---

## 4. IDE Enforcement Rules

**Rule 1 — No Orphan Engines**  
Every engine in `PX_Engine/`, `PX_Engine/operations/`, and `PX_Laboratory/` must be imported by at least one orchestrator or by a test in `PX_Validation/tests/` or `tests/`.

**Rule 2 — No Dead Orchestrator Paths**  
Every canonical orchestrator must call at least one engine from each required stage for its role (Feed / PRV / Trial / Finalization).

**Rule 3 — Schema Enforcement**  
`GradingSchema_Discovery.json` must be loaded by GradingEngine. Finalization (run_finalize_dossiers.py → Finalization_Pipeline) uses GradingEngine.

**Rule 4 — Governance Enforcement**  
Every orchestrator that gates on constitution must use:

- `from PX_System.foundation.ZeusLaws import run_zeus_gate` (Finalization_Pipeline), or  
- `from PX_System.foundation.ZeusLaws import check_constitutional` (PRV_24H_Orchestrator),

or equivalent governance hook. No bypass of Zeus gate.

**Rule 5 — Data Flow Integrity**  
Outputs from Feed → PRV → Trial → Finalization must be passed along the pipeline, not regenerated or bypassed.

**Rule 6 — No Legacy Imports**  
The IDE must block imports from `99_LEGACY_CODE/`. Do not import or run from there.

**Rule 7 — Validation Hooks**  
Every canonical orchestrator must be reachable or covered by:

- `PX_Validation/tests/run_all_tests.py`
- `run_e2e_layers.py`

---

## 5. IDE Connectivity Test (Auto-Generated)

Run:

```bash
python PX_Validation/tests/run_all_tests.py
python run_e2e_layers.py
```

If both pass:

- All engines are connected
- All orchestrators are functional
- Governance is enforced
- No dead modules exist
- No bypasses exist
- Full throughput is active

---

## 6. Status

When the test suite and E2E layers both pass, the system is at **100% throughput**. See `STATUS.md` and `JUSTFILE.md` for commands (`just test`, `just govern`).
