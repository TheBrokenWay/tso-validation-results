# PX_WAREHOUSE_REDUNDANCY_SCANNER_SPEC

**Safe-Mode Redundancy, Legacy Detection & Pre-Action Testing Engine**

This specification defines the **governed, step-wise process** the IDE must follow to detect redundant, obsolete, or unused files inside every subfolder of `PX_Warehouse/` **without modifying anything** until each step is validated.

Aligned with: **PX_WAREHOUSE_ENTERPRISE_MODEL.md**, **PX_WAREHOUSE_MIGRATION_INSTRUCTION.md**, **PX_WAREHOUSE_IDE_ENFORCEMENT_SPEC.md**.

---

## §1. Safe-Mode Execution Model

The scanner must operate in **three phases**:

| Phase | Mode | Purpose |
|-------|------|---------|
| **Discovery** | Read-only | Scan and classify; build dependency graph |
| **Simulation** | No changes | Propose actions; produce Simulation Report |
| **Execution** | Apply only validated actions | Archive, delete, or keep only after validation |

**No file or folder may be altered until:**

1. It is **classified**  
2. It is **validated**  
3. It is **simulated**  
4. It is **approved** by the migration directive  

---

## §2. Discovery Phase (Read-Only)

The IDE scans **every subfolder** under:

```
PX_Warehouse/
```

including:

- CommercialAssets  
- WorldLines  
- TrialSimulations  
- Operations  
- Archive  
- Any legacy or unexpected folder  

**For each file**, the IDE collects:

- **Path**  
- **Size**  
- **Modified timestamp**  
- **File type**  
- Whether **referenced by any script**  
- Whether **referenced by any dossier**  
- Whether **referenced by any world-line record**  
- Whether **referenced by any consolidation script**  
- Whether **referenced by any trial simulation**  
- Whether **referenced by any commercial asset**  
- Whether **referenced by any PK/PD run**  
- Whether **referenced by any evidence package**  

This produces a **complete dependency graph**. No files or folders are modified in this phase.

---

## §3. Redundancy Classification Rules

Each file is classified into **exactly one** of the following:

### 3.1 ACTIVE

Used by the current enterprise pipeline (scripts, orchestrator, consolidation, trial engine, evidence package, etc.).

### 3.2 REFERENCED

Not active in the pipeline, but **referenced by**:

- a dossier  
- a world-line  
- a trial simulation  
- a consolidation script  
- a commercial asset  
- a PK/PD run  
- an evidence package  

### 3.3 LEGACY

Matches known legacy patterns:

- v3 analytics  
- old dashboards  
- old benchmark reports  
- monetization / stratification  
- deduplication logs  
- pre-enterprise snapshots  
- obsolete world-line formats  

### 3.4 DUPLICATE

Same **content hash** as another file in the warehouse (redundant copy).

### 3.5 UNUSED

**Not** referenced by anything in the dependency graph.

### 3.6 UNKNOWN

Cannot be classified with the rules above; **requires manual review**.

---

## §4. Simulation Phase (No Changes)

For each file, the IDE generates a **proposed action**:

| Classification | Proposed Action |
|----------------|------------------|
| ACTIVE | **KEEP** |
| REFERENCED | **KEEP** |
| LEGACY | **ARCHIVE** |
| DUPLICATE | **DELETE** (after confirming which copy to keep) |
| UNUSED | **DELETE** |
| UNKNOWN | **REVIEW** |

The IDE then produces a **Simulation Report** (per file or aggregated):

```
SIMULATION REPORT
-----------------
File: <path>
Classification: <ACTIVE|REFERENCED|LEGACY|DUPLICATE|UNUSED|UNKNOWN>
Proposed Action: <KEEP|ARCHIVE|DELETE|REVIEW>
Reason: <detailed reason>
Dependencies: <list of referencing files or NONE>
```

**No changes are made in this phase.**

---

## §5. Execution Phase (Validated Actions Only)

The IDE applies actions **only after** simulation passes and each item is validated.

### 5.1 KEEP

Do nothing.

### 5.2 ARCHIVE

Move file to:

```
PX_Warehouse/Archive/Legacy/<timestamp>/
```

(or `PX_LOGS/archive/warehouse_legacy/` if the migration directive specifies logs/analytics there.)

### 5.3 DELETE

Delete **only if**:

- File is **UNUSED** or **DUPLICATE**  
- Simulation confirms **no dependencies**  
- Migration directive **allows** deletion  
- Step-wise tests (§6) pass  

### 5.4 REVIEW

Leave **untouched**; flag for manual inspection. Do not archive or delete.

---

## §6. Step-Wise Testing Before Alteration

Before **altering** any file (archive or delete), the IDE must:

1. **Re-scan dependencies** — confirm no new references appeared  
2. **Confirm** the file is not part of a current run  
3. **Confirm** the file is not required by any script  
4. **Confirm** the file is not required by any dossier  
5. **Confirm** the file is not required by any world-line  
6. **Confirm** the file is not required by any trial simulation  
7. **Confirm** the file is not required by any commercial asset, PK/PD run, or evidence package  

**Only then** may the action (ARCHIVE or DELETE) proceed.

---

## §7. Drift Protection

After execution, the IDE must:

1. **Re-scan** the warehouse  
2. **Re-classify** all remaining files  
3. **Confirm** no new drift  
4. **Confirm** no missing dependencies  
5. **Confirm** canonical structure remains intact  

If any issue is detected, the IDE must generate a **Warehouse Drift Report** (using **WAREHOUSE_DRIFT_REPORT_TEMPLATE.md**).

---

## §8. Lock Clause

Once redundancy cleanup is complete:

- **Warehouse is locked** to the canonical structure and classification rules  
- **No new legacy files** may be introduced without going through classification  
- **Any new file** must be classified as ACTIVE or REFERENCED (or REVIEW if unknown)  
- **Any drift** (e.g. unclassified file, broken reference, non-canonical path) triggers **governance failure** and a drift report  

---

**Governed by:** PX_WAREHOUSE_ENTERPRISE_MODEL.md, PX_WAREHOUSE_MIGRATION_INSTRUCTION.md, PX_WAREHOUSE_IDE_ENFORCEMENT_SPEC.md

**Outputs:** Simulation Report (per file or summary); optional Redundancy Cleanup Report; Warehouse Drift Report if §7 detects issues.
