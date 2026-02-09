# PX Warehouse Migration Instruction

**Authoritative Migration Directive for Converting the Current Warehouse Into the Canonical Enterprise Warehouse**

This document defines the **mandatory migration steps** required to bring the existing `PX_Warehouse/` directory into full alignment with:

- **PX_WAREHOUSE_ENTERPRISE_MODEL.md**
- **PX_ENTERPRISE_LAYOUT.md**
- **PX_ENTERPRISE_EDITION.md**

The IDE must treat this as a **governed, enforceable** migration plan.

---

## §1. Canonical Warehouse Structure (Target State)

After migration, `PX_Warehouse/` must contain **only** the following top-level directories:

```
PX_Warehouse/
    CommercialAssets/
        Diamond/
        Gold/
        Silver/
        Bronze/
        Dossier_Final/
        Executive_Summary/
        Audit_Trail/
        Learning_Material/

    WorldLines/
        Diamond/
        Gold/
        Silver/
        Bronze/

    TrialSimulations/
        LiveRuns/
        (trial dossiers)

    Operations/
        Control_Scripts/
        Inputs/
        manifests/
        reports/
        scripts/
        System_Health/

    Archive/
        Snapshot_<timestamp>/
        Dossiers/
        TrialSimulations/
        WorldLines/
        CommercialAssets/

    ENTERPRISE_MODEL_README.md
```

**Everything else** must be moved, renamed, or archived.

---

## §2. Mandatory Folder Consolidation

The IDE must perform the following consolidations.

### 2.1 Remove legacy world-line structures

The following must be removed or migrated:

| Item | Action |
|------|--------|
| `Operational_Data/WorldLine_Records/` | Move contents into `WorldLines/<tier>/` |
| `WorldLine_Database.py` | Move to `PX_Warehouse/Operations/scripts/` or archive if deprecated |
| Any folder named `WorldLine_*` outside `WorldLines/` | Migrate or archive |

### 2.2 Remove legacy commercial asset structures

The following must be removed or migrated:

- Any **CommercialAssets/** subfolder **not** in the canonical list (§1) → rename, merge, or archive.
- Any **Learning_Material** folder outside `CommercialAssets/Learning_Material/` → move to canonical location or archive.
- Any **Executive_Summary** or **Audit_Trail** folder outside canonical location → move or archive.

### 2.3 Remove legacy logs

The folder:

```
PX_Warehouse/logs/
```

contains v3-era analytics and must be moved to:

```
PX_LOGS/archive/warehouse_legacy/
```

### 2.4 Remove duplicate or obsolete folders

The IDE must archive (or merge then remove):

- **Operational_Data/** (after extracting world-line records into WorldLines)
- **Learning_Material/** (if duplicate of CommercialAssets/Learning_Material)
- **__pycache__/**
- Any folder **not** listed in §1

---

## §3. Mandatory File Migration

### 3.1 Trial Simulation Dossiers

All trial simulation dossiers must be moved to:

```
PX_Warehouse/TrialSimulations/
```

(Subfolders such as `LiveRuns/` as defined in §1.)

### 3.2 Commercial Dossiers

All commercial dossiers must be moved to:

```
PX_Warehouse/CommercialAssets/Dossier_Final/
```

(Tier subfolders under Dossier_Final if needed; otherwise flat or by timestamp.)

### 3.3 World-Line Snapshots

All snapshot folders must be moved to:

```
PX_Warehouse/Archive/Snapshot_<timestamp>/
```

(Preserve existing timestamp in folder name or consolidate under a single Snapshot_<timestamp>.)

### 3.4 Legacy analytics

Files such as:

- predator_x_v3_audit
- stratification
- monetization
- deduplication
- portfolio_analytics.json
- stability_metrics.json
- drift_metrics.json
- EXECUTIVE_BENCHMARK_REPORT.*
- PORTFOLIO_DASHBOARD.*

must be moved to:

```
PX_LOGS/archive/warehouse_legacy/
```

---

## §4. Tier Assignment Rules

The IDE must classify assets into:

- **Diamond**
- **Gold**
- **Silver**
- **Bronze**

based on (in order of precedence):

1. Existing folder name (if asset is already in a tier-named folder)
2. Metadata in dossiers (e.g. `grade`, `tier`, `predator_x_v3.grade`)
3. Metadata in world-line records (e.g. `physical_realization.toxicity_tier`)

**Default to Bronze** if no tier is known.

---

## §5. Lifecycle Enforcement

The IDE must enforce:

### 5.1 Created → CommercialAssets

Any asset **without** validation or dossier must be placed in:

```
PX_Warehouse/CommercialAssets/<tier>/
```

### 5.2 Validated → WorldLines

Any asset **with** validation metadata must be placed in:

```
PX_Warehouse/WorldLines/<tier>/
```

### 5.3 Dossier Created → Archive

Any asset **with** a dossier must be moved to:

```
PX_Warehouse/Archive/<appropriate_subfolder>/
```

(e.g. `Archive/Dossiers/`, `Archive/TrialSimulations/`, `Archive/CommercialAssets/`, or `Archive/Snapshot_<timestamp>/` as appropriate.)

---

## §6. Cleanup Rules

The IDE must:

- **Delete** all `__pycache__/` folders under PX_Warehouse.
- **Delete** all `.tmp`, `.cache`, `.log` files under PX_Warehouse that are **not** part of canonical logs (canonical logs live in PX_LOGS).
- **Remove** any folder not listed in §1 (by moving to Archive or PX_LOGS/archive/warehouse_legacy, then removing the duplicate).
- **Remove** any file not referenced by the enterprise model (or archive first, then remove).

---

## §7. Script Updates

The IDE must update:

- **unified_warehouse_consolidation.py** (root and/or PX_Warehouse)
- **PX_Warehouse/Operations/scripts/** (all scripts referencing warehouse paths)
- Any other script referencing old paths (e.g. WorldLine_Records, Operational_Data, legacy CommercialAssets layout, PX_Warehouse/logs)

…to use the **canonical structure** defined in §1 and the vocabulary in **PX_WAREHOUSE_ENTERPRISE_MODEL.md**.

---

## §8. Verification Step

After migration, the IDE must:

1. **Compare** the resulting structure to §1 (exact top-level list; allowed subfolders as specified).
2. **Confirm** no legacy folders remain at PX_Warehouse top level (e.g. no `logs/`, no `Operational_Data/`, no `WorldLine_Records/`).
3. **Confirm** all assets are tiered (in a Diamond/Gold/Silver/Bronze path where applicable).
4. **Confirm** all dossiers are in CommercialAssets/Dossier_Final or Archive (as per lifecycle).
5. **Confirm** all world-lines are in the canonical WorldLines tree (WorldLines/<tier>/).
6. **Confirm** all warehouse-origin logs are in PX_LOGS (e.g. PX_LOGS/archive/warehouse_legacy).
7. **Confirm** no duplicate folder names exist (e.g. two "Gold" at same level in different organs is allowed; two "CommercialAssets" is not).
8. **Confirm** **PX_WAREHOUSE_ENTERPRISE_MODEL.md** is satisfied.

If any mismatch is found, the IDE must produce a **Warehouse Drift Report** using **WAREHOUSE_DRIFT_REPORT_TEMPLATE.md** (e.g. `PX_Warehouse/WAREHOUSE_DRIFT_REPORT_<timestamp>.md`). Verification steps are defined in **WAREHOUSE_DRIFT_REPORT_EXECUTION_LOGIC.md**; enforcement steps are defined in **PX_WAREHOUSE_IDE_ENFORCEMENT_SPEC.md**.

---

## §9. Lock Clause

Once migration is complete:

- The warehouse is **locked** to the canonical structure.
- Any **new folder** requires amendment to **PX_WAREHOUSE_ENTERPRISE_MODEL.md** (and this instruction if the migration target changes).
- Any **new asset** must follow the lifecycle rules (§5).
- Any **drift** from §1 or the enterprise model triggers **governance failure** and must be reported (Warehouse Drift Report) and remediated.

---

**Related:** **PX_WAREHOUSE_FULL_REBUILD_DIRECTIVE.md** defines the full rebuild command (recursive scan, reclassification, renaming, legacy removal, **system-wide reference update**). Use it when a complete warehouse rebuild and cross-system path update is required.

**Document version:** 1.0  
**Governed by:** PX_WAREHOUSE_ENTERPRISE_MODEL.md, PX_ENTERPRISE_LAYOUT.md, PX_ENTERPRISE_EDITION.md
