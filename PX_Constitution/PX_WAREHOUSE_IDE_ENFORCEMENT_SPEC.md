# PX_WAREHOUSE_IDE_ENFORCEMENT_SPEC

**Authoritative IDE Execution Script for Warehouse Migration & Drift Elimination**  
**Enterprise Edition — Constitutional Enforcement Logic**

This specification defines the **exact operational steps** the IDE must execute to:

1. Migrate the current warehouse into the canonical structure  
2. Apply lifecycle and tiering rules  
3. Archive legacy content  
4. Update scripts  
5. Verify zero drift  
6. Produce a **Warehouse Drift Report** if needed  

The IDE must treat this specification as **mandatory**.

---

## §1. Initialization

The IDE loads:

- **PX_WAREHOUSE_ENTERPRISE_MODEL.md**
- **PX_WAREHOUSE_MIGRATION_INSTRUCTION.md**
- **PX_ENTERPRISE_LAYOUT.md**
- **DOSSIER_COMMERCIALASSETS_FILING_RULES.md** (Dossier_Final admission + CommercialAssets biopath routing; fail-closed)

It constructs the **canonical warehouse tree** from these documents.

It then scans:

```
PX_Warehouse/
```

and builds the **actual warehouse tree**.

---

## §2. Folder-Level Enforcement

The IDE must enforce the canonical folder structure exactly.

### 2.1 Create missing canonical folders

If any canonical folder is missing, create it:

```
CommercialAssets/
WorldLines/
TrialSimulations/
Operations/
Archive/
```

And their required subfolders:

```
CommercialAssets/{Diamond,Gold,Silver,Bronze,Dossier_Final,Executive_Summary,Audit_Trail,Learning_Material}
WorldLines/{Diamond,Gold,Silver,Bronze}
TrialSimulations/LiveRuns
Archive/{Dossiers,TrialSimulations,WorldLines,CommercialAssets}
```

*(Snapshot_<timestamp> under Archive is created when migrating snapshots.)*

### 2.2 Remove or migrate non-canonical folders

Any folder **not** in the canonical list must be:

- **moved** (if it contains assets)  
- **archived** (if legacy)  
- **deleted** (if cache/temp)  

Examples:

- **Operational_Data/** → extract world-line records → migrate → delete  
- **WorldLine_Records/** → migrate to WorldLines/<tier>/  
- **logs/** → move to PX_LOGS/archive/warehouse_legacy/  
- **Learning_Material/** (duplicate) → merge into canonical location  
- **__pycache__/** → delete  

---

## §3. File-Level Enforcement

### 3.1 Trial Simulation Dossiers

Move all trial simulation dossiers to:

```
PX_Warehouse/TrialSimulations/
```

### 3.2 Commercial Dossiers

Move all commercial dossiers to:

```
PX_Warehouse/CommercialAssets/Dossier_Final/
```

### 3.3 World-Line Snapshots

Move all snapshot folders to:

```
PX_Warehouse/Archive/Snapshot_<timestamp>/
```

### 3.4 Legacy Analytics

Move all v3-era analytics to:

```
PX_LOGS/archive/warehouse_legacy/
```

This includes: monetization; stratification; deduplication; portfolio analytics; stability metrics; drift metrics; executive benchmark reports; dashboards.

---

## §4. Tier Assignment Logic

For each asset, the IDE determines its tier using the following priority:

1. Explicit metadata (in dossier or world-line record)  
2. Folder name (if already tiered)  
3. Heuristic (if metadata indicates quality)  
4. **Default: Bronze**  

The IDE then moves the asset to:

```
CommercialAssets/<tier>/
```

or

```
WorldLines/<tier>/
```

depending on lifecycle state.

---

## §5. Lifecycle Enforcement Logic

The IDE determines lifecycle state:

### 5.1 Created

No validation metadata → must be in:

```
CommercialAssets/<tier>/
```

### 5.2 Validated

Has validation metadata → must be in:

```
WorldLines/<tier>/
```

### 5.3 Dossier Created

Has dossier → must be in:

```
Archive/<appropriate_subfolder>/
```

Any mismatch triggers a **drift entry**.

---

## §6. Cleanup Logic

The IDE must remove:

- **__pycache__/**
- **.tmp**, **.cache**, **.log** (non-canonical; canonical logs live in PX_LOGS)
- Any **folder** not in canonical structure (after migrating contents)
- Any **file** not referenced by enterprise model (after archiving if needed)

---

## §7. Script Update Logic

The IDE must update all scripts referencing old paths:

- **unified_warehouse_consolidation.py**
- **PX_Warehouse/Operations/scripts/***

and any script referencing:

- Operational_Data  
- WorldLine_Records  
- logs/  
- CommercialAssets duplicates  
- old snapshot paths  

All paths must be rewritten to the **canonical structure**.

---

## §8. Verification Logic

After migration, the IDE:

1. Re-scans the warehouse  
2. Re-builds the actual tree  
3. Compares against canonical tree  
4. Confirms: no legacy folders; no duplicate structures; no unclassified assets; no lifecycle violations; no tiering violations; no unexpected files; **no drift**  

If any mismatch remains, the IDE generates:

```
PX_Warehouse/WAREHOUSE_DRIFT_REPORT_<timestamp>.md
```

using **WAREHOUSE_DRIFT_REPORT_TEMPLATE.md**.

---

## §9. Lock Enforcement

Once verification passes:

- **Warehouse is locked**  
- Any **new folder** requires amendment to the enterprise model  
- Any **drift** triggers governance failure  
- **IDE must run drift verification on every commit** (or on every warehouse-touching change, as configured)  

---

**Related:** **PX_WAREHOUSE_REDUNDANCY_SCANNER_SPEC.md** defines the safe-mode redundancy scanner (Discovery → Simulation → Execution) for detecting redundant, legacy, or unused files without modifying anything until validated. **PX_WAREHOUSE_FULL_REBUILD_DIRECTIVE.md** defines the full rebuild command, including **system-wide reference update** (all scripts and components that connect to the warehouse must be updated to canonical paths).

**Governed by:** PX_WAREHOUSE_ENTERPRISE_MODEL.md, PX_WAREHOUSE_MIGRATION_INSTRUCTION.md, PX_ENTERPRISE_LAYOUT.md
