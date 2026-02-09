# PX_WAREHOUSE_FULL_REBUILD_DIRECTIVE

**Authoritative Command for Complete Warehouse Rebuild, Reclassification, Renaming, Legacy Removal, and System-Wide Reference Update**  
**Enterprise Edition — Final, Locked Instruction for the IDE**

This directive instructs the IDE to **fully rebuild** the PX_Warehouse, **reclassify** all assets, **remove** all legacy content, and **update every system that connects to the warehouse** so that the entire platform uses the new canonical structure.

This directive **supersedes** all legacy warehouse structures and all legacy warehouse references.

Aligned with: **PX_WAREHOUSE_ENTERPRISE_MODEL.md**, **PX_WAREHOUSE_MIGRATION_INSTRUCTION.md**, **PX_WAREHOUSE_IDE_ENFORCEMENT_SPEC.md**, **PX_WAREHOUSE_REDUNDANCY_SCANNER_SPEC.md**.

---

## §1. Canonical Warehouse Structure (Target State)

After rebuild, the warehouse **must contain only**:

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

**No other folders** may remain at PX_Warehouse top level (or outside the canonical tree). Any other folder must be migrated, archived, or removed per this directive.

---

## §2. Full Recursive Scan (All Files, Folders, Subfolders)

The IDE must **recursively scan every item** under:

```
PX_Warehouse/
```

and **classify each file and folder** using:

- **ACTIVE**
- **REFERENCED**
- **VALIDATED**
- **DOSSIER**
- **LEGACY**
- **DUPLICATE**
- **UNUSED**
- **UNKNOWN**

This includes:

- CommercialAssets (all subfolders and files)
- WorldLines (all subfolders and files)
- TrialSimulations (all subfolders and files)
- Operations (all subfolders and files)
- Archive (all subfolders and files)
- v3 analytics
- snapshots
- world-line records
- dossiers (e.g. 16,510+ or current count)
- GOLD_CORE, DIAMOND_PRV, BRONZE inside SILVER, etc.
- **any folder not in the canonical model**

---

## §3. Tier Assignment Rules

Tier is determined by (in order of precedence):

1. **Metadata** (in dossier or world-line record)
2. **Folder name** (if valid: Diamond, Gold, Silver, Bronze)
3. **Heuristic** (if metadata indicates quality)
4. **Default: Bronze**

**Valid tiers only:**

- Diamond  
- Gold  
- Silver  
- Bronze  

Any other tier name must be **renamed** to one of the above (or mapped via a deterministic rule).

---

## §4. Lifecycle Rules

The IDE must enforce:

| State | Location |
|-------|----------|
| **Created** | CommercialAssets/<tier>/ |
| **Validated** | WorldLines/<tier>/ |
| **Dossier** | Archive/<appropriate_subfolder>/ |

**Any mismatch must be corrected** (move asset to correct location per lifecycle state).

---

## §5. Renaming Rules

The IDE must **rename** any folder or file that uses:

- non-canonical vocabulary  
- mixed-era naming  
- inconsistent tier names  
- legacy prefixes  

**Examples:**

| Current | Canonical |
|---------|-----------|
| GOLD_CORE | Gold/ |
| DIAMOND_PRV | Diamond/ |
| BRONZE inside Silver/ | move to Bronze/ |
| WorldLine_Records | WorldLines/<tier>/ |
| Commercial_Assets | CommercialAssets/ |
| WorldLine_Database (if folder) | per §1 or archive |

---

## §6. Relocation Rules

The IDE must **move** files to their correct canonical location:

| Content Type | Destination |
|--------------|-------------|
| Commercial assets | CommercialAssets/<tier>/ |
| World-lines | WorldLines/<tier>/ |
| Trial sims | TrialSimulations/ (e.g. LiveRuns/) |
| Dossiers | CommercialAssets/Dossier_Final/ (or Archive/Dossiers/ when finalized) |
| Snapshots | Archive/Snapshot_<timestamp>/ |
| Legacy analytics | PX_LOGS/archive/warehouse_legacy/ |

---

## §7. Legacy Removal Rules

The IDE must **remove** (after migrating contents where required):

- **__pycache__/**
- **Operational_Data/** (after extracting world-line records into WorldLines/<tier>/)
- **WorldLine_Records/** (after migrating to WorldLines/<tier>/)
- **logs/** under PX_Warehouse (move to PX_LOGS/archive/warehouse_legacy/)
- v3 analytics (move to PX_LOGS/archive/warehouse_legacy/)
- old dashboards (move or delete per classification)
- old benchmark reports (move or delete per classification)
- **duplicate** Learning_Material folders (merge into CommercialAssets/Learning_Material/)
- **any folder not in §1** (migrate, archive, or delete per classification)

---

## §8. Safe-Mode Testing (Before Any Change)

**Before altering any file or folder**, the IDE must:

1. **Re-scan dependencies**  
2. **Confirm** no script references (or update script first)  
3. **Confirm** no dossier references  
4. **Confirm** no world-line references  
5. **Confirm** no trial simulation references  
6. **Confirm** no pipeline references  
7. **Confirm** no active run references  

**If any dependency exists**, the action is **blocked** until the dependency is updated or the dependency graph is re-validated.

---

## §9. Simulation Phase (Dry Run)

The IDE must generate a **Simulation Report** listing:

- every **proposed move**  
- every **proposed rename**  
- every **proposed archive**  
- every **proposed deletion**  
- every **tier assignment**  
- every **lifecycle correction**  
- every **script update**  

**No changes occur in this phase.** The report is for validation and approval before execution.

---

## §10. Execution Phase (After Simulation Approval)

The IDE **applies** (only after simulation is approved):

- **MOVE** — relocate files/folders to canonical paths  
- **RENAME** — apply canonical vocabulary and tier names  
- **ARCHIVE** — move legacy to Archive/Legacy/<timestamp>/ or PX_LOGS/archive/warehouse_legacy/  
- **DELETE** — remove __pycache__, duplicates, or unused only when safe  
- **TIER** — assign and place assets in correct tier folder  
- **LIFECYCLE FIXES** — move assets to correct lifecycle location  
- **SCRIPT UPDATES** — see §11  

**All actions must match the simulation.** No ad-hoc changes.

---

## §11. System-Wide Reference Update (MANDATORY)

**Any system that connects to the warehouse must be updated so it has the updated information.**

The IDE must update **every component** that references the warehouse, including but not limited to:

- **unified_warehouse_consolidation.py** (root and/or PX_Warehouse)
- **PX_Warehouse/Operations/scripts/** *
- **PX_Warehouse/Operations/Control_Scripts/** *
- **PX_Warehouse/Operations/manifests/** * (if they contain paths)
- **PX_Warehouse/Operations/reports/** * (if they contain paths)
- **PX_Warehouse/Operations/System_Health/** * (if present and referencing paths)
- **PX_Executive/Acquisition/** *
- **PX_Engine/** * (if referencing warehouse paths)
- **PX_System/** * (if referencing warehouse paths)
- **PX_Executive/** * (orchestrators, pipelines, miners)
- **PredatorXOrchestratorV2** (PX_Executive/orchestrators/PX_Live_Orchestrator_v2.py and related)
- **Any script** referencing: CommercialAssets, WorldLines, TrialSimulations, Operations, Archive, snapshots, world-line records, legacy folders, old analytics paths

**Requirement:**

- **All references** must be **rewritten** to the **canonical structure** defined in §1.  
- **No script** may reference: legacy folder names, legacy paths, legacy tiers, legacy analytics, legacy snapshots, legacy world-line structures.

---

## §12. Post-Rebuild Verification

After execution, the IDE must:

1. **Re-scan** the warehouse  
2. **Re-classify** all files  
3. **Confirm** canonical structure (§1)  
4. **Confirm** no legacy folders  
5. **Confirm** no duplicates (or only intentional copies in Archive)  
6. **Confirm** no unexpected files  
7. **Confirm** no outdated script references (system-wide)  
8. **Confirm** no drift  

If any mismatch remains, generate:

```
PX_Warehouse/WAREHOUSE_DRIFT_REPORT_<timestamp>.md
```

using **WAREHOUSE_DRIFT_REPORT_TEMPLATE.md**.

---

## §13. Lock Clause

After **successful** rebuild:

- **Warehouse is locked** to the canonical structure.  
- **Any new folder** requires amendment to the enterprise model (PX_WAREHOUSE_ENTERPRISE_MODEL.md and migration instruction).  
- **Any drift** triggers **governance failure** and a drift report.  
- **IDE must verify warehouse on every commit** (or on every warehouse-touching change, as configured).  

---

**Governed by:** PX_WAREHOUSE_ENTERPRISE_MODEL.md, PX_WAREHOUSE_MIGRATION_INSTRUCTION.md, PX_WAREHOUSE_IDE_ENFORCEMENT_SPEC.md, PX_WAREHOUSE_REDUNDANCY_SCANNER_SPEC.md

**Outputs:** Simulation Report (dry run); post-rebuild verification; WAREHOUSE_DRIFT_REPORT_<timestamp>.md if §12 finds mismatch.
