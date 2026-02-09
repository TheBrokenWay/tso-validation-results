# PX Warehouse — Enterprise Edition Model

**Single source of truth for the governed, unified warehouse.**

This document defines the **one warehouse**, **one CommercialAssets**, **one WorldLines**, **one Archive**, and **one terminology** model. No duplicates, no parallel eras, no legacy folder names. Aligns with **PX_ENTERPRISE_LAYOUT.md** and **PX_ENTERPRISE_EDITION.md**.

---

## 1. ONE Warehouse

All warehouse content lives under:

```
PX_Warehouse/
```

- **No duplicates.** No parallel "Commercial_Assets" or "WorldLine_Records" at root.
- **No parallel eras.** No mixed legacy and new structures.
- **No legacy folder names** that conflict with this model.

Everything warehouse-related is under `PX_Warehouse/`.

---

## 2. ONE CommercialAssets Structure

**Purpose:** Active, in-progress, or pre-dossier commercial assets only. **Not** long-term storage.

**Canonical structure:**

```
PX_Warehouse/CommercialAssets/
    Gold/
    Silver/
    Bronze/
    Diamond/           (optional — exceptional assets)
    Dossier_Final/      (finalized dossiers only, before move to Archive)
    Executive_Summary/
    Audit_Trail/
    Learning_Material/
```

**What belongs here:**

- Assets being validated, graded, or prepared for dossier.
- Learning material and executive summaries for active campaigns.
- Audit trail for in-flight operations.
- Dossier_Final holds dossiers **until** they are moved to Archive. **Admission:** Only assets that have passed all checks (scientific, ADMET, legal/FTO, commercial, provenance, timestamp, reproducibility) may enter Dossier_Final. See **DOSSIER_COMMERCIALASSETS_FILING_RULES.md** and **Filing_Rules.py**.

**What does NOT stay here:**

- Once an asset is **validated**, **verified**, and **given a dossier**, it moves to **Archive**.
- CommercialAssets is **not** long-term storage. Completed work goes to Archive.

---

## 3. ONE WorldLines Structure

**Purpose:** The reference library of world-line data, sorted by quality tier.

**Canonical structure:**

```
PX_Warehouse/WorldLines/
    Gold/
    Silver/
    Bronze/
    Diamond/
```

- **One** WorldLines tree. No separate "WorldLine_Records", "Operational_Data", or snapshot folders as alternate world-line stores.
- Tiered subfolders (Gold, Silver, Bronze, Diamond) hold the **reference** copies.

**Lifecycle:**

- After **dossier creation**, the world-line snapshot moves to **Archive** (e.g. `Archive/Snapshot_<timestamp>/`).
- The reference copy may remain in WorldLines for lookup; the **authoritative completed artifact** is in Archive.

---

## 4. ONE Archive

**Purpose:** Final resting place for completed work. **Append-only and immutable.**

**Canonical location:**

```
PX_Warehouse/Archive/
    Snapshot_<timestamp>/   (world-line snapshots post-dossier)
    Dossiers/               (completed commercial dossiers)
    TrialSimulations/       (completed trial runs)
    EvidencePackages/       (completed evidence packages)
    PKPD_Runs/              (completed PK/PD runs)
    …
```

**What belongs here:**

- Completed dossiers  
- Completed trial simulations  
- Completed world-line snapshots  
- Completed commercial assets (post-dossier)  
- Completed PK/PD runs  
- Completed evidence packages  

**Governance:** Archive is append-only. No in-place edits. New completions create new entries (e.g. new snapshot or dated subfolder).

---

## 5. ONE Terminology Set

**Eliminate:**

- "CommercialAssets" vs "Commercial Assets"
- "WorldLines" vs "WorldLine_Records"
- "Archive" vs "Archive_Snapshots"
- "Gold/Silver/Bronze" vs "Diamond/Gold/Silver/Bronze"
- "Operational_Data" vs "WorldLine_Records"

**Use exactly:**

| Term | Meaning |
|------|---------|
| **Asset Tiers** | Diamond, Gold, Silver, Bronze |
| **CommercialAssets** | Active/pre-dossier commercial assets (see §2) |
| **WorldLines** | Reference library of world-lines, tiered (see §3) |
| **TrialSimulations** | Trial simulation outputs (e.g. under Warehouse or Archive) |
| **Operations** | Scripts, inputs, control (e.g. `PX_Warehouse/Operations/`) |
| **Archive** | Append-only, immutable completed artifacts (see §4) |

**Lifecycle (single vocabulary):**

1. **Created** → CommercialAssets (or ingestion path into CommercialAssets).
2. **Validated** → WorldLines (tiered: Gold/Silver/Bronze/Diamond).
3. **Dossier created** → Asset/snapshot moves to **Archive**.

---

## 6. Summary Diagram

```
PX_Warehouse/
├── CommercialAssets/     # Active only; completed → Archive
│   ├── Gold/
│   ├── Silver/
│   ├── Bronze/
│   ├── Diamond/
│   ├── Dossier_Final/
│   ├── Executive_Summary/
│   ├── Audit_Trail/
│   └── Learning_Material/
├── WorldLines/          # Reference library; tiered
│   ├── Gold/
│   ├── Silver/
│   ├── Bronze/
│   └── Diamond/
├── Archive/             # Append-only; completed dossiers, snapshots, trials
│   ├── Snapshot_<timestamp>/
│   ├── Dossiers/
│   └── …
├── TrialSimulations/    # Trial run outputs (or under Archive when completed)
├── Operations/          # Scripts, Inputs, Control_Scripts
└── (no duplicate or legacy top-level folders)
```

---

## 7. Alignment with Enterprise Edition

- **PX_ENTERPRISE_LAYOUT.md:** PX_Warehouse is the single warehouse organ; this document defines its internal model.
- **PX_ENTERPRISE_EDITION.md:** Governance and determinism apply; Archive is immutable, lifecycle is deterministic.
- **Constitutional tests:** No references to legacy names (WorldLine_Records, Operational_Data as world-line store, etc.); scripts and docs should migrate to this vocabulary and structure.

---

## 8. Migration

Migration of the current warehouse into this model is governed by:

**PX_WAREHOUSE_MIGRATION_INSTRUCTION.md**

That document defines mandatory folder consolidation, file migration, tier assignment, lifecycle enforcement, cleanup rules, script updates, and verification. The IDE must treat it as enforceable.

---

## 9. Amendment

Changes to this model require:

- Update to this file  
- Update to **PX_ENTERPRISE_LAYOUT.md** if top-level structure changes  
- Passing governance tests  
- Zero architecture drift  

This is the **locked** Enterprise Edition warehouse model.
