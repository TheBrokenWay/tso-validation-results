# WAREHOUSE DRIFT REPORT — IDE EXECUTION LOGIC

**Authoritative Verification Algorithm for PX_Warehouse**

The IDE must run the following steps in order, producing a **Warehouse Drift Report** if any mismatch is detected.

---

## 1. Load Canonical Model

The IDE loads the canonical structure from:

- **PX_WAREHOUSE_ENTERPRISE_MODEL.md**
- **PX_WAREHOUSE_MIGRATION_INSTRUCTION.md**
- **PX_ENTERPRISE_LAYOUT.md**

It extracts:

- Allowed top-level directories  
- Allowed subdirectories  
- Allowed vocabulary  
- Lifecycle rules  
- Tiering rules  
- Archive rules  
- Cleanup rules  

This becomes the **canonical warehouse tree**.

---

## 2. Scan Actual Warehouse

The IDE recursively scans:

```
PX_Warehouse/
```

Collecting:

- All folders  
- All files  
- All metadata (timestamps, names, dossier presence, validation markers)  
- All world-line records  
- All trial simulation outputs  
- All commercial assets  
- All logs  
- All snapshots  

This produces the **actual warehouse tree**.

---

## 3. Compare Canonical vs Actual

The IDE performs a node-by-node comparison.

### 3.1 Structural Drift

Any folder **not** in the canonical tree is flagged.

### 3.2 Naming Drift

Any folder or file whose name does **not** match the canonical vocabulary is flagged.

### 3.3 Lifecycle Drift

The IDE checks each asset:

- If it has **no** validation → must be in **CommercialAssets**
- If it has **validation metadata** → must be in **WorldLines**
- If it has a **dossier** → must be in **Archive**

Any mismatch is flagged.

### 3.4 Tiering Drift

The IDE checks: folder tier vs metadata tier; missing tier metadata; unknown tier names. Any mismatch is flagged.

### 3.5 Legacy Artifacts

The IDE identifies: v3 analytics; old dashboards; old benchmark reports; old monetization/stratification; deduplication logs; any file older than the enterprise edition cutoff. These are flagged for archival.

### 3.6 Unclassified Assets

Any file that **cannot** be assigned to CommercialAssets, WorldLines, TrialSimulations, Operations, or Archive is flagged.

### 3.7 Unexpected Files

Any file **not** referenced by the enterprise model is flagged.

---

## 4. Generate Required Actions

For each drift item, the IDE generates a deterministic action:

- **MOVE**
- **RENAME**
- **ARCHIVE**
- **DELETE**
- **TIER**
- **UPDATE_SCRIPT**

Each action is listed with: source path, destination path, reason, rule violated.

---

## 5. Re-Verification

After applying all actions, the IDE:

1. Re-scans the warehouse  
2. Re-compares canonical vs actual  
3. Confirms zero drift  

If drift remains, a **new report** is generated.

---

## 6. Certification Block

If the warehouse matches the canonical model:

- **Warehouse Verified:** YES  
- **Drift Eliminated:** YES  
- **Canonical Structure Achieved:** YES  

If not:

- **Warehouse Verified:** NO  
- **Drift Eliminated:** NO  
- **Canonical Structure Achieved:** NO  

---

**Report output:** Use **WAREHOUSE_DRIFT_REPORT_TEMPLATE.md** to produce  
`PX_Warehouse/WAREHOUSE_DRIFT_REPORT_<timestamp>.md` (or `PX_LOGS/warehouse_drift_<timestamp>.md`).
