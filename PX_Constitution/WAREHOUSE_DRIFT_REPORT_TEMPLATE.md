# WAREHOUSE DRIFT REPORT

**Template for IDE-Generated Warehouse Drift Analysis**  
**Enterprise Edition — Governance Enforcement Artifact**

This report is generated when the IDE performs a warehouse verification pass and detects any deviation from the canonical structure defined in:

- **PX_WAREHOUSE_ENTERPRISE_MODEL.md**
- **PX_WAREHOUSE_MIGRATION_INSTRUCTION.md**
- **PX_ENTERPRISE_LAYOUT.md**

The purpose of this report is to provide a complete, auditable record of all structural, naming, lifecycle, and classification drift within `PX_Warehouse/`.

---

## 1. Report Metadata

| Field | Value |
|-------|--------|
| Report ID | *auto-generated* |
| Timestamp (UTC) | |
| IDE Version | |
| Warehouse Path | PX_Warehouse/ |
| Migration Directive Version | |
| Enterprise Model Version | |
| Drift Status | PASS / FAIL |

---

## 2. Summary of Drift

A high-level summary of all drift categories.

```
Total Issues Detected: <N>
Structural Drift: <N>
Naming Drift: <N>
Lifecycle Drift: <N>
Tiering Drift: <N>
Legacy Artifacts: <N>
Unclassified Assets: <N>
Unexpected Folders: <N>
Unexpected Files: <N>
```

---

## 3. Structural Drift (Folder-Level)

List any folder that exists but is not part of the canonical structure defined in §1 of the migration directive.

```
- <folder_path> — Not part of canonical warehouse structure
- <folder_path> — Duplicate of <canonical_path>
- <folder_path> — Legacy folder (should be archived)
```

---

## 4. Naming Drift

List any folder or file whose name does not match the canonical vocabulary: **Diamond**, **Gold**, **Silver**, **Bronze**, **CommercialAssets**, **WorldLines**, **TrialSimulations**, **Operations**, **Archive**.

```
- <path> — Non-canonical name; expected <canonical_name>
- <path> — Mixed terminology (e.g., WorldLine_Records)
```

---

## 5. Lifecycle Drift

List any asset that is in the wrong lifecycle location:

- Created assets not in CommercialAssets
- Validated assets not in WorldLines
- Dossier-created assets not in Archive

```
- <path> — Validated asset found outside WorldLines
- <path> — Dossier found outside Archive
- <path> — Created asset found outside CommercialAssets
```

---

## 6. Tiering Drift

List any asset that: has no tier; has conflicting tier metadata; is placed in the wrong tier folder.

```
- <path> — Missing tier metadata; defaulted to Bronze
- <path> — Tier mismatch: metadata=Gold, folder=Silver
- <path> — Unrecognized tier
```

---

## 7. Legacy Artifacts

List all v2/v3-era or deprecated files that must be archived.

```
- <path> — v3 analytics (monetization, stratification, deduplication)
- <path> — legacy benchmark reports
- <path> — old dashboards
- <path> — obsolete logs
```

---

## 8. Unclassified Assets

List any file the IDE could not classify into: CommercialAssets, WorldLines, TrialSimulations, Operations, Archive.

```
- <path> — Unclassified; requires manual review
```

---

## 9. Unexpected Files

List any file not referenced by the enterprise model.

```
- <path> — Unexpected file; not part of canonical warehouse
```

---

## 10. Required Actions

A deterministic list of actions the IDE must take to resolve drift.

```
MOVE: <path> → <canonical_path>
RENAME: <path> → <canonical_name>
ARCHIVE: <path> → PX_Warehouse/Archive/<timestamp>/
DELETE: <path>
TIER: <path> → <Diamond/Gold/Silver/Bronze>
UPDATE_SCRIPT: <script_path> — Update references to canonical structure
```

---

## 11. Post-Migration Verification

After applying required actions, the IDE must: Re-scan the warehouse; Re-generate this report; Confirm zero drift.

```
Final Drift Status: PASS / FAIL
```

---

## 12. Certification Block

```
Warehouse Verified: YES / NO
Drift Eliminated: YES / NO
Canonical Structure Achieved: YES / NO
Signed By: <IDE / Governance Engine>
Timestamp: <timestamp>
```
