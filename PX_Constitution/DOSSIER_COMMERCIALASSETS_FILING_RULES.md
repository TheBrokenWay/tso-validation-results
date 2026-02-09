# DIRECTIVE: ENFORCE_FINAL_PRODUCT_FILING

**Dossier + CommercialAssets Filing Rules**  
**Deterministic • Enforceable • Drift-proof**  
**Authority:** PX_Constitution — IDE and all writers must comply.

---

## 1. Dossier_Final Admission Rule

- **Only** assets that have **PASSED ALL CHECKS** may enter  
  `PX_Warehouse/CommercialAssets/Dossier_Final`.
- **Required checks** include:
  - Scientific filters (Lipinski, Veber, Ghose, PAINS, structural alerts)
  - ADMET and exposure modeling
  - Legal/FTO clearance
  - Commercial viability
  - Provenance, timestamp, reproducibility
- Any asset that **fails ANY check** is automatically **excluded**.
- **No** partial, draft, or intermediate files may enter Dossier_Final.

---

## 2. CommercialAssets Biopath Routing

- All finished commercial-grade outputs must be routed into the **correct tier**:
  - **Diamond** → highest-value, fully validated, multi-pathway assets
  - **Gold** → validated, high-confidence assets
  - **Silver** → validated, moderate-confidence assets
  - **Bronze** → validated, low-confidence assets
  - **Executive_Summary** → human-readable summaries for each asset
  - **Audit_Trail** → provenance, logs, fingerprints, and validation traces
- Routing is determined by the **final commercial score** and **governance metadata**.
- **No** asset may bypass tier assignment.

---

## 3. Enforcement

- The IDE (and any script writing into the warehouse) must **reject** any write operation that attempts to place:
  - an **unvalidated** asset into Dossier_Final
  - a **non-final** asset into any CommercialAssets tier
  - any file into a tier that **does not match** its assigned classification
- Violations trigger **FAIL-CLOSED** behavior and **must** be logged.

---

## 4. Canonical Filing Structure

- All final products follow this path:

```
PX_Warehouse/
  CommercialAssets/
    <Tier>/
      <AssetID>/
        dossier.json
        evidence/
        summary.md
        audit/
```

- **No** deviations or alternate paths permitted.
- `<Tier>` is one of: **Diamond**, **Gold**, **Silver**, **Bronze**.
- **Executive_Summary** and **Audit_Trail** are sibling folders under CommercialAssets; content is placed per asset or aggregated as defined by the pipeline.

---

## 5. Reference Implementation

- **PX_Constitution/Filing_Rules.py** provides deterministic, callable enforcement:
  - `may_enter_dossier_final(asset_metadata)` → bool
  - `get_tier_for_asset(asset_metadata)` → str
  - `validate_filing_path(path, asset_metadata)` → (ok: bool, reason: str)
- All writers (scripts, IDE, pipelines) should call these before writing into CommercialAssets or Dossier_Final.
