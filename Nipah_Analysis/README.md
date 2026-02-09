# Nipah Analysis

Dedicated workspace for **Nipah virus** (NiV) analysis — **strain-aware** and **validation-gated**.

## Scope

- **Pathogen:** Nipah virus (Henipavirus, Paramyxoviridae) — **NiV_Malaysia** and **NiV_Bangladesh** (no monolithic NiV).
- **Use:** Outbreak/epidemiology data, antiviral discovery constraints, strain-specific CFR and temporal bounds.
- **Validation:** Blocking. Raw data must pass checksum, schema, strain classification, CFR envelope, and temporal consistency before normalization or PX analysis.

## Layout

```
Nipah_Analysis/
├── README.md
├── RUNBOOK.md
├── config/
│   ├── nipah_manifest.json   # Strain-split (NiV_Malaysia, NiV_Bangladesh)
│   ├── nipah_schema.json     # Required fields, strain enum, record shape
│   └── nipah_ontology.json   # Disease-specific bounds (CFR, first outbreak year)
├── validation/               # Blocking gate (Python)
│   ├── validate_manifest.py
│   └── validate_raw_data.py  # Checksum, schema, strain, CFR, temporal
├── data/
│   ├── raw/                  # Top-level JSON/CSV only (discovered for validation)
│   │   ├── valid_*.json      # Valid strain-annotated data
│   │   └── edge_cases/      # Intentional failure cases for testing
│   └── normalized/
├── results/                  # summary.json, validation_report.json
├── run_nipah_analysis.py     # Entry: manifest → raw validation → summary (exit 1 on fail)
└── run_nipah_analysis.ps1    # Calls Python; fails if validation fails
```

## Control flow

1. **Load manifest** → validate manifest schema (strain-split) → **FAIL** if invalid.
2. **Discover raw files** → **FAIL** if count == 0.
3. **For each file:** checksum, schema, forbidden family, strain, CFR envelope, temporal → **FAIL** if any reject.
4. **Strain isolation (default):** >1 strain present → **FAIL** unless `--allow-multi-strain`. **Coverage:** required only when `--allow-multi-strain` (all manifest strains must have data). **Single-strain mode:** the one strain in data/raw must be fully valid; other manifest strains need not be present (definition ≠ requirement).
5. **Normalization** — units (e.g. viral_load → copies/mL), PII anonymized → `data/normalized/`. If normalization fails → **STOP** (exit 1).
6. **Optional:** `--tso-validate` — run TSO_Validator on normalized data as a second gate; **STOP** if TSO status is TSO_FAILED.
7. **Analytics** — cross-strain comparative, constraint divergence, evolutionary drift, intervention sensitivity → `results/*.json`.

No partial runs. CFR and first-outbreak-year bounds are derived from the manifest (single source of truth).

## Data

- **Raw:** Place **strain-annotated** JSON in `data/raw/` (top-level only). Each file must have `strain` (`NiV_Malaysia` or `NiV_Bangladesh`) and `records` with `outbreak_year` and `cfr` per record. CFR must fall within strain envelope (Malaysia 0.35–0.45, Bangladesh 0.70–1.00). Outbreak year must be ≥ first outbreak year for that strain (Malaysia 1998, Bangladesh 2001).
- **Edge cases:** `data/raw/edge_cases/` contains intentional failure cases (mixed strain, schema violation, CFR anomaly, temporal drift) for testing; they are not discovered by default (discovery is top-level only).

## Repo context

- **TSO_Validator** — Schema/ontology validation is the **gatekeeper**: run TSO (or this validation gate) first; only on PASS proceed to normalization and PX.
- **PX_System/foundation/Disease_Constraint_Model.py** — Nipah DCMs are wired: `get_disease_constraints("NiV_Malaysia")` or `("NiV_Bangladesh")`, or `create_nipah_malaysia_dcm()`, `create_nipah_bangladesh_dcm()`, `get_nipah_dcm(strain)`. Use strain-specific DCM for constraint fitting and intervention modeling.
