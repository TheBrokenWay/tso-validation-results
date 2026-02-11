# Nipah Analysis Runbook

## Prerequisites

- Repo root (e.g. `e:\foundation`) or Nipah_Analysis as cwd.
- Python 3.x (no extra deps for validation gate).

## 1. Prepare data

- Place **strain-annotated** Nipah raw files in `Nipah_Analysis/data/raw/` (top-level only).
- Each JSON must include:
  - `strain`: `"NiV_Malaysia"` or `"NiV_Bangladesh"` (no generic NiV).
  - `records`: array of objects with `outbreak_year` (integer) and `cfr` (number in [0,1]).
- CFR must lie within strain envelope (see `config/nipah_ontology.json`).
- Outbreak year must be ≥ first outbreak year for that strain (Malaysia 1998, Bangladesh 2001).

## 2. Config

- **Manifest:** `config/nipah_manifest.json` — strain-split (NiV_Malaysia, NiV_Bangladesh), CFR ranges, vectors, human-to-human, spillover pattern.
- **Schema:** `config/nipah_schema.json` — required fields, strain enum, record shape.
- **Ontology:** `config/nipah_ontology.json` — disease-specific CFR bounds and temporal bounds by strain.

## 3. Run analysis (validation as hard gate)

From Nipah_Analysis:

```powershell
cd e:\foundation\Nipah_Analysis
.\run_nipah_analysis.ps1
```

Or run Python directly:

```powershell
python run_nipah_analysis.py
```

Optional second gate (TSO_Validator on normalized data):

```powershell
python run_nipah_analysis.py --tso-validate
```

- **Definition vs requirement:** The manifest **defines** what strains *can* be analyzed; it does **not** require every strain in every batch.
- **Strain isolation (default):** If raw data contains **more than one strain**, the run **FAILs** unless you pass `--allow-multi-strain`. Single-strain analysis (e.g. only NiV-M files) **passes** without needing Bangladesh data; coverage is **not** required for the other manifest strain.
- **Full coverage (multi-strain only):** When you pass `--allow-multi-strain`, **all** manifest strains must have at least one valid file in the batch (no “forgot to include Bangladesh data”).
- **If validation fails:** script exits with code 1. No partial runs.
- **If validation passes:** Normalization (units + PII anonymization) runs, then analytics; outputs: `results/summary.json`, `results/validation_report.json`, `validation/checksums.json`, `data/normalized/*.json`, and analytics (cross_strain_comparative, constraint_divergence, evolutionary_drift, intervention_sensitivity).

## 4. Pipeline order (high-integrity)

1. **RAW DATA** → **Validation gate (blocking)** — manifest schema + raw file: checksum, file-type, schema, ontology (CFR envelope, temporal), strain classification, forbidden family. If **any** file fails or count == 0 → **STOP** (exit 1). No partial runs.
2. **Only if PASS** → **Normalization** (`normalize_data.py`: units → copies/mL, PII anonymized). If normalization fails → **STOP** (exit 1).
3. **Optional:** `--tso-validate` — run TSO_Validator on normalized data as a **second gate**; if TSO status is TSO_FAILED → **STOP** (exit 1).
4. **Only if PASS** → **Analytics** (cross-strain comparative, constraint divergence, evolutionary drift, intervention sensitivity).
5. **Only if PASS** → PX_System analysis (constraint fitting, DCMs).

**Schema and ontology validation is the primary gatekeeper.** Never run normalization or PX on unvalidated raw data. TSO_Validator is an optional second gate on normalized output when `--tso-validate` is set.

## 5. Edge case testing

- **Location:** `data/raw/edge_cases/` (not discovered by default; discovery is top-level `data/raw/*.json` only).
- **Use:** Copy an edge-case file to `data/raw/` to test failure behavior:
  - `edge_mixed_strain.json` — no strain; expect **FAIL** (Unclassifiable / schema).
  - `edge_schema_violation.json` — wrong field name (CFR vs cfr); expect **FAIL** (Schema violation).
  - `edge_cfr_anomaly.json` — Malaysia CFR 18%; expect **FAIL** (CFR outside envelope).
  - `edge_temporal_drift.json` — Bangladesh outbreak_year 1998; expect **FAIL** (Temporal inconsistency).

## 6. PII hashing (normalization)

- **NIPAH_SALT** (optional): Cryptographic salt for PII anonymization. If set, hashing is deterministic (sha256(salt + value)). If unset and **NIPAH_SALT_STRICT** is not true, a random ephemeral salt is used (non-reproducible).
- **NIPAH_SALT_STRICT** (optional): If `true` and NIPAH_SALT is not set, normalization fails (exit 1). Use for production when reproducible anonymization is required.

## 7. Outputs

- **results/summary.json** — analysis_id, pathogen, validation PASS, **gates** (raw_validation, normalization, tso_validation: PASS | SKIPPED | FAIL), accepted_files, strains_seen, normalized_count, normalized_files.
- **results/validation_report.json** — manifest_valid, raw_validation PASS, checksums, **provenance**, **constraint_drift**.
- **validation/checksums.json** — checksum registry (last-seen SHA256, passed, last_seen per file).
- **data/normalized/** — one JSON per raw file: units standardized (e.g. viral_load → copies_per_mL), PII anonymized.
- **results/cross_strain_comparative.json** — CFR mean/min/max and years per strain (when multi-strain or comparative).
- **results/constraint_divergence.json** — flags when observed CFR drifts outside manifest bounds.
- **results/evolutionary_drift.json** — early-warning: current vs prior CFR by strain; threshold **max_mean_cfr_drift** from manifest `analytics_config` (default 0.05).
- **results/intervention_sensitivity.json** — strain-specific baseline and intervention-sensitivity placeholders.

**Manifest analytics_config:** `config/nipah_manifest.json` includes `analytics_config.max_mean_cfr_drift` (e.g. 0.05) for evolutionary drift; override there instead of hardcoding.
