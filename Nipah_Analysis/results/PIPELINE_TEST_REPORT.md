# Full Pipeline Test Report

**Date:** 2026-01-30  
**Scope:** Nipah_Analysis, TSO_Validator, PX_System (Nipah DCM wiring).

---

## 1. Nipah_Analysis Pipeline

| Step | Result | Notes |
|------|--------|--------|
| Manifest load + validate | PASS | Strain-split manifest valid |
| Raw validation (multi-strain) | PASS | 2 files accepted (NiV_Malaysia, NiV_Bangladesh) |
| Raw validation (single-strain) | PASS | 1 file (NiV_Malaysia only); mode=single_strain |
| Normalization | PASS | 2 → 2 normalized (units + PII anonymized) |
| Checksum registry | PASS | validation/checksums.json updated |
| Summary + mode | PASS | summary.json has "mode": "single_strain" | "multi_strain" |
| Analytics (cross_strain, divergence, drift, intervention) | PASS | All 4 outputs with analytics_version "1.0" |

**Tests:** `pytest tests/` — **9 passed** (manifest, edge cases, forbidden family, strain coverage, single-strain no-coverage).

---

## 2. TSO_Validator Pipeline

| Step | Result | Notes |
|------|--------|--------|
| run_validation.py | PASS | Exit 0; run_YYYYMMDD_HHMMSS created |
| run_summary.json, provenance.json, historical_context.json | PASS | Present in run dir |
| Deprecation fix | DONE | datetime.utcnow() → datetime.now(timezone.utc) |

**Tests:** `pytest tests/` — **17 passed** (ethical_guard, history, provenance, run_summary_schema).  
**Fix applied:** test_run_summary_schema.py had invalid syntax (`when`); corrected to `assert not run_summary["blind_mode"] or run_summary["blind_mode_note"] is not None`.

---

## 3. Nipah → TSO Gate (--tso-validate)

| Step | Result | Notes |
|------|--------|--------|
| Nipah run with --tso-validate | FAIL (expected) | TSO runs on Nipah normalized data; schema mismatch → TSO_FAILED; gate correctly fails |

Nipah normalized output (strain, records, cfr, outbreak_year) is not TSO-shaped (delta_energy, delta_entropy, control_parameter, time). Use --tso-validate only when normalized data is TSO-compatible.

---

## 4. PX_System Nipah DCM Wiring

| Check | Result |
|-------|--------|
| get_nipah_dcm("NiV_Malaysia") | OK (18 keys) |
| get_nipah_dcm("NiV_Bangladesh") | OK (18 keys) |
| get_disease_constraints("NiV_Malaysia") | OK (returns Malaysia DCM) |

---

## 5. Resolved Issues

1. **TSO_Validator:** Replaced `datetime.utcnow()` with `datetime.now(timezone.utc)` to remove deprecation warnings.
2. **TSO_Validator tests:** Fixed syntax in test_run_summary_schema.py (invalid `when` in assert).
3. **Nipah analytics:** All analytics JSON outputs include `analytics_version` (from analytics/_version.py).
4. **Nipah summary:** Explicit `"mode": "single_strain" | "multi_strain"` in summary.json.

---

## 6. Layer Verification Summary

- **Nipah:** Validation → Normalization → Analytics — all layers working; mode and analytics_version present.
- **TSO_Validator:** Config → Ingestion → Normalization → Physics → Audit → Reports → run_summary/provenance/history — working; deprecations fixed.
- **PX_System:** Nipah DCMs (create_nipah_malaysia_dcm, create_nipah_bangladesh_dcm, get_nipah_dcm, get_disease_constraints) — working.
