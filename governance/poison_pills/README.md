# Governance Poison Pills — CI-Wide Registry

Convention for IDE and CI to scan: **governance as a testable surface, not a philosophy.**

## Layout

```
governance/poison_pills/
├── registry.json          # Master list: payload, sidecar, test_runner, expected_failure_class
├── README.md              # This file
├── nipah_cfr_envelope/    # Active pill: CFR envelope violation
│   └── manifest.json
├── temporal_paradox/       # Reserved
├── duplicate_domain_truth/ # Reserved
└── shadow_validation/      # Reserved
```

## Per-pill contract

Each pill (active or reserved) has:

| Field | Description |
|-------|-------------|
| **payload** | Path to the poison-pill data file (e.g. raw JSON that must be rejected). |
| **sidecar** | Path to classification/sidecar (e.g. `.px_classification.json`). |
| **test_runner** | Script that runs the gate (e.g. `run_poison_pill_test.py`). |
| **expected_failure_class** | e.g. `ONTOLOGY_VIOLATION`, `TEMPORAL_VIOLATION`. |
| **expected_failure_layer** | Layer that must reject (e.g. `PX_Validation`). Locks ontology to the correct layer. |

## IDE / CI usage

1. **Scan** `registry.json` (or each `*/manifest.json`) for `payload`, `sidecar`, `test_runner`, `expected_failure_class`.
2. **Run** the `test_runner` from repo root; it must exit 0 and the adapter must set `last_failure_layer` to the expected layer.
3. **Assert** layer pinning: failure must occur at the documented layer (stops Engine/Audit from becoming de facto validators).

## Active pills

- **nipah_cfr_envelope**: 0% CFR for NiV_Malaysia → rejected at PX_Validation (CFR envelope). Test: `python run_poison_pill_test.py`.
- **temporal_paradox**: Outbreak year 1998 for NiV_Bangladesh (first outbreak 2001) → rejected at PX_Validation (TEMPORAL_DRIFT). Test: `python run_temporal_paradox_test.py`.

## Reserved (placeholders)

- **duplicate_domain_truth**: Conflicting domain truth (e.g. dual strain / duplicate source).
- **shadow_validation**: Validation logic duplicated outside PX_Validation (architectural violation).
- **energy_paradox**: Thermodynamic/energy paradox (e.g. non-zero energy delta). Expected at PX_Engine.
- **shadow_ontology**: Ontology logic duplicated outside PX_Validation (shadow ontology).
- **dual_truth_collision**: Conflicting domain truth (same entity, incompatible strain/source).

Add payload, sidecar, and test_runner when implementing.

## Governance coverage

`governance/governance_coverage.json` is auto-generated from the registry (on successful poison-pill gate run). Purely informational, but powerful for audits: `total_pills`, `active`, `reserved`, `coverage` (active/total).

## Signed gate stamp (future)

The gate stamp writes `hash(GOVERNANCE_SIGNATURE)` and, when configured, **HMAC-SHA256(payload, secret)**. Set `PX_GATE_STAMP_SECRET` (hex or raw) for production; or `PX_GATE_STAMP_DEV_ALLOW=1` to use the dev secret so the stamp file includes `hmac_sha256=<hex>`. Verification: `PX_Constitution.gate_stamp.verify_gate_stamp(payload_bytes, expected_hmac_hex)`. Not required for internal audit; use for hostile-modification resistance.
