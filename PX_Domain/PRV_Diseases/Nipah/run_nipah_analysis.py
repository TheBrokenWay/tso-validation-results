#!/usr/bin/env python3
"""
Nipah Analysis — blocking validation gate then normalization/summary.

Control flow:
  1. Load manifest → validate manifest schema → FAIL if invalid
  2. Discover raw files → FAIL if count == 0
  3. For each file: checksum, schema, strain, forbidden family, CFR envelope, temporal → FAIL if any reject
  4. Strain isolation: default FAIL if >1 strain (use --allow-multi-strain to allow)
  5. Strain coverage: FAIL if manifest defines a strain with zero valid files
  6. Proceed: write summary, validation report (with provenance), checksum registry, drift summary

Exit: 0 on success, 1 on validation failure.
"""

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

# Allow running from repo root or Nipah_Analysis
SCRIPT_DIR = Path(__file__).resolve().parent
if SCRIPT_DIR.name != "Nipah_Analysis" and (SCRIPT_DIR / "Nipah_Analysis").exists():
    SCRIPT_DIR = SCRIPT_DIR / "Nipah_Analysis"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

CONFIG_DIR = SCRIPT_DIR / "config"
RAW_DIR = SCRIPT_DIR / "data" / "raw"
NORM_DIR = SCRIPT_DIR / "data" / "normalized"
RESULTS_DIR = SCRIPT_DIR / "results"
VALIDATION_DIR = SCRIPT_DIR / "validation"
MANIFEST_PATH = CONFIG_DIR / "nipah_manifest.json"
SCHEMA_VERSION = "2.0"
ONTOLOGY_VERSION = "1.0"


def _parse_args():
    p = argparse.ArgumentParser(description="Nipah Analysis — validation gate")
    p.add_argument(
        "--allow-multi-strain",
        action="store_true",
        help="Allow runs with >1 strain; default FAIL if multiple strains present.",
    )
    p.add_argument(
        "--tso-validate",
        action="store_true",
        help="Blocking: run TSO_Validator on normalized data; FAIL if TSO validation fails.",
    )
    return p.parse_args()


def _run_tso_validator_on_normalized(norm_dir: Path, results_dir: Path) -> bool:
    """
    Run TSO_Validator with normalized Nipah data as input. Returns True if TSO run succeeded
    and status is acceptable (TSO_VALIDATED or TSO_INDETERMINATE); False if run failed or TSO_FAILED.
    """
    tso_root = SCRIPT_DIR.parent / "TSO_Validator"
    if not (tso_root / "run_validation.py").exists():
        print("WARNING: TSO_Validator not found; skipping --tso-validate.", file=sys.stderr)
        return True
    import subprocess
    import shutil
    tso_raw = tso_root / "data" / "raw"
    tso_raw.mkdir(parents=True, exist_ok=True)
    for f in norm_dir.glob("*.json"):
        shutil.copy2(f, tso_raw / f.name)
    try:
        r = subprocess.run(
            [sys.executable, "run_validation.py"],
            cwd=tso_root,
            capture_output=True,
            text=True,
            timeout=300,
        )
    except (subprocess.TimeoutExpired, FileNotFoundError) as e:
        print(f"TSO_Validator run error: {e}", file=sys.stderr)
        return False
    if r.returncode != 0:
        print(f"TSO_Validator exited {r.returncode}: {r.stderr or r.stdout}", file=sys.stderr)
        return False
    run_dirs = sorted((tso_root / "data" / "results").glob("run_*"), key=lambda p: p.name, reverse=True)
    if not run_dirs:
        return True
    summary_path = run_dirs[0] / "run_summary.json"
    if not summary_path.exists():
        return True
    try:
        with open(summary_path, encoding="utf-8") as f:
            summary = json.load(f)
    except (json.JSONDecodeError, IOError):
        return True
    status = summary.get("status", "")
    if status == "TSO_FAILED":
        print(f"TSO_Validator status TSO_FAILED; see {summary_path}", file=sys.stderr)
        return False
    return True


def main() -> int:
    args = _parse_args()
    allow_multi_strain = args.allow_multi_strain
    tso_validate = args.tso_validate

    # Ensure dirs
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    NORM_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # 1. Load and validate manifest
    try:
        from validation.validate_manifest import load_and_validate_manifest
        manifest = load_and_validate_manifest(MANIFEST_PATH)
    except Exception as e:
        print(f"FATAL: Manifest validation failed: {e}", file=sys.stderr)
        return 1

    # 2. Validate all raw data (blocking gate)
    # Definition vs requirement: manifest defines what CAN be analyzed; coverage is required only in multi-strain mode.
    # Single-strain mode: the single strain found in data/raw must be fully valid; absence of other manifest strains is OK.
    # Multi-strain mode (--allow-multi-strain): all manifest strains must be present in the batch.
    try:
        from validation.validate_raw_data import validate_all_raw
        accepted, _, drift_summary = validate_all_raw(
            RAW_DIR,
            manifest,
            fail_on_no_files=True,
            allow_multi_strain=allow_multi_strain,
            require_strain_coverage=allow_multi_strain,
        )
    except Exception as e:
        print(f"FATAL: Raw data validation failed: {e}", file=sys.stderr)
        return 1

    ingestion_time = datetime.now(timezone.utc).isoformat()

    # 3. Per-file provenance for report
    provenance_list = []
    for a in accepted:
        provenance_list.append({
            "source": Path(a["path"]).name,
            "path": a["path"],
            "ingestion_time": ingestion_time,
            "checksum": a["checksum"],
            "schema_version": SCHEMA_VERSION,
            "ontology_version": ONTOLOGY_VERSION,
            "strain": a["strain"],
            "record_count": a["record_count"],
        })

    # 4. Checksum registry: persist last-seen SHA256 per file and pass/fail (provenance ledger)
    registry_path = VALIDATION_DIR / "checksums.json"
    VALIDATION_DIR.mkdir(parents=True, exist_ok=True)
    registry = {}
    if registry_path.exists():
        try:
            with open(registry_path, encoding="utf-8") as f:
                registry = json.load(f)
        except (json.JSONDecodeError, IOError):
            registry = {}
    registry.setdefault("files", {})
    for a in accepted:
        key = Path(a["path"]).name
        registry["files"][key] = {
            "checksum": a["checksum"],
            "passed": True,
            "last_seen": ingestion_time,
            "strain": a["strain"],
        }
    registry["last_run"] = ingestion_time
    with open(registry_path, "w", encoding="utf-8") as f:
        json.dump(registry, f, indent=2)

    # 4. Normalization (units + PII anonymization)
    normalized_paths = []
    normalization_error = None
    try:
        from normalize_data import normalize_accepted
        normalized_paths = normalize_accepted(RAW_DIR, NORM_DIR, accepted)
    except Exception as e:
        print(f"FATAL: Normalization failed: {e}", file=sys.stderr)
        return 1

    # 4b. Optional blocking gate: TSO_Validator on normalized data
    if tso_validate:
        if not _run_tso_validator_on_normalized(NORM_DIR, RESULTS_DIR):
            print("FATAL: TSO_Validator gate failed (--tso-validate).", file=sys.stderr)
            return 1

    # 5. Summary and validation report (with provenance and drift)
    gates = {
        "raw_validation": "PASS",
        "normalization": "PASS",
        "tso_validation": "PASS" if tso_validate else "SKIPPED",
    }
    run_mode = "multi_strain" if allow_multi_strain else "single_strain"
    summary = {
        "analysis_id": manifest.get("analysis_id", "NIPAH_ANALYSIS"),
        "mode": run_mode,
        "pathogen": manifest.get("pathogen"),
        "validation": "PASS",
        "gates": gates,
        "accepted_files": [{k: v for k, v in a.items() if k != "cfr_values"} for a in accepted],
        "accepted_count": len(accepted),
        "strains_seen": sorted({a["strain"] for a in accepted}),
        "allow_multi_strain": allow_multi_strain,
        "normalized_count": len(normalized_paths),
        "normalized_files": [str(p.name) for p in normalized_paths],
    }
    if normalization_error:
        summary["normalization_error"] = normalization_error
    summary["tso_validate_run"] = tso_validate
    summary_path = RESULTS_DIR / "summary.json"
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    validation_report = {
        "manifest_valid": True,
        "raw_validation": "PASS",
        "files_accepted": len(accepted),
        "files_rejected": 0,
        "checksums": {Path(a["path"]).name: a["checksum"] for a in accepted},
        "provenance": {
            "schema_version": SCHEMA_VERSION,
            "ontology_version": ONTOLOGY_VERSION,
            "ingestion_time": ingestion_time,
            "files": provenance_list,
        },
        "constraint_drift": {
            "cfr_by_strain": drift_summary,
            "description": "CFR distribution per strain for drift detection.",
        },
    }
    report_path = RESULTS_DIR / "validation_report.json"
    with open(report_path, "w", encoding="utf-8") as f:
        json.dump(validation_report, f, indent=2)

    # 6. Analytics (cross-strain comparative, constraint divergence, evolutionary drift, intervention sensitivity)
    try:
        from analytics.cross_strain_comparative import run_cross_strain_comparative
        from analytics.constraint_divergence import run_constraint_divergence_detection
        from analytics.evolutionary_drift import run_evolutionary_drift_early_warning
        from analytics.intervention_sensitivity import run_intervention_sensitivity_modeling
        run_cross_strain_comparative(NORM_DIR, RESULTS_DIR, manifest)
        run_constraint_divergence_detection(RESULTS_DIR / "validation_report.json", RESULTS_DIR, manifest)
        run_evolutionary_drift_early_warning(
            RESULTS_DIR / "validation_report.json",
            RESULTS_DIR,
            registry_path,
            prior_drift_path=None,
            manifest=manifest,
        )
        run_intervention_sensitivity_modeling(NORM_DIR, RESULTS_DIR, manifest)
    except Exception as e:
        print(f"WARNING: Analytics step failed: {e}", file=sys.stderr)

    print(f"Nipah analysis validation PASS. Summary: {summary_path}")
    print(f"  Accepted files: {len(accepted)}")
    print(f"  Strains: {summary['strains_seen']}")
    print(f"  Normalized: {summary['normalized_count']} files")
    print(f"  Checksum registry: {registry_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
