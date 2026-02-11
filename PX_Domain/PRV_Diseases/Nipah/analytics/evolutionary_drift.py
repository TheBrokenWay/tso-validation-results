"""
Evolutionary drift early-warning.

Temporal and distribution shifts: compare current run CFR/temporal to prior run (checksum registry or last results).
Input: validation/checksums.json + results/validation_report.json; optional prior run. Output: results/evolutionary_drift.json.
"""

import json
from pathlib import Path
from typing import Any

from ._version import ANALYTICS_OUTPUT_VERSION


DEFAULT_MAX_MEAN_CFR_DRIFT = 0.05


def run_evolutionary_drift_early_warning(
    current_validation_report_path: Path,
    results_dir: Path,
    checksum_registry_path: Path | None = None,
    prior_drift_path: Path | None = None,
    manifest: dict | None = None,
) -> dict[str, Any]:
    """
    Compare current cfr_by_strain to prior (if available); flag significant mean shift or range expansion.
    Threshold max_mean_cfr_drift from manifest.analytics_config (default 0.05).
    """
    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)
    max_drift = DEFAULT_MAX_MEAN_CFR_DRIFT
    if manifest and isinstance(manifest.get("analytics_config"), dict):
        mc = manifest["analytics_config"].get("max_mean_cfr_drift")
        if isinstance(mc, (int, float)) and 0 < mc <= 1:
            max_drift = float(mc)

    out = {"analytics_version": ANALYTICS_OUTPUT_VERSION, "warnings": [], "current_cfr_by_strain": {}, "prior_cfr_by_strain": None, "status": "ok", "max_mean_cfr_drift": max_drift}

    if not Path(current_validation_report_path).exists():
        out["status"] = "no_report"
        _write(out, results_dir)
        return out

    with open(current_validation_report_path, encoding="utf-8") as f:
        report = json.load(f)
    current = report.get("constraint_drift", {}).get("cfr_by_strain", {})
    out["current_cfr_by_strain"] = current

    prior = None
    if prior_drift_path and Path(prior_drift_path).exists():
        try:
            with open(prior_drift_path, encoding="utf-8") as f:
                prior_data = json.load(f)
            prior = prior_data.get("current_cfr_by_strain") or prior_data.get("cfr_by_strain")
        except (json.JSONDecodeError, IOError) as e:
            import sys
            print(f"    WARN: prior drift data unreadable, skipping comparison: {e}", file=sys.stderr)
    if prior:
        out["prior_cfr_by_strain"] = prior
        for strain in set(current) | set(prior):
            c_cur = current.get(strain, {}).get("mean")
            c_pri = prior.get(strain, {}).get("mean")
            if c_cur is not None and c_pri is not None:
                delta = abs(c_cur - c_pri)
                if delta > max_drift:
                    out["warnings"].append({
                        "strain": strain,
                        "message": f"CFR mean shift: {c_pri:.4f} → {c_cur:.4f} (delta={delta:.4f})",
                        "severity": "high" if delta > 0.10 else "medium",
                    })
        if out["warnings"]:
            out["status"] = "drift_detected"

    _write(out, results_dir)
    return out


def _write(out: dict[str, Any], results_dir: Path) -> None:
    with open(results_dir / "evolutionary_drift.json", "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)

