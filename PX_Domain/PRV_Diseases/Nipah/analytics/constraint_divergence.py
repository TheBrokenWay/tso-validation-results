"""
Constraint divergence detection.

Detect when observed CFR (or other constraints) drift outside manifest/strain bounds.
Input: validation_report constraint_drift + manifest. Output: results/constraint_divergence.json.
"""

import json
from pathlib import Path
from typing import Any

from ._version import ANALYTICS_OUTPUT_VERSION

CFR_BOUNDS = {"NiV_Malaysia": (0.35, 0.45), "NiV_Bangladesh": (0.70, 1.00)}


def run_constraint_divergence_detection(
    validation_report_path: Path,
    results_dir: Path,
    manifest: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """
    Compare cfr_by_strain from validation report to manifest CFR bounds; flag divergence.
    """
    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)

    if not Path(validation_report_path).exists():
        return {"divergences": [], "status": "no_report"}

    with open(validation_report_path, encoding="utf-8") as f:
        report = json.load(f)
    drift = report.get("constraint_drift", {}).get("cfr_by_strain", {})

    divs = []
    for strain, stats in drift.items():
        if strain not in CFR_BOUNDS:
            continue
        lo, hi = CFR_BOUNDS[strain]
        mean = stats.get("mean")
        mn = stats.get("min")
        mx = stats.get("max")
        if mean is not None and (mean < lo or mean > hi):
            divs.append({"strain": strain, "metric": "cfr_mean", "value": mean, "bounds": [lo, hi], "severity": "high"})
        if mn is not None and mn < lo:
            divs.append({"strain": strain, "metric": "cfr_min", "value": mn, "bounds": [lo, hi], "severity": "medium"})
        if mx is not None and mx > hi:
            divs.append({"strain": strain, "metric": "cfr_max", "value": mx, "bounds": [lo, hi], "severity": "medium"})

    out = {"analytics_version": ANALYTICS_OUTPUT_VERSION, "divergences": divs, "status": "divergent" if divs else "within_bounds"}
    out_path = results_dir / "constraint_divergence.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    return out

