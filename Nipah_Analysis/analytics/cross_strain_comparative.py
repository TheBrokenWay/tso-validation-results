"""
Cross-strain comparative analytics.

Compare CFR, transmission patterns, outcomes across NiV_Malaysia vs NiV_Bangladesh.
Input: normalized data (or accepted + raw paths). Output: results/cross_strain_comparative.json.
"""

import json
from pathlib import Path
from typing import Any

from ._version import ANALYTICS_OUTPUT_VERSION


def run_cross_strain_comparative(
    normalized_dir: Path,
    results_dir: Path,
    manifest: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """
    Load normalized JSON per strain; compute comparative metrics (CFR mean/range, case counts, temporal span).
    """
    normalized_dir = Path(normalized_dir)
    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)

    by_strain: dict[str, list[dict]] = {}
    for p in sorted(normalized_dir.glob("*.json")):
        try:
            with open(p, encoding="utf-8") as f:
                data = json.load(f)
        except (json.JSONDecodeError, IOError):
            continue
        strain = data.get("strain")
        records = data.get("records", [])
        if strain not in by_strain:
            by_strain[strain] = []
        by_strain[strain].extend(records)

    comparative = {}
    for strain, records in by_strain.items():
        if not records:
            comparative[strain] = {"n": 0, "cfr_mean": None, "cfr_min": None, "cfr_max": None, "years": []}
            continue
        cfr_vals = [float(r["cfr"]) for r in records if r.get("cfr") is not None]
        years = [int(r["outbreak_year"]) for r in records if r.get("outbreak_year") is not None]
        comparative[strain] = {
            "n": len(records),
            "cfr_mean": round(sum(cfr_vals) / len(cfr_vals), 4) if cfr_vals else None,
            "cfr_min": round(min(cfr_vals), 4) if cfr_vals else None,
            "cfr_max": round(max(cfr_vals), 4) if cfr_vals else None,
            "years": sorted(set(years)) if years else [],
        }

    out = {"analytics_version": ANALYTICS_OUTPUT_VERSION, "strains": list(by_strain.keys()), "by_strain": comparative}
    out_path = results_dir / "cross_strain_comparative.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    return out

