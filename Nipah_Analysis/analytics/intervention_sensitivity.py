"""
Intervention sensitivity modeling.

Model impact of interventions (e.g. vaccination coverage, contact reduction) by strain.
Input: normalized data + manifest (CFR, human_to_human). Output: results/intervention_sensitivity.json.
"""

import json
from pathlib import Path
from typing import Any

from ._version import ANALYTICS_OUTPUT_VERSION


def run_intervention_sensitivity_modeling(
    normalized_dir: Path,
    results_dir: Path,
    manifest: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """
    Placeholder: strain-specific baseline CFR and human_to_human flag;
    sensitivity outputs (e.g. relative reduction in effective R per strain) can be computed here.
    """
    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)

    pathogen = (manifest or {}).get("pathogen", {})
    strains_def = pathogen.get("strains", {})

    by_strain = {}
    for strain_id, defn in strains_def.items():
        cfr_range = defn.get("cfr_range", [0, 1])
        h2h = defn.get("human_to_human", False)
        by_strain[strain_id] = {
            "cfr_baseline_range": cfr_range,
            "human_to_human": h2h,
            "intervention_sensitivity_note": "Vaccination/contact-reduction impact model: plug-in strain-specific R0 and efficacy here.",
        }

    out = {"analytics_version": ANALYTICS_OUTPUT_VERSION, "by_strain": by_strain, "status": "computed"}
    out_path = results_dir / "intervention_sensitivity.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    return out

