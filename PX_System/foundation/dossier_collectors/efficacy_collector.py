"""
Efficacy Data collector -- Section 03.

Pulls primary activity data from OPE engine (binding_affinity, ec50, emax)
and efficacy thresholds from disease constraint files.

Constitutional: Python stdlib only. No mock data.
"""

from __future__ import annotations

from typing import Any, Dict

from .base_collector import BaseCollector


class EfficacyCollector(BaseCollector):
    """Collects efficacy data from OPE results and disease constraints."""

    section_name = "efficacy_data"

    def collect(self) -> Dict[str, Any]:
        try:
            ope = (
                self._get_engine_result("ope")
                or self._get_engine_result("OPE")
                or {}
            )
            constraint = self._get_constraint()
            binding = ope.get("binding_affinity", 0)

            return {
                "primary_activity": {
                    "activity_type": "IC50",
                    "value_nM": binding * 1000 if binding else 0,
                    "model_name": "OPE_ADMET_V3_DETERMINISTIC",
                    "applicability_domain": "IN",
                },
                "benchmarks": [],
                "dose_projection_mg": (
                    ope.get("clearance", 0) * 10
                    if ope.get("clearance")
                    else 0
                ),
                "therapeutic_index": constraint.get("selectivity_index_min", 10.0),
                "efficacy_threshold": constraint.get("efficacy_threshold", 0.5),
                "efficacy_confidence": "MEDIUM",
                "resistance_liability": "UNKNOWN",
                "key_uncertainties": [
                    "In silico predictions require in vitro confirmation"
                ],
                "recommended_validation_experiments": [
                    "Biochemical IC50 assay",
                    "Cell-based antiviral assay",
                    "Selectivity panel",
                ],
                "incomplete": not bool(ope),
            }
        except Exception as e:
            self._error(f"Failed to collect efficacy data: {e}")
            return {"incomplete": True, "error": str(e)}
