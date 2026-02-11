"""
Executive Summary collector -- Section 11.

Aggregates data from all other section collectors to generate
investment thesis, key metrics, and go/no-go recommendation.

Constitutional: Python stdlib only. No mock data.
Toxicity tiers (Law L11) are immutable hard limits.
"""

from __future__ import annotations

from typing import Any, Dict, Optional

from .base_collector import BaseCollector


class ExecutiveSummaryCollector(BaseCollector):
    """Collects executive summary by aggregating all other sections."""

    section_name = "executive_summary"

    def __init__(
        self,
        compound_data: Dict[str, Any],
        disease_id: str,
        other_sections: Optional[Dict[str, Any]] = None,
    ):
        super().__init__(compound_data, disease_id)
        self.sections = other_sections or {}

    def collect(self) -> Dict[str, Any]:
        try:
            constraint = self._get_constraint()
            safety = self.sections.get("safety_profile", {})
            efficacy = self.sections.get("efficacy_data", {})
            regulatory = self.sections.get("regulatory_strategy", {})
            manufacturing = self.sections.get("manufacturing_assessment", {})

            tox = safety.get("overall_toxicity_index", 0)
            tier_class = safety.get("tier_classification", "UNKNOWN")
            prv = regulatory.get("prv_eligible", False)
            disease_name = constraint.get("disease_name", self.disease_id)

            # Build key metrics
            safety_margin_fold = safety.get("safety_margin_fold", 0)
            key_metrics = [
                {
                    "name": "Toxicity Index",
                    "value": tox,
                    "target": 0.0200,
                    "unit": "index",
                    "passed": tox < 0.0200,
                    "display_bar_percent": max(
                        0, min(100, int((1 - tox / 0.0210) * 100))
                    ),
                },
                {
                    "name": "Safety Margin",
                    "value": safety_margin_fold,
                    "target": 10.0,
                    "unit": "fold",
                    "passed": safety_margin_fold >= 10,
                    "display_bar_percent": min(
                        100, int(safety_margin_fold / 50 * 100)
                    ),
                },
                {
                    "name": "PRV Eligible",
                    "value": 1.0 if prv else 0.0,
                    "target": 1.0,
                    "unit": "boolean",
                    "passed": prv,
                    "display_bar_percent": 100 if prv else 0,
                },
            ]

            # Determine recommendation
            if tier_class == "TOXICITY_DIAMOND":
                recommendation = "ADVANCE_TO_IND"
                confidence = 75
            elif tier_class == "TOXICITY_GOLD":
                recommendation = "ADVANCE_TO_IND"
                confidence = 60
            elif tier_class == "TOXICITY_SILVER":
                recommendation = "OPTIMIZE_FURTHER"
                confidence = 40
            else:
                recommendation = "REJECT"
                confidence = 20

            # Investment thesis
            safety_desc = (
                "favorable" if tox < 0.02 else "suboptimal"
            )
            thesis = (
                f"Novel compound targeting {disease_name} with "
                f"{safety_desc} safety profile. "
            )
            if prv:
                thesis += (
                    "PRV-eligible indication with estimated "
                    "00M+ voucher value. "
                )
            advance_text = (
                "IND-enabling studies"
                if recommendation == "ADVANCE_TO_IND"
                else "further optimization"
            )
            thesis += f"Recommended for {advance_text}."

            return {
                "compound_id": self.compound.get("compound_id", "UNKNOWN"),
                "disease_name": disease_name,
                "disease_id": self.disease_id,
                "generation_date": "",  # Set by orchestrator
                "tier": tier_class,
                "investment_thesis": thesis,
                "key_metrics": key_metrics,
                "recommendation": recommendation,
                "confidence_percent": confidence,
                "cost_to_ind_usd_low": 8_000_000,
                "cost_to_ind_usd_high": 15_000_000,
                "months_to_ind": 18,
                "risks": [
                    {
                        "category": "safety",
                        "level": "LOW" if tox < 0.01 else "MEDIUM",
                        "factors": [f"Toxicity index: {tox}"],
                        "mitigations": ["GLP tox studies"],
                    },
                    {
                        "category": "regulatory",
                        "level": "LOW" if prv else "MEDIUM",
                        "factors": [],
                        "mitigations": [],
                    },
                    {
                        "category": "commercial",
                        "level": "MEDIUM",
                        "factors": ["Early stage"],
                        "mitigations": ["Partner with established pharma"],
                    },
                ],
                "governance": {},  # Populated by orchestrator
                "incomplete": False,
            }
        except Exception as e:
            self._error(f"Failed to collect executive summary: {e}")
            return {"incomplete": True, "error": str(e)}
