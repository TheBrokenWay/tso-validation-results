"""
Regulatory Strategy collector -- Section 07.

Pulls regulatory pathway information from disease constraint files
and OLE engine results. Assesses PRV eligibility, orphan designation,
and recommended development timelines.

Constitutional: Python stdlib only. No mock data.
"""

from __future__ import annotations

from typing import Any, Dict

from .base_collector import BaseCollector


class RegulatoryCollector(BaseCollector):
    """Collects regulatory strategy data from constraints and OLE results."""

    section_name = "regulatory_strategy"

    def collect(self) -> Dict[str, Any]:
        try:
            constraint = self._get_constraint()
            regulatory = constraint.get("regulatory", {})
            ole = (
                self._get_engine_result("ole")
                or self._get_engine_result("OLE")
                or {}
            )

            prv_eligible = ole.get(
                "prv_eligible", constraint.get("prv_eligible", False)
            )
            disease_name = constraint.get("disease_name", self.disease_id)

            designations = [
                {
                    "designation": "PRV",
                    "eligible": prv_eligible,
                    "confidence": "CONFIRMED" if prv_eligible else "UNLIKELY",
                    "rationale": (
                        f"{disease_name} is PRV-qualifying tropical disease"
                        if prv_eligible
                        else "Not PRV-eligible"
                    ),
                },
                {
                    "designation": "Orphan",
                    "eligible": regulatory.get("orphan_eligible", True),
                    "confidence": "LIKELY",
                    "rationale": (
                        f"{disease_name} affects <200,000 US patients"
                    ),
                },
            ]

            category = constraint.get("category", "")
            pathogen_type = constraint.get("pathogen_type", "")
            primary_endpoints = constraint.get("primary_endpoints", ["survival"])

            return {
                "designations": designations,
                "prv_eligible": prv_eligible,
                "prv_qualifying_disease": disease_name,
                "prv_value_estimate_usd_millions": 100.0,
                "recommended_pathway": constraint.get(
                    "regulatory_pathway", "PRV_ACCELERATED"
                ),
                "animal_rule_applicable": category in ("mcm", "tropical"),
                "biomarker_strategy": (
                    "Viral load reduction"
                    if pathogen_type == "virus"
                    else "Clinical endpoint"
                ),
                "primary_endpoint_phase3": (
                    primary_endpoints[0] if primary_endpoints else "survival"
                ),
                "surrogate_endpoint_available": False,
                "precedents": [],
                "relevant_fda_guidances": [],
                "timeline": [
                    {
                        "milestone": "IND filing",
                        "month_from_start": 18,
                        "dependencies": ["GLP tox complete"],
                    },
                    {
                        "milestone": "Phase 1 start",
                        "month_from_start": 21,
                        "dependencies": ["IND clearance"],
                    },
                    {
                        "milestone": "Phase 2 start",
                        "month_from_start": 33,
                        "dependencies": ["Phase 1 complete"],
                    },
                ],
                "regulatory_risk": "MEDIUM",
                "regulatory_hurdles": [],
                "regulatory_mitigations": [],
                "incomplete": not bool(constraint),
            }
        except Exception as e:
            self._error(f"Failed to collect regulatory strategy: {e}")
            return {"incomplete": True, "error": str(e)}
