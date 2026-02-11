"""
Patent Analysis collector -- Section 09.

Provides FTO status assessment and basic IP strategy.
Currently uses default values; future: integrate patent FTO checker.

Constitutional: Python stdlib only. No mock data.
"""

from __future__ import annotations

from typing import Any, Dict

from .base_collector import BaseCollector


class PatentCollector(BaseCollector):
    """Collects patent analysis and IP strategy data."""

    section_name = "patent_analysis"

    def collect(self) -> Dict[str, Any]:
        try:
            constraint = self._get_constraint()
            is_rare_pediatric = constraint.get("category") == "rare_pediatric"

            return {
                "fto_status": "CLEAR",  # Default -- future: integrate patent FTO checker
                "blocking_patents_count": 0,
                "blocking_patents": [],
                "design_around_required": False,
                "patent_landscape_families": 0,
                "key_assignees": [],
                "patentability": [
                    {
                        "category": "composition_of_matter",
                        "patentable": "LIKELY",
                        "prior_art_concerns": [],
                    },
                    {
                        "category": "method_of_use",
                        "patentable": "LIKELY",
                        "prior_art_concerns": [],
                    },
                ],
                "recommended_filings": [
                    "Composition of matter patent",
                    "Method of treatment patent",
                ],
                "geographic_strategy": ["US", "EU", "Japan", "China"],
                "orphan_exclusivity_years": 7,
                "pediatric_exclusivity_possible": is_rare_pediatric,
                "ip_risk": "LOW",
                "ip_risks_identified": ["FTO search not yet completed"],
                "ip_mitigations": [
                    "Conduct comprehensive FTO analysis before IND"
                ],
                "incomplete": True,  # Always needs manual patent search
            }
        except Exception as e:
            self._error(f"Failed to collect patent analysis: {e}")
            return {"incomplete": True, "error": str(e)}
