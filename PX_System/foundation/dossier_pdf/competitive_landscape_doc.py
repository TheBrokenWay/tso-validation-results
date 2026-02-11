"""Competitive Landscape document generator (Section 08)."""
from __future__ import annotations

from .base_generator import BaseDocumentGenerator


class CompetitiveLandscapeDocument(BaseDocumentGenerator):
    """Generates the 08_COMPETITIVE_LANDSCAPE report."""

    section_name = "competitive_landscape"
    section_number = "08"

    def generate(self) -> str:
        data = self.sections_data.get("competitive_landscape", {})
        market = data.get("market_overview", {})
        pipeline = data.get("pipeline_analysis", [])
        approved = data.get("approved_therapies", [])
        opportunity = data.get("opportunity", {})
        differentiation = data.get("differentiation", {})

        out = self._header("COMPETITIVE LANDSCAPE")

        out += self._section_header("MARKET OVERVIEW")
        out += self._kv("Disease Area", market.get("disease_area", "N/A")) + "\n"
        out += self._kv("Unmet Need", market.get("unmet_need", "N/A")) + "\n"
        out += self._kv("Patient Population", f"{market.get('patient_population', 0):,}") + "\n"
        out += self._kv("Current Standard of Care", market.get("standard_of_care", "N/A")) + "\n"
        out += self._kv("Market Size (USD)", f"${market.get('market_size_usd', 0):,}") + "\n"
        out += self._kv("CAGR (%)", f"{market.get('cagr_pct', 0):.1f}") + "\n"

        out += self._section_header("APPROVED THERAPIES")
        if approved:
            rows = [[a.get("drug_name", "N/A"), a.get("company", "N/A"),
                      a.get("mechanism", "N/A"), str(a.get("approval_year", "N/A")),
                      a.get("limitations", "N/A")]
                     for a in approved]
            out += self._table(["Drug", "Company", "Mechanism", "Approved", "Limitations"],
                               rows, [14, 14, 14, 10, 20])
        else:
            out += "  No approved therapies in this indication.\n"

        out += self._section_header("PIPELINE ANALYSIS")
        if pipeline:
            rows = [[p.get("compound", "N/A"), p.get("company", "N/A"),
                      p.get("phase", "N/A"), p.get("mechanism", "N/A"),
                      p.get("expected_approval", "N/A")]
                     for p in pipeline]
            out += self._table(["Compound", "Company", "Phase", "Mechanism", "Exp. Approval"],
                               rows, [14, 14, 8, 16, 14])
        else:
            out += "  No pipeline competitors identified.\n"

        out += self._section_header("OPPORTUNITY ASSESSMENT")
        out += self._kv("Market Opportunity", opportunity.get("assessment", "N/A")) + "\n"
        out += self._kv("First-in-Class", opportunity.get("first_in_class", "N/A")) + "\n"
        out += self._kv("Best-in-Class Potential", opportunity.get("best_in_class", "N/A")) + "\n"
        out += self._kv("Time to Market Advantage", opportunity.get("time_advantage", "N/A")) + "\n"
        out += self._kv("Revenue Potential", f"${opportunity.get('revenue_potential_usd', 0):,}") + "\n"

        out += self._section_header("DIFFERENTIATION STRATEGY")
        advantages = differentiation.get("advantages", [])
        if advantages:
            out += "  Key Advantages:\n"
            for adv in advantages:
                out += f"    + {adv}\n"
        challenges = differentiation.get("challenges", [])
        if challenges:
            out += "\n  Key Challenges:\n"
            for ch in challenges:
                out += f"    - {ch}\n"
        out += "\n"
        out += self._kv("Overall Position", differentiation.get("overall_position", "N/A")) + "\n"

        out += self._footer()
        return out
