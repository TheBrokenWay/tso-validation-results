"""Manufacturing Assessment document generator (Section 06)."""
from __future__ import annotations

from .base_generator import BaseDocumentGenerator


class ManufacturingAssessmentDocument(BaseDocumentGenerator):
    """Generates the 06_MANUFACTURING_ASSESSMENT report."""

    section_name = "manufacturing_assessment"
    section_number = "06"

    def generate(self) -> str:
        data = self.sections_data.get("manufacturing_assessment", {})
        complexity = data.get("complexity", {})
        synthesis = data.get("synthesis_route", {})
        cogs = data.get("cogs", {})
        cmc_risks = data.get("cmc_risks", [])
        scalability = data.get("scalability", {})

        out = self._header("MANUFACTURING ASSESSMENT")

        out += self._section_header("SYNTHETIC COMPLEXITY")
        out += self._kv("Complexity Score", f"{complexity.get('score', 0):.2f} / 10") + "\n"
        out += self._kv("Number of Steps", complexity.get("num_steps", "N/A")) + "\n"
        out += self._kv("Longest Linear Sequence", complexity.get("longest_linear_sequence", "N/A")) + "\n"
        out += self._kv("Chiral Centers", complexity.get("chiral_centers", 0)) + "\n"
        out += self._kv("Protection Steps", complexity.get("protection_steps", 0)) + "\n"
        out += self._kv("Overall Yield (%)", f"{complexity.get('overall_yield_pct', 0):.1f}") + "\n"

        out += self._section_header("SYNTHESIS ROUTE")
        steps = synthesis.get("steps", [])
        if steps:
            for i, step in enumerate(steps, 1):
                out += f"\n  Step {i}: {step.get('description', 'N/A')}\n"
                reagents = ", ".join(step.get("reagents", []))
                out += self._kv("Reagents", reagents, indent=4) + "\n"
                out += self._kv("Conditions", step.get("conditions", "N/A"), indent=4) + "\n"
                out += self._kv("Yield (%)", f"{step.get('yield_pct', 0):.1f}", indent=4) + "\n"
                out += self._kv("Hazard Class", step.get("hazard_class", "N/A"), indent=4) + "\n"
        else:
            out += "  No synthesis route data available.\n"
        out += self._kv("Starting Materials", synthesis.get("starting_materials", "N/A")) + "\n"
        out += self._kv("Key Intermediates", synthesis.get("key_intermediates", "N/A")) + "\n"

        out += self._section_header("COST OF GOODS (COGS)")
        out += self._kv("API Cost ($/kg)", f"${cogs.get('api_cost_per_kg', 0):,.0f}") + "\n"
        out += self._kv("Formulation Cost ($/unit)", f"${cogs.get('formulation_cost_per_unit', 0):.2f}") + "\n"
        out += self._kv("Total COGS ($/unit)", f"${cogs.get('total_cogs_per_unit', 0):.2f}") + "\n"
        out += self._kv("Scale", cogs.get("scale", "N/A")) + "\n"
        out += self._kv("Supply Chain Risk", cogs.get("supply_chain_risk", "N/A")) + "\n"

        out += self._section_header("CMC RISKS")
        if cmc_risks:
            rows = [[r.get("category", "N/A"), r.get("risk_level", "N/A"),
                      r.get("description", "N/A"), r.get("mitigation", "N/A")]
                     for r in cmc_risks]
            out += self._table(["Category", "Risk", "Description", "Mitigation"],
                               rows, [14, 8, 26, 22])
        else:
            out += "  No significant CMC risks identified.\n"

        out += self._section_header("SCALABILITY ASSESSMENT")
        out += self._kv("Lab Scale", scalability.get("lab_scale", "N/A")) + "\n"
        out += self._kv("Pilot Scale", scalability.get("pilot_scale", "N/A")) + "\n"
        out += self._kv("Commercial Scale", scalability.get("commercial_scale", "N/A")) + "\n"
        out += self._kv("GMP Readiness", scalability.get("gmp_readiness", "N/A")) + "\n"

        out += self._footer()
        return out
