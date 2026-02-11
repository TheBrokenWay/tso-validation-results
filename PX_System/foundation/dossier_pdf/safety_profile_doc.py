"""Safety Profile document generator (Section 04)."""
from __future__ import annotations

from .base_generator import BaseDocumentGenerator


class SafetyProfileDocument(BaseDocumentGenerator):
    """Generates the 04_SAFETY_PROFILE report."""

    section_name = "safety_profile"
    section_number = "04"

    def generate(self) -> str:
        data = self.sections_data.get("safety_profile", {})
        tox = data.get("toxicity", {})
        organ_tox = data.get("organ_toxicity", [])
        cardiac = data.get("cardiac_safety", {})
        cyp = data.get("cyp_interactions", [])
        genotox = data.get("genotoxicity", {})

        out = self._header("SAFETY PROFILE")

        # Overall toxicity
        out += self._section_header("TOXICITY OVERVIEW")
        tox_index = tox.get("toxicity_index", 0)
        tox_tier = tox.get("tier", "N/A")
        out += self._kv("Toxicity Index", f"{tox_index:.4f}") + "\n"
        out += self._kv("Toxicity Tier", tox_tier) + "\n"
        out += self._kv("Safety Margin", f"{tox.get('safety_margin', 0):.1f}x") + "\n"
        out += self._kv("LD50 (mg/kg)", f"{tox.get('ld50_mg_kg', 0):.1f}") + "\n"
        out += self._kv("NOAEL (mg/kg)", f"{tox.get('noael_mg_kg', 0):.1f}") + "\n"

        # Toxicity bar (against hard limit 0.0210)
        out += "\n  Toxicity Index vs Hard Limit (Law L11):\n"
        out += f"    {self._metric_bar(tox_index, 0.0210)}\n"
        out += "    Hard limit: 0.0210 (TOXICITY_FAILURE threshold)\n"

        # Organ toxicity table
        out += self._section_header("ORGAN TOXICITY PREDICTIONS")
        if organ_tox:
            rows = []
            for o in organ_tox:
                risk = o.get("risk_level", "N/A")
                flag = " **" if risk in ("HIGH", "CRITICAL") else ""
                rows.append([
                    o.get("organ", "N/A"),
                    f"{o.get('probability', 0):.3f}",
                    risk + flag,
                    o.get("mechanism", "N/A"),
                ])
            out += self._table(
                ["Organ", "Probability", "Risk", "Mechanism"],
                rows,
                [14, 12, 12, 30],
            )
        else:
            out += "  No organ toxicity data available.\n"

        # Cardiac safety
        out += self._section_header("CARDIAC SAFETY")
        out += self._kv("hERG IC50 (uM)", f"{cardiac.get('herg_ic50_um', 0):.2f}") + "\n"
        out += self._kv("hERG Margin", f"{cardiac.get('herg_margin', 0):.1f}x") + "\n"
        out += self._kv("QT Prolongation Risk", cardiac.get("qt_risk", "N/A")) + "\n"
        out += self._kv("Cardiac Ion Channels", cardiac.get("ion_channel_flags", "N/A")) + "\n"
        herg_status = "PASS" if cardiac.get("herg_margin", 0) >= 30 else "CAUTION"
        out += self._kv("Assessment", herg_status) + "\n"
        out += "  (hERG margin >= 30x required for DIAMOND tier)\n"

        # CYP interactions
        out += self._section_header("CYP450 INTERACTIONS")
        if cyp:
            rows = []
            for c in cyp:
                rows.append([
                    c.get("enzyme", "N/A"),
                    c.get("interaction_type", "N/A"),
                    f"{c.get('ic50_um', 0):.2f}",
                    c.get("clinical_significance", "N/A"),
                ])
            out += self._table(
                ["CYP Enzyme", "Interaction", "IC50 (uM)", "Clinical Impact"],
                rows,
                [12, 14, 12, 22],
            )
        else:
            out += "  No significant CYP interactions predicted.\n"

        # Genotoxicity
        out += self._section_header("GENOTOXICITY ASSESSMENT")
        out += self._kv("Ames Test", genotox.get("ames", "N/A")) + "\n"
        out += self._kv("Micronucleus", genotox.get("micronucleus", "N/A")) + "\n"
        out += self._kv("Chromosomal Aberration", genotox.get("chromosomal_aberration", "N/A")) + "\n"
        out += self._kv("Overall Genotoxicity", genotox.get("overall_risk", "N/A")) + "\n"

        out += self._footer()
        return out
