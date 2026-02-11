"""Patent Analysis document generator (Section 09)."""
from __future__ import annotations

from .base_generator import BaseDocumentGenerator


class PatentAnalysisDocument(BaseDocumentGenerator):
    """Generates the 09_PATENT_ANALYSIS report."""

    section_name = "patent_analysis"
    section_number = "09"

    def generate(self) -> str:
        data = self.sections_data.get("patent_analysis", {})
        fto = data.get("fto_status", {})
        blocking = data.get("blocking_patents", [])
        patentability = data.get("patentability", {})
        ip_strategy = data.get("ip_strategy", {})
        landscape = data.get("patent_landscape", {})

        out = self._header("PATENT & IP ANALYSIS")

        out += self._section_header("FREEDOM TO OPERATE (FTO)")
        out += self._kv("FTO Status", fto.get("status", "N/A")) + "\n"
        out += self._kv("Risk Level", fto.get("risk_level", "N/A")) + "\n"
        out += self._kv("Jurisdiction", fto.get("jurisdiction", "N/A")) + "\n"
        out += self._kv("Analysis Date", fto.get("analysis_date", "N/A")) + "\n"
        out += self._kv("Confidence", f"{fto.get('confidence_pct', 0):.0f}%") + "\n"
        summary = fto.get("summary", "")
        if summary:
            out += f"\n  Summary:\n    {summary}\n"

        out += self._section_header("BLOCKING PATENTS")
        if blocking:
            for i, p in enumerate(blocking, 1):
                out += f"\n  [{i}] {p.get('patent_number', 'N/A')}\n"
                out += self._kv("Title", p.get("title", "N/A"), indent=4) + "\n"
                out += self._kv("Assignee", p.get("assignee", "N/A"), indent=4) + "\n"
                out += self._kv("Expiry", p.get("expiry_date", "N/A"), indent=4) + "\n"
                out += self._kv("Claims Overlap", p.get("claims_overlap", "N/A"), indent=4) + "\n"
                out += self._kv("Risk", p.get("risk_level", "N/A"), indent=4) + "\n"
                workaround = p.get("design_around", "")
                if workaround:
                    out += self._kv("Design-Around", workaround, indent=4) + "\n"
        else:
            out += "  No blocking patents identified.\n"

        out += self._section_header("PATENTABILITY ASSESSMENT")
        out += self._kv("Novelty", patentability.get("novelty", "N/A")) + "\n"
        out += self._kv("Non-Obviousness", patentability.get("non_obviousness", "N/A")) + "\n"
        out += self._kv("Enablement", patentability.get("enablement", "N/A")) + "\n"
        out += self._kv("Written Description", patentability.get("written_description", "N/A")) + "\n"
        out += self._kv("Overall Patentability", patentability.get("overall", "N/A")) + "\n"
        claims = patentability.get("potential_claims", [])
        if claims:
            out += "\n  Potential Patent Claims:\n"
            for c in claims:
                out += f"    - {c}\n"

        out += self._section_header("IP STRATEGY")
        out += self._kv("Filing Strategy", ip_strategy.get("filing_strategy", "N/A")) + "\n"
        out += self._kv("Priority Date Target", ip_strategy.get("priority_date", "N/A")) + "\n"
        out += self._kv("PCT Filing", ip_strategy.get("pct_filing", "N/A")) + "\n"
        out += self._kv("Key Jurisdictions", ip_strategy.get("jurisdictions", "N/A")) + "\n"
        out += self._kv("Patent Term", ip_strategy.get("patent_term", "N/A")) + "\n"
        out += self._kv("Extension Eligible", ip_strategy.get("extension_eligible", "N/A")) + "\n"

        out += self._section_header("PATENT LANDSCAPE SUMMARY")
        out += self._kv("Total Patents in Space", landscape.get("total_patents", "N/A")) + "\n"
        out += self._kv("Active Patent Families", landscape.get("active_families", "N/A")) + "\n"
        out += self._kv("Key Assignees", landscape.get("key_assignees", "N/A")) + "\n"
        out += self._kv("Trend", landscape.get("filing_trend", "N/A")) + "\n"
        out += self._kv("White Space Identified", landscape.get("white_space", "N/A")) + "\n"

        out += self._footer()
        return out
