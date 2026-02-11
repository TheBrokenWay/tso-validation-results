"""Executive Summary document generator (Section 00)."""
from __future__ import annotations

from .base_generator import BaseDocumentGenerator


class ExecutiveSummaryDocument(BaseDocumentGenerator):
    """Generates the 00_EXECUTIVE_SUMMARY report."""

    section_name = "executive_summary"
    section_number = "00"

    def generate(self) -> str:
        data = self.sections_data.get("executive_summary", {})

        out = self._header(
            f"EXECUTIVE SUMMARY -- {data.get('compound_id', 'UNKNOWN')}"
        )

        # Compound overview
        out += self._section_header("COMPOUND OVERVIEW")
        out += self._kv("Compound ID", data.get("compound_id", "N/A")) + "\n"
        out += self._kv("Disease", data.get("disease_name", "N/A")) + "\n"
        out += self._kv("Tier", data.get("tier", "N/A")) + "\n"
        out += self._kv("Recommendation", data.get("recommendation", "N/A")) + "\n"
        out += self._kv("Confidence", f"{data.get('confidence_percent', 0)}%") + "\n"

        # Investment thesis
        out += self._section_header("INVESTMENT THESIS")
        out += f"  {data.get('investment_thesis', 'N/A')}\n"

        # Key metrics
        out += self._section_header("KEY METRICS")
        for m in data.get("key_metrics", []):
            passed = "PASS" if m.get("passed") else "FAIL"
            out += (
                f"  {m.get('name', ''):25s} "
                f"{m.get('value', 0):>10.4f} / "
                f"{m.get('target', 0):>10.4f} [{passed}]\n"
            )

        # Cost and timeline
        out += self._section_header("COST & TIMELINE")
        out += (
            self._kv(
                "Cost to IND",
                f"${data.get('cost_to_ind_usd_low', 0):,} - "
                f"${data.get('cost_to_ind_usd_high', 0):,}",
            )
            + "\n"
        )
        out += self._kv("Months to IND", data.get("months_to_ind", 0)) + "\n"

        # Risk assessment
        out += self._section_header("RISK ASSESSMENT")
        for r in data.get("risks", []):
            factors = ", ".join(r.get("factors", []))
            out += (
                f"  {r.get('category', ''):15s} "
                f"[{r.get('level', 'N/A'):6s}]  {factors}\n"
            )

        out += self._footer()
        return out
