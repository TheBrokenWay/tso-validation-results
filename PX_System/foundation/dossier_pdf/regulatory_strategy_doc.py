"""Regulatory Strategy document generator (Section 07)."""
from __future__ import annotations

from .base_generator import BaseDocumentGenerator


class RegulatoryStrategyDocument(BaseDocumentGenerator):
    """Generates the 07_REGULATORY_STRATEGY report."""

    section_name = "regulatory_strategy"
    section_number = "07"

    def generate(self) -> str:
        data = self.sections_data.get("regulatory_strategy", {})
        designations = data.get("designations", [])
        pathway = data.get("pathway", {})
        timeline = data.get("timeline", {})
        prv = data.get("priority_review_voucher", {})
        precedent = data.get("regulatory_precedent", [])

        out = self._header("REGULATORY STRATEGY")

        out += self._section_header("REGULATORY DESIGNATIONS")
        if designations:
            rows = [[d.get("designation", "N/A"), d.get("agency", "N/A"),
                      d.get("status", "N/A"), d.get("benefit", "N/A")]
                     for d in designations]
            out += self._table(["Designation", "Agency", "Status", "Key Benefit"],
                               rows, [22, 8, 14, 28])
        else:
            out += "  No regulatory designations identified.\n"

        out += self._section_header("APPROVAL PATHWAY")
        out += self._kv("Primary Pathway", pathway.get("primary", "N/A")) + "\n"
        out += self._kv("Accelerated Approval", pathway.get("accelerated", "N/A")) + "\n"
        out += self._kv("Breakthrough Therapy", pathway.get("breakthrough", "N/A")) + "\n"
        out += self._kv("Fast Track", pathway.get("fast_track", "N/A")) + "\n"
        out += self._kv("Orphan Drug", pathway.get("orphan_drug", "N/A")) + "\n"
        out += self._kv("Animal Rule (21 CFR 314.600)", pathway.get("animal_rule", "N/A")) + "\n"
        out += self._kv("Recommended Approach", pathway.get("recommended_approach", "N/A")) + "\n"

        out += self._section_header("REGULATORY TIMELINE")
        milestones = timeline.get("milestones", [])
        if milestones:
            rows = [[m.get("milestone", "N/A"), m.get("target_date", "N/A"),
                      f"{m.get('duration_months', 0)} months", m.get("status", "N/A")]
                     for m in milestones]
            out += self._table(["Milestone", "Target Date", "Duration", "Status"],
                               rows, [24, 14, 12, 14])
        else:
            out += self._kv("Pre-IND Meeting", timeline.get("pre_ind_meeting", "N/A")) + "\n"
            out += self._kv("IND Filing", timeline.get("ind_filing", "N/A")) + "\n"
            out += self._kv("Phase 1 Start", timeline.get("phase1_start", "N/A")) + "\n"
            out += self._kv("NDA/BLA Submission", timeline.get("nda_submission", "N/A")) + "\n"

        out += self._section_header("PRIORITY REVIEW VOUCHER (PRV)")
        out += self._kv("PRV Eligible", prv.get("eligible", "N/A")) + "\n"
        out += self._kv("Qualifying Disease", prv.get("qualifying_disease", "N/A")) + "\n"
        out += self._kv("PRV Estimated Value", f"${prv.get('estimated_value_usd', 0):,}") + "\n"
        out += self._kv("Precedent Sales", prv.get("precedent_sales", "N/A")) + "\n"
        out += self._kv("Strategic Value", prv.get("strategic_value", "N/A")) + "\n"

        out += self._section_header("REGULATORY PRECEDENT")
        if precedent:
            for p in precedent:
                out += f"\n  {p.get('drug_name', 'N/A')} ({p.get('year', 'N/A')})\n"
                out += self._kv("Indication", p.get("indication", "N/A"), indent=4) + "\n"
                out += self._kv("Pathway", p.get("pathway", "N/A"), indent=4) + "\n"
                out += self._kv("Relevance", p.get("relevance", "N/A"), indent=4) + "\n"
        else:
            out += "  No regulatory precedent data available.\n"

        out += self._footer()
        return out
