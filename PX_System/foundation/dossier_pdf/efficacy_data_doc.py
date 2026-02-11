"""Efficacy Data document generator (Section 03)."""
from __future__ import annotations

from .base_generator import BaseDocumentGenerator


class EfficacyDataDocument(BaseDocumentGenerator):
    """Generates the 03_EFFICACY_DATA report."""

    section_name = "efficacy_data"
    section_number = "03"

    def generate(self) -> str:
        data = self.sections_data.get("efficacy_data", {})
        activity = data.get("activity_predictions", {})
        benchmarks = data.get("benchmarks", [])
        dose_proj = data.get("dose_projection", {})
        uncertainties = data.get("uncertainties", [])

        out = self._header("EFFICACY DATA")

        # Activity predictions
        out += self._section_header("ACTIVITY PREDICTIONS")
        out += self._kv("IC50 (nM)", f"{activity.get('ic50_nm', 0):.2f}") + "\n"
        out += self._kv("EC50 (nM)", f"{activity.get('ec50_nm', 0):.2f}") + "\n"
        out += self._kv("Ki (nM)", f"{activity.get('ki_nm', 0):.2f}") + "\n"
        out += self._kv("Emax (%)", f"{activity.get('emax_pct', 0):.1f}") + "\n"
        out += self._kv("Binding Affinity", f"{activity.get('binding_affinity', 0):.4f}") + "\n"
        out += self._kv("Selectivity Index", f"{activity.get('selectivity_index', 0):.2f}") + "\n"
        out += self._kv("pChEMBL", f"{activity.get('pchembl', 0):.2f}") + "\n"

        # Activity bar
        ic50 = activity.get("ic50_nm", 0)
        if ic50 > 0:
            out += "\n  IC50 Potency:\n"
            out += f"    {self._metric_bar(ic50, 100.0)}\n"
            out += "    (Target: < 100 nM for DIAMOND tier)\n"

        # Benchmarks
        out += self._section_header("BENCHMARK COMPARISONS")
        if benchmarks:
            rows = []
            for b in benchmarks:
                rows.append([
                    b.get("compound", "N/A"),
                    f"{b.get('ic50_nm', 0):.2f}",
                    b.get("status", "N/A"),
                    f"{b.get('fold_improvement', 0):.1f}x",
                ])
            out += self._table(
                ["Benchmark", "IC50 (nM)", "Status", "Fold Change"],
                rows,
                [22, 12, 14, 12],
            )
        else:
            out += "  No benchmark data available.\n"

        # Dose projection
        out += self._section_header("DOSE PROJECTION")
        out += self._kv("Projected Human Dose", f"{dose_proj.get('human_dose_mg', 0):.1f} mg") + "\n"
        out += self._kv("Dosing Frequency", dose_proj.get("frequency", "N/A")) + "\n"
        out += self._kv("Route", dose_proj.get("route", "N/A")) + "\n"
        out += self._kv("Therapeutic Index", f"{dose_proj.get('therapeutic_index', 0):.1f}") + "\n"
        out += self._kv("Confidence", dose_proj.get("confidence", "N/A")) + "\n"

        # Uncertainties
        out += self._section_header("KEY UNCERTAINTIES")
        if uncertainties:
            for i, u in enumerate(uncertainties, 1):
                out += f"  {i}. [{u.get('severity', 'N/A'):6s}] {u.get('description', 'N/A')}\n"
                mitigation = u.get("mitigation", "")
                if mitigation:
                    out += f"     Mitigation: {mitigation}\n"
        else:
            out += "  No significant uncertainties identified.\n"

        out += self._footer()
        return out
