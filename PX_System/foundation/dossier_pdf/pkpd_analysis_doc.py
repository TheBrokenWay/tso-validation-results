"""PKPD Analysis document generator (Section 05)."""
from __future__ import annotations

from .base_generator import BaseDocumentGenerator


class PKPDAnalysisDocument(BaseDocumentGenerator):
    """Generates the 05_PKPD_ANALYSIS report."""

    section_name = "pkpd_analysis"
    section_number = "05"

    def generate(self) -> str:
        data = self.sections_data.get("pkpd_analysis", {})
        adme = data.get("adme_summary", {})
        pk_params = data.get("pk_parameters", {})
        compartment = data.get("compartmental_model", {})
        dose_proj = data.get("dose_projection", {})

        out = self._header("PKPD ANALYSIS")

        # ADME summary
        out += self._section_header("ADME SUMMARY")
        absorption = adme.get("absorption", {})
        distribution = adme.get("distribution", {})
        metabolism = adme.get("metabolism", {})
        excretion = adme.get("excretion", {})

        out += "  Absorption:\n"
        out += self._kv("Bioavailability (%)", f"{absorption.get('bioavailability_pct', 0):.1f}", indent=4) + "\n"
        out += self._kv("Caco-2 Permeability", f"{absorption.get('caco2_permeability', 0):.2e}", indent=4) + "\n"
        out += self._kv("Oral Absorption", absorption.get("oral_absorption", "N/A"), indent=4) + "\n"

        out += "\n  Distribution:\n"
        out += self._kv("Vd (L/kg)", f"{distribution.get('vd_l_kg', 0):.2f}", indent=4) + "\n"
        out += self._kv("Protein Binding (%)", f"{distribution.get('protein_binding_pct', 0):.1f}", indent=4) + "\n"
        out += self._kv("BBB Penetration", distribution.get("bbb_penetration", "N/A"), indent=4) + "\n"

        out += "\n  Metabolism:\n"
        out += self._kv("Primary CYP", metabolism.get("primary_cyp", "N/A"), indent=4) + "\n"
        out += self._kv("Clearance (mL/min/kg)", f"{metabolism.get('clearance_ml_min_kg', 0):.2f}", indent=4) + "\n"
        out += self._kv("Metabolic Stability", metabolism.get("stability", "N/A"), indent=4) + "\n"

        out += "\n  Excretion:\n"
        out += self._kv("Half-life (h)", f"{excretion.get('half_life_h', 0):.1f}", indent=4) + "\n"
        out += self._kv("Route", excretion.get("route", "N/A"), indent=4) + "\n"
        out += self._kv("Renal Clearance (%)", f"{excretion.get('renal_clearance_pct', 0):.1f}", indent=4) + "\n"

        # PK parameters table
        out += self._section_header("PHARMACOKINETIC PARAMETERS")
        pk_rows = []
        for param, details in pk_params.items():
            if isinstance(details, dict):
                pk_rows.append([
                    param,
                    f"{details.get('value', 0):.4f}",
                    details.get("unit", ""),
                    details.get("range", "N/A"),
                ])
            else:
                pk_rows.append([param, str(details), "", ""])
        if pk_rows:
            out += self._table(
                ["Parameter", "Value", "Unit", "Range"],
                pk_rows,
                [22, 12, 14, 20],
            )
        else:
            out += "  No PK parameters available.\n"

        # Compartmental model
        out += self._section_header("COMPARTMENTAL MODEL")
        out += self._kv("Model Type", compartment.get("model_type", "1-compartment")) + "\n"
        out += self._kv("Ka (1/h)", f"{compartment.get('ka', 0):.4f}") + "\n"
        out += self._kv("Ke (1/h)", f"{compartment.get('ke', 0):.4f}") + "\n"
        out += self._kv("AUC (ng*h/mL)", f"{compartment.get('auc', 0):.2f}") + "\n"
        out += self._kv("Cmax (ng/mL)", f"{compartment.get('cmax', 0):.2f}") + "\n"
        out += self._kv("Tmax (h)", f"{compartment.get('tmax', 0):.1f}") + "\n"

        # Dose projection
        out += self._section_header("DOSE PROJECTION")
        out += self._kv("Recommended Dose", f"{dose_proj.get('dose_mg', 0):.1f} mg") + "\n"
        out += self._kv("Dosing Interval", dose_proj.get("interval", "N/A")) + "\n"
        out += self._kv("Steady-State Cmin", f"{dose_proj.get('ss_cmin', 0):.2f} ng/mL") + "\n"
        out += self._kv("Steady-State Cmax", f"{dose_proj.get('ss_cmax', 0):.2f} ng/mL") + "\n"
        out += self._kv("Accumulation Ratio", f"{dose_proj.get('accumulation_ratio', 0):.2f}") + "\n"

        out += self._footer()
        return out
