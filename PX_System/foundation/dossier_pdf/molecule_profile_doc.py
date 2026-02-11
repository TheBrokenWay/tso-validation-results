"""Molecule Profile document generator (Section 01)."""
from __future__ import annotations

from .base_generator import BaseDocumentGenerator


class MoleculeProfileDocument(BaseDocumentGenerator):
    """Generates the 01_MOLECULE_PROFILE report."""

    section_name = "molecule_profile"
    section_number = "01"

    def generate(self) -> str:
        data = self.sections_data.get("molecule_profile", {})
        structure = data.get("structure", {})
        physchem = data.get("physicochemical", {})
        druglike = data.get("druglikeness", {})
        alerts = data.get("structural_alerts", [])

        out = self._header("MOLECULE PROFILE")

        # Structure information
        out += self._section_header("STRUCTURE INFORMATION")
        out += self._kv("SMILES", structure.get("smiles", "N/A")) + "\n"
        out += self._kv("InChI Key", structure.get("inchi_key", "N/A")) + "\n"
        out += self._kv("Molecular Formula", structure.get("molecular_formula", "N/A")) + "\n"
        out += self._kv("Canonical Name", structure.get("canonical_name", "N/A")) + "\n"
        out += self._kv("Stereochemistry", structure.get("stereochemistry", "N/A")) + "\n"

        # Physicochemical properties
        out += self._section_header("PHYSICOCHEMICAL PROPERTIES")
        props = [
            ("Molecular Weight", physchem.get("molecular_weight", 0), "Da", 500),
            ("LogP", physchem.get("logp", 0), "", 5.0),
            ("TPSA", physchem.get("tpsa", 0), "A^2", 140),
            ("HBD", physchem.get("hbd", 0), "", 5),
            ("HBA", physchem.get("hba", 0), "", 10),
            ("Rotatable Bonds", physchem.get("rotatable_bonds", 0), "", 10),
            ("Aromatic Rings", physchem.get("aromatic_rings", 0), "", 4),
            ("Heavy Atoms", physchem.get("heavy_atoms", 0), "", "N/A"),
        ]
        rows = []
        for name, val, unit, limit in props:
            val_str = f"{val:.2f}" if isinstance(val, float) else str(val)
            limit_str = f"<= {limit}" if isinstance(limit, (int, float)) else str(limit)
            rows.append([name, val_str, unit, limit_str])
        out += self._table(
            ["Property", "Value", "Unit", "Ro5 Limit"],
            rows,
            [22, 10, 6, 12],
        )

        # Drug-likeness assessment
        out += self._section_header("DRUG-LIKENESS ASSESSMENT")
        out += self._kv("Lipinski Violations", druglike.get("lipinski_violations", "N/A")) + "\n"
        out += self._kv("Veber Compliance", druglike.get("veber_compliant", "N/A")) + "\n"
        out += self._kv("QED Score", f"{druglike.get('qed', 0):.4f}") + "\n"
        out += self._kv("Lead-likeness", druglike.get("lead_like", "N/A")) + "\n"
        out += self._kv("PAINS Alerts", druglike.get("pains_count", 0)) + "\n"

        # Structural alerts
        out += self._section_header("STRUCTURAL ALERTS")
        if alerts:
            rows = []
            for a in alerts:
                rows.append([
                    a.get("type", "N/A"),
                    a.get("substructure", "N/A"),
                    a.get("severity", "N/A"),
                    a.get("recommendation", "N/A"),
                ])
            out += self._table(
                ["Alert Type", "Substructure", "Severity", "Recommendation"],
                rows,
                [18, 20, 10, 26],
            )
        else:
            out += "  No structural alerts detected.\n"

        out += self._footer()
        return out
