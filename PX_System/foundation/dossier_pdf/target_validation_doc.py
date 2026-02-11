"""Target Validation document generator (Section 02)."""
from __future__ import annotations

from .base_generator import BaseDocumentGenerator


class TargetValidationDocument(BaseDocumentGenerator):
    """Generates the 02_TARGET_VALIDATION report."""

    section_name = "target_validation"
    section_number = "02"

    def generate(self) -> str:
        data = self.sections_data.get("target_validation", {})
        target = data.get("target", {})
        mechanism = data.get("mechanism", {})
        tissues = data.get("tissue_targets", [])
        validation = data.get("validation_evidence", {})

        out = self._header("TARGET VALIDATION")

        # Target information
        out += self._section_header("PRIMARY TARGET")
        out += self._kv("Target Name", target.get("name", "N/A")) + "\n"
        out += self._kv("Gene Symbol", target.get("gene_symbol", "N/A")) + "\n"
        out += self._kv("UniProt ID", target.get("uniprot_id", "N/A")) + "\n"
        out += self._kv("Target Class", target.get("target_class", "N/A")) + "\n"
        out += self._kv("Organism", target.get("organism", "N/A")) + "\n"
        out += self._kv("Druggability Score", f"{target.get('druggability', 0):.2f}") + "\n"

        # Mechanism of action
        out += self._section_header("MECHANISM OF ACTION")
        out += self._kv("Action Type", mechanism.get("action_type", "N/A")) + "\n"
        out += self._kv("Binding Site", mechanism.get("binding_site", "N/A")) + "\n"
        out += self._kv("Selectivity", mechanism.get("selectivity", "N/A")) + "\n"
        out += self._kv("Reversibility", mechanism.get("reversibility", "N/A")) + "\n"
        out += "\n  Description:\n"
        desc = mechanism.get("description", "N/A")
        # Word-wrap description at ~74 chars
        words = str(desc).split()
        line = "    "
        for w in words:
            if len(line) + len(w) + 1 > 76:
                out += line + "\n"
                line = "    " + w
            else:
                line += " " + w if len(line) > 4 else w
        out += line + "\n"

        # Tissue targets
        out += self._section_header("TISSUE EXPRESSION & TARGETS")
        if tissues:
            rows = []
            for t in tissues:
                rows.append([
                    t.get("tissue", "N/A"),
                    f"{t.get('expression_level', 0):.2f}",
                    t.get("relevance", "N/A"),
                    t.get("off_target_risk", "N/A"),
                ])
            out += self._table(
                ["Tissue", "Expression", "Relevance", "Off-Target Risk"],
                rows,
                [18, 12, 20, 18],
            )
        else:
            out += "  No tissue expression data available.\n"

        # Validation confidence
        out += self._section_header("VALIDATION CONFIDENCE")
        out += self._kv("Overall Confidence", f"{validation.get('confidence_score', 0):.2f}") + "\n"
        out += self._kv("Genetic Evidence", validation.get("genetic_evidence", "N/A")) + "\n"
        out += self._kv("Clinical Precedent", validation.get("clinical_precedent", "N/A")) + "\n"
        out += self._kv("Literature Support", validation.get("literature_support", "N/A")) + "\n"
        out += self._kv("Biomarker Available", validation.get("biomarker_available", "N/A")) + "\n"

        pathways = validation.get("pathways", [])
        if pathways:
            out += "\n  Associated Pathways:\n"
            for p in pathways:
                out += f"    - {p}\n"

        out += self._footer()
        return out
