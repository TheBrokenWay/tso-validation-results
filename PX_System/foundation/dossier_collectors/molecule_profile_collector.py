"""
Molecule Profile collector -- Section 01.

Gathers physicochemical properties, drug-likeness assessment,
and structural information from OPE engine results.

Constitutional: Python stdlib only. No mock data.
"""

from __future__ import annotations

from typing import Any, Dict

from .base_collector import BaseCollector


class MoleculeProfileCollector(BaseCollector):
    """Collects molecular profile data from OPE engine results."""

    section_name = "molecule_profile"

    def collect(self) -> Dict[str, Any]:
        try:
            ope = self._get_engine_result("ope") or self._get_engine_result("OPE") or {}
            smiles = self.compound.get("smiles", "")

            return {
                "smiles": smiles,
                "canonical_smiles": smiles,
                "molecular_weight": ope.get("molecular_weight", 0),
                "logp": ope.get("logp", 0),
                "hbd": ope.get("hbd", 0),
                "hba": ope.get("hba", 0),
                "tpsa": ope.get("tpsa", 0),
                "rotatable_bonds": ope.get("rotatable_bonds", 0),
                "physicochemical": {
                    "logp": ope.get("logp", 0),
                    "psa": ope.get("tpsa", 0),
                    "hbd": ope.get("hbd", 0),
                    "hba": ope.get("hba", 0),
                    "rotatable_bonds": ope.get("rotatable_bonds", 0),
                    "aromatic_rings": ope.get("aromatic_rings", 0),
                    "heavy_atoms": ope.get("heavy_atoms", 0),
                    "fraction_sp3": ope.get("fraction_sp3", 0.0),
                },
                "druglikeness": {
                    "lipinski_violations": self._count_lipinski_violations(ope),
                    "lipinski_pass": self._count_lipinski_violations(ope) <= 1,
                    "qed_score": ope.get("qed", 0.0),
                },
                "structural_alerts": [],
                "nearest_known_drugs": [],
                "incomplete": not bool(smiles),
            }
        except Exception as e:
            self._error(f"Failed to collect molecule profile: {e}")
            return {"incomplete": True, "error": str(e)}

    def _count_lipinski_violations(self, ope: Dict[str, Any]) -> int:
        """Count Lipinski Rule-of-Five violations."""
        violations = 0
        if ope.get("molecular_weight", 0) > 500:
            violations += 1
        if ope.get("logp", 0) > 5:
            violations += 1
        if ope.get("hbd", 0) > 5:
            violations += 1
        if ope.get("hba", 0) > 10:
            violations += 1
        return violations
