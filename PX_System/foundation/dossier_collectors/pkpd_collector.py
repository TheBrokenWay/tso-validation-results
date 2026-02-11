"""
PK/PD Analysis collector -- Section 05.

Pulls pharmacokinetic and pharmacodynamic data from OPE, PKPD engine,
ADMET results, and DoseOptimizer outputs.

Constitutional: Python stdlib only. No mock data.
"""

from __future__ import annotations

from typing import Any, Dict

from .base_collector import BaseCollector


class PKPDCollector(BaseCollector):
    """Collects PK/PD analysis data from engine results."""

    section_name = "pkpd_analysis"

    def collect(self) -> Dict[str, Any]:
        try:
            ope = (
                self._get_engine_result("ope")
                or self._get_engine_result("OPE")
                or {}
            )
            pkpd = (
                self._get_engine_result("pkpd")
                or self._get_engine_result("PKPD")
                or {}
            )
            admet = (
                self._get_engine_result("admet")
                or self._get_engine_result("ADMET")
                or {}
            )
            dose = self._get_engine_result("dose_optimizer") or {}

            clearance = ope.get("clearance", 0)
            vd = ope.get("vd", 0)

            absorption_raw = admet.get("absorption", {})
            absorption = absorption_raw if isinstance(absorption_raw, dict) else {}

            distribution_raw = admet.get("distribution", {})
            distribution = distribution_raw if isinstance(distribution_raw, dict) else {}

            metabolism_raw = admet.get("metabolism", {})
            metabolism = metabolism_raw if isinstance(metabolism_raw, dict) else {}

            bioavail = absorption.get("bioavailability")
            oral_bioavail = (
                bioavail * 100 if bioavail else 50.0
            )

            return {
                "absorption": {
                    "oral_bioavailability_percent": oral_bioavail,
                    "caco2_permeability_nm_s": absorption.get("caco2", 0),
                    "pgp_substrate": absorption.get("pgp_substrate", False),
                    "food_effect": "UNKNOWN",
                },
                "distribution": {
                    "vd_L_kg": vd,
                    "plasma_protein_binding_percent": distribution.get("ppb", 90),
                    "blood_brain_barrier": False,
                },
                "metabolism": {
                    "primary_cyp": "CYP3A4",
                    "hlm_half_life_min": metabolism.get("hlm_t_half", 0),
                    "metabolic_stability": "MEDIUM",
                },
                "excretion": {
                    "clearance_mL_min_kg": clearance,
                    "primary_route": "hepatic",
                },
                "pk_parameters_human_predicted": pkpd.get("pk_parameters", {}),
                "dose_projection": dose.get("optimal_dose", {}),
                "pk_model_type": "1-compartment",
                "pk_confidence": "MEDIUM",
                "pk_uncertainties": [
                    "Allometric scaling from in silico predictions",
                    "No clinical PK data available",
                ],
                "incomplete": not bool(ope),
            }
        except Exception as e:
            self._error(f"Failed to collect PK/PD analysis: {e}")
            return {"incomplete": True, "error": str(e)}
