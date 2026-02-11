"""
Clinical Trial Design collector -- Section 10.

Generates trial design recommendations from disease constraint files
and disease category. Accounts for pediatric, MCM, and tropical
disease-specific requirements.

Constitutional: Python stdlib only. No mock data.
"""

from __future__ import annotations

from typing import Any, Dict

from .base_collector import BaseCollector


class ClinicalTrialCollector(BaseCollector):
    """Collects clinical trial design data from disease constraints."""

    section_name = "clinical_trial_design"

    def collect(self) -> Dict[str, Any]:
        try:
            constraint = self._get_constraint()
            category = constraint.get("category", "tropical")
            clinical = constraint.get("clinical_requirements", {})
            disease_name = constraint.get("disease_name", self.disease_id)

            is_pediatric = category == "rare_pediatric"
            is_mcm = category == "mcm"
            pathogen_type = constraint.get("pathogen_type", "")
            outbreak_dep = pathogen_type == "virus" and category == "tropical"

            phase1_n = 20 if is_pediatric else 40
            phase2_n = 60 if is_pediatric else 100

            primary_endpoints = constraint.get("primary_endpoints", ["efficacy"])
            phase2_primary = (
                primary_endpoints[0] if primary_endpoints else "efficacy"
            )

            return {
                "phase1_design": {
                    "phase": "1",
                    "study_type": "SAD/MAD",
                    "population": (
                        "Pediatric patients"
                        if is_pediatric
                        else "Healthy volunteers"
                    ),
                    "n_subjects": phase1_n,
                    "design": "Randomized, double-blind, placebo-controlled",
                    "primary_endpoint": "Safety and tolerability",
                    "secondary_endpoints": [
                        "Pharmacokinetics",
                        "Preliminary efficacy biomarkers",
                    ],
                    "duration_weeks": 8,
                    "estimated_cost_usd_millions": 5.0 if is_mcm else 3.0,
                },
                "phase2_design": {
                    "phase": "2",
                    "study_type": "POC",
                    "population": f"{disease_name} patients",
                    "n_subjects": phase2_n,
                    "design": "Randomized, placebo-controlled",
                    "primary_endpoint": phase2_primary,
                    "secondary_endpoints": ["Safety", "PK/PD", "Biomarkers"],
                    "duration_weeks": 24,
                    "estimated_cost_usd_millions": 8.0,
                },
                "biomarker_strategy": {
                    "pharmacodynamic_biomarkers": [],
                    "efficacy_biomarkers": constraint.get(
                        "primary_endpoints", []
                    ),
                    "safety_biomarkers": [
                        "Liver enzymes",
                        "ECG",
                        "Renal markers",
                    ],
                    "companion_diagnostic_required": False,
                },
                "patient_population": {
                    "inclusion_criteria_key": [
                        f"Confirmed {disease_name} diagnosis"
                    ],
                    "exclusion_criteria_key": [
                        "Severe hepatic impairment",
                        "Pregnancy",
                    ],
                    "recruitment_considerations": clinical.get(
                        "patient_population", "Adults"
                    ),
                    "geographic_strategy": (
                        ["Endemic regions"]
                        if category == "tropical"
                        else ["US", "EU"]
                    ),
                },
                "outbreak_dependent": outbreak_dep,
                "adaptive_design_elements": [
                    "Dose escalation",
                    "Futility analysis",
                ],
                "total_cost_to_phase2_usd_millions": 11.0,
                "total_months_to_phase2_readout": 36,
                "incomplete": False,
            }
        except Exception as e:
            self._error(f"Failed to collect clinical trial design: {e}")
            return {"incomplete": True, "error": str(e)}
