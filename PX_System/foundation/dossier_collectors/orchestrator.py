"""
Dossier Collector Orchestrator -- runs all section collectors.

Assembles a complete dossier package by executing each section collector
in sequence. Executive summary runs last since it aggregates other sections.

Constitutional: Python stdlib only. No mock data.
"""

from __future__ import annotations

import sys
from datetime import datetime, timezone
from typing import Any, Dict

from .base_collector import BaseCollector
from .molecule_profile_collector import MoleculeProfileCollector
from .target_validation_collector import TargetValidationCollector
from .efficacy_collector import EfficacyCollector
from .safety_collector import SafetyCollector
from .pkpd_collector import PKPDCollector
from .manufacturing_collector import ManufacturingCollector
from .regulatory_collector import RegulatoryCollector
from .competitive_collector import CompetitiveCollector
from .patent_collector import PatentCollector
from .clinical_trial_collector import ClinicalTrialCollector
from .executive_summary_collector import ExecutiveSummaryCollector


class DossierCollectorOrchestrator:
    """Runs all section collectors and assembles a DossierPackage."""

    COLLECTORS = [
        ("molecule_profile", MoleculeProfileCollector),
        ("target_validation", TargetValidationCollector),
        ("efficacy_data", EfficacyCollector),
        ("safety_profile", SafetyCollector),
        ("pkpd_analysis", PKPDCollector),
        ("manufacturing_assessment", ManufacturingCollector),
        ("regulatory_strategy", RegulatoryCollector),
        ("competitive_landscape", CompetitiveCollector),
        ("patent_analysis", PatentCollector),
        ("clinical_trial_design", ClinicalTrialCollector),
    ]

    def __init__(
        self,
        compound_data: Dict[str, Any],
        disease_id: str,
        tier: str = "DIAMOND",
    ):
        self.compound = compound_data
        self.disease_id = disease_id
        self.tier = tier
        self.errors: list[str] = []
        self.warnings: list[str] = []

    def collect_all(self) -> Dict[str, Any]:
        """Run all collectors and return assembled dossier package."""
        sections: Dict[str, Any] = {}

        # Run the 10 primary section collectors
        for section_name, collector_class in self.COLLECTORS:
            try:
                collector = collector_class(self.compound, self.disease_id)
                sections[section_name] = collector.collect()
                self.errors.extend(collector.errors)
                self.warnings.extend(collector.warnings)
            except Exception as e:
                self.errors.append(f"{section_name}: {e}")
                sections[section_name] = {"error": str(e), "incomplete": True}
                print(
                    f"    ERROR [{section_name}]: {e}", file=sys.stderr
                )

        # Executive summary runs last -- it aggregates other sections
        try:
            exec_collector = ExecutiveSummaryCollector(
                self.compound, self.disease_id, sections
            )
            exec_data = exec_collector.collect()
            exec_data["generation_date"] = datetime.now(
                timezone.utc
            ).isoformat()
            sections["executive_summary"] = exec_data
            self.errors.extend(exec_collector.errors)
            self.warnings.extend(exec_collector.warnings)
        except Exception as e:
            self.errors.append(f"executive_summary: {e}")
            sections["executive_summary"] = {
                "error": str(e),
                "incomplete": True,
            }

        return {
            "schema_version": "1.0.0",
            "tier": self.tier,
            "disease_id": self.disease_id,
            "generation_timestamp": datetime.now(timezone.utc).isoformat(),
            "sections": sections,
            "collection_errors": self.errors,
            "collection_warnings": self.warnings,
            "complete": len(self.errors) == 0,
        }
