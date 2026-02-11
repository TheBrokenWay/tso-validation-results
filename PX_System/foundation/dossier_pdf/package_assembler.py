"""Package assembler -- creates complete DIAMOND dossier folder."""
from __future__ import annotations

import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict

from .executive_summary_doc import ExecutiveSummaryDocument
from .molecule_profile_doc import MoleculeProfileDocument
from .target_validation_doc import TargetValidationDocument
from .efficacy_data_doc import EfficacyDataDocument
from .safety_profile_doc import SafetyProfileDocument
from .pkpd_analysis_doc import PKPDAnalysisDocument
from .manufacturing_assessment_doc import ManufacturingAssessmentDocument
from .regulatory_strategy_doc import RegulatoryStrategyDocument
from .competitive_landscape_doc import CompetitiveLandscapeDocument
from .patent_analysis_doc import PatentAnalysisDocument
from .clinical_trial_design_doc import ClinicalTrialDesignDocument

DOCUMENT_GENERATORS = [
    ("00_EXECUTIVE_SUMMARY", ExecutiveSummaryDocument),
    ("01_MOLECULE_PROFILE", MoleculeProfileDocument),
    ("02_TARGET_VALIDATION", TargetValidationDocument),
    ("03_EFFICACY_DATA", EfficacyDataDocument),
    ("04_SAFETY_PROFILE", SafetyProfileDocument),
    ("05_PKPD_ANALYSIS", PKPDAnalysisDocument),
    ("06_MANUFACTURING_ASSESSMENT", ManufacturingAssessmentDocument),
    ("07_REGULATORY_STRATEGY", RegulatoryStrategyDocument),
    ("08_COMPETITIVE_LANDSCAPE", CompetitiveLandscapeDocument),
    ("09_PATENT_ANALYSIS", PatentAnalysisDocument),
    ("10_CLINICAL_TRIAL_DESIGN", ClinicalTrialDesignDocument),
]


class DossierPackageAssembler:
    """Assembles a complete DIAMOND dossier document package.

    Creates a directory structure with formatted text reports for each
    dossier section, plus raw data and governance certificates.
    """

    def __init__(self, dossier_data: Dict[str, Any] | None = None):
        self._data = dossier_data

    def generate_all(self) -> Dict[str, str]:
        """Generate all section documents as in-memory text strings.

        Returns dict mapping section name to generated report text.
        Does NOT write to disk.
        """
        if self._data is None:
            return {}
        documents: Dict[str, str] = {}
        for filename, generator_class in DOCUMENT_GENERATORS:
            try:
                generator = generator_class(self._data)
                documents[filename] = generator.generate()
            except Exception as e:
                documents[filename] = f"[Generation error: {e}]"
        return documents

    def assemble(self, dossier_data: Dict[str, Any], output_dir: str) -> str:
        """Assemble a full dossier package into output_dir.

        Returns the path to the created package directory.
        """
        exec_data = dossier_data.get("sections", {}).get("executive_summary", {})
        compound_id = exec_data.get("compound_id", "UNKNOWN")
        disease_id = dossier_data.get("disease_id", "unknown")
        date_str = datetime.now(timezone.utc).strftime("%Y%m%d")
        package_name = f"DOSSIER_{compound_id}_{disease_id}_{date_str}"
        package_dir = Path(output_dir) / package_name

        # Create subdirectories
        for subdir in ["DATA", "CERTIFICATES", "APPENDICES"]:
            (package_dir / subdir).mkdir(parents=True, exist_ok=True)

        # Generate all section documents
        for filename, generator_class in DOCUMENT_GENERATORS:
            generator = generator_class(dossier_data)
            text = generator.generate()
            (package_dir / f"{filename}.txt").write_text(text, encoding="utf-8")

        # Write raw JSON data
        with open(package_dir / "DATA" / "dossier_data.json", "w", encoding="utf-8") as f:
            json.dump(dossier_data, f, indent=2, default=str)

        # Write governance certificates
        governance = exec_data.get("governance", {})
        cert = {
            "governance": governance,
            "generated": date_str,
            "constitutional_compliance": True,
        }
        with open(package_dir / "CERTIFICATES" / "governance.json", "w", encoding="utf-8") as f:
            json.dump(cert, f, indent=2)

        return str(package_dir)
