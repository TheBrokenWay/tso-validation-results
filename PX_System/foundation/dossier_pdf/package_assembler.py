"""Package assembler -- creates complete DIAMOND dossier folder."""
from __future__ import annotations

import hashlib
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
from .base_generator import extract_compound_id, extract_disease_id

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

    def assemble(
        self,
        dossier_data: Dict[str, Any],
        output_dir: str,
        compound_id: str | None = None,
        disease_id: str | None = None,
    ) -> str:
        """Assemble a full dossier package into output_dir.

        Returns the path to the created package directory.
        """
        exec_data = dossier_data.get("sections", {}).get("executive_summary", {})
        if not compound_id:
            compound_id = extract_compound_id(dossier_data)
        if not disease_id:
            disease_id = extract_disease_id(dossier_data)
        date_str = datetime.now(timezone.utc).strftime("%Y%m%d")
        package_name = f"DOSSIER_{compound_id}_{disease_id}_{date_str}"
        package_dir = Path(output_dir) / package_name

        # Enrich dossier copy so document generators see resolved IDs in headers
        gen_data = {**dossier_data, "compound_id": compound_id, "disease_id": disease_id}

        # Create subdirectories
        for subdir in ["DATA", "CERTIFICATES", "APPENDICES"]:
            (package_dir / subdir).mkdir(parents=True, exist_ok=True)

        # Track files for manifest
        manifest_files: list[Dict[str, str]] = []

        # Generate all section documents
        for filename, generator_class in DOCUMENT_GENERATORS:
            generator = generator_class(gen_data)
            text = generator.generate()
            rel = f"{filename}.txt"
            (package_dir / rel).write_text(text, encoding="utf-8")
            manifest_files.append({
                "file": rel,
                "sha256": hashlib.sha256(text.encode("utf-8")).hexdigest(),
            })

        # Write raw JSON data
        data_json = json.dumps(dossier_data, indent=2, default=str)
        (package_dir / "DATA" / "full_dossier.json").write_text(data_json, encoding="utf-8")
        manifest_files.append({
            "file": "DATA/full_dossier.json",
            "sha256": hashlib.sha256(data_json.encode("utf-8")).hexdigest(),
        })

        # Write governance certificates
        governance = exec_data.get("governance", {})
        cert = {
            "governance": governance,
            "generated": date_str,
            "constitutional_compliance": True,
        }
        cert_json = json.dumps(cert, indent=2)
        (package_dir / "CERTIFICATES" / "governance.json").write_text(cert_json, encoding="utf-8")
        manifest_files.append({
            "file": "CERTIFICATES/governance.json",
            "sha256": hashlib.sha256(cert_json.encode("utf-8")).hexdigest(),
        })

        # Write Zeus gate approval certificate
        fin = dossier_data.get("finalization", {})
        zeus_cert = {
            "zeus_verdict": fin.get("zeus_verdict", {}),
            "authorization_chain": fin.get("authorization_chain", dossier_data.get("authorization_chain", [])),
            "constitutional_seal": fin.get("constitutional_seal", ""),
            "compound_id": compound_id,
            "disease_id": disease_id,
            "generated": date_str,
        }
        zeus_json = json.dumps(zeus_cert, indent=2, default=str)
        (package_dir / "CERTIFICATES" / "zeus_gate_approval.json").write_text(zeus_json, encoding="utf-8")
        manifest_files.append({
            "file": "CERTIFICATES/zeus_gate_approval.json",
            "sha256": hashlib.sha256(zeus_json.encode("utf-8")).hexdigest(),
        })

        # Write dossier manifest
        manifest = {
            "package_name": package_name,
            "compound_id": compound_id,
            "disease_id": disease_id,
            "tier": dossier_data.get("finalization", {}).get("tier", "DIAMOND"),
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "files": manifest_files,
        }
        (package_dir / "dossier_manifest.json").write_text(
            json.dumps(manifest, indent=2), encoding="utf-8"
        )

        return str(package_dir)
