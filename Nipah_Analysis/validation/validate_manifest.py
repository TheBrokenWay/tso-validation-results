"""
Manifest schema validation — load and validate nipah_manifest.json.

Ensures strain-split architecture (NiV_Malaysia, NiV_Bangladesh) and required keys.
"""

import json
from pathlib import Path
from typing import Any

REQUIRED_TOP_KEYS = {"analysis_id", "pathogen", "validation"}
REQUIRED_PATHOGEN_KEYS = {"species", "family", "genus", "strains"}
REQUIRED_STRAIN_KEYS = {"cfr_range", "primary_vectors", "human_to_human", "spillover_pattern"}


class ManifestValidationError(Exception):
    """Raised when manifest fails schema or ontology check."""

    def __init__(self, message: str, field: str | None = None):
        self.field = field
        super().__init__(message)


def validate_manifest_schema(manifest: dict[str, Any]) -> None:
    """
    Validate manifest has strain-split pathogen and required structure.
    Raises ManifestValidationError on failure.
    """
    if not isinstance(manifest, dict):
        raise ManifestValidationError("Manifest must be a JSON object", field="manifest")

    for key in REQUIRED_TOP_KEYS:
        if key not in manifest:
            raise ManifestValidationError(f"Missing required key: {key}", field=key)

    pathogen = manifest.get("pathogen")
    if not isinstance(pathogen, dict):
        raise ManifestValidationError("pathogen must be an object", field="pathogen")

    for key in REQUIRED_PATHOGEN_KEYS:
        if key not in pathogen:
            raise ManifestValidationError(f"pathogen missing key: {key}", field=f"pathogen.{key}")

    strains = pathogen.get("strains")
    if not isinstance(strains, dict):
        raise ManifestValidationError("pathogen.strains must be an object", field="pathogen.strains")

    if "NiV_Malaysia" not in strains or "NiV_Bangladesh" not in strains:
        raise ManifestValidationError(
            "pathogen.strains must contain NiV_Malaysia and NiV_Bangladesh",
            field="pathogen.strains",
        )

    for strain_id, strain_def in strains.items():
        if not isinstance(strain_def, dict):
            raise ManifestValidationError(f"strain {strain_id} must be an object", field=f"strains.{strain_id}")
        for key in REQUIRED_STRAIN_KEYS:
            if key not in strain_def:
                raise ManifestValidationError(f"strain {strain_id} missing key: {key}", field=f"strains.{strain_id}.{key}")
        cfr = strain_def.get("cfr_range")
        if not isinstance(cfr, list) or len(cfr) != 2 or not (0 <= cfr[0] <= 1 and 0 <= cfr[1] <= 1):
            raise ManifestValidationError(
                f"strain {strain_id} cfr_range must be [min, max] in [0,1]",
                field=f"strains.{strain_id}.cfr_range",
            )


def load_and_validate_manifest(manifest_path: Path) -> dict[str, Any]:
    """Load manifest JSON and validate; return manifest dict. Raises on failure."""
    manifest_path = Path(manifest_path)
    if not manifest_path.exists():
        raise ManifestValidationError(f"Manifest not found: {manifest_path}", field="path")
    try:
        with open(manifest_path, encoding="utf-8") as f:
            manifest = json.load(f)
    except json.JSONDecodeError as e:
        raise ManifestValidationError(f"Manifest JSON invalid: {e}", field="json")
    validate_manifest_schema(manifest)
    return manifest
