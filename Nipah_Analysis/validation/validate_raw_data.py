"""
Raw data validation gate — checksum, file type, schema, strain, CFR envelope, temporal, forbidden family.

If ANY check fails → reject file; pipeline must halt if any rejection.
"""

import hashlib
import json
from pathlib import Path
from typing import Any

ALLOWED_EXTENSIONS = {".json", ".csv"}
BLOCKED_EXTENSIONS = {".exe", ".bat", ".sh", ".py", ".md"}
STRAIN_ENUM = frozenset({"NiV_Malaysia", "NiV_Bangladesh"})
# Ontology: CFR bounds and first outbreak year by strain (must match manifest)
CFR_BOUNDS = {"NiV_Malaysia": (0.35, 0.45), "NiV_Bangladesh": (0.70, 1.00)}
FIRST_OUTBREAK_YEAR = {"NiV_Malaysia": 1998, "NiV_Bangladesh": 2001}
# Forbidden families: prevents cross-family data entering Nipah pipeline
FORBIDDEN_FAMILIES = frozenset({"Filoviridae", "Coronaviridae", "Flaviviridae", "Orthomyxoviridae"})
SCHEMA_VERSION = "2.0"
ONTOLOGY_VERSION = "1.0"


class RawValidationError(Exception):
    """Raised when a raw file fails validation."""

    def __init__(self, message: str, file_path: str | None = None, code: str | None = None):
        self.file_path = file_path
        self.code = code
        super().__init__(message)


def discover_raw_files(raw_dir: Path) -> list[Path]:
    """Discover raw data files; exclude .gitkeep and blocked types. Returns list of Paths."""
    raw_dir = Path(raw_dir)
    if not raw_dir.exists():
        return []
    out = []
    for p in raw_dir.iterdir():
        if not p.is_file():
            continue
        if p.name == ".gitkeep":
            continue
        suf = p.suffix.lower()
        if suf in BLOCKED_EXTENSIONS:
            continue
        if suf not in ALLOWED_EXTENSIONS:
            continue
        out.append(p)
    return sorted(out)


def validate_file_checksum(path: Path) -> str:
    """Compute SHA256 of file. Returns hex digest. Does not raise."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _load_json_for_validation(path: Path) -> dict[str, Any]:
    """Load JSON; raise RawValidationError if invalid or not object."""
    try:
        with open(path, encoding="utf-8") as f:
            data = json.load(f)
    except json.JSONDecodeError as e:
        raise RawValidationError(f"Invalid JSON: {e}", file_path=str(path), code="INVALID_JSON")
    if not isinstance(data, dict):
        raise RawValidationError("Root must be a JSON object", file_path=str(path), code="SCHEMA_ROOT")
    return data


def validate_forbidden_family(data: dict[str, Any], path: Path) -> None:
    """
    If data declares pathogen family (family or pathogen_family), it must not be in forbidden_families.
    Prevents Ebola/SARS/influenza etc. from entering Nipah pipeline.
    """
    for key in ("family", "pathogen_family"):
        val = data.get(key)
        if val is None or not isinstance(val, str):
            continue
        fam = val.strip()
        if fam in FORBIDDEN_FAMILIES:
            raise RawValidationError(
                f"Forbidden family {fam!r}; Nipah pipeline accepts Paramyxoviridae/Henipavirus only.",
                file_path=str(path),
                code="FORBIDDEN_FAMILY",
            )


def validate_file_schema(data: dict[str, Any], path: Path) -> None:
    """
    Validate structure: required strain (enum), records array, each record outbreak_year + cfr.
    Raises RawValidationError on schema violation.
    """
    if "strain" not in data:
        raise RawValidationError("Missing required field: strain", file_path=str(path), code="SCHEMA_VIOLATION")
    strain = data.get("strain")
    if strain not in STRAIN_ENUM:
        raise RawValidationError(
            f"strain must be one of {sorted(STRAIN_ENUM)}; got {strain!r}",
            file_path=str(path),
            code="SCHEMA_VIOLATION",
        )
    if "records" not in data:
        raise RawValidationError("Missing required field: records", file_path=str(path), code="SCHEMA_VIOLATION")
    records = data.get("records")
    if not isinstance(records, list) or len(records) == 0:
        raise RawValidationError("records must be a non-empty array", file_path=str(path), code="SCHEMA_VIOLATION")
    for i, rec in enumerate(records):
        if not isinstance(rec, dict):
            raise RawValidationError(f"records[{i}] must be an object", file_path=str(path), code="SCHEMA_VIOLATION")
        if "outbreak_year" not in rec or "cfr" not in rec:
            raise RawValidationError(
                f"records[{i}] must have outbreak_year and cfr",
                file_path=str(path),
                code="SCHEMA_VIOLATION",
            )
        try:
            y = int(rec["outbreak_year"])
            if y < 1998 or y > 2030:
                raise RawValidationError(
                    f"records[{i}].outbreak_year must be 1998–2030; got {y}",
                    file_path=str(path),
                    code="SCHEMA_VIOLATION",
                )
        except (TypeError, ValueError):
            raise RawValidationError(
                f"records[{i}].outbreak_year must be integer",
                file_path=str(path),
                code="SCHEMA_VIOLATION",
            )
        try:
            c = float(rec["cfr"])
            if c < 0 or c > 1:
                raise RawValidationError(
                    f"records[{i}].cfr must be in [0,1]; got {c}",
                    file_path=str(path),
                    code="SCHEMA_VIOLATION",
                )
        except (TypeError, ValueError):
            raise RawValidationError(
                f"records[{i}].cfr must be number",
                file_path=str(path),
                code="SCHEMA_VIOLATION",
            )


def classify_strain(data: dict[str, Any]) -> str | None:
    """Return strain from data (NiV_Malaysia or NiV_Bangladesh) or None if ambiguous/missing."""
    s = data.get("strain")
    if s in STRAIN_ENUM:
        return s
    return None


def _ontology_from_manifest(manifest: dict[str, Any]) -> tuple[dict[str, tuple[float, float]], dict[str, int]]:
    """Derive CFR bounds and first_outbreak_year from manifest (single source of truth)."""
    cfr_bounds: dict[str, tuple[float, float]] = {}
    first_outbreak: dict[str, int] = {}
    strains = (manifest.get("pathogen") or {}).get("strains") or {}
    for sid, sdef in strains.items():
        if not isinstance(sdef, dict):
            continue
        cfr = sdef.get("cfr_range")
        if isinstance(cfr, list) and len(cfr) == 2:
            cfr_bounds[sid] = (float(cfr[0]), float(cfr[1]))
        yr = sdef.get("first_outbreak_year")
        if isinstance(yr, int):
            first_outbreak[sid] = yr
    return cfr_bounds or CFR_BOUNDS, first_outbreak or FIRST_OUTBREAK_YEAR


def validate_cfr_envelope(
    data: dict[str, Any],
    path: Path,
    *,
    cfr_bounds: dict[str, tuple[float, float]] | None = None,
) -> None:
    """
    All record CFRs must fall within strain CFR envelope (ontology).
    Raises RawValidationError if any CFR outside envelope (e.g. 18% for NiV).
    Use cfr_bounds from manifest when called from validate_all_raw.
    """
    bounds = cfr_bounds if cfr_bounds is not None else CFR_BOUNDS
    strain = data.get("strain")
    if strain not in bounds:
        return
    lo, hi = bounds[strain]
    for i, rec in enumerate(data.get("records") or []):
        cfr = rec.get("cfr")
        if cfr is None:
            continue
        try:
            c = float(cfr)
            if c < lo or c > hi:
                raise RawValidationError(
                    f"records[{i}].cfr={c} outside strain envelope [{lo},{hi}] for {strain}",
                    file_path=str(path),
                    code="CFR_ANOMALY",
                )
        except (TypeError, ValueError):
            pass


def validate_temporal_consistency(
    data: dict[str, Any],
    path: Path,
    *,
    first_outbreak_year: dict[str, int] | None = None,
) -> None:
    """
    Outbreak years must be >= first outbreak year for that strain (historical consistency).
    E.g. Bangladesh first 2001 → 1998 labeled Bangladesh fails.
    Use first_outbreak_year from manifest when called from validate_all_raw.
    """
    foy = first_outbreak_year if first_outbreak_year is not None else FIRST_OUTBREAK_YEAR
    strain = data.get("strain")
    if strain not in foy:
        return
    first_year = foy[strain]
    for i, rec in enumerate(data.get("records") or []):
        y = rec.get("outbreak_year")
        if y is None:
            continue
        try:
            yr = int(y)
            if yr < first_year:
                raise RawValidationError(
                    f"records[{i}].outbreak_year={yr} before first outbreak year {first_year} for {strain}",
                    file_path=str(path),
                    code="TEMPORAL_DRIFT",
                )
        except (TypeError, ValueError):
            pass


def _cfr_drift_summary(accepted: list[dict[str, Any]]) -> dict[str, Any]:
    """Build CFR distribution by strain for drift detection."""
    by_strain: dict[str, list[float]] = {}
    for a in accepted:
        s = a.get("strain")
        vals = a.get("cfr_values") or []
        if s not in by_strain:
            by_strain[s] = []
        by_strain[s].extend(vals)
    out = {}
    for s, vals in by_strain.items():
        if not vals:
            out[s] = {"n": 0, "mean": None, "min": None, "max": None}
        else:
            import statistics
            out[s] = {
                "n": len(vals),
                "mean": round(statistics.mean(vals), 4),
                "min": round(min(vals), 4),
                "max": round(max(vals), 4),
            }
    return out


def validate_all_raw(
    raw_dir: Path,
    manifest: dict[str, Any],
    *,
    fail_on_no_files: bool = True,
    allow_multi_strain: bool = False,
    require_strain_coverage: bool = True,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    """
    Discover raw files; for each: checksum, schema, strain, CFR envelope, temporal.
    Returns (accepted_list, rejected_list, drift_summary). accepted_list items include path, checksum, strain, record_count, cfr_values.
    If fail_on_no_files and no files discovered → raise RawValidationError.
    If any file rejected → raise RawValidationError with first rejection reason.
    If not allow_multi_strain and >1 strain present → raise RawValidationError (strain isolation).
    If require_strain_coverage (typically True only when allow_multi_strain): every manifest strain must have >=1 accepted file.
    When allow_multi_strain is False (single-strain mode), require_strain_coverage should be False: manifest defines what CAN be analyzed, not what MUST be in this batch.
    """
    files = discover_raw_files(raw_dir)
    if fail_on_no_files and len(files) == 0:
        raise RawValidationError(
            "No raw data files found; pipeline cannot proceed. Add JSON/CSV to data/raw/.",
            code="NO_RAW_FILES",
        )

    cfr_bounds, first_outbreak_year = _ontology_from_manifest(manifest)
    accepted = []
    rejected = []
    for path in files:
        if path.suffix.lower() != ".json":
            # CSV: minimal validation for now — require strain column or reject
            rejected.append({"path": str(path), "reason": "CSV strain validation not implemented; use JSON", "code": "CSV_STRAIN"})
            continue

        try:
            data = _load_json_for_validation(path)
            validate_forbidden_family(data, path)
            validate_file_schema(data, path)
            strain = classify_strain(data)
            if strain is None:
                rejected.append({"path": str(path), "reason": "Unclassifiable or ambiguous strain", "code": "UNCLASSIFIABLE_STRAIN"})
                continue
            validate_cfr_envelope(data, path, cfr_bounds=cfr_bounds)
            validate_temporal_consistency(data, path, first_outbreak_year=first_outbreak_year)
            checksum = validate_file_checksum(path)
            records = data.get("records", [])
            cfr_values = [float(r["cfr"]) for r in records if r.get("cfr") is not None]
            accepted.append({
                "path": str(path),
                "checksum": checksum,
                "strain": strain,
                "record_count": len(records),
                "cfr_values": cfr_values,
            })
        except RawValidationError as e:
            rejected.append({"path": str(path), "reason": str(e), "code": e.code or "VALIDATION_FAILED"})

    if rejected:
        first = rejected[0]
        raise RawValidationError(
            f"Validation failed: {first['path']} — {first['reason']}",
            file_path=first["path"],
            code=first.get("code", "VALIDATION_FAILED"),
        )

    strains_seen = {a["strain"] for a in accepted}
    if not allow_multi_strain and len(strains_seen) > 1:
        raise RawValidationError(
            f"Strain isolation: multiple strains present {sorted(strains_seen)}. Use --allow-multi-strain to allow.",
            code="STRAIN_ISOLATION",
        )

    if require_strain_coverage:
        defined_strains = set((manifest.get("pathogen") or {}).get("strains") or {})
        missing = defined_strains - strains_seen
        if defined_strains and missing:
            raise RawValidationError(
                f"Strain coverage: no valid raw files for strain(s) {sorted(missing)}. Include data for all manifest strains.",
                code="STRAIN_COVERAGE",
            )
    drift_summary = _cfr_drift_summary(accepted)
    return accepted, rejected, drift_summary


def validate_file_content(
    raw_data_path: str | Path,
    manifest: dict[str, Any],
) -> tuple[bool, dict[str, Any] | str]:
    """
    Validate a single raw Nipah data file against manifest-derived ontology.
    Returns (True, meta_dict) on success; (False, error_message_str) on failure.
    meta_dict includes: strain, checksum, record_count, path.
    """
    path = Path(raw_data_path)
    if not path.exists():
        return False, f"File not found: {path}"
    if path.suffix.lower() != ".json":
        return False, "Only JSON raw files are supported"
    try:
        data = _load_json_for_validation(path)
    except RawValidationError as e:
        return False, str(e)
    cfr_bounds, first_outbreak_year = _ontology_from_manifest(manifest)
    try:
        validate_forbidden_family(data, path)
        validate_file_schema(data, path)
        strain = classify_strain(data)
        if strain is None:
            return False, "Unclassifiable or ambiguous strain"
        validate_cfr_envelope(data, path, cfr_bounds=cfr_bounds)
        validate_temporal_consistency(data, path, first_outbreak_year=first_outbreak_year)
    except RawValidationError as e:
        return False, str(e)
    checksum = validate_file_checksum(path)
    records = data.get("records", [])
    return True, {
        "strain": strain,
        "checksum": checksum,
        "record_count": len(records),
        "path": str(path),
    }
