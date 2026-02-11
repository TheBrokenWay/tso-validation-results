#!/usr/bin/env python3
"""
Nipah Analysis — Normalization (runs after validation gate).

- Standardize units: viral load → copies/mL; CFR already 0–1; cases integer.
- Anonymize PII: location, patient_id, site_id, reporter_id → sha256(salt + value); salt from NIPAH_SALT or ephemeral.
- Output: data/normalized/{strain}_{basename}.json with schema_version and normalization_notes.
"""

import hashlib
import json
import os
import secrets
from pathlib import Path
from typing import Any

# PII fields: anonymize by hashing with salt (deterministic per run) or REDACTED
PII_FIELDS = frozenset({"location", "patient_id", "site_id", "reporter_id", "institution", "contact"})
# Unit conversions: key = field name, value = (target_unit, conversion_factor or callable)
VIRAL_LOAD_UNIT = "copies_per_mL"
NORMALIZATION_SCHEMA_VERSION = "1.0"


def get_pii_salt() -> str:
    """
    Salt for PII hashing: NIPAH_SALT env if set; else NIPAH_SALT_STRICT=true → fail; else random ephemeral.
    Hashing uses sha256(salt + value) for reproducibility when salt is set.
    """
    salt = os.environ.get("NIPAH_SALT", "").strip()
    strict = os.environ.get("NIPAH_SALT_STRICT", "false").lower() in ("1", "true", "yes")
    if salt:
        return salt
    if strict:
        raise RuntimeError(
            "NIPAH_SALT must be set when NIPAH_SALT_STRICT=true. "
            "Set NIPAH_SALT to a fixed value for reproducible anonymization."
        )
    return secrets.token_hex(16)


def _anonymize(value: str | None, *, salt: str) -> str:
    """Deterministic anonymization: sha256(salt + value); no reversible PII."""
    if value is None or (isinstance(value, str) and not value.strip()):
        return "REDACTED"
    raw = value.strip()
    payload = salt + raw
    return "anon_" + hashlib.sha256(payload.encode("utf-8")).hexdigest()[:16]


def _to_copies_per_ml(raw_value: Any, unit_hint: str | None) -> float | None:
    """Convert viral load to copies/mL. Accepts numeric (assume already copies/mL) or log10 copies/mL."""
    if raw_value is None:
        return None
    try:
        v = float(raw_value)
    except (TypeError, ValueError):
        return None
    if unit_hint and "log" in unit_hint.lower():
        import math
        return math.pow(10, v) if v > 0 else None
    return v


def normalize_record(rec: dict[str, Any], *, pii_salt: str) -> dict[str, Any]:
    """One record: standardize units, anonymize PII. Returns new dict."""
    out = dict(rec)
    for key in list(out.keys()):
        if key in PII_FIELDS and out[key] is not None:
            out[key] = _anonymize(str(out[key]), salt=pii_salt)
    if "viral_load" in out and out["viral_load"] is not None:
        unit_hint = out.get("viral_load_unit") or ""
        converted = _to_copies_per_ml(out["viral_load"], unit_hint)
        if converted is not None:
            out["viral_load_copies_per_mL"] = round(converted, 4)
        out.pop("viral_load_unit", None)
    return out


def normalize_file(raw_path: Path, normalized_dir: Path, *, pii_salt: str) -> Path:
    """
    Read raw JSON, normalize units and PII, write to normalized_dir.
    Returns path to written file. Expects validated raw (strain + records).
    """
    normalized_dir = Path(normalized_dir)
    normalized_dir.mkdir(parents=True, exist_ok=True)
    with open(raw_path, encoding="utf-8") as f:
        data = json.load(f)
    strain = data.get("strain", "unknown")
    records = data.get("records", [])
    normalized_records = [normalize_record(r, pii_salt=pii_salt) for r in records]
    out = {
        "strain": strain,
        "records": normalized_records,
        "schema_version": NORMALIZATION_SCHEMA_VERSION,
        "normalization_notes": ["units: viral_load → copies_per_mL; PII anonymized"],
        "source_file": raw_path.name,
    }
    for k in ("source", "year_collected", "ontology_context"):
        if k in data:
            out[k] = data[k]
    out_path = normalized_dir / f"{strain}_{raw_path.stem}.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    return out_path


def normalize_accepted(raw_dir: Path, normalized_dir: Path, accepted: list[dict[str, Any]]) -> list[Path]:
    """
    Normalize every file in accepted (list from validate_all_raw). Returns list of output paths.
    Uses get_pii_salt() for PII hashing (NIPAH_SALT env or ephemeral if not strict).
    """
    raw_dir = Path(raw_dir)
    normalized_dir = Path(normalized_dir)
    pii_salt = get_pii_salt()
    written = []
    for a in accepted:
        path_str = a.get("path")
        if not path_str:
            continue
        p = Path(path_str)
        if not p.is_absolute():
            p = raw_dir / p.name
        if not p.exists():
            continue
        out_path = normalize_file(p, normalized_dir, pii_salt=pii_salt)
        written.append(out_path)
    return written
