"""
DIRECTIVE: ENFORCE_FINAL_PRODUCT_FILING
Dossier + CommercialAssets filing rules — deterministic, enforceable, drift-proof.

Authority: DOSSIER_COMMERCIALASSETS_FILING_RULES.md
- Only assets that PASS ALL CHECKS may enter Dossier_Final.
- All commercial-grade outputs must be routed to the correct tier (Diamond/Gold/Silver/Bronze).
- Violations trigger FAIL-CLOSED and must be logged.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Tuple

# Canonical tiers under CommercialAssets (biopath routing)
COMMERCIAL_TIERS = frozenset({"Diamond", "Gold", "Silver", "Bronze"})

# Company-specific filing zones for discontinued/withdrawn/failed assets (supplements tiers)
COMPANY_FOLDERS = frozenset({
    "Pfizer", "Janssen_JnJ", "Eli_Lilly",
    "Merck", "Novartis", "Roche", "GSK", "Sanofi", "AstraZeneca",
    "Takeda", "Bayer", "BMS", "Amgen",
})
SPONSOR_TO_COMPANY = {
    "pfizer": "Pfizer",
    "janssen": "Janssen_JnJ",
    "johnson & johnson": "Janssen_JnJ",
    "jnj": "Janssen_JnJ",
    "eli lilly": "Eli_Lilly",
    "lilly": "Eli_Lilly",
    "merck": "Merck",
    "novartis": "Novartis",
    "roche": "Roche",
    "gsk": "GSK",
    "glaxosmithkline": "GSK",
    "sanofi": "Sanofi",
    "astrazeneca": "AstraZeneca",
    "astra zeneca": "AstraZeneca",
    "takeda": "Takeda",
    "bayer": "Bayer",
    "bms": "BMS",
    "bristol-myers squibb": "BMS",
    "bristol myers squibb": "BMS",
    "amgen": "Amgen",
}
DISCONTINUED_STATUSES = frozenset({"discontinued", "withdrawn", "terminated", "failed_phase"})

COMMERCIAL_SUBFOLDERS = frozenset({
    "Diamond", "Gold", "Silver", "Bronze",
    "Dossier_Final", "Executive_Summary", "Audit_Trail", "Learning_Material",
    "Pfizer", "Janssen_JnJ", "Eli_Lilly",
    "Merck", "Novartis", "Roche", "GSK", "Sanofi", "AstraZeneca",
    "Takeda", "Bayer", "BMS", "Amgen",
})

# Required check keys that must be present and True for Dossier_Final admission
REQUIRED_CHECKS = frozenset({
    "scientific_filters_pass",
    "admet_exposure_pass",
    "legal_fto_pass",
    "commercial_viability_pass",
    "provenance_required",
    "timestamp_required",
    "reproducibility_required",
})

# Canonical structure under <Tier>/<AssetID>/
CANONICAL_ASSET_SUBDIRS = frozenset({"evidence", "audit"})
CANONICAL_ASSET_FILES = frozenset({"dossier.json", "summary.md"})


def may_enter_dossier_final(asset_metadata: Dict[str, Any]) -> bool:
    """
    True only if the asset has PASSED ALL CHECKS. Fail-closed: any missing or False check → False.
    Required checks: scientific_filters_pass, admet_exposure_pass, legal_fto_pass,
    commercial_viability_pass, provenance_required, timestamp_required, reproducibility_required.
    """
    if not asset_metadata or not isinstance(asset_metadata, dict):
        return False
    checks = asset_metadata.get("checks") or asset_metadata.get("validation") or {}
    if isinstance(checks, dict):
        for key in REQUIRED_CHECKS:
            if not checks.get(key):
                return False
        return True
    return False


def get_filing_location(asset_metadata: Dict[str, Any]) -> Tuple[str, str]:
    """
    Data Engine directive: if sponsor is Pfizer/Janssen/J&J/Eli Lilly AND status is
    discontinued/withdrawn/terminated/failed_phase → file under company folder.
    Else → file under tier by commercial score.
    Returns ("company", "Pfizer") or ("tier", "Gold").
    """
    if not asset_metadata or not isinstance(asset_metadata, dict):
        return "tier", "Bronze"
    status = (asset_metadata.get("status") or "").strip().lower()
    sponsor = (asset_metadata.get("sponsor") or "").strip().lower()
    if status in DISCONTINUED_STATUSES and sponsor:
        company = SPONSOR_TO_COMPANY.get(sponsor)
        if company:
            return "company", company
    return "tier", get_tier_for_asset(asset_metadata)


def get_tier_for_asset(asset_metadata: Dict[str, Any]) -> str:
    """
    Return the canonical tier (Diamond, Gold, Silver, Bronze) for biopath routing.
    Determined by final commercial score and governance metadata. Default: Bronze.
    """
    if not asset_metadata or not isinstance(asset_metadata, dict):
        return "Bronze"
    tier = (asset_metadata.get("tier") or asset_metadata.get("commercial_tier") or "").strip()
    if tier in COMMERCIAL_TIERS:
        return tier
    score = asset_metadata.get("commercial_score")
    if score is not None:
        try:
            s = float(score)
            if s >= 0.9:
                return "Diamond"
            if s >= 0.75:
                return "Gold"
            if s >= 0.5:
                return "Silver"
        except (TypeError, ValueError):
            pass
    return "Bronze"


def validate_filing_path(
    path: Path | str, asset_metadata: Dict[str, Any] | None
) -> Tuple[bool, str]:
    """
    Validate that a write path complies with the filing rules.
    Returns (ok: bool, reason: str). Fail-closed: violations return (False, reason).
    """
    path = Path(path).resolve() if path else Path()
    path_str = str(path)
    parts = path.parts

    # Must be under CommercialAssets
    if "CommercialAssets" not in parts:
        return False, "Path is not under PX_Warehouse/CommercialAssets"

    # Dossier_Final: only allowed if asset passed all checks
    if "Dossier_Final" in parts:
        if not asset_metadata:
            return False, "Dossier_Final requires asset_metadata for validation"
        if not may_enter_dossier_final(asset_metadata):
            return False, "Asset has not passed all required checks; may not enter Dossier_Final"
        return True, "OK"

    # Tier folder (Diamond, Gold, Silver, Bronze): must match assigned tier
    for tier in COMMERCIAL_TIERS:
        if tier in parts:
            if asset_metadata:
                loc_type, loc_name = get_filing_location(asset_metadata)
                if loc_type == "company":
                    return False, f"Asset is assigned to company folder {loc_name}, not tier {tier}"
                if loc_name != tier:
                    return False, f"Tier mismatch: path has {tier} but asset is assigned {loc_name}"
            return True, "OK"

    # Company folder (Pfizer, Janssen_JnJ, Eli_Lilly): only for discontinued/withdrawn/failed + matching sponsor
    for company in COMPANY_FOLDERS:
        if company in parts:
            if asset_metadata:
                loc_type, loc_name = get_filing_location(asset_metadata)
                if loc_type != "company" or loc_name != company:
                    return False, f"Company folder {company} requires discontinued/withdrawn/failed asset from matching sponsor; got {loc_type}={loc_name}"
            return True, "OK"

    # Executive_Summary, Audit_Trail: allowed (routing by pipeline)
    if "Executive_Summary" in parts or "Audit_Trail" in parts:
        return True, "OK"

    # Learning_Material: allowed
    if "Learning_Material" in parts:
        return True, "OK"

    # Unknown CommercialAssets subfolder
    try:
        idx = parts.index("CommercialAssets")
        if idx + 1 < len(parts):
            sub = parts[idx + 1]
            if sub not in COMMERCIAL_SUBFOLDERS:
                return False, f"Non-canonical CommercialAssets subfolder: {sub}"
    except ValueError:
        pass
    return True, "OK"


def canonical_asset_path(
    repo_root: Path | str, tier_or_company: str, asset_id: str
) -> Path:
    """
    Return the canonical path for an asset: CommercialAssets/<TierOrCompany>/<AssetID>/.
    tier_or_company may be a tier (Diamond/Gold/Silver/Bronze) or company folder (Pfizer/Janssen_JnJ/Eli_Lilly).
    """
    root = Path(repo_root)
    if tier_or_company in COMPANY_FOLDERS:
        folder = tier_or_company
    elif tier_or_company in COMMERCIAL_TIERS:
        folder = tier_or_company
    else:
        folder = "Bronze"
    safe_id = "".join(c for c in asset_id if c.isalnum() or c in "-_") or "unknown"
    return root / "PX_Warehouse" / "CommercialAssets" / folder / safe_id


def log_violation(message: str, path: str | None = None, asset_id: str | None = None) -> None:
    """Log a filing-rule violation (caller may also write to PX_LOGS or sovereign log)."""
    parts = [f"[FILING_VIOLATION] {message}"]
    if path:
        parts.append(f" path={path}")
    if asset_id:
        parts.append(f" asset_id={asset_id}")
    print("".join(parts))


__all__ = [
    "COMMERCIAL_TIERS",
    "COMPANY_FOLDERS",
    "COMMERCIAL_SUBFOLDERS",
    "REQUIRED_CHECKS",
    "SPONSOR_TO_COMPANY",
    "DISCONTINUED_STATUSES",
    "may_enter_dossier_final",
    "get_filing_location",
    "get_tier_for_asset",
    "validate_filing_path",
    "canonical_asset_path",
    "log_violation",
]
