#!/usr/bin/env python3
"""
Olympus Bridge v2: NCT raw API → PX_Warehouse/Learning_Material.

Input: Scan Olympus_Research for raw_api_results.json.
Mapping: adverseEventsModule → safety_margin (sm); outcomeMeasuresModule → response_ratio (rr).
OPE: Flag missing binding_affinity_nM (ba) and oral_bioavailability_percent (bio) as PENDING_OPE_CALCULATION.
Output: JSONs in PX_Warehouse/Learning_Material with audit_note for OPE.
Rules: No stochastic libs, strict type hints, Zero Interpolation. TraceabilityError if toxicity_index missing.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, TypedDict

# -----------------------------------------------------------------------------
# Constants & Paths (canonical: PX_Warehouse/Learning_Material)
# -----------------------------------------------------------------------------
OLYMPUS_RESEARCH = Path("Olympus_Research")


def _learning_material_dir() -> Path:
    try:
        from PX_Warehouse.warehouse_layout import get_learning_material_dir
        repo = Path(__file__).resolve().parents[2]  # Acquisition -> PX_Executive -> foundation
        return get_learning_material_dir(repo)
    except Exception:
        return Path("PX_Warehouse") / "Learning_Material"


LEARNING_MATERIAL_DIR = _learning_material_dir()
PENDING_OPE = "PENDING_OPE_CALCULATION"
BRIDGE_VERSION = "2.0.0"


class TraceabilityError(Exception):
    """Raised when a critical value is missing from raw data and cannot be calculated (Zero Interpolation)."""


# -----------------------------------------------------------------------------
# Typed structures (.cursorrules: strict type hinting)
# -----------------------------------------------------------------------------
class DoseOptimizationOut(TypedDict):
    safety_margin: float


class ResponderRateOut(TypedDict):
    response_ratio: float


class VirtualEfficacyOut(TypedDict):
    responder_rate: ResponderRateOut


class OpeOut(TypedDict):
    binding_affinity_nM: str | float


class AbsorptionOut(TypedDict):
    oral_bioavailability_percent: str | float


class AdmetOut(TypedDict):
    absorption: AbsorptionOut
    toxicity: dict[str, Any] | None


class MetadataOut(TypedDict, total=False):
    audit_note: str
    nct_id: str
    source: str
    bridge_version: str


class DerivedMetricOut(TypedDict):
    value: float
    confidence: str


class DatasetFingerprintOut(TypedDict):
    source: str
    api_version: str
    retrieved_at: str
    hash: str


class LearningMaterialDoc(TypedDict, total=False):
    dose_optimization: DoseOptimizationOut
    virtual_efficacy: VirtualEfficacyOut
    ope: OpeOut
    admet: AdmetOut
    metadata: MetadataOut
    derived_metrics: dict[str, DerivedMetricOut]
    dataset_fingerprint: DatasetFingerprintOut
    protocolSection: dict[str, Any]
    resultsSection: dict[str, Any]


def _derive_safety_margin(raw: dict[str, Any]) -> float:
    """
    Map adverseEventsModule to safety_margin (sm). Zero interpolation.
    Uses eventGroups: seriousNumAffected, seriousNumAtRisk (and deaths).
    sm = 100 - 100 * (total_affected / total_at_risk) when total_at_risk > 0, else 0.
    """
    results: dict[str, Any] = raw.get("resultsSection") or {}
    aem: dict[str, Any] = results.get("adverseEventsModule") or {}
    groups: list[dict[str, Any]] = aem.get("eventGroups") or []
    total_affected = 0
    total_at_risk = 0
    for g in groups:
        serious_aff = g.get("seriousNumAffected")
        serious_risk = g.get("seriousNumAtRisk")
        if serious_risk is not None and serious_risk > 0:
            total_at_risk += int(serious_risk)
            total_affected += int(serious_aff or 0)
    if total_at_risk <= 0:
        return 0.0
    rate = total_affected / total_at_risk
    return float(100.0 - 100.0 * rate)


def _derive_response_ratio(raw: dict[str, Any]) -> float:
    """
    Map outcomeMeasuresModule (and primary outcome counts) to response_ratio (rr). Zero interpolation.
    Uses primary outcome denoms/counts and enrollment; no p-value interpolation.
    """
    results: dict[str, Any] = raw.get("resultsSection") or {}
    omm: dict[str, Any] = results.get("outcomeMeasuresModule") or {}
    measures: list[dict[str, Any]] = omm.get("outcomeMeasures") or []
    enrollment = 0
    proto: dict[str, Any] = raw.get("protocolSection") or {}
    design: dict[str, Any] = proto.get("designModule") or {}
    enroll_info: dict[str, Any] = design.get("enrollmentInfo") or {}
    enroll_count = enroll_info.get("count")
    if enroll_count is not None:
        enrollment = int(enroll_count)
    if enrollment <= 0:
        return 0.0
    total_responders = 0
    for m in measures:
        if m.get("type") != "PRIMARY":
            continue
        denoms: list[dict[str, Any]] = m.get("denoms") or []
        for d in denoms:
            counts: list[dict[str, Any]] = d.get("counts") or []
            for c in counts:
                val = c.get("value")
                if val is not None:
                    try:
                        total_responders += int(val)
                    except (ValueError, TypeError):
                        pass
    return float(total_responders) / float(enrollment)


def _check_toxicity_index(raw: dict[str, Any]) -> None:
    """
    toxicity_index is critical and not present in ClinicalTrials.gov API v2.
    Cannot be calculated in the bridge (no interpolation). Raise TraceabilityError.
    """
    results: dict[str, Any] = raw.get("resultsSection") or {}
    # ADMET/toxicology not in API v2; toxicity_index would come from OPE/ADMET
    if results.get("toxicity_index") is not None:
        return
    proto: dict[str, Any] = raw.get("protocolSection") or {}
    if proto.get("toxicity_index") is not None:
        return
    raise TraceabilityError(
        "toxicity_index missing from raw data and cannot be calculated; OPE/ADMET must provide."
    )


def _build_audit_note(
    sm: float,
    rr: float,
    missing_properties: list[str],
    toxicity_missing: bool = False,
) -> str:
    """
    Format audit_note for predator_x_v3_stratification fallback: SM=, RR=, BA=, and MISSING list.
    OPE uses MISSING to know what to calculate first.
    """
    parts: list[str] = [
        f"SM={sm}",
        f"RR={rr}",
        f"BA={PENDING_OPE}",
        f"MISSING: {', '.join(missing_properties)}. OPE: calculate these first.",
    ]
    if toxicity_missing:
        parts.append("toxicity_index missing; OPE/ADMET must provide.")
    return ", ".join(parts)


def _build_derived_metrics(sm: float, rr: float) -> dict[str, dict[str, Any]]:
    """Confidence-annotated derived metrics for interpretability and legal defensibility."""
    return {
        "safety_margin": {
            "value": sm,
            "confidence": "clinical_adverse_events_proxy",
        },
        "response_ratio": {
            "value": rr,
            "confidence": "primary_outcome_denominator_proxy",
        },
    }


def _build_dataset_fingerprint(raw_path: Path, raw_bytes: bytes) -> DatasetFingerprintOut:
    """Provenance fingerprint for government/defense/legal audit and compliance."""
    retrieved_at: str
    try:
        mtime = raw_path.stat().st_mtime
        retrieved_at = datetime.fromtimestamp(mtime, tz=timezone.utc).isoformat()
    except OSError:
        retrieved_at = datetime.now(timezone.utc).isoformat()
    content_hash = hashlib.sha256(raw_bytes).hexdigest()
    return {
        "source": "ClinicalTrials.gov",
        "api_version": "v2",
        "retrieved_at": retrieved_at,
        "hash": content_hash,
    }


def _build_doc(
    raw: dict[str, Any],
    raw_bytes: bytes,
    raw_path: Path,
    nct_id: str,
    sm: float,
    rr: float,
    audit_note: str,
) -> LearningMaterialDoc:
    """Build learning material doc with versioning, derived_metrics, and dataset_fingerprint."""
    return {
        "dose_optimization": {"safety_margin": sm},
        "virtual_efficacy": {"responder_rate": {"response_ratio": rr}},
        "ope": {"binding_affinity_nM": PENDING_OPE},
        "admet": {
            "absorption": {"oral_bioavailability_percent": PENDING_OPE},
            "toxicity": None,
        },
        "metadata": {
            "audit_note": audit_note,
            "nct_id": nct_id,
            "source": "olympus_bridge_v2",
            "bridge_version": BRIDGE_VERSION,
        },
        "derived_metrics": _build_derived_metrics(sm, rr),
        "dataset_fingerprint": _build_dataset_fingerprint(raw_path, raw_bytes),
        "protocolSection": raw.get("protocolSection") or {},
        "resultsSection": raw.get("resultsSection") or {},
    }


def process_one_raw(raw_path: Path, nct_id: str) -> LearningMaterialDoc:
    """
    Read raw_api_results.json, derive sm/rr, flag ba/bio as PENDING_OPE_CALCULATION,
    enforce toxicity_index (raise TraceabilityError if missing). Return learning material doc.
    """
    raw_bytes = raw_path.read_bytes()
    raw: dict[str, Any] = json.loads(raw_bytes.decode("utf-8"))

    _check_toxicity_index(raw)

    sm = _derive_safety_margin(raw)
    rr = _derive_response_ratio(raw)

    missing: list[str] = [
        "binding_affinity_nM (ba)",
        "oral_bioavailability_percent (bio)",
    ]
    audit_note = _build_audit_note(sm, rr, missing, toxicity_missing=False)

    return _build_doc(raw, raw_bytes, raw_path, nct_id, sm, rr, audit_note)


def run_bridge(
    olympus_root: Path = OLYMPUS_RESEARCH,
    learning_dir: Path = LEARNING_MATERIAL_DIR,
    strict_toxicity: bool = True,
) -> int:
    """
    Scan olympus_root for raw_api_results.json; write Learning_Material JSONs.
    If strict_toxicity=True (default), raises TraceabilityError when toxicity_index is missing.
    Returns count of files written.
    """
    if not olympus_root.exists():
        print(f"Olympus_Research not found: {olympus_root}")
        return 0
    learning_dir.mkdir(parents=True, exist_ok=True)
    count = 0
    for raw_path in olympus_root.rglob("raw_api_results.json"):
        nct_id = raw_path.parent.name
        if not nct_id.startswith("NCT"):
            continue
        try:
            if strict_toxicity:
                doc = process_one_raw(raw_path, nct_id)
            else:
                raw_bytes = raw_path.read_bytes()
                raw = json.loads(raw_bytes.decode("utf-8"))
                sm = _derive_safety_margin(raw)
                rr = _derive_response_ratio(raw)
                missing = ["binding_affinity_nM (ba)", "oral_bioavailability_percent (bio)"]
                audit_note = _build_audit_note(sm, rr, missing, toxicity_missing=True)
                doc = _build_doc(raw, raw_bytes, raw_path, nct_id, sm, rr, audit_note)
            out_name = f"{nct_id}_Learning_Material.json"
            out_path = learning_dir / out_name
            with open(out_path, "w", encoding="utf-8") as f:
                json.dump(doc, f, indent=2)
            count += 1
            print(f"  Bridge: {nct_id} -> {out_name}")
        except TraceabilityError as e:
            print(f"  TraceabilityError ({nct_id}): {e}")
            if strict_toxicity:
                raise
        except Exception as e:
            print(f"  Error processing {nct_id}: {e}")
    return count


def main() -> None:
    import argparse
    p = argparse.ArgumentParser(description="Olympus Bridge v2: NCT raw -> Learning_Material")
    p.add_argument("--olympus", type=Path, default=OLYMPUS_RESEARCH, help="Olympus_Research root")
    p.add_argument("--output", type=Path, default=LEARNING_MATERIAL_DIR, help="Learning_Material dir")
    p.add_argument("--allow-missing-toxicity", action="store_true", help="Do not raise TraceabilityError for missing toxicity_index")
    args = p.parse_args()
    print("Olympus Bridge v2: scanning for raw_api_results.json...")
    n = run_bridge(
        olympus_root=args.olympus,
        learning_dir=args.output,
        strict_toxicity=not args.allow_missing_toxicity,
    )
    print(f"Done. Wrote {n} Learning_Material JSON(s).")


if __name__ == "__main__":
    main()
