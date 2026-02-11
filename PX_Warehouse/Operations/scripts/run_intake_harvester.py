#!/usr/bin/env python3
"""
Intake harvester: seed libraries → normalized, filtered, governed queue.

Flow:
  1. Load intake_policy.json; abort if missing or no enabled sources.
  2. Harvest from each enabled source (ChEMBL, DrugBank, ClinicalTrials, FDA_OrangeBook, PubChem, NCATS_Repurposing).
  3. Normalize and dedupe (canonical SMILES, InChIKey dedup).
  4. Apply scientific filters (Lipinski, Veber, Ghose, PAINS, etc. when enabled).
  5. Apply governance (provenance, timestamp, reproducibility).
  6. Write queue file (prv_24h_queue.json or policy.output.queue_file).

Run from repo root or ensure PYTHONPATH includes repo root.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import random
import sys
from pathlib import Path
from datetime import datetime, timezone
from typing import Any, Dict, List

# Repo root: this script lives at PX_Warehouse/Operations/scripts/run_intake_harvester.py
_SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = _SCRIPT_DIR.parents[1]  # PX_Warehouse -> foundation
PX_LOGS = REPO_ROOT / "PX_LOGS"
PX_LOGS.mkdir(parents=True, exist_ok=True)


def log_info(message: str) -> None:
    print("[INFO]", message)
    _append_log("INFO", message)


def log_warning(message: str) -> None:
    print("[WARN]", message)
    _append_log("WARN", message)


def log_error(message: str, error: Exception | None = None) -> None:
    print("[ERROR]", message)
    if error is not None:
        print("        ", str(error))
    _append_log("ERROR", message, error)


def _append_log(level: str, message: str, error: Exception | None = None) -> None:
    try:
        ts = datetime.now(timezone.utc).strftime("%Y%m%d")
        log_path = PX_LOGS / f"intake_harvester_{ts}.log"
        with open(log_path, "a", encoding="utf-8") as f:
            f.write(f"{datetime.now(timezone.utc).isoformat()} [{level}] {message}\n")
            if error is not None:
                f.write(f"        {error}\n")
    except Exception:
        pass


def detect_repo_root() -> Path:
    """Mirror what PRV_24H_Orchestrator uses: foundation as repo root."""
    return Path(__file__).resolve().parents[1]


def get_intake_policy(repo_root: Path) -> dict[str, Any]:
    """Load intake_policy.json; return {} if missing, log and return {} on parse error."""
    policy_path = repo_root / "PX_Warehouse" / "Operations" / "Inputs" / "intake_policy.json"
    if not policy_path.exists():
        return {}
    try:
        return json.loads(policy_path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError) as e:
        log_error("Failed to parse intake_policy.json", e)
        return {}


def get_enabled_sources_from_policy(policy: dict[str, Any]) -> List[str]:
    """Return list of source names where sources[name].enabled is True."""
    enabled: List[str] = []
    sources = policy.get("sources") or {}
    for source_name, cfg in sources.items():
        if isinstance(cfg, dict) and cfg.get("enabled") is True:
            enabled.append(source_name)
    return enabled


# --- Helpers for harvesters (candidate shape: id, name, smiles, source, metadata) ---


def _load_csv(path: Path) -> List[dict[str, Any]]:
    with path.open("r", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def _load_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def _hash_str(s: str) -> str:
    return hashlib.sha256(s.encode("utf-8")).hexdigest()


def _make_candidate(
    id_: str, name: str, smiles: str, source: str, metadata: dict[str, Any] | None = None
) -> dict[str, Any]:
    return {
        "id": id_ or "",
        "name": name or "",
        "smiles": smiles or "",
        "source": source,
        "metadata": metadata or {},
    }


def _resolve_path(path_str: str, repo_root: Path) -> Path | None:
    """Resolve path_str; if relative, under repo_root. Return None if path_str empty."""
    if not path_str or not path_str.strip():
        return None
    p = Path(path_str.strip())
    if not p.is_absolute():
        p = repo_root / p
    return p if p.exists() else None


def harvest_from_source(
    source_name: str, source_cfg: dict[str, Any], repo_root: Path
) -> List[dict[str, Any]]:
    """Dispatch to source-specific harvester. Returns list of candidate dicts (id, name, smiles, source, metadata)."""
    if source_name == "ChEMBL":
        return harvest_chembl(source_cfg, repo_root)
    if source_name == "DrugBank":
        return harvest_drugbank(source_cfg, repo_root)
    if source_name == "ClinicalTrials":
        return harvest_clinicaltrials(source_cfg, repo_root)
    if source_name == "FDA_OrangeBook":
        return harvest_orange_book(source_cfg, repo_root)
    if source_name == "PubChem":
        return harvest_pubchem(source_cfg, repo_root)
    if source_name == "NCATS_Repurposing":
        return harvest_ncats(source_cfg, repo_root)
    if source_name == "ZINC":
        return harvest_zinc(source_cfg, repo_root)
    log_warning(f"Unknown source: {source_name}")
    return []


def harvest_chembl(source_cfg: dict[str, Any], repo_root: Path) -> List[dict[str, Any]]:
    """
    Harvest candidates from ChEMBL.
    source_cfg: include_approved, include_investigational, include_withdrawn; optional local_cache_path.
    """
    candidates: List[dict[str, Any]] = []
    path_str = source_cfg.get("local_cache_path")
    path = _resolve_path(path_str, repo_root) if path_str else None
    if not path:
        return candidates
    try:
        rows = _load_csv(path)
    except Exception as e:
        log_warning(f"ChEMBL cache unreadable: {path}: {e}")
        return []
    for row in rows:
        status = (row.get("status") or "").lower()
        if status == "approved" and not source_cfg.get("include_approved", True):
            continue
        if status == "investigational" and not source_cfg.get("include_investigational", True):
            continue
        if status == "withdrawn" and not source_cfg.get("include_withdrawn", True):
            continue
        chembl_id = row.get("chembl_id") or row.get("id", "")
        name = row.get("pref_name") or row.get("name", "")
        smiles = row.get("canonical_smiles") or row.get("smiles", "")
        if not smiles:
            continue
        metadata = {
            "status": status,
            "targets": row.get("targets"),
            "mechanism": row.get("mechanism_of_action"),
            "source_record": "ChEMBL",
        }
        candidates.append(_make_candidate(chembl_id, name, smiles, "ChEMBL", metadata))
    return candidates


def harvest_drugbank(source_cfg: dict[str, Any], repo_root: Path) -> List[dict[str, Any]]:
    """
    Harvest candidates from DrugBank public data.
    source_cfg: include_approved, include_investigational; local_file path to CSV/TSV.
    """
    candidates: List[dict[str, Any]] = []
    path = _resolve_path(source_cfg.get("local_file") or "", repo_root)
    if not path:
        return []
    try:
        rows = _load_csv(path)
    except Exception as e:
        log_warning(f"DrugBank file unreadable: {path}: {e}")
        return []
    for row in rows:
        status = (row.get("status") or "").lower()
        if "approved" in status and not source_cfg.get("include_approved", True):
            continue
        if "investigational" in status and not source_cfg.get("include_investigational", True):
            continue
        dbid = row.get("drugbank_id") or row.get("id", "")
        name = row.get("name", "")
        smiles = row.get("smiles", "")
        if not smiles:
            continue
        metadata = {
            "status": status,
            "atc_codes": row.get("atc_codes"),
            "targets": row.get("targets"),
            "source_record": "DrugBank",
        }
        candidates.append(_make_candidate(dbid, name, smiles, "DrugBank", metadata))
    return candidates


def harvest_clinicaltrials(source_cfg: dict[str, Any], repo_root: Path) -> List[dict[str, Any]]:
    """
    Harvest candidates from ClinicalTrials.gov.
    source_cfg: include_phase_1/2/3/4; optional local_file (CSV).
    """
    candidates: List[dict[str, Any]] = []
    path = _resolve_path(source_cfg.get("local_file") or "", repo_root)
    if not path:
        return []
    try:
        rows = _load_csv(path)
    except Exception as e:
        log_warning(f"ClinicalTrials file unreadable: {path}: {e}")
        return []
    for row in rows:
        phase = (row.get("phase") or "").lower()
        if "phase 1" in phase and not source_cfg.get("include_phase_1", True):
            continue
        if "phase 2" in phase and not source_cfg.get("include_phase_2", True):
            continue
        if "phase 3" in phase and not source_cfg.get("include_phase_3", True):
            continue
        if "phase 4" in phase and not source_cfg.get("include_phase_4", True):
            continue
        nct_id = row.get("nct_id") or row.get("id", "")
        drug_name = row.get("intervention_name") or row.get("name", "")
        smiles = row.get("smiles", "")
        if not smiles:
            continue
        metadata = {
            "phase": phase,
            "condition": row.get("condition"),
            "sponsor": row.get("sponsor"),
            "source_record": "ClinicalTrials",
        }
        candidates.append(_make_candidate(nct_id, drug_name, smiles, "ClinicalTrials", metadata))
    return candidates


def harvest_orange_book(source_cfg: dict[str, Any], repo_root: Path) -> List[dict[str, Any]]:
    """
    Harvest candidates from FDA Orange Book.
    source_cfg: include_small_molecules, include_patent_expired, include_exclusivity_expired; local_file CSV.
    """
    candidates: List[dict[str, Any]] = []
    path = _resolve_path(source_cfg.get("local_file") or "", repo_root)
    if not path:
        return []
    try:
        rows = _load_csv(path)
    except Exception as e:
        log_warning(f"FDA Orange Book file unreadable: {path}: {e}")
        return []
    for row in rows:
        if not source_cfg.get("include_small_molecules", True):
            if (row.get("is_small_molecule") or "Y") != "Y":
                continue
        patent_expiry = row.get("patent_expiry")
        exclusivity_expiry = row.get("exclusivity_expiry")
        appl_no = row.get("appl_no", "")
        prod_no = row.get("product_no", "")
        name = row.get("drug_name") or row.get("name", "")
        smiles = row.get("smiles", "")
        if not smiles:
            continue
        metadata = {
            "dosage_form": row.get("dosage_form"),
            "route": row.get("route"),
            "patent_expiry": patent_expiry,
            "exclusivity_expiry": exclusivity_expiry,
            "source_record": "FDA_OrangeBook",
        }
        id_ = f"OB_{appl_no}_{prod_no}"
        candidates.append(_make_candidate(id_, name, smiles, "FDA_OrangeBook", metadata))
    return candidates


def harvest_pubchem(source_cfg: dict[str, Any], repo_root: Path) -> List[dict[str, Any]]:
    """
    Harvest candidates from PubChem.
    source_cfg: include_bioactive, include_assay_hits; optional local_file (CSV).
    """
    candidates: List[dict[str, Any]] = []
    path = _resolve_path(source_cfg.get("local_file") or "", repo_root)
    if not path:
        return []
    try:
        rows = _load_csv(path)
    except Exception as e:
        log_warning(f"PubChem file unreadable: {path}: {e}")
        return []
    for row in rows:
        cid = row.get("cid") or row.get("id", "")
        name = row.get("name", "")
        smiles = row.get("smiles", "")
        if not smiles:
            continue
        bioactivity_flag = (row.get("bioactive") or "N") == "Y"
        assay_hit_flag = (row.get("assay_hit") or "N") == "Y"
        if bioactivity_flag and not source_cfg.get("include_bioactive", True):
            continue
        if assay_hit_flag and not source_cfg.get("include_assay_hits", True):
            continue
        metadata = {
            "bioactive": bioactivity_flag,
            "assay_hit": assay_hit_flag,
            "source_record": "PubChem",
        }
        candidates.append(_make_candidate(f"CID_{cid}", name, smiles, "PubChem", metadata))
    return candidates


def harvest_ncats(source_cfg: dict[str, Any], repo_root: Path) -> List[dict[str, Any]]:
    """
    Harvest candidates from NCATS Repurposing Collection.
    source_cfg: include_all; local_file CSV path.
    """
    candidates: List[dict[str, Any]] = []
    path = _resolve_path(source_cfg.get("local_file") or "", repo_root)
    if not path:
        return []
    try:
        rows = _load_csv(path)
    except Exception as e:
        log_warning(f"NCATS file unreadable: {path}: {e}")
        return []
    for row in rows:
        ncats_id = row.get("ncats_id") or row.get("id", "")
        name = row.get("name", "")
        smiles = row.get("smiles", "")
        if not smiles:
            continue
        metadata = {
            "status": row.get("status"),
            "atc": row.get("atc"),
            "mechanism": row.get("mechanism"),
            "source_record": "NCATS_Repurposing",
        }
        candidates.append(_make_candidate(ncats_id, name, smiles, "NCATS_Repurposing", metadata))
    return candidates


def harvest_zinc(source_cfg: dict[str, Any], repo_root: Path) -> List[dict[str, Any]]:
    """
    Harvest candidates from ZINC (ZINC15/20 commercially purchasable).
    source_cfg: include_in_stock; local_file CSV with zinc_id, name, smiles, catalog.
    """
    candidates: List[dict[str, Any]] = []
    path = _resolve_path(source_cfg.get("local_file") or "", repo_root)
    if not path:
        return []
    try:
        rows = _load_csv(path)
    except Exception as e:
        log_warning(f"ZINC file unreadable: {path}: {e}")
        return []
    for row in rows:
        zinc_id = row.get("zinc_id") or row.get("id", "")
        name = row.get("name", "")
        smiles = row.get("smiles", "")
        if not smiles:
            continue
        metadata = {
            "catalog": row.get("catalog", "ZINC20"),
            "source_record": "ZINC",
        }
        candidates.append(_make_candidate(zinc_id, name, smiles, "ZINC", metadata))
    return candidates


def _rdkit_available() -> bool:
    try:
        from rdkit import Chem
        return True
    except ImportError:
        return False


def _canonicalize_smiles(smiles: str) -> str:
    if _rdkit_available():
        try:
            from rdkit import Chem
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                return Chem.MolToSmiles(mol)
        except Exception:
            pass
    return smiles.strip()


def _is_valid_smiles(smiles: str) -> bool:
    if not smiles or not isinstance(smiles, str):
        return False
    if _rdkit_available():
        try:
            from rdkit import Chem
            return Chem.MolFromSmiles(smiles) is not None
        except Exception:
            return False
    return len(smiles) > 0 and " " not in smiles


def _compute_inchikey_or_hash(smiles: str) -> str:
    if _rdkit_available():
        try:
            from rdkit import Chem
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                key = Chem.MolToInchiKey(mol)
                if key:
                    return key
        except Exception:
            pass
    return hashlib.sha256(smiles.encode()).hexdigest()[:32]


def normalize_and_dedupe(
    raw_candidates: List[dict[str, Any]], filter_cfg: dict[str, Any]
) -> List[dict[str, Any]]:
    """Normalize SMILES, optionally drop invalid; dedupe by InChIKey or SMILES hash."""
    remove_invalid = filter_cfg.get("remove_invalid_smiles", True)
    normalized: List[dict[str, Any]] = []

    for c in raw_candidates:
        smiles = (c.get("smiles") or "").strip()
        if not smiles:
            continue
        if remove_invalid and not _is_valid_smiles(smiles):
            continue
        c = dict(c)
        c["smiles"] = _canonicalize_smiles(smiles)
        c.setdefault("metadata", {})
        normalized.append(c)

    seen: set[str] = set()
    deduped: List[dict[str, Any]] = []
    for c in normalized:
        key = _compute_inchikey_or_hash(c["smiles"])
        if key in seen:
            continue
        seen.add(key)
        deduped.append(c)
    return deduped


def _compute_physchem_properties(smiles: str) -> dict[str, Any] | None:
    if not _rdkit_available():
        return None
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, rdMolDescriptors
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return {
            "mw": Descriptors.MolWt(mol),
            "logp": Descriptors.MolLogP(mol),
            "hbd": rdMolDescriptors.CalcNumHBD(mol),
            "hba": rdMolDescriptors.CalcNumHBA(mol),
            "tpsa": Descriptors.TPSA(mol),
            "rot_bonds": Descriptors.NumRotatableBonds(mol),
        }
    except Exception:
        return None


def _passes_lipinski(props: dict[str, Any]) -> bool:
    if not props:
        return True
    violations = 0
    if (props.get("mw") or 0) > 500:
        violations += 1
    if (props.get("logp") or 0) > 5:
        violations += 1
    if (props.get("hbd") or 0) > 5:
        violations += 1
    if (props.get("hba") or 0) > 10:
        violations += 1
    return violations <= 1


def _passes_veber(props: dict[str, Any]) -> bool:
    if not props:
        return True
    tpsa = props.get("tpsa") or 0
    rot = props.get("rot_bonds") or 0
    return tpsa <= 140 and rot <= 10


def _passes_ghose(props: dict[str, Any]) -> bool:
    """Ghose drug-likeness. Stub: pass if we have props."""
    if not props:
        return True
    mw = props.get("mw") or 0
    logp = props.get("logp") or 0
    return 160 <= mw <= 480 and -0.4 <= logp <= 5.6


def _is_pains_hit(smiles: str) -> bool:
    """PAINS filter. Stub: no PAINS library; assume pass."""
    return False


def _has_structural_alert(smiles: str) -> bool:
    """Structural alerts. Stub: assume pass."""
    return False


def _passes_admet_prefilter(smiles: str, props: dict[str, Any] | None) -> bool:
    """ADMET prefilter. Stub: pass."""
    return True


def _passes_bbb_prediction(smiles: str, props: dict[str, Any] | None) -> bool:
    """BBB prediction. Stub: pass."""
    return True


def _is_efflux_liability(smiles: str, props: dict[str, Any] | None) -> bool:
    """Efflux liability. Stub: assume no liability."""
    return False


def _passes_sas(smiles: str) -> bool:
    """Synthetic accessibility. Stub: pass."""
    return True


def apply_scientific_filters(
    candidates: List[dict[str, Any]], sci_cfg: dict[str, Any]
) -> List[dict[str, Any]]:
    """Apply Lipinski, Veber, Ghose, PAINS, structural alerts, ADMET, BBB, efflux, SAS when enabled."""
    filtered: List[dict[str, Any]] = []
    for c in candidates:
        smiles = c.get("smiles") or ""
        props = _compute_physchem_properties(smiles)

        if sci_cfg.get("lipinski") and not _passes_lipinski(props):
            continue
        if sci_cfg.get("veber") and not _passes_veber(props):
            continue
        if sci_cfg.get("ghose") and not _passes_ghose(props):
            continue
        if sci_cfg.get("pains") and _is_pains_hit(smiles):
            continue
        if sci_cfg.get("structural_alerts") and _has_structural_alert(smiles):
            continue
        if sci_cfg.get("admet_prefilter") and not _passes_admet_prefilter(smiles, props):
            continue
        if sci_cfg.get("bbb_prediction") and not _passes_bbb_prediction(smiles, props):
            continue
        if sci_cfg.get("efflux_prediction") and _is_efflux_liability(smiles, props):
            continue
        if sci_cfg.get("synthetic_accessibility") and not _passes_sas(smiles):
            continue
        filtered.append(c)
    return filtered


def _has_provenance(c: dict[str, Any]) -> bool:
    """Candidate has minimal provenance (source + id or smiles)."""
    return bool(c.get("source")) and (bool(c.get("id")) or bool(c.get("smiles")))


def _hash_of_current_policy(repo_root: Path) -> str:
    """Hash of intake_policy.json for reproducibility."""
    p = repo_root / "PX_Warehouse" / "Operations" / "Inputs" / "intake_policy.json"
    if not p.exists():
        return ""
    try:
        return hashlib.sha256(p.read_bytes()).hexdigest()[:16]
    except Exception:
        return ""


def _compute_source_fingerprint(c: dict[str, Any]) -> str:
    """Fingerprint for reproducibility: source + id + smiles hash."""
    parts = [c.get("source", ""), c.get("id", ""), c.get("smiles", "")]
    return hashlib.sha256(json.dumps(parts, sort_keys=True).encode()).hexdigest()[:16]


def apply_governance(
    candidates: List[dict[str, Any]], gov_cfg: dict[str, Any], repo_root: Path
) -> List[dict[str, Any]]:
    """Attach timestamp; require provenance when governance.provenance_required; add policy hash and source fingerprint when reproducibility_required."""
    governed: List[dict[str, Any]] = []
    now_ts = datetime.now(timezone.utc).isoformat()
    policy_hash = _hash_of_current_policy(repo_root) if gov_cfg.get("reproducibility_required") else ""

    for c in candidates:
        c = dict(c)
        c.setdefault("metadata", {})

        if gov_cfg.get("provenance_required") and not _has_provenance(c):
            continue
        if gov_cfg.get("timestamp_required"):
            c["metadata"]["intake_timestamp"] = now_ts
        if gov_cfg.get("reproducibility_required"):
            if policy_hash:
                c["metadata"]["intake_policy_hash"] = policy_hash
            c["metadata"]["source_fingerprint"] = _compute_source_fingerprint(c)
        governed.append(c)
    return governed


def build_queue_path(repo_root: Path, queue_filename: str) -> Path:
    from PX_Warehouse.warehouse_layout import get_feeder_dir
    return get_feeder_dir(repo_root) / queue_filename


def _generate_queue_id(c: dict[str, Any], idx: int) -> str:
    """Generate id like PRV_REP_<short_hash> or PRV_NOV_<short_hash>."""
    source = (c.get("source") or "").upper()
    is_repurposed = source in ("CHEMBL", "DRUGBANK", "CLINICALTRIALS", "FDA_ORANGEBOOK", "PUBCHEM", "NCATS_REPURPOSING", "ZINC")
    prefix = "PRV_REP" if is_repurposed else "PRV_NOV"
    raw = f"{c.get('id', '')}_{c.get('smiles', '')[:20]}_{idx}"
    h = hashlib.sha256(raw.encode()).hexdigest()[:8]
    return f"{prefix}_{h}"


def _infer_type(c: dict[str, Any]) -> str:
    """Infer R (repurposed) or N (novel) from source or metadata."""
    source = (c.get("source") or "").upper()
    if source in ("CHEMBL", "DRUGBANK", "CLINICALTRIALS", "FDA_ORANGEBOOK", "PUBCHEM", "NCATS_REPURPOSING", "ZINC"):
        return "R"
    return "N"


def build_queue_items(
    candidates: List[dict[str, Any]], output_cfg: dict[str, Any]
) -> dict[str, Any]:
    """Build { "candidates": [ { id, type, smiles, name, source, metadata }, ... ] }; apply max_items and shuffle."""
    items: List[dict[str, Any]] = []
    for idx, c in enumerate(candidates):
        items.append({
            "id": _generate_queue_id(c, idx),
            "type": _infer_type(c),
            "smiles": c.get("smiles", ""),
            "name": c.get("name") or c.get("id", ""),
            "source": c.get("source", ""),
            "metadata": c.get("metadata", {}),
        })

    max_items = output_cfg.get("max_items", "unbounded")
    if max_items != "unbounded":
        try:
            n = int(max_items)
            items = items[:n]
        except (TypeError, ValueError):
            pass
    if output_cfg.get("shuffle") is True:
        random.shuffle(items)
    return {"candidates": items}


def write_queue_file(queue_path: Path, queue_items: dict[str, Any]) -> None:
    queue_path.parent.mkdir(parents=True, exist_ok=True)
    with open(queue_path, "w", encoding="utf-8") as f:
        json.dump(queue_items, f, indent=2, ensure_ascii=False)


def _load_existing_queue(queue_path: Path) -> List[dict[str, Any]]:
    """Load existing candidates from queue file; return [] if missing or invalid."""
    if not queue_path.exists():
        return []
    try:
        data = json.loads(queue_path.read_text(encoding="utf-8"))
        cands = data.get("candidates", []) if isinstance(data, dict) else []
        return cands if isinstance(cands, list) else []
    except Exception:
        return []


def _merge_with_existing(
    existing: List[dict[str, Any]], new_items: List[dict[str, Any]]
) -> List[dict[str, Any]]:
    """Merge new_items into existing; dedupe by normalized SMILES (keep existing first)."""
    seen_smiles: set[str] = set()
    result: List[dict[str, Any]] = []
    for c in existing:
        s = (c.get("smiles") or "").strip()
        if not s or s in seen_smiles:
            continue
        seen_smiles.add(s)
        result.append(c)
    for c in new_items:
        s = (c.get("smiles") or "").strip()
        if not s or s in seen_smiles:
            continue
        seen_smiles.add(s)
        result.append(c)
    return result


def _dedupe_queue_by_id(candidates: List[dict[str, Any]]) -> List[dict[str, Any]]:
    """Return list with first occurrence per id so queue has unique ids."""
    seen: set[str] = set()
    out: List[dict[str, Any]] = []
    for c in candidates:
        kid = (c.get("id") or "").strip()
        if not kid or kid in seen:
            continue
        seen.add(kid)
        out.append(c)
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description="Intake harvester: seed libraries → PRV queue.")
    parser.add_argument(
        "--merge",
        action="store_true",
        help="Merge harvested repurposed candidates into existing queue (Feeder); do not overwrite.",
    )
    args = parser.parse_args()

    repo_root = detect_repo_root()
    policy = get_intake_policy(repo_root)

    if not policy:
        log_error("No intake_policy.json found; aborting intake.")
        return 1

    enabled_sources = get_enabled_sources_from_policy(policy)
    if not enabled_sources:
        log_error("No enabled sources in intake policy; aborting intake.")
        return 1

    raw_candidates: List[dict[str, Any]] = []
    for source_name in enabled_sources:
        sources = policy.get("sources") or {}
        source_cfg = sources.get(source_name) or {}
        source_candidates = harvest_from_source(source_name, source_cfg, repo_root)
        log_info(f"{source_name}: harvested {len(source_candidates)} candidates")
        raw_candidates.extend(source_candidates)

    log_info(f"Total raw candidates before normalization: {len(raw_candidates)}")

    filter_cfg = policy.get("filters") or {}
    normalized_candidates = normalize_and_dedupe(raw_candidates, filter_cfg)
    log_info(f"Total candidates after normalization/deduplication: {len(normalized_candidates)}")

    sci_cfg = policy.get("scientific_filters") or {}
    filtered_candidates = apply_scientific_filters(normalized_candidates, sci_cfg)
    log_info(f"Total candidates after scientific filters: {len(filtered_candidates)}")

    gov_cfg = policy.get("governance") or {}
    final_candidates = apply_governance(filtered_candidates, gov_cfg, repo_root)
    log_info(f"Total candidates after governance checks: {len(final_candidates)}")

    output_cfg = policy.get("output") or {}
    queue_filename = output_cfg.get("queue_file", "prv_24h_queue.json")
    queue_path = build_queue_path(repo_root, queue_filename)
    queue_items = build_queue_items(final_candidates, output_cfg)
    new_candidates = queue_items["candidates"]

    if args.merge and queue_path.exists():
        existing = _load_existing_queue(queue_path)
        merged = _merge_with_existing(existing, new_candidates)
        queue_items = {"candidates": merged}
        log_info(f"Merge: kept {len(existing)} existing, added {len(new_candidates)} harvested → {len(merged)} total")

    # Ensure written queue has unique ids (no duplicate queue entries)
    final = _dedupe_queue_by_id(queue_items["candidates"])
    if len(final) < len(queue_items["candidates"]):
        log_info(f"Dedupe by id: {len(queue_items['candidates'])} → {len(final)} unique")
    queue_items = {"candidates": final}

    write_queue_file(queue_path, queue_items)
    log_info(f"Wrote queue file: {queue_path} with {len(queue_items['candidates'])} items")

    return 0


if __name__ == "__main__":
    sys.exit(main())
