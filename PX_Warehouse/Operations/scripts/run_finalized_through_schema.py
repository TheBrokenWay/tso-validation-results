#!/usr/bin/env python3
"""
Run every dossier in Finalized_Dossiers/Diamond and Finalized_Dossiers/Gold
through the new schema (DOSSIER_SCHEMA_v2.json): normalize keys and ensure
schema-aligned shape, then overwrite in place.

- Renames predicted_therapeutic_window -> therapeutic_window anywhere in the tree.
- Ensures dossier_version is "2.0.0".
- Ensures top-level schema keys exist (add null if missing).
- Writes back to the same file.

Usage (from repo root):
  python PX_Warehouse/run_finalized_through_schema.py
  python PX_Warehouse/run_finalized_through_schema.py --dry-run
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

WAREHOUSE = _REPO_ROOT / "PX_Warehouse"
FINALIZED = WAREHOUSE / "Finalized_Dossiers"
SCHEMA_PATH = WAREHOUSE / "DOSSIER_SCHEMA_v2.json"

# Top-level keys from DOSSIER_SCHEMA_v2 that should exist (add as null if missing)
SCHEMA_TOP_LEVEL = [
    "dossier_version",
    "generated",
    "alcoa_metadata",
    "causal_trace_log",
    "candidate",
    "engines",
    "harm_energy",
    "context",
    "fda_compliance",
    "constitutional_seal",
    "trial_simulation_id",
    "trial_simulation_hash",
    "trial_outcome_summary",
    "finalization",
    "worldline_generated",
    "worldline_path",
    "refinery_status",
]


def _normalize_node(obj):
    """Recursively normalize: rename predicted_therapeutic_window -> therapeutic_window."""
    if isinstance(obj, dict):
        out = {}
        for k, v in obj.items():
            if k == "predicted_therapeutic_window":
                out["therapeutic_window"] = _normalize_node(v)
            else:
                out[k] = _normalize_node(v)
        return out
    if isinstance(obj, list):
        return [_normalize_node(x) for x in obj]
    return obj


def normalize_dossier(d: dict) -> dict:
    """Ensure schema-aligned shape: version, top-level keys, therapeutic_window name."""
    d = _normalize_node(d)
    d["dossier_version"] = "2.0.0"
    fin = d.get("finalization") or {}
    # Map top-level key <- finalization key when missing at top level
    copy_from_fin = {
        "trial_simulation_id": "trial_simulation_id",
        "trial_simulation_hash": "trial_simulation_hash",
        "trial_outcome_summary": "trial_outcome_summary",
        "worldline_generated": "worldline_generated",
        "worldline_path": "worldline_path",
    }
    for key in SCHEMA_TOP_LEVEL:
        if key not in d:
            d[key] = None
        if d[key] is None and key in copy_from_fin and copy_from_fin[key] in fin:
            d[key] = fin[copy_from_fin[key]]
    return d


def main() -> int:
    dry_run = "--dry-run" in sys.argv
    tiers = ["Diamond", "Gold"]
    total = 0
    ok = 0
    errors = []

    for tier in tiers:
        tier_dir = FINALIZED / tier
        if not tier_dir.exists():
            print(f"Skip {tier}: directory not found")
            continue
        for path in sorted(tier_dir.glob("*.json")):
            total += 1
            try:
                text = path.read_text(encoding="utf-8")
                dossier = json.loads(text)
                dossier = normalize_dossier(dossier)
                if not dry_run:
                    path.write_text(json.dumps(dossier, indent=2, default=str), encoding="utf-8")
                ok += 1
                print(f"  {tier}/{path.name}")
            except Exception as e:
                errors.append((str(path), str(e)))
                print(f"  {tier}/{path.name} ERROR: {e}")

    print()
    print(f"Processed: {ok}/{total} dossiers" + (" (dry-run, no writes)" if dry_run else ""))
    if errors:
        print(f"Errors: {len(errors)}")
        for p, e in errors:
            print(f"  {p}: {e}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
