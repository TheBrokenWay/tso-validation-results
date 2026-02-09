#!/usr/bin/env python3
"""
Restore data that was written to false paths back into canonical warehouse folders.

False paths → Correct locations:
- TrialSimulations/*.json (root)           → TrialSimulations/LiveRuns/run_recovery_<ts>/
- CommercialAssets/Gold/<subdir>/*         → CommercialAssets/Gold/ (flat)
- CommercialAssets/Diamond/TRIAL_SIM_*     → TrialSimulations/LiveRuns/run_recovery_<ts>/
- CommercialAssets/Gold/TRIAL_SIM_*        → TrialSimulations/LiveRuns/run_recovery_<ts>/
- CommercialAssets/Learning_Material/TRIAL_SIM_* → TrialSimulations/LiveRuns/run_recovery_<ts>/ (optional)
- CommercialAssets/Silver/TRIAL_SIM_*      → TrialSimulations/LiveRuns/run_recovery_<ts>/
- CommercialAssets/Bronze/TRIAL_SIM_*     → TrialSimulations/LiveRuns/run_recovery_<ts>/

Run from repo root or script dir. Writes a manifest of all moves.
"""
from __future__ import annotations

import json
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

_SCRIPT_DIR = Path(__file__).resolve().parent
_REPO_ROOT = _SCRIPT_DIR.parents[1]
_WAREHOUSE = _REPO_ROOT / "PX_Warehouse"
_TS = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
_RECOVERY_RUN = _WAREHOUSE / "TrialSimulations" / "LiveRuns" / f"run_recovery_{_TS}"

# Canonical targets
COMMERCIAL = _WAREHOUSE / "CommercialAssets"
GOLD_FLAT = COMMERCIAL / "Gold"
LIVE_RUNS = _WAREHOUSE / "TrialSimulations" / "LiveRuns"
TRIAL_SIM_ROOT = _WAREHOUSE / "TrialSimulations"


def _ensure_recovery_run() -> Path:
    _RECOVERY_RUN.mkdir(parents=True, exist_ok=True)
    return _RECOVERY_RUN


def _unique_dest(parent: Path, base_name: str) -> Path:
    dest = parent / base_name
    if not dest.exists():
        return dest
    stem = base_name.rsplit(".", 1)[0] if "." in base_name else base_name
    ext = "." + base_name.rsplit(".", 1)[-1] if "." in base_name else ""
    n = 1
    while True:
        dest = parent / f"{stem}_recovered_{n}{ext}"
        if not dest.exists():
            return dest
        n += 1


def move_trial_sim_dossiers_from_root(manifest: list, dry_run: bool) -> int:
    """TrialSimulations/*.json -> LiveRuns/run_recovery_<ts>/"""
    count = 0
    for f in list(TRIAL_SIM_ROOT.glob("TRIAL_SIMULATION_DOSSIER-*.json")):
        dest = _unique_dest(_RECOVERY_RUN, f.name) if not dry_run else _RECOVERY_RUN / f.name
        manifest.append({"from": "TrialSimulations_root", "src": str(f), "dest": str(dest)})
        if not dry_run:
            shutil.move(str(f), str(dest))
        count += 1
    return count


def move_trial_sim_from_tier(tier_name: str, manifest: list, dry_run: bool) -> int:
    """CommercialAssets/<tier>/TRIAL_SIMULATION_DOSSIER-*.json -> LiveRuns/run_recovery_<ts>/"""
    tier_dir = COMMERCIAL / tier_name
    if not tier_dir.exists():
        return 0
    count = 0
    for f in list(tier_dir.glob("TRIAL_SIMULATION_DOSSIER-*.json")):
        dest = _unique_dest(_RECOVERY_RUN, f.name) if not dry_run else _RECOVERY_RUN / f.name
        manifest.append({"from": f"CommercialAssets/{tier_name}", "src": str(f), "dest": str(dest)})
        if not dry_run:
            shutil.move(str(f), str(dest))
        count += 1
    return count


def move_trial_sim_from_learning_material(manifest: list, dry_run: bool, limit: int | None = 5000) -> int:
    """CommercialAssets/Learning_Material/TRIAL_SIMULATION_DOSSIER-*.json -> LiveRuns/run_recovery_<ts>/"""
    learn_dir = COMMERCIAL / "Learning_Material"
    if not learn_dir.exists():
        return 0
    count = 0
    for f in list(learn_dir.glob("TRIAL_SIMULATION_DOSSIER-*.json")):
        if limit is not None and count >= limit:
            break
        dest = _unique_dest(_RECOVERY_RUN, f.name) if not dry_run else _RECOVERY_RUN / f.name
        manifest.append({"from": "CommercialAssets/Learning_Material", "src": str(f), "dest": str(dest)})
        if not dry_run:
            shutil.move(str(f), str(dest))
        count += 1
    return count


def flatten_gold_subfolders(manifest: list, dry_run: bool) -> int:
    """Gold/<subdir>/* -> Gold/<unique_name> (flat)."""
    if not GOLD_FLAT.exists():
        return 0
    count = 0
    for subdir in list(GOLD_FLAT.iterdir()):
        if not subdir.is_dir():
            continue
        for f in list(subdir.glob("*")):
            if not f.is_file():
                continue
            # Avoid moving PATH_TEST into flat Gold as "recovery" - keep as-is or skip
            if f.name.startswith("PATH_TEST"):
                continue
            base = f"{subdir.name}_{f.name}"
            dest = _unique_dest(GOLD_FLAT, base) if not dry_run else GOLD_FLAT / base
            manifest.append({"from": "CommercialAssets/Gold/subdir", "src": str(f), "dest": str(dest)})
            if not dry_run:
                shutil.move(str(f), str(dest))
            count += 1
        if not dry_run and subdir.exists() and not any(subdir.iterdir()):
            subdir.rmdir()
    return count


def run_recovery(dry_run: bool = False, include_learning_material: bool = True, learning_limit: int = 5000) -> dict:
    if not dry_run:
        _ensure_recovery_run()
    manifest = []
    stats = {
        "run_ts": datetime.now(timezone.utc).isoformat(),
        "recovery_run_dir": str(_RECOVERY_RUN),
        "dry_run": dry_run,
        "moves": manifest,
        "counts": {
            "trial_sim_from_root": 0,
            "trial_sim_from_gold": 0,
            "trial_sim_from_diamond": 0,
            "trial_sim_from_silver": 0,
            "trial_sim_from_bronze": 0,
            "trial_sim_from_learning": 0,
            "gold_flatten": 0,
        },
    }

    stats["counts"]["trial_sim_from_root"] = move_trial_sim_dossiers_from_root(manifest, dry_run)
    stats["counts"]["trial_sim_from_gold"] = move_trial_sim_from_tier("Gold", manifest, dry_run)
    stats["counts"]["trial_sim_from_diamond"] = move_trial_sim_from_tier("Diamond", manifest, dry_run)
    stats["counts"]["trial_sim_from_silver"] = move_trial_sim_from_tier("Silver", manifest, dry_run)
    stats["counts"]["trial_sim_from_bronze"] = move_trial_sim_from_tier("Bronze", manifest, dry_run)
    if include_learning_material:
        stats["counts"]["trial_sim_from_learning"] = move_trial_sim_from_learning_material(
            manifest, dry_run, limit=learning_limit
        )
    stats["counts"]["gold_flatten"] = flatten_gold_subfolders(manifest, dry_run)

    total = sum(stats["counts"].values())
    stats["total_moves"] = total
    return stats


def main() -> int:
    import argparse
    p = argparse.ArgumentParser(description="Restore misplaced warehouse data to canonical folders.")
    p.add_argument("--dry-run", action="store_true", help="Only report moves; do not move files.")
    p.add_argument("--no-learning-material", action="store_true", help="Do not move TRIAL_SIM_* from Learning_Material.")
    p.add_argument("--learning-limit", type=int, default=5000, help="Max TRIAL_SIM files to move from Learning_Material (default 5000).")
    p.add_argument("--manifest", type=str, default="", help="Write manifest JSON to this path (e.g. PX_LOGS/recovery_manifest.json).")
    args = p.parse_args()

    print("[RECOVER] Restoring misplaced data to canonical paths...")
    print(f"  Recovery run dir: {_RECOVERY_RUN}")
    print(f"  Dry run: {args.dry_run}")

    stats = run_recovery(
        dry_run=args.dry_run,
        include_learning_material=not args.no_learning_material,
        learning_limit=args.learning_limit,
    )

    print(f"\n[RECOVER] Counts: {stats['counts']}")
    print(f"[RECOVER] Total moves: {stats['total_moves']}")

    if args.manifest and not args.dry_run:
        out = Path(args.manifest)
        if not out.is_absolute():
            out = _REPO_ROOT / out
        out.parent.mkdir(parents=True, exist_ok=True)
        with open(out, "w", encoding="utf-8") as f:
            json.dump(stats, f, indent=2)
        print(f"[RECOVER] Manifest written: {out}")

    # Write manifest into recovery run dir too
    if not args.dry_run and stats["total_moves"]:
        manifest_path = _RECOVERY_RUN / "recovery_manifest.json"
        with open(manifest_path, "w", encoding="utf-8") as f:
            json.dump(stats, f, indent=2)
        print(f"[RECOVER] Manifest in run dir: {manifest_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
