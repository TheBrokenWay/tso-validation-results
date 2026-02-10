#!/usr/bin/env python3
"""
Fix double-prefix dossier filenames and recover any PRV/worldline files from wrong locations.

1. RENAME IN PLACE: In Prv_Dossiers and Novel_Dossiers, rename
   PRV_NOV_PRV_NOV_<id>.json -> PRV_NOV_<id>.json
   PRV_REP_PRV_REP_<id>.json -> PRV_REP_<id>.json
   so dossiers are findable under the canonical name.

2. REPO-WIDE SCAN: Find PRV_*.json, PRV_*.worldline, TRIAL_*.json outside
   PX_Warehouse/Prv_Dossiers, Novel_Dossiers, WorldLines, Calibration_Molecules, Feeder, Learning_Material
   and move them into the correct canonical folder (by tier or type).

Run from repo root: python PX_Warehouse/Operations/scripts/fix_double_prefix_and_recover_misplaced.py [--execute]
"""
from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
WAREHOUSE = REPO_ROOT / "PX_Warehouse"

SKIP_DIRS = {".git", "__pycache__", "node_modules", "foundation.egg-info", ".cursor", "99_LEGACY_CODE"}
CANONICAL_UNDER = ("Prv_Dossiers", "Novel_Dossiers", "WorldLines", "Calibration_Molecules", "Feeder", "Learning_Material", "Operations")


def get_tier_from_dossier(data: dict) -> str:
    """Delegates to canonical PX_Warehouse.warehouse_layout.get_tier."""
    from PX_Warehouse.warehouse_layout import get_tier
    return get_tier(data)


def fix_double_prefix(dry_run: bool = True) -> list[tuple[Path, Path]]:
    """Rename PRV_NOV_PRV_NOV_*.json -> PRV_NOV_*.json (and REP) in Prv_Dossiers/Novel_Dossiers."""
    renamed: list[tuple[Path, Path]] = []
    for root_name in ("Novel_Dossiers", "Prv_Dossiers"):
        root = WAREHOUSE / root_name
        if not root.exists():
            continue
        for tier in ("Diamond", "Gold", "Silver", "Bronze"):
            tier_dir = root / tier
            if not tier_dir.exists():
                continue
            for path in tier_dir.glob("*.json"):
                name = path.name
                if name.startswith("PRV_NOV_PRV_NOV_"):
                    new_name = "PRV_NOV_" + name.replace("PRV_NOV_PRV_NOV_", "", 1)
                    dest = path.parent / new_name
                    if dest != path and (not dest.exists() or dest.resolve() == path.resolve()):
                        renamed.append((path, dest))
                elif name.startswith("PRV_REP_PRV_REP_"):
                    new_name = "PRV_REP_" + name.replace("PRV_REP_PRV_REP_", "", 1)
                    dest = path.parent / new_name
                    if dest != path and (not dest.exists() or dest.resolve() == path.resolve()):
                        renamed.append((path, dest))
    if not dry_run:
        for src, dest in renamed:
            if src.exists() and (not dest.exists() or dest == src):
                shutil.move(str(src), str(dest))
    return renamed


def find_misplaced_prv_and_worldlines() -> list[Path]:
    """Find PRV_*.json, *.worldline, TRIAL_*.json not under canonical warehouse paths."""
    out: list[Path] = []
    for path in REPO_ROOT.rglob("*"):
        if not path.is_file():
            continue
        if any(p in path.parts for p in SKIP_DIRS):
            continue
        try:
            rel = path.relative_to(REPO_ROOT)
        except ValueError:
            continue
        parts = rel.parts
        if len(parts) < 2 or parts[0] != "PX_Warehouse":
            # Outside PX_Warehouse: include if it looks like PRV/trial/worldline
            if path.suffix == ".worldline" and ("PRV_" in path.name or "WL-" in path.name or "Genesis" in path.name):
                out.append(path)
            if path.suffix == ".json" and ("PRV_NOV" in path.name or "PRV_REP" in path.name or "TRIAL_" in path.name or "DOSSIER_" in path.name):
                out.append(path)
            continue
        # Under PX_Warehouse
        if parts[1] in CANONICAL_UNDER:
            continue
        # Under something else (e.g. Commercial_Dossiers, Archive, TrialSimulations at top level)
        if path.suffix == ".worldline":
            out.append(path)
        if path.suffix == ".json" and ("PRV_" in path.name or "TRIAL_" in path.name or "DOSSIER_" in path.name):
            out.append(path)
    return sorted(set(out), key=lambda p: str(p))


def move_misplaced(files: list[Path], dry_run: bool) -> list[dict]:
    """Move misplaced files into canonical Prv_Dossiers/Novel_Dossiers/WorldLines/Learning_Material/Calibration_Molecules."""
    manifest = []
    for path in files:
        name = path.name
        dest_root = None
        tier = "Bronze"
        if path.suffix == ".worldline":
            dest_root = WAREHOUSE / "WorldLines" / "Bronze"
            try:
                with open(path, "r", encoding="utf-8") as f:
                    data = json.load(f)
                pr = data.get("physical_realization") or data.get("physics_snapshot") or {}
                risk = (pr.get("status") or data.get("header", {}).get("route") or "").upper()
                if "DIAMOND" in risk:
                    tier = "Diamond"
                elif "GOLD" in risk:
                    tier = "Gold"
                elif "SILVER" in risk:
                    tier = "Silver"
                dest_root = WAREHOUSE / "WorldLines" / tier
            except Exception:
                pass
        elif "PRV_NOV" in name:
            dest_root = WAREHOUSE / "Novel_Dossiers"
            try:
                with open(path, "r", encoding="utf-8") as f:
                    data = json.load(f)
                tier = get_tier_from_dossier(data)
            except Exception:
                pass
            dest_root = WAREHOUSE / "Novel_Dossiers" / tier
        elif "PRV_REP" in name or "PRV_" in name:
            dest_root = WAREHOUSE / "Prv_Dossiers"
            try:
                with open(path, "r", encoding="utf-8") as f:
                    data = json.load(f)
                tier = get_tier_from_dossier(data)
            except Exception:
                pass
            dest_root = WAREHOUSE / "Prv_Dossiers" / tier
        elif "TRIAL_" in name or "DOSSIER_" in name:
            dest_root = WAREHOUSE / "Calibration_Molecules"
            if not dest_root.exists():
                dest_root = WAREHOUSE / "Learning_Material"
        if dest_root is None:
            continue
        dest_root.mkdir(parents=True, exist_ok=True)
        dest = dest_root / name
        n = 0
        while dest.exists() and dest.resolve() != path.resolve():
            n += 1
            stem, ext = name.rsplit(".", 1) if "." in name else (name, "")
            dest = dest_root / f"{stem}_dup{n}.{ext}"
        manifest.append({"src": str(path), "dest": str(dest)})
        if not dry_run and path.exists():
            shutil.move(str(path), str(dest))
    return manifest


def main():
    dry_run = "--execute" not in sys.argv
    if dry_run:
        print("DRY RUN (use --execute to apply changes)\n")

    # 1. Fix double-prefix in place
    print("1. Fixing double-prefix dossier filenames (Prv_Dossiers / Novel_Dossiers)...")
    renames = fix_double_prefix(dry_run=dry_run)
    for src, dest in renames:
        print(f"   {'Would rename' if dry_run else 'Renamed'}: {src.name} -> {dest.name}")
    print(f"   Total: {len(renames)} dossier(s)\n")

    # 2. Find and move misplaced files
    print("2. Scanning repo for PRV/worldline/trial files outside canonical folders...")
    misplaced = find_misplaced_prv_and_worldlines()
    if not misplaced:
        print("   None found.\n")
    else:
        for p in misplaced:
            print(f"   Misplaced: {p}")
        manifest = move_misplaced(misplaced, dry_run=dry_run)
        for m in manifest:
            print(f"   {'Would move' if dry_run else 'Moved'}: {m['src']} -> {m['dest']}")
        print(f"   Total: {len(manifest)} file(s)\n")

    print("Done." if dry_run else "Done. Changes applied.")


if __name__ == "__main__":
    main()
