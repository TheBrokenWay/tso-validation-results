#!/usr/bin/env python3
"""
Scan the entire repo for PRV_ and TRIAL_ JSON files; move them into the canonical
warehouse structure (Prv_Dossiers / Novel_Dossiers / Learning_Material) with tiering.

Use after warehouse architect renovation or anytime dossiers may have been written
outside PX_Warehouse (e.g. PX_LOGS, Nipah_Analysis, tests).

Run from repo root: python PX_Warehouse/Operations/scripts/recover_prv_repo_wide.py [--execute]
"""
from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
WAREHOUSE = REPO_ROOT / "PX_Warehouse"

# Skip these when scanning
SKIP_DIRS = {".git", "__pycache__", "node_modules", "foundation.egg-info", ".cursor"}
# Already-sorted warehouse paths (skip so we don't re-move)
CANONICAL_PARTS = ("Prv_Dossiers", "Novel_Dossiers", "Learning_Material")


def get_tier(data: dict) -> str:
    """Delegates to canonical PX_Warehouse.warehouse_layout.get_tier."""
    from PX_Warehouse.warehouse_layout import get_tier as _canonical
    return _canonical(data)


def ensure_warehouse_dirs() -> None:
    for folder, tiers in [("Prv_Dossiers", ["Diamond", "Gold", "Silver", "Bronze"]),
                          ("Novel_Dossiers", ["Diamond", "Gold", "Silver", "Bronze"]),
                          ("Learning_Material", [])]:
        (WAREHOUSE / folder).mkdir(parents=True, exist_ok=True)
        for t in tiers:
            (WAREHOUSE / folder / t).mkdir(exist_ok=True)


def find_prv_trial_jsons() -> list[Path]:
    out: list[Path] = []
    for path in REPO_ROOT.rglob("*.json"):
        if not path.is_file():
            continue
        if any(part in path.parts for part in SKIP_DIRS):
            continue
        if any(part in path.parts for part in CANONICAL_PARTS):
            # Already in canonical warehouse
            continue
        name = path.name
        if "PRV_" in name or "TRIAL" in name or "trial" in name.lower():
            out.append(path)
    return sorted(out, key=lambda p: str(p))


def run(dry_run: bool = True) -> dict:
    ensure_warehouse_dirs()
    candidates = find_prv_trial_jsons()
    manifest = []
    moved = 0
    for file_path in candidates:
        fname = file_path.name
        if "PRV_NOV" in fname:
            dest_root = "Novel_Dossiers"
        elif "PRV_REP" in fname or ("PRV_" in fname and "PRV_NOV" not in fname):
            dest_root = "Prv_Dossiers"
        elif "TRIAL" in fname or "trial" in fname.lower():
            dest_root = "Learning_Material"
        else:
            continue

        if dest_root == "Learning_Material":
            tier = ""
            dest_path = WAREHOUSE / dest_root / fname
        else:
            try:
                with open(file_path, "r", encoding="utf-8") as f:
                    data = json.load(f)
                tier = get_tier(data)
            except Exception:
                tier = "Bronze"
            dest_path = WAREHOUSE / dest_root / tier / fname

        n = 0
        while dest_path.exists() and dest_path.resolve() != file_path.resolve():
            n += 1
            stem, ext = (fname.rsplit(".", 1)[0], "." + fname.rsplit(".", 1)[1]) if "." in fname else (fname, "")
            dest_path = (WAREHOUSE / dest_root / tier / f"{stem}_dup{n}{ext}") if tier else (WAREHOUSE / dest_root / f"{stem}_dup{n}{ext}")

        rel_src = file_path.relative_to(REPO_ROOT) if REPO_ROOT in file_path.parents else file_path
        manifest.append({"src": str(rel_src), "dest": str(dest_path.relative_to(REPO_ROOT)), "tier": tier or "."})
        if not dry_run:
            shutil.move(str(file_path), str(dest_path))
            moved += 1
        else:
            moved += 1
    return {"dry_run": dry_run, "candidates": len(candidates), "moved": moved, "manifest": manifest}


def main() -> int:
    import argparse
    p = argparse.ArgumentParser(description="Recover PRV_ and TRIAL_ JSON from entire repo into warehouse.")
    p.add_argument("--execute", action="store_true", help="Perform moves (default is dry-run).")
    args = p.parse_args()
    dry_run = not args.execute
    if dry_run:
        print("DRY RUN (use --execute to move files)\n")
    print(f"Scanning repo: {REPO_ROOT}")
    result = run(dry_run=dry_run)
    print(f"Found {result['candidates']} PRV_/TRIAL_ JSON file(s) outside canonical warehouse.")
    for m in result["manifest"]:
        print(f"   -> {m['src']} -> {m['dest']}")
    print(f"\nTotal to move: {result['moved']}")
    if not dry_run and result["manifest"]:
        out = WAREHOUSE / "Operations" / "Inputs" / "recover_prv_manifest.json"
        out.parent.mkdir(parents=True, exist_ok=True)
        with open(out, "w", encoding="utf-8") as f:
            json.dump(result, f, indent=2)
        print(f"Manifest written: {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
