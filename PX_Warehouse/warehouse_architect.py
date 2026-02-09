#!/usr/bin/env python3
"""
Warehouse Architect — One-time renovation to strict tiered structure.

Creates: Prv_Dossiers, Novel_Dossiers, Learning_Material (each with Diamond/Gold/Silver/Bronze
         where applicable). Moves PRV_REP → Prv_Dossiers/<tier>, PRV_NOV → Novel_Dossiers/<tier>,
         TRIAL_* → Learning_Material. Then removes old folders (Commercial_Dossiers,
         CommercialAssets, TrialSimulations, Archive).

Run from repo root: python PX_Warehouse/Operations/scripts/warehouse_architect.py
"""
from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

# Repo root (foundation)
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
ROOT = REPO_ROOT / "PX_Warehouse"

NEW_STRUCTURE = {
    "Prv_Dossiers": ["Diamond", "Gold", "Silver", "Bronze"],
    "Novel_Dossiers": ["Diamond", "Gold", "Silver", "Bronze"],
    "Learning_Material": [],
}


def get_tier(data: dict) -> str:
    try:
        physics = data.get("physics", {}) or data.get("engines", {}).get("ope", {}) or {}
        tox = physics.get("toxicity_index", 1.0)
        if isinstance(data.get("engines"), dict) and isinstance(data["engines"].get("admet"), dict):
            admet_tox = (data["engines"]["admet"].get("toxicity") or {}).get("toxicity_index")
            if admet_tox is not None:
                tox = float(admet_tox)
        if "candidate_profile" in data and data["candidate_profile"]:
            cp_tox = data["candidate_profile"].get("toxicity_index")
            if cp_tox is not None:
                tox = float(cp_tox)
        safety_margin = float(data.get("safety_margin", 0))
        if tox < 0.01 or (tox > 0.02 and safety_margin > 50.0):
            return "Diamond"
        if tox < 0.021:
            return "Gold"
        if tox < 0.10:
            return "Silver"
        return "Bronze"
    except Exception:
        return "Bronze"


def main() -> int:
    print("CONSTRUCTING NEW WAREHOUSE LAYOUT...")

    for folder, subfolders in NEW_STRUCTURE.items():
        main_path = ROOT / folder
        main_path.mkdir(parents=True, exist_ok=True)
        for sub in subfolders:
            (main_path / sub).mkdir(exist_ok=True)

    print("TRANSFERRING INVENTORY...")
    count = 0
    skip_parts = ("Prv_Dossiers", "Novel_Dossiers", "Learning_Material", "Operations", "Inputs", "WorldLines")
    for file_path in ROOT.rglob("*.json"):
        if any(x in file_path.parts for x in skip_parts):
            continue
        fname = file_path.name
        dest_root = ""
        if "PRV_NOV" in fname:
            dest_root = "Novel_Dossiers"
        elif "PRV_REP" in fname or ("PRV_" in fname and "PRV_NOV" not in fname):
            dest_root = "Prv_Dossiers"
        elif "TRIAL" in fname or "trial" in fname.lower():
            dest_root = "Learning_Material"
        else:
            continue

        tier = ""
        if dest_root == "Learning_Material":
            dest_path = ROOT / dest_root / fname
        elif dest_root != "Learning_Material":
            try:
                with open(file_path, "r", encoding="utf-8") as f:
                    data = json.load(f)
                tier = get_tier(data)
            except Exception:
                tier = "Bronze"
            dest_path = ROOT / dest_root / tier / fname
        # Dedup: avoid overwriting
        n = 0
        while dest_path.exists() and dest_path.resolve() != file_path.resolve():
            n += 1
            stem, ext = (fname.rsplit(".", 1)[0], "." + fname.rsplit(".", 1)[1]) if "." in fname else (fname, "")
            dest_path = (ROOT / dest_root / tier / f"{stem}_dup{n}{ext}") if tier else (ROOT / dest_root / f"{stem}_dup{n}{ext}")
        try:
            shutil.move(str(file_path), str(dest_path))
            count += 1
            print(f"   -> Moved {fname} to {dest_root}/{tier or '.'}")
        except Exception as e:
            print(f"   Failed to move {fname}: {e}")

    print(f"MOVED {count} ASSETS.")

    print("DEMOLISHING OLD STRUCTURES...")
    # Do NOT add "Prv_Dossiers" / "PRV_Dossiers": that is our new folder (case-collides on Windows).
    folders_to_erase = [
        "Commercial_Dossiers",
        "CommercialAssets",
        "TrialSimulations",
        "Archive",
    ]
    for folder in folders_to_erase:
        target = ROOT / folder
        if target.exists():
            try:
                shutil.rmtree(target)
                print(f"   -> Demolished: {folder}")
            except Exception as e:
                print(f"   Could not fully remove {folder}: {e}")

    print("\nWAREHOUSE RENOVATION COMPLETE. Use Prv_Dossiers, Novel_Dossiers, Learning_Material.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
