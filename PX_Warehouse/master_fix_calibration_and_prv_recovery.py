#!/usr/bin/env python3
"""
Master Fix Script:
  1. Calibration Correction — Move TRIAL_* files from Learning_Material to Calibration_Molecules.
  2. Deep Search — Scan entire warehouse (including Archive, Legacy, hidden) for PRV_REP files;
     move them to Prv_Dossiers/<tier>.

Run from repo root: python PX_Warehouse/Operations/scripts/master_fix_calibration_and_prv_recovery.py
"""
from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

_SCRIPT_DIR = Path(__file__).resolve().parent
_REPO_ROOT = _SCRIPT_DIR.parents[1]
ROOT = _REPO_ROOT / "PX_Warehouse"
CALIBRATION_VAULT = ROOT / "Calibration_Molecules"
PRV_VAULT = ROOT / "Prv_Dossiers"
LEARNING_BIN = ROOT / "Learning_Material"


def get_prv_tier(file_path: Path) -> str:
    """Read dossier to determine tier (Diamond/Gold/Silver/Bronze)."""
    try:
        with open(file_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        physics = data.get("physics") or data.get("engines", {}).get("ope", {}) or {}
        tox = physics.get("toxicity_index", 1.0) if isinstance(physics, dict) else 1.0
        if tox == 1.0 and isinstance(data.get("engines"), dict):
            admet = (data["engines"].get("admet") or {}).get("toxicity") or {}
            if isinstance(admet, dict) and admet.get("toxicity_index") is not None:
                tox = float(admet["toxicity_index"])
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
    print("STARTING EMERGENCY FIX & RECOVERY...")

    # 1. CALIBRATION: Move TRIAL_* from Learning_Material to Calibration_Molecules
    CALIBRATION_VAULT.mkdir(parents=True, exist_ok=True)
    moved_cal = 0
    if LEARNING_BIN.exists():
        for file_path in list(LEARNING_BIN.glob("*.json")):
            if "TRIAL" in file_path.name or "trial" in file_path.name.lower():
                dest = CALIBRATION_VAULT / file_path.name
                if dest.exists() and dest.resolve() != file_path.resolve():
                    stem, ext = file_path.stem, file_path.suffix
                    n = 1
                    while dest.exists():
                        dest = CALIBRATION_VAULT / f"{stem}_dup{n}{ext}"
                        n += 1
                shutil.move(str(file_path), str(dest))
                moved_cal += 1
                print(f"   -> FIXED: Moved {file_path.name} to Calibration_Molecules")

    # 2. DEEP SEARCH: Find PRV_REP anywhere in warehouse (Archive, Legacy, etc.)
    print("\nHUNTING FOR LOST PRV_REP FILES...")
    recovered_prv = 0
    for file_path in list(ROOT.rglob("*.json")):
        if not file_path.is_file():
            continue
        if "Prv_Dossiers" in file_path.parts:
            continue
        if "PRV_REP" not in file_path.name:
            continue
        tier = get_prv_tier(file_path)
        dest_dir = PRV_VAULT / tier
        dest_dir.mkdir(parents=True, exist_ok=True)
        dest = dest_dir / file_path.name
        if dest.exists() and dest.resolve() != file_path.resolve():
            stem, ext = file_path.stem, file_path.suffix
            n = 1
            while dest.exists():
                dest = dest_dir / f"{stem}_dup{n}{ext}"
                n += 1
        try:
            shutil.move(str(file_path), str(dest))
            recovered_prv += 1
            print(f"   -> RECOVERED: {file_path.name} to Prv_Dossiers/{tier}")
        except Exception as e:
            print(f"   Failed to move {file_path.name}: {e}")

    print("\nOPERATION COMPLETE.")
    print(f"   - Calibration Fixed: {moved_cal} files relocated.")
    print(f"   - Repurposed Recovered: {recovered_prv} files found.")
    if recovered_prv == 0:
        print("\nNO LOST PRV_REP FILES FOUND. Run 'run_prv_repurposed.py' to regenerate.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
