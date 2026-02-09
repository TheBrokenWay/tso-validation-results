#!/usr/bin/env python3
"""
Warehouse Architect — Phase 2: Feeder & Calibration_Molecules.

Adds two zones:
  Feeder              — Loading dock for raw queues (prv_24h_queue.json, reprocess_candidates.json).
  Calibration_Molecules — QA vault for reference standards (TRIAL_SIMULATION*.json).

Flow: Feeder (In) → Factory (Process) → Dossiers (Out) / Calibration_Molecules (Test).

Run from repo root: python PX_Warehouse/Operations/scripts/warehouse_architect_phase2.py
"""
from __future__ import annotations

import shutil
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
ROOT = REPO_ROOT / "PX_Warehouse"

NEW_ZONES = {
    "Feeder": [],
    "Calibration_Molecules": [],
}


def main() -> int:
    print("EXPANDING WAREHOUSE: PHASE 2...")

    # 1. Build new zones
    for folder in NEW_ZONES:
        (ROOT / folder).mkdir(parents=True, exist_ok=True)
        print(f"   -> Built Zone: {folder}")

    # 2. Relocate calibration standards (TRIAL_SIMULATION)
    print("SECURING CALIBRATION STANDARDS...")
    calibration_count = 0
    for source in [ROOT, ROOT / "Learning_Material"]:
        if not source.exists():
            continue
        for file_path in source.glob("*.json"):
            if "TRIAL_SIMULATION" in file_path.name:
                dest = ROOT / "Calibration_Molecules" / file_path.name
                try:
                    shutil.move(str(file_path), str(dest))
                    calibration_count += 1
                    print(f"   -> Calibrated: {file_path.name}")
                except Exception as e:
                    print(f"   Error moving {file_path.name}: {e}")

    # 3. Organize feeder queues
    print("ORGANIZING FEEDER QUEUES...")
    feeder_count = 0
    for file_path in ROOT.rglob("*queue*.json"):
        if "Feeder" in str(file_path):
            continue
        if "Inputs" in str(file_path.parent):
            # Copy active queue into Feeder (do not move so engine still sees it until code cutover)
            dest = ROOT / "Feeder" / file_path.name
            try:
                shutil.copy(str(file_path), str(dest))
                feeder_count += 1
                print(f"   -> Copied to Feeder: {file_path.name}")
            except Exception as e:
                print(f"   Error copying {file_path.name}: {e}")
        else:
            dest = ROOT / "Feeder" / file_path.name
            try:
                shutil.move(str(file_path), str(dest))
                feeder_count += 1
                print(f"   -> Stored in Feeder: {file_path.name}")
            except Exception:
                pass

    print("\nPHASE 2 COMPLETE.")
    print(f"   - Calibration Standards: {calibration_count} items secured.")
    print(f"   - Feeder Lists: {feeder_count} items organized.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
