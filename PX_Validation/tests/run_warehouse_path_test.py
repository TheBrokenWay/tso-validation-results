"""
Drop one test document into each warehouse classification folder to verify paths.
Run from repo root: python tests/run_warehouse_path_test.py
"""
import os
import json
from pathlib import Path
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[2]
WAREHOUSE = ROOT / "PX_Warehouse"

# Canonical classification folders (align with warehouse_layout)
CLASSIFICATION_FOLDERS = [
    "Feeder",
    "Calibration_Molecules",
    "Prv_Dossiers",
    "Prv_Dossiers/Diamond",
    "Prv_Dossiers/Gold",
    "Prv_Dossiers/Silver",
    "Prv_Dossiers/Bronze",
    "Novel_Dossiers",
    "Novel_Dossiers/Diamond",
    "Novel_Dossiers/Gold",
    "Novel_Dossiers/Silver",
    "Novel_Dossiers/Bronze",
    "Learning_Material",
    "WorldLines",
    "WorldLines/Diamond",
    "WorldLines/Gold",
    "WorldLines/Silver",
    "WorldLines/Bronze",
    "Operations",
    "Operations/Inputs",
]

TEST_DOC = {
    "path_test": True,
    "purpose": "Visual verification that this classification path is writable",
    "timestamp_utc": datetime.now(timezone.utc).isoformat(),
    "script": "tests/run_warehouse_path_test.py",
    "governance": "PATH_VERIFICATION_ONLY",
}


def main():
    print("PX_Warehouse path test — dropping one test document per classification folder\n")
    print(f"Warehouse root: {WAREHOUSE}\n")
    results = []
    for rel in CLASSIFICATION_FOLDERS:
        folder = WAREHOUSE / rel
        folder.mkdir(parents=True, exist_ok=True)
        doc = {**TEST_DOC, "folder": rel}
        filename = "PATH_TEST_2026-02-05.json"
        filepath = folder / filename
        with open(filepath, "w", encoding="utf-8") as f:
            json.dump(doc, f, indent=2)
        results.append((rel, str(filepath)))
        print(f"  OK  {rel}")
    print(f"\nDone. Wrote {len(results)} test documents.")
    print("\nYou can open each folder and look for PATH_TEST_2026-02-05.json to confirm paths work.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
