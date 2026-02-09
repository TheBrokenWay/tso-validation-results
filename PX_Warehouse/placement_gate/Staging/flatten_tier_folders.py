"""
Flatten tier folders: Diamond, Gold, Silver, Bronze must have no subfolders.
Moves contents of Gold/DIAMOND_PRV -> Diamond/, Gold/GOLD_CORE -> Gold/,
Silver/SILVER -> Silver/, Silver/BRONZE -> Bronze/, then removes empty subdirs.
Run from repo root: python PX_Warehouse/Operations/scripts/flatten_tier_folders.py
"""
import shutil
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
COMMERCIAL = REPO_ROOT / "PX_Warehouse" / "CommercialAssets"

# (parent_tier_folder, subfolder_name) -> canonical flat tier folder for contents
MOVES = [
    (COMMERCIAL / "Gold", "DIAMOND_PRV", COMMERCIAL / "Diamond"),
    (COMMERCIAL / "Gold", "GOLD_CORE", COMMERCIAL / "Gold"),
    (COMMERCIAL / "Silver", "SILVER", COMMERCIAL / "Silver"),
    (COMMERCIAL / "Silver", "BRONZE", COMMERCIAL / "Bronze"),
]


def main():
    print("[flatten_tier_folders] Ensuring Diamond, Gold, Silver, Bronze are flat (no subfolders).\n")
    for parent, subname, dest_dir in MOVES:
        sub = parent / subname
        if not sub.exists() or not sub.is_dir():
            continue
        dest_dir.mkdir(parents=True, exist_ok=True)
        for f in list(sub.glob("*")):
            if f.is_file():
                dest = dest_dir / f.name
                if dest.exists():
                    dest = dest_dir / f"{f.stem}_flattened{f.suffix}"
                shutil.move(str(f), str(dest))
                print(f"  {sub.relative_to(COMMERCIAL)}/{f.name} -> {dest_dir.name}/")
        try:
            sub.rmdir()
            print(f"  Removed empty subfolder: {sub.relative_to(COMMERCIAL)}")
        except OSError:
            print(f"  (subfolder not empty or already removed: {sub})")
    print("\n[flatten_tier_folders] Done.")


if __name__ == "__main__":
    main()
