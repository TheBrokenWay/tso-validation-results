#!/usr/bin/env python
"""
Reclassify non-PRV-eligible dossiers as Learning Material.

Moves dossiers with prv_eligible=false (or missing) from active warehouse dirs
into PX_Warehouse/Learning_Material/<source>/<tier>/, preserving tier structure.

Usage:
    python PX_Executive/reclassify_learning_material.py --dry-run   # preview
    python PX_Executive/reclassify_learning_material.py              # execute
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO_ROOT))

WAREHOUSE = _REPO_ROOT / "PX_Warehouse"
LEARNING = WAREHOUSE / "Learning_Material"

# Source dirs to scan: (label, path, is_finalized)
SCAN_SOURCES = [
    ("Novel", WAREHOUSE / "Novel_Dossiers", False),
    ("Prv", WAREHOUSE / "Prv_Dossiers", False),
    ("Finalized", WAREHOUSE / "Finalized_Dossiers", True),
]

TIERS = ("Diamond", "Gold", "Silver", "Bronze", "Discovery_Accepted")


def is_prv_eligible(filepath: Path) -> bool:
    """Return True only if dossier explicitly has prv_eligible=true."""
    try:
        with open(filepath, "r", encoding="utf-8") as f:
            data = json.load(f)
        return data.get("prv_eligible", False) is True
    except (json.JSONDecodeError, OSError):
        return False


def ensure_learning_structure() -> None:
    """Create Learning_Material subdirs for all sources and tiers."""
    for label in ("Novel", "Prv", "Finalized"):
        for tier in TIERS:
            (LEARNING / label / tier).mkdir(parents=True, exist_ok=True)
            if label == "Finalized":
                (LEARNING / label / tier / "packages").mkdir(exist_ok=True)


def migrate(dry_run: bool = True) -> dict:
    """Scan and move non-PRV-eligible dossiers. Returns summary dict."""
    ensure_learning_structure()

    stats = {"moved": 0, "skipped_eligible": 0, "skipped_non_json": 0, "errors": 0}
    details: list[str] = []

    for label, source_root, is_finalized in SCAN_SOURCES:
        for tier in TIERS:
            tier_dir = source_root / tier
            if not tier_dir.is_dir():
                continue

            dest_dir = LEARNING / label / tier

            # Move dossier JSON files
            for filepath in sorted(tier_dir.glob("*.json")):
                if is_prv_eligible(filepath):
                    stats["skipped_eligible"] += 1
                    details.append(f"  KEEP  {filepath.relative_to(WAREHOUSE)}")
                    continue

                dest = dest_dir / filepath.name
                if dry_run:
                    details.append(f"  MOVE  {filepath.relative_to(WAREHOUSE)} -> {dest.relative_to(WAREHOUSE)}")
                else:
                    try:
                        shutil.move(str(filepath), str(dest))
                        details.append(f"  MOVED {filepath.relative_to(WAREHOUSE)} -> {dest.relative_to(WAREHOUSE)}")
                    except OSError as e:
                        details.append(f"  ERROR {filepath.relative_to(WAREHOUSE)}: {e}")
                        stats["errors"] += 1
                        continue
                stats["moved"] += 1

                # Move matching pharma package if finalized
                if is_finalized:
                    pkg_name = filepath.stem + "_pharma_package.json"
                    pkg_src = tier_dir / "packages" / pkg_name
                    if pkg_src.exists():
                        pkg_dest = dest_dir / "packages" / pkg_name
                        if dry_run:
                            details.append(f"  MOVE  {pkg_src.relative_to(WAREHOUSE)} -> {pkg_dest.relative_to(WAREHOUSE)}")
                        else:
                            try:
                                shutil.move(str(pkg_src), str(pkg_dest))
                            except OSError as e:
                                details.append(f"  ERROR pkg {pkg_src.relative_to(WAREHOUSE)}: {e}")
                                stats["errors"] += 1

    return {"stats": stats, "details": details}


def main() -> None:
    parser = argparse.ArgumentParser(description="Reclassify non-PRV-eligible dossiers as Learning Material")
    parser.add_argument("--dry-run", action="store_true", help="Preview moves without executing")
    args = parser.parse_args()

    mode = "DRY RUN" if args.dry_run else "EXECUTING"
    print(f"\n=== Reclassify Learning Material ({mode}) ===\n")

    result = migrate(dry_run=args.dry_run)
    stats = result["stats"]

    # Print details (abbreviated for large runs)
    details = result["details"]
    if len(details) <= 50:
        for line in details:
            print(line)
    else:
        for line in details[:20]:
            print(line)
        print(f"  ... ({len(details) - 40} more lines) ...")
        for line in details[-20:]:
            print(line)

    print(f"\n--- Summary ---")
    print(f"  Moved:             {stats['moved']}")
    print(f"  Kept (eligible):   {stats['skipped_eligible']}")
    print(f"  Errors:            {stats['errors']}")

    if args.dry_run:
        print(f"\n  (dry run — no files were moved)")
    else:
        print(f"\n  Migration complete.")

    if stats["errors"] > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()
