#!/usr/bin/env python
"""
bootstrap_canonical_warehouse.py

Creates and verifies the canonical PX_Warehouse structure for the Enterprise Edition.
Idempotent: safe to run multiple times.

Aligned with:
- PX_WAREHOUSE_ENTERPRISE_MODEL.md
- PX_WAREHOUSE_MIGRATION_INSTRUCTION.md
- PX_ENTERPRISE_LAYOUT.md
"""

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]  # PX_Warehouse -> foundation
PX_WAREHOUSE = REPO_ROOT / "PX_Warehouse"

# Tier folders (Diamond, Gold, Silver, Bronze) are FLAT only — no subfolders (no GOLD_CORE, DIAMOND_PRV, Silver/SILVER, etc.)
CANONICAL_DIRS = [
    # Top-level warehouse organs
    PX_WAREHOUSE / "CommercialAssets",
    PX_WAREHOUSE / "WorldLines",
    PX_WAREHOUSE / "TrialSimulations",
    PX_WAREHOUSE / "Operations",
    PX_WAREHOUSE / "Archive",

    # CommercialAssets tiers + surfaces
    PX_WAREHOUSE / "CommercialAssets" / "Diamond",
    PX_WAREHOUSE / "CommercialAssets" / "Gold",
    PX_WAREHOUSE / "CommercialAssets" / "Silver",
    PX_WAREHOUSE / "CommercialAssets" / "Bronze",
    PX_WAREHOUSE / "CommercialAssets" / "Dossier_Final",
    PX_WAREHOUSE / "CommercialAssets" / "Executive_Summary",
    PX_WAREHOUSE / "CommercialAssets" / "Audit_Trail",
    PX_WAREHOUSE / "CommercialAssets" / "Learning_Material",

    # WorldLines tiers
    PX_WAREHOUSE / "WorldLines" / "Diamond",
    PX_WAREHOUSE / "WorldLines" / "Gold",
    PX_WAREHOUSE / "WorldLines" / "Silver",
    PX_WAREHOUSE / "WorldLines" / "Bronze",

    # TrialSimulations
    PX_WAREHOUSE / "TrialSimulations" / "LiveRuns",

    # Operations
    PX_WAREHOUSE / "Operations" / "Control_Scripts",
    PX_WAREHOUSE / "Operations" / "Inputs",
    PX_WAREHOUSE / "Operations" / "manifests",
    PX_WAREHOUSE / "Operations" / "reports",
    PX_WAREHOUSE / "Operations" / "scripts",
    PX_WAREHOUSE / "Operations" / "System_Health",

    # Archive subtrees
    PX_WAREHOUSE / "Archive" / "Dossiers",
    PX_WAREHOUSE / "Archive" / "TrialSimulations",
    PX_WAREHOUSE / "Archive" / "WorldLines",
    PX_WAREHOUSE / "Archive" / "CommercialAssets",
]


def ensure_dirs():
    """Create all canonical directories that do not exist. Returns list of created paths."""
    created = []
    for d in CANONICAL_DIRS:
        if not d.exists():
            d.mkdir(parents=True, exist_ok=True)
            created.append(d)
    return created


def main():
    print(f"[bootstrap] Repo root: {REPO_ROOT}")
    print(f"[bootstrap] Ensuring canonical PX_Warehouse structure...")

    created = ensure_dirs()

    if created:
        print("[bootstrap] Created directories:")
        for d in created:
            print(f"  - {d.relative_to(REPO_ROOT)}")
    else:
        print("[bootstrap] All canonical directories already exist. No changes made.")

    print("[bootstrap] Canonical warehouse structure is now present.")
    print("[bootstrap] Next steps: run the migration/enforcement logic to move/rename/archive legacy content.")


if __name__ == "__main__":
    main()
