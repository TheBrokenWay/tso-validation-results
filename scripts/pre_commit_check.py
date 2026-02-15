"""Pre-commit: lightweight warehouse structure check.

Only scans top-level PX_Warehouse folders (~instant).
Full warehouse simulation available via: just warehouse-audit
"""
import os
import sys
from pathlib import Path

root = Path(__file__).resolve().parent.parent
warehouse = root / "PX_Warehouse"

# Canonical top-level folders — what SHOULD exist per §1 directive
CANONICAL_TOP = {
    "CommercialAssets", "WorldLines", "TrialSimulations", "Operations", "Archive",
    "Feeder", "Calibration_Molecules", "Prv_Dossiers", "Novel_Dossiers",
    "Finalized_Dossiers", "Learning_Material", "placement_gate",
    "Discovery_Accepted",
}
# Known legacy dirs — flagged but don't block commits
KNOWN_LEGACY_TOP = {
    "Archive_Novel", "Archive_Primary", "Backup_Pre_Refinery",
    "00_COMMERCIAL_DOSSIERS",
}
# Infrastructure dirs — auto-generated, not a governance concern
INFRA_DIRS = {"__pycache__", "PX_LOGS", "PX_Warehouse"}

if not warehouse.exists():
    sys.exit(0)

top_dirs = sorted(d.name for d in warehouse.iterdir() if d.is_dir())
unknown = [d for d in top_dirs if d not in CANONICAL_TOP and d not in KNOWN_LEGACY_TOP and d not in INFRA_DIRS]

if unknown:
    print("[PRE-COMMIT] FAIL: Unknown top-level warehouse folders detected:")
    for d in unknown:
        print(f"  - {d}")
    print("Fix before committing. Run 'just warehouse-audit' for full analysis.")
    sys.exit(1)

legacy = [d for d in top_dirs if d in KNOWN_LEGACY_TOP]
if legacy:
    print(f"[PRE-COMMIT] WARN: {len(legacy)} known-legacy folders pending consolidation")

sys.exit(0)
