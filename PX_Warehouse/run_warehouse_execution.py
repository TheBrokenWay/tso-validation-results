#!/usr/bin/env python
"""
PHASE 2 — EXECUTION MODE (after simulation approval)
Applies MOVE, RENAME, ARCHIVE, DELETE per simulation; updates script references;
runs verification and generates WAREHOUSE_DRIFT_REPORT_<timestamp>.md.
"""
from pathlib import Path
from datetime import datetime, timezone
import shutil
import sys

REPO_ROOT = Path(__file__).resolve().parents[1]  # PX_Warehouse -> foundation
PX_WAREHOUSE = REPO_ROOT / "PX_Warehouse"
PX_LOGS = REPO_ROOT / "PX_LOGS"

# Root-level items to apply (from Phase 1 simulation)
def apply_phase2():
    # 1. Ensure canonical structure
    try:
        from bootstrap_canonical_warehouse import ensure_dirs
        ensure_dirs()
    except Exception as e:
        print(f"[Phase 2] Bootstrap warning: {e}")

    # 2. Create Archive/Legacy and PX_LOGS/archive/warehouse_legacy
    (PX_WAREHOUSE / "Archive" / "Legacy").mkdir(parents=True, exist_ok=True)
    (PX_LOGS / "archive" / "warehouse_legacy").mkdir(parents=True, exist_ok=True)

    # 3. MOVE root-level Learning_Material → CommercialAssets/Learning_Material
    src_learn = PX_WAREHOUSE / "Learning_Material"
    dest_learn = PX_WAREHOUSE / "CommercialAssets" / "Learning_Material"
    if src_learn.exists() and src_learn.is_dir():
        if dest_learn.exists():
            # Merge: move contents into existing
            for item in src_learn.iterdir():
                dest_item = dest_learn / item.name
                if dest_item.exists():
                    dest_item = dest_learn / f"merged_{item.name}"
                shutil.move(str(item), str(dest_item))
            shutil.rmtree(src_learn)
            print("[Phase 2] Merged Learning_Material -> CommercialAssets/Learning_Material")
        else:
            shutil.move(str(src_learn), str(dest_learn))
            print("[Phase 2] Moved Learning_Material -> CommercialAssets/Learning_Material")

    # 4. ARCHIVE root-level logs → PX_LOGS/archive/warehouse_legacy
    src_logs = PX_WAREHOUSE / "logs"
    if src_logs.exists():
        dest_logs = PX_LOGS / "archive" / "warehouse_legacy" / "logs"
        dest_logs.mkdir(parents=True, exist_ok=True)
        for item in src_logs.iterdir():
            shutil.move(str(item), str(dest_logs / item.name))
        src_logs.rmdir()
        print("[Phase 2] Archived logs -> PX_LOGS/archive/warehouse_legacy/")

    # 5. DELETE root-level __pycache__
    pycache = PX_WAREHOUSE / "__pycache__"
    if pycache.exists():
        shutil.rmtree(pycache)
        print("[Phase 2] Deleted PX_Warehouse/__pycache__")

    # 6. ARCHIVE root-level legacy files to Archive/Legacy (do NOT archive WorldLine_Database.py — canonical)
    legacy_root_files = [
        "FINAL_WAREHOUSE_INVENTORY.md",
        "WAREHOUSE_CONSOLIDATION_REPORT.md",
    ]
    archive_legacy = PX_WAREHOUSE / "Archive" / "Legacy"
    for name in legacy_root_files:
        p = PX_WAREHOUSE / name
        if p.exists():
            shutil.move(str(p), str(archive_legacy / name))
            print(f"[Phase 2] Archived {name} -> Archive/Legacy/")

    # 7. Update script references (canonical paths)
    update_script_references()
    return 0


def update_script_references():
    """Update every script that references legacy warehouse paths to canonical structure."""
    # olympus_bridge_v2.py: PX_Warehouse/Learning_Material -> PX_Warehouse/CommercialAssets/Learning_Material
    p = REPO_ROOT / "PX_Executive" / "Acquisition" / "olympus_bridge_v2.py"
    if p.exists():
        text = p.read_text(encoding="utf-8")
        old = 'LEARNING_MATERIAL_DIR = Path("PX_Warehouse") / "Learning_Material"'
        new = 'LEARNING_MATERIAL_DIR = Path("PX_Warehouse") / "CommercialAssets" / "Learning_Material"'
        if old in text and new not in text:
            text = text.replace(old, new)
            p.write_text(text, encoding="utf-8")
            print("[Phase 2] Updated PX_Executive/Acquisition/olympus_bridge_v2.py")

    # PX_Warehouse/Operations/scripts/unified_warehouse_consolidation.py
    p = PX_WAREHOUSE / "Operations" / "scripts" / "unified_warehouse_consolidation.py"
    if p.exists():
        text = p.read_text(encoding="utf-8")
        replacements = [
            ('LEARNING_MASTER = WAREHOUSE / "Learning_Material"', 'LEARNING_MASTER = COMMERCIAL / "Learning_Material"'),
            ('ARCHIVE_LEGACY = WAREHOUSE / "Archive_Legacy_2026"', 'ARCHIVE_LEGACY = WAREHOUSE / "Archive" / "Legacy"'),
            ('OPERATIONAL = WAREHOUSE / "Operational_Data"', 'OPERATIONAL = WAREHOUSE / "Operations"'),
        ]
        for old, new in replacements:
            if old in text:
                text = text.replace(old, new)
        # ensure_dirs: Archive_Legacy_2026 -> Archive/Legacy, Operational_Data -> Operations
        text = text.replace('("Archive_Legacy_2026", ARCHIVE_LEGACY)', '("Archive/Legacy", ARCHIVE_LEGACY)')
        text = text.replace('("Operational_Data", OPERATIONAL)', '("Operations", OPERATIONAL)')
        text = text.replace('log(f"  -> Archive_Legacy_2026:', 'log(f"  -> Archive/Legacy:')
        p.write_text(text, encoding="utf-8")
        print("[Phase 2] Updated PX_Warehouse/Operations/scripts/unified_warehouse_consolidation.py")

    # Root-level unified_warehouse_consolidation.py was moved to Archive/Legacy; if a copy exists at repo root, update it
    p_root = REPO_ROOT / "PX_Warehouse" / "unified_warehouse_consolidation.py"
    if not p_root.exists():
        p_ops = PX_WAREHOUSE / "Operations" / "scripts" / "unified_warehouse_consolidation.py"
        if p_ops.exists():
            pass  # already updated above


def run_verification_and_drift_report():
    """Generate WAREHOUSE_DRIFT_REPORT_<timestamp>.md per WAREHOUSE_DRIFT_REPORT_EXECUTION_LOGIC.md."""
    ts = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
    drift_dir = PX_LOGS / "drift_reports"
    drift_dir.mkdir(parents=True, exist_ok=True)
    report_path = drift_dir / f"WAREHOUSE_DRIFT_REPORT_{ts}.md"
    canonical_top = {"CommercialAssets", "WorldLines", "TrialSimulations", "Operations", "Archive", "ENTERPRISE_MODEL_README.md"}
    structural = []
    for item in PX_WAREHOUSE.iterdir():
        if item.name.startswith("."):
            continue
        if item.is_dir() and item.name not in canonical_top and item.name != "Archive":
            structural.append(f"- {item.name} — Not part of canonical warehouse structure")
        elif item.is_file() and item.name not in ("ENTERPRISE_MODEL_README.md",):
            structural.append(f"- {item.name} — Unexpected file at warehouse root")

    drift_status = "PASS" if not structural else "FAIL"
    lines = [
        "# WAREHOUSE DRIFT REPORT (Post–Phase 2 Verification)",
        "",
        f"**Report ID:** drift_report_{ts}",
        f"**Timestamp (UTC):** {datetime.now(timezone.utc).isoformat()}",
        f"**Warehouse Path:** PX_Warehouse/",
        f"**Drift Status:** {drift_status}",
        "",
        "## 1. Summary of Drift",
        "",
        f"Total structural issues at warehouse root: **{len(structural)}**",
        "",
        "## 2. Structural Drift (Folder-Level)",
        "",
    ]
    if structural:
        lines.extend(structural)
    else:
        lines.append("None. Canonical top-level structure satisfied.")
    lines.extend([
        "",
        "## 3. Certification Block",
        "",
        f"Warehouse Verified: {'YES' if drift_status == 'PASS' else 'NO'}",
        f"Drift Eliminated: {'YES' if drift_status == 'PASS' else 'NO'}",
        "Signed By: Phase 2 Execution Script",
        f"Timestamp: {ts}",
        "",
    ])
    report_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"[Phase 2] Drift report written: {report_path}")
    return drift_status


if __name__ == "__main__":
    print("[Phase 2] EXECUTION MODE — applying simulation actions and script updates...")
    rc = apply_phase2()
    if rc == 0:
        print("[Phase 2] Running post-execution verification...")
        status = run_verification_and_drift_report()
        print(f"[Phase 2] Verification: {status}")
    sys.exit(rc)
