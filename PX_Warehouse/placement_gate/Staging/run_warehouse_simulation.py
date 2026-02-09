#!/usr/bin/env python
"""
PHASE 1 — SIMULATION MODE (NO CHANGES)
Full recursive scan of PX_Warehouse; classify; tier; lifecycle; propose actions;
identify script references; generate Simulation Report. Does not modify any files.
"""
from pathlib import Path
from datetime import datetime, timezone
import re
import json
import sys

REPO_ROOT = Path(__file__).resolve().parents[1]  # PX_Warehouse -> foundation
PX_WAREHOUSE = REPO_ROOT / "PX_Warehouse"

# Canonical top-level folders — what SHOULD exist per §1 directive
CANONICAL_TOP = {
    # Business zones
    "CommercialAssets", "WorldLines", "TrialSimulations", "Operations", "Archive",
    # Placement-gate canonical zones
    "Feeder", "Calibration_Molecules", "Prv_Dossiers", "Novel_Dossiers",
    "Finalized_Dossiers", "Learning_Material", "placement_gate",
}
# Known legacy dirs — exist, need consolidation into Archive/ or deletion.
# Flagged in simulation reports but do NOT block --enforce commits.
KNOWN_LEGACY_TOP = {
    "Archive_Novel", "Archive_Primary", "Backup_Pre_Refinery",
    "00_COMMERCIAL_DOSSIERS",
}
# Infrastructure dirs — auto-generated, not a governance concern.
# __pycache__ proposed for deletion; PX_LOGS and PX_Warehouse (self-ref) are misplaced.
INFRA_DIRS = {"__pycache__", "PX_LOGS", "PX_Warehouse"}
# Files (not dirs) that are expected at root
CANONICAL_TOP_FILES = {"ENTERPRISE_MODEL_README.md"}
CANONICAL_COMMERCIAL = {"Diamond", "Gold", "Silver", "Bronze", "Dossier_Final", "Executive_Summary", "Audit_Trail", "Learning_Material"}
CANONICAL_WORLDLINES = {"Diamond", "Gold", "Silver", "Bronze"}
CANONICAL_TRIAL = {"LiveRuns"}
CANONICAL_OPS = {"Control_Scripts", "Inputs", "manifests", "reports", "scripts", "System_Health"}
CANONICAL_ARCHIVE_SUB = {"Dossiers", "TrialSimulations", "WorldLines", "CommercialAssets"}  # + Snapshot_*

# Legacy path segments (classify as LEGACY only when path contains these as folder/file names)
# Do NOT match canonical tier names: Diamond, Gold, Silver, Bronze as single folder names
LEGACY_SEGMENTS = {
    "__pycache__", "logs", "Operational_Data", "WorldLine_Records", "Archive_Legacy_2026",
    "v3_audit", "monetization", "stratification", "deduplication", "benchmark", "dashboard",
    "EXECUTIVE_BENCHMARK", "PORTFOLIO_DASHBOARD", "GOLD_CORE", "DIAMOND_PRV",
    "FINAL_WAREHOUSE_INVENTORY", "WAREHOUSE_CONSOLIDATION_REPORT", "WorldLine_Database.py",
    "unified_warehouse_consolidation.py",
}
# Combined set for classification (canonical + known legacy + infra + files)
CANONICAL_TOP_ALL = CANONICAL_TOP | KNOWN_LEGACY_TOP | INFRA_DIRS | CANONICAL_TOP_FILES

# Tier from folder name
TIER_NAMES = {"Diamond", "Gold", "Silver", "Bronze"}
RENAME_MAP = {"GOLD_CORE": "Gold", "DIAMOND_PRV": "Diamond", "SILVER": "Silver", "BRONZE": "Bronze"}


def scan_warehouse():
    """Recursively scan PX_Warehouse; return list of (path, is_dir, rel_path)."""
    if not PX_WAREHOUSE.exists():
        return []
    out = []
    for p in PX_WAREHOUSE.rglob("*"):
        try:
            rel = p.relative_to(PX_WAREHOUSE)
        except ValueError:
            continue
        if ".git" in rel.parts:
            continue
        out.append((p, p.is_dir(), str(rel)))
    return out


def classify_path(rel_path: str, is_dir: bool) -> str:
    """Classify as ACTIVE, REFERENCED, VALIDATED, DOSSIER, LEGACY, DUPLICATE, UNUSED, UNKNOWN."""
    parts = rel_path.replace("\\", "/").split("/")
    top = parts[0] if parts else ""

    if "__pycache__" in rel_path or rel_path.endswith(".pyc"):
        return "LEGACY"
    if top == "logs":
        return "LEGACY"
    # Any path segment that is a known legacy segment (not canonical tier name)
    for seg in parts:
        if seg in LEGACY_SEGMENTS:
            return "LEGACY"
    if top == "Operations" and "scripts" in parts:
        return "ACTIVE"
    if top == "Operations":
        return "ACTIVE"
    if top == "CommercialAssets":
        if "Dossier_Final" in parts:
            return "DOSSIER"
        return "REFERENCED"
    if top == "WorldLines":
        return "VALIDATED"
    if top == "Archive":
        return "DOSSIER"
    if top == "TrialSimulations":
        return "REFERENCED"
    if top in CANONICAL_TOP:
        return "REFERENCED"
    if top in KNOWN_LEGACY_TOP:
        return "LEGACY"
    if top in INFRA_DIRS:
        return "LEGACY"
    # Non-canonical top-level folder — truly unknown
    if top and top not in CANONICAL_TOP_FILES:
        return "LEGACY"
    return "UNKNOWN"


def tier_from_path(rel_path: str) -> str:
    """Extract tier from path; default Bronze."""
    parts = rel_path.replace("\\", "/").split("/")
    for p in parts:
        if p in TIER_NAMES:
            return p
    if "GOLD" in rel_path.upper() or "GOLD_CORE" in rel_path:
        return "Gold"
    if "DIAMOND" in rel_path.upper() or "DIAMOND_PRV" in rel_path:
        return "Diamond"
    if "SILVER" in rel_path.upper():
        return "Silver"
    if "BRONZE" in rel_path.upper():
        return "Bronze"
    return "Bronze"


def lifecycle_from_path(rel_path: str) -> str:
    """Created | Validated | Dossier."""
    if rel_path.startswith("Archive"):
        return "Dossier"
    if rel_path.startswith("WorldLines"):
        return "Validated"
    if rel_path.startswith("CommercialAssets"):
        if "Dossier_Final" in rel_path:
            return "Dossier"
        return "Created"
    return "Unknown"


def proposed_action(classification: str, rel_path: str, is_dir: bool) -> tuple:
    """(action, reason, destination_or_None)."""
    parts = rel_path.replace("\\", "/").split("/")
    top = parts[0] if parts else ""
    # Top-level Learning_Material → MOVE to canonical CommercialAssets/Learning_Material
    if top == "Learning_Material":
        return ("MOVE", "Canonical location per §1", "PX_Warehouse/CommercialAssets/Learning_Material/")
    if classification == "LEGACY":
        if "__pycache__" in rel_path or rel_path.endswith(".pyc"):
            return ("DELETE", "Cache", None)
        if top == "logs" or "logs" in rel_path or "benchmark" in rel_path or "dashboard" in rel_path:
            return ("ARCHIVE", "Legacy analytics/logs", "PX_LOGS/archive/warehouse_legacy/")
        # Root-level legacy files (WorldLine_Database.py, etc.) and __pycache__ folder
        return ("ARCHIVE", "Legacy folder/file", "PX_Warehouse/Archive/Legacy/")
    if classification == "UNKNOWN" and rel_path not in ("", "ENTERPRISE_MODEL_README.md"):
        return ("REVIEW", "Unclassified", None)
    if "__pycache__" in rel_path:
        return ("DELETE", "Cache", None)
    # Non-canonical folder name → RENAME
    if parts and parts[0] == "CommercialAssets" and len(parts) >= 2:
        sub = parts[1]
        if sub in RENAME_MAP:
            return ("RENAME", "Canonical tier name", f"CommercialAssets/{RENAME_MAP[sub]}")
    return ("KEEP", "Canonical or active", None)


def find_script_references():
    """Find all .py files that reference warehouse paths."""
    refs = []
    for py in REPO_ROOT.rglob("*.py"):
        if "__pycache__" in str(py) or ".git" in str(py):
            continue
        try:
            text = py.read_text(encoding="utf-8", errors="ignore")
        except Exception:
            continue
        rel = py.relative_to(REPO_ROOT)
        if "PX_Warehouse" in text or "CommercialAssets" in text or "WorldLines" in text:
            refs.append(str(rel))
        elif "TrialSimulations" in text or "Operations/Inputs" in text:
            refs.append(str(rel))
    return sorted(set(refs))


def main():
    ts = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
    sim_dir = REPO_ROOT / "PX_LOGS" / "warehouse_simulation"
    sim_dir.mkdir(parents=True, exist_ok=True)
    report_path = sim_dir / f"PX_Warehouse_SIMULATION_REPORT_{ts}.md"

    print("[Phase 1] Scanning PX_Warehouse (read-only)...")
    items = scan_warehouse()
    print(f"  Total items: {len(items)}")

    print("[Phase 1] Classifying and proposing actions...")
    by_class = {}
    actions = []
    folders_non_canonical = []
    script_refs = find_script_references()

    for p, is_dir, rel in items:
        classification = classify_path(rel, is_dir)
        by_class[classification] = by_class.get(classification, 0) + 1
        tier = tier_from_path(rel)
        lifecycle = lifecycle_from_path(rel)
        action, reason, dest = proposed_action(classification, rel, is_dir)
        if action != "KEEP":
            actions.append({"path": rel, "is_dir": is_dir, "classification": classification, "tier": tier, "lifecycle": lifecycle, "action": action, "reason": reason, "dest": dest})
        # Non-canonical top-level folder
        top = rel.split("/")[0].split("\\")[0]
        if is_dir and rel.count("/") + rel.count("\\") == 0 and top not in CANONICAL_TOP:
            folders_non_canonical.append(rel)

    # Action breakdown by type
    by_action = {}
    for a in actions:
        by_action[a["action"]] = by_action.get(a["action"], 0) + 1

    # Build report
    lines = [
        "# WAREHOUSE SIMULATION REPORT (PHASE 1 — NO CHANGES)",
        "",
        f"**Generated:** {datetime.now(timezone.utc).isoformat()}",
        f"**Warehouse:** {PX_WAREHOUSE}",
        "",
        "## 1. Summary",
        "",
        f"- Total items scanned: **{len(items)}**",
        f"- Classifications: {dict(by_class)}",
        f"- Proposed non-KEEP actions: **{len(actions)}**",
        f"- Action breakdown: {dict(by_action)}",
        f"- Non-canonical top-level folders: **{len(folders_non_canonical)}**",
        f"- Scripts referencing warehouse: **{len(script_refs)}**",
        "",
        "## 2. Non-Canonical Top-Level Folders",
        "",
    ]
    for f in sorted(folders_non_canonical):
        lines.append(f"- `{f}` — not in §1 canonical structure")
    lines.extend(["", "## 3. Proposed Actions (MOVE / RENAME / ARCHIVE / DELETE / REVIEW)", ""])
    for a in actions[:500]:  # Cap at 500
        lines.append(f"- **{a['action']}** `{a['path']}` — {a['reason']}" + (f" → `{a['dest']}`" if a.get("dest") else ""))
    if len(actions) > 500:
        lines.append(f"- ... and {len(actions) - 500} more (see full list in script output)")
    lines.extend(["", "## 4. Scripts / Systems to Update (System-Wide Reference Update)", ""])
    for r in script_refs:
        lines.append(f"- `{r}`")
    lines.extend([
        "",
        "## 5. Tier Assignment (from path)",
        "",
        "Tiers used: Diamond, Gold, Silver, Bronze (default Bronze when unknown).",
        "",
        "## 6. Lifecycle",
        "",
        "Created → CommercialAssets; Validated → WorldLines; Dossier → Archive.",
        "",
        "## 7. Certification",
        "",
        "**Phase 1 complete. No files or folders were modified.**",
        "After approval, run Phase 2 (execution) to apply MOVE, RENAME, ARCHIVE, DELETE, TIER, LIFECYCLE, and script updates.",
        "",
    ])

    report_path.write_text("\n".join(lines), encoding="utf-8")
    # Write machine-readable manifest for Phase 2 (summary + scripts; full actions recomputed in Phase 2)
    manifest_path = sim_dir / f"PX_Warehouse_SIMULATION_MANIFEST_{ts}.json"
    manifest = {
        "timestamp": ts,
        "warehouse": str(PX_WAREHOUSE),
        "total_items": len(items),
        "action_breakdown": by_action,
        "non_canonical_top_folders": folders_non_canonical,
        "scripts_to_update": script_refs,
        "action_sample": [{"path": a["path"], "action": a["action"], "dest": a.get("dest")} for a in actions[:200]],
    }
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"[Phase 1] Report written: {report_path}")
    print(f"[Phase 1] Manifest written: {manifest_path}")
    print(f"  Classifications: {by_class}")
    print(f"  Proposed actions: {len(actions)}")
    print(f"  Scripts to update: {len(script_refs)}")

    # --enforce: fail only on truly unknown folders (used by pre-commit hook)
    # Known legacy and infra dirs are flagged but don't block commits.
    enforce = "--enforce" in sys.argv or "-e" in sys.argv
    if enforce:
        unknown_folders = [
            f for f in folders_non_canonical
            if f not in KNOWN_LEGACY_TOP and f not in INFRA_DIRS
        ]
        legacy_folders = [
            f for f in folders_non_canonical
            if f in KNOWN_LEGACY_TOP
        ]
        infra_folders = [
            f for f in folders_non_canonical
            if f in INFRA_DIRS
        ]
        if unknown_folders:
            print("[ENFORCE] FAIL: Unknown top-level folders detected. Fix before committing.")
            for f in unknown_folders:
                print(f"  - {f}")
            return 1
        if legacy_folders:
            print(f"[ENFORCE] WARN: {len(legacy_folders)} known-legacy folders pending consolidation:")
            for f in legacy_folders:
                print(f"  - {f} (queued for Archive/)")
        if infra_folders:
            print(f"[ENFORCE] WARN: {len(infra_folders)} infra folders to clean up:")
            for f in infra_folders:
                print(f"  - {f}")
        print("[ENFORCE] PASS: No unknown structural drift at warehouse root.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
