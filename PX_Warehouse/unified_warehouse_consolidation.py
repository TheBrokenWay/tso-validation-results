"""
Unified Warehouse Consolidation & Purge
- WorldLine recovery & purge (verification-first; Archive/Legacy safety)
- Tier folders (Diamond, Gold, Silver, Bronze) are FLAT only — no subfolders.
- Dossier_Final is flat; grade routing goes to CommercialAssets/ Diamond|Gold|Silver|Bronze.
- Learning_Material consolidation (CommercialAssets/Learning_Material)
- Operations canonical paths; benchmark/OPE move, zero-randomness audit
"""
import os
import sys
import json
import shutil
import subprocess
from pathlib import Path

# Workspace root (syz)
WORKSPACE = Path(__file__).resolve().parent.parent.parent.parent
WAREHOUSE = WORKSPACE / "PX_Warehouse"
COMMERCIAL = WAREHOUSE / "CommercialAssets"
DOSSIER_FINAL = COMMERCIAL / "Dossier_Final"
LEARNING_MASTER = COMMERCIAL / "Learning_Material"
ARCHIVE_LEGACY = WAREHOUSE / "Archive" / "Legacy"
OPERATIONAL = WAREHOUSE / "Operations"
WORLDLINES = WAREHOUSE / "WorldLines"


def log(msg: str) -> None:
    print(f"[CONSOLIDATE] {msg}")


def ensure_dirs() -> None:
    # Tier folders (Diamond, Gold, Silver, Bronze) are flat only — no subfolders under them
    for path in [
        ARCHIVE_LEGACY,
        OPERATIONAL,
        DOSSIER_FINAL,
        COMMERCIAL / "Diamond",
        COMMERCIAL / "Gold",
        COMMERCIAL / "Silver",
        COMMERCIAL / "Bronze",
    ]:
        path.mkdir(parents=True, exist_ok=True)
        try:
            log(f"Ensured: {path.relative_to(WAREHOUSE)}")
        except ValueError:
            log(f"Ensured: {path}")


# --- Step 1: WorldLine Recovery & Purge (Verification-First) ---
def _has_monetizable_fields(data: dict) -> bool:
    if isinstance(data, dict):
        if data.get("dose_optimization") or data.get("virtual_efficacy"):
            return True
        # Nested: physical_realization may hold results
        pr = data.get("physical_realization") or data.get("physics_snapshot") or {}
        if isinstance(pr, dict) and (pr.get("dose_optimization") or pr.get("virtual_efficacy")):
            return True
    return False


def step1_worldline_recovery_and_purge() -> None:
    log("Step 1: WorldLine recovery & purge (verification-first)")
    if not WORLDLINES.exists():
        log("WorldLines folder does not exist; skip.")
        return
    to_learning = []
    to_archive = []
    for p in WORLDLINES.rglob("*"):
        if not p.is_file():
            continue
        if p.suffix not in (".json", ".worldline"):
            continue
        try:
            with open(p, "r", encoding="utf-8") as f:
                data = json.load(f)
        except Exception as e:
            log(f"  Skip (unreadable): {p.name} - {e}")
            to_archive.append(p)
            continue
        if _has_monetizable_fields(data):
            to_learning.append(p)
        else:
            to_archive.append(p)
    for src in to_learning:
        dest = LEARNING_MASTER / src.name
        if dest.exists():
            dest = LEARNING_MASTER / f"WL_{src.stem}_{src.name}"
        shutil.copy2(src, dest)
        log(f"  -> Learning_Material (re-grade): {src.name}")
    for src in to_archive:
        dest = ARCHIVE_LEGACY / src.name
        if dest.exists():
            dest = ARCHIVE_LEGACY / f"WL_{src.stem}_{src.name}"
        shutil.copy2(src, dest)
        log(f"  -> Archive/Legacy: {src.name}")
    # Delete WorldLines folder
    if WORLDLINES.exists():
        shutil.rmtree(WORLDLINES)
        log("  Deleted root WorldLines folder.")
    # Do not touch 99_WAREHOUSE_ARCHIVE/WorldLines (cold storage)


# --- Step 2: Dossier_Final flat; route by prefix to flat tier folders (Diamond, Gold, Silver, Bronze) ---
def step2_dossier_final_tiers() -> None:
    log("Step 2: Dossier_Final flat; route to CommercialAssets tier folders (flat)")
    if not DOSSIER_FINAL.exists():
        DOSSIER_FINAL.mkdir(parents=True, exist_ok=True)
    # Tier folders are flat: CommercialAssets/Diamond, Gold, Silver, Bronze (no subfolders)
    tier_dest = {
        "GOLD": COMMERCIAL / "Gold",
        "SILVER": COMMERCIAL / "Silver",
        "BRONZE": COMMERCIAL / "Bronze",
        "DIAMOND": COMMERCIAL / "Diamond",
    }
    for d in tier_dest.values():
        d.mkdir(parents=True, exist_ok=True)
    # Move flat DOSSIER_* from Dossier_Final into canonical flat tier folders
    for f in list(DOSSIER_FINAL.glob("*.json")):
        name = f.name.upper()
        dest_dir = None
        if name.startswith("DOSSIER_GOLD_"):
            dest_dir = tier_dest["GOLD"]
        elif name.startswith("DOSSIER_SILVER_"):
            dest_dir = tier_dest["SILVER"]
        elif name.startswith("DOSSIER_BRONZE_"):
            dest_dir = tier_dest["BRONZE"]
        elif name.startswith("DOSSIER_DIAMOND_"):
            dest_dir = tier_dest["DIAMOND"]
        if dest_dir is not None:
            dest = dest_dir / f.name
            if dest.exists():
                dest = dest_dir / f"{f.stem}_moved.{f.suffix}"
            shutil.move(str(f), str(dest))
            log(f"  Dossier_Final -> {dest_dir.name}/{dest.name}")
    # Move any files inside legacy Dossier_Final/GOLD, SILVER, BRONZE subdirs up to flat tier folders
    for grade, dest_dir in tier_dest.items():
        sub = DOSSIER_FINAL / grade
        if sub.exists() and sub.is_dir():
            for f in list(sub.glob("*.json")):
                dest = dest_dir / f.name
                if dest.exists():
                    dest = dest_dir / f"{f.stem}_from_sub.{f.suffix}"
                shutil.move(str(f), str(dest))
                log(f"  Dossier_Final/{grade} -> {dest_dir.name}/{dest.name}")
            try:
                sub.rmdir()
            except OSError:
                pass


# --- Step 3: Learning consolidation ---
def step3_learning_consolidation() -> None:
    log("Step 3: Learning_Material consolidation")
    LEARNING_MASTER.mkdir(parents=True, exist_ok=True)
    # Merge CommercialAssets/Learning_Material into master
    comm_learn = COMMERCIAL / "Learning_Material"
    if comm_learn.exists():
        for f in comm_learn.glob("*"):
            if f.is_file():
                dest = LEARNING_MASTER / f.name
                if dest.exists():
                    dest = LEARNING_MASTER / f"{f.stem}_from_Commercial.{f.suffix}"
                shutil.copy2(str(f), str(dest))
                log(f"  CommercialAssets/Learning_Material -> Learning_Material/{dest.name}")
        for f in comm_learn.glob("*"):
            if f.is_file():
                f.unlink()
        try:
            comm_learn.rmdir()
        except OSError:
            pass
    # Optionally merge Archive snapshot Learning_Material (only copy unique; Archive stays as cold storage)
    archive_dir = WAREHOUSE / "Archive"
    if archive_dir.exists():
        for snap in archive_dir.rglob("Learning_Material"):
            if not snap.is_dir():
                continue
            for f in snap.glob("*.json"):
                dest = LEARNING_MASTER / f.name
                if not dest.exists():
                    shutil.copy2(str(f), str(dest))
                    log(f"  Archive snapshot -> Learning_Material/{f.name}")


# --- Step 4: Operational_Data ---
def step4_operational_data() -> None:
    log("Step 4: Operational_Data (benchmarks + OPE + root .json/.py)")
    OPERATIONAL.mkdir(parents=True, exist_ok=True)
    benchmarks_src = WORKSPACE / "PX_Validation" / "benchmarks"
    if benchmarks_src.exists():
        for py in benchmarks_src.glob("*.py"):
            dest = OPERATIONAL / py.name
            shutil.copy2(str(py), str(dest))
            log(f"  benchmarks/{py.name} -> Operational_Data/")
    ope_src = WORKSPACE / "PX_Engine" / "operations" / "OPE.py"
    if ope_src.exists():
        shutil.copy2(str(ope_src), str(OPERATIONAL / "OPE.py"))
        log("  PX_Engine/operations/OPE.py -> Operational_Data/")
    # Clean Desk: harvest/bridge in Acquisition (enterprise only; no archive/tests/manual)
    root_operational = [
        ("reprocess_candidates.json", WORKSPACE / "PX_Warehouse" / "Feeder"),
        ("harvest_failed_trials.py", WORKSPACE / "PX_Executive" / "Acquisition"),
        ("olympus_bridge_v2.py", WORKSPACE / "PX_Executive" / "Acquisition"),
    ]
    for name, base in root_operational:
        p = base / name
        if p.exists() and p.is_file():
            dest = OPERATIONAL / name
            shutil.copy2(str(p), str(dest))
            log(f"  {base.relative_to(WORKSPACE) if WORKSPACE in p.parents else 'root'}/{name} -> Operational_Data/")


# --- Step 5: Zero-randomness audit ---
def step5_zero_randomness_audit() -> None:
    log("Step 5: Zero-randomness audit (Operational_Data)")
    flagged = []
    for py in OPERATIONAL.glob("*.py"):
        try:
            for line in py.read_text(encoding="utf-8").splitlines():
                stripped = line.strip()
                if stripped.startswith("#"):
                    continue
                if "import random" in stripped or "from random import" in stripped:
                    flagged.append(str(py))
                    break
                if "from numpy" in stripped and "random" in stripped:
                    flagged.append(str(py))
                    break
        except Exception:
            pass
    if flagged:
        for f in flagged:
            log(f"  FLAG (stochastic): {f}")
    else:
        log("  No stochastic imports found in Operational_Data (LAW 10 compliant).")


def main() -> None:
    os.chdir(WORKSPACE)
    ensure_dirs()
    step1_worldline_recovery_and_purge()
    step2_dossier_final_tiers()
    step3_learning_consolidation()
    step4_operational_data()
    step5_zero_randomness_audit()
    log("Consolidation complete. Run stability.py from Operational_Data to verify.")


if __name__ == "__main__":
    main()
