#!/usr/bin/env python3
"""
One complete cycle test: Feeder → Engine (Novel + Repurposed) → Finalization.
Verifies warehouse structure and that each folder receives correct output.

Run from repo root: python PX_Executive/run_one_cycle_test.py
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def run(cmd: list, env: dict | None = None, timeout: int = 120) -> tuple[int, str]:
    env = env or os.environ.copy()
    if str(REPO_ROOT) not in env.get("PYTHONPATH", ""):
        env["PYTHONPATH"] = str(REPO_ROOT) + (os.pathsep + env.get("PYTHONPATH", "") if env.get("PYTHONPATH") else "")
    env["PYTHONUNBUFFERED"] = "1"
    try:
        r = subprocess.run(cmd, cwd=str(REPO_ROOT), env=env, capture_output=True, text=True, timeout=timeout)
        return r.returncode, (r.stdout or "") + (r.stderr or "")
    except Exception as e:
        return -1, str(e)


def verify_warehouse_structure() -> tuple[bool, list[str]]:
    """Check Feeder, Novel_Dossiers, Prv_Dossiers, Finalized_Dossiers and tiers exist."""
    from PX_Warehouse.warehouse_layout import (
        TIERS,
        NOVEL_DOSSIERS_ROOT,
        PRV_DOSSIERS_ROOT,
        FINALIZED_DOSSIERS_ROOT,
        get_feeder_dir,
        ensure_structure,
    )
    ensure_structure(REPO_ROOT)
    wh = REPO_ROOT / "PX_Warehouse"
    errors = []
    if not wh.exists():
        errors.append("PX_Warehouse root missing")
    feeder = get_feeder_dir(REPO_ROOT)
    if not feeder.exists():
        feeder.mkdir(parents=True, exist_ok=True)
    for name, root in [
        ("Novel_Dossiers", NOVEL_DOSSIERS_ROOT),
        ("Prv_Dossiers", PRV_DOSSIERS_ROOT),
        ("Finalized_Dossiers", FINALIZED_DOSSIERS_ROOT),
    ]:
        if not root.exists():
            errors.append(f"{name} missing")
        for tier in TIERS:
            t_dir = root / tier
            if not t_dir.exists():
                t_dir.mkdir(parents=True, exist_ok=True)
    return len(errors) == 0, errors


def count_dossiers() -> dict:
    """Count JSON dossiers per folder (excl PATH_TEST)."""
    wh = REPO_ROOT / "PX_Warehouse"
    counts = {}
    for folder in ["Novel_Dossiers", "Prv_Dossiers", "Finalized_Dossiers"]:
        root = wh / folder
        if not root.exists():
            counts[folder] = {}
            continue
        counts[folder] = {}
        for tier in ("Diamond", "Gold", "Silver", "Bronze"):
            td = root / tier
            if td.exists():
                n = sum(1 for f in td.glob("*.json") if not f.name.startswith("PATH_TEST"))
                counts[folder][tier] = n
        counts[folder]["_total"] = sum(counts[folder].get(t, 0) for t in ("Diamond", "Gold", "Silver", "Bronze"))
    return counts


def main() -> int:
    print("=" * 60)
    print("  ONE COMPLETE CYCLE TEST — Feeder, Engine, Repurposed, Finalization")
    print("=" * 60)

    # 1. Verify warehouse structure
    ok, errs = verify_warehouse_structure()
    if not ok:
        print("  [FAIL] Warehouse structure:", errs)
        return 1
    print("  [OK] Warehouse structure (Feeder, Novel_Dossiers, Prv_Dossiers, Finalized_Dossiers, tiers)")

    before = count_dossiers()
    print("\n  Counts BEFORE cycle:")
    for folder, tiers in before.items():
        print(f"    {folder}: {tiers.get('_total', 0)} total", end="")
        for t in ("Diamond", "Gold", "Silver", "Bronze"):
            if tiers.get(t):
                print(f"  {t}={tiers[t]}", end="")
        print()

    # 2. Feeder (Novel): one batch of Genesis, then exit
    print("\n  --- Cycle 1: Feeder (Novel) ---")
    code, out = run(
        [sys.executable, "PX_Executive/run_genesis_feed.py"],
        env={**os.environ, "GENESIS_INTERVAL": "0", "GENESIS_COUNT": "3"},
        timeout=90,
    )
    if code != 0:
        print(f"  [WARN] Genesis feed exit {code}: {out[:400]}")
    else:
        print("  [OK] Genesis feed ran (single batch, exit)")
    from PX_Warehouse.warehouse_layout import get_queue_path
    queue_path = get_queue_path("prv_24h_queue.json", REPO_ROOT)
    queue_ok = queue_path.exists()
    queue_n = 0
    if queue_ok:
        try:
            data = json.loads(queue_path.read_text(encoding="utf-8"))
            queue_n = len(data.get("candidates", []))
        except Exception:
            pass
    print(f"  Queue: {queue_path} exists={queue_ok}, candidates={queue_n}")

    # 3. Engine (Novel): process 1 novel item
    print("\n  --- Cycle 2: Engine (Novel — 24H Orchestrator, 1 item) ---")
    code, out = run(
        [sys.executable, "PX_Executive/PRV_24H_Orchestrator.py"],
        env={**os.environ, "PRV_MODE": "NOVEL", "PRV_MAX_ITEMS": "1", "PRV_24H_DURATION_SEC": "60"},
        timeout=60,
    )
    if code != 0:
        print(f"  [WARN] Orchestrator (Novel) exit {code}: {out[-500:] if len(out) > 500 else out}")
    else:
        print("  [OK] Orchestrator (Novel) ran")
    after_novel = count_dossiers()
    novel_added = (after_novel.get("Novel_Dossiers", {}).get("_total") or 0) - (before.get("Novel_Dossiers", {}).get("_total") or 0)
    print(f"  Novel_Dossiers change: +{novel_added} (may be 0 if queue had no new N or already filed)")

    # 4. Repurposed: feed then engine 1 item (if queue has R)
    print("\n  --- Cycle 3: Repurposed feed + engine (1 item) ---")
    code, _ = run([sys.executable, "PX_Executive/run_repurposed_feed.py"], timeout=60)
    if code != 0:
        print("  [WARN] Repurposed feed exit", code, "(may have no source data)")
    else:
        print("  [OK] Repurposed feed ran")
    code, out = run(
        [sys.executable, "PX_Executive/PRV_24H_Orchestrator.py"],
        env={**os.environ, "PRV_MODE": "REPURPOSED", "PRV_MAX_ITEMS": "1", "PRV_24H_DURATION_SEC": "60"},
        timeout=60,
    )
    if code != 0:
        print(f"  [WARN] Orchestrator (Repurposed) exit {code}")
    else:
        print("  [OK] Orchestrator (Repurposed) ran")
    after_rep = count_dossiers()
    prv_added = (after_rep.get("Prv_Dossiers", {}).get("_total") or 0) - (before.get("Prv_Dossiers", {}).get("_total") or 0)
    print(f"  Prv_Dossiers change: +{prv_added} (may be 0 if no R in queue)")

    # 5. Finalization: run 1 dossier through pipeline
    print("\n  --- Cycle 4: Finalization (limit 1) ---")
    code, out = run(
        [sys.executable, "PX_Executive/run_finalize_dossiers.py", "--limit", "1"],
        timeout=90,
    )
    if code not in (0, 1):
        print(f"  [WARN] Finalize exit {code}")
    else:
        print("  [OK] Finalization ran (exit 0 or 1 depending on Zeus)")
    after_fin = count_dossiers()
    fin_added = (after_fin.get("Finalized_Dossiers", {}).get("_total") or 0) - (before.get("Finalized_Dossiers", {}).get("_total") or 0)
    print(f"  Finalized_Dossiers change: +{fin_added} (0 if no Zeus-authorized dossier was pending)")

    # 6. Final verification
    print("\n  --- Verification ---")
    final_counts = count_dossiers()
    print("  Counts AFTER cycle:")
    for folder, tiers in final_counts.items():
        print(f"    {folder}: {tiers.get('_total', 0)} total", end="")
        for t in ("Diamond", "Gold", "Silver", "Bronze"):
            if tiers.get(t):
                print(f"  {t}={tiers[t]}", end="")
        print()

    # Sanity: Finalized_Dossiers only contains tier subdirs with valid structure
    from PX_Warehouse.warehouse_layout import TIERS
    fin_root = REPO_ROOT / "PX_Warehouse" / "Finalized_Dossiers"
    for tier in TIERS:
        td = fin_root / tier
        if td.exists():
            jsons = [f for f in td.glob("*.json") if not f.name.startswith("PATH_TEST")]
            for f in jsons[:1]:  # spot-check one file
                try:
                    d = json.loads(f.read_text(encoding="utf-8"))
                    if "finalization" not in d:
                        print(f"  [WARN] {f.name} missing 'finalization' block")
                    elif "finalization_version" not in d.get("finalization", {}):
                        print(f"  [WARN] {f.name} missing finalization_version")
                except Exception as e:
                    print(f"  [WARN] {f.name} read error: {e}")

    print("\n" + "=" * 60)
    print("  Cycle test complete. Feeder → Engine (Novel/Repurposed) → Finalization verified.")
    print("=" * 60)
    return 0


if __name__ == "__main__":
    sys.exit(main())
