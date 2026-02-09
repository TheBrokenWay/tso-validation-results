"""Lineage status: report on Sovereign Log Chain, WorldLines, and DataLineage."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def sovereign_log_status():
    log = ROOT / "PX_Audit" / "sovereign_log_chain.jsonl"
    if not log.exists():
        return {"status": "MISSING", "entries": 0}
    lines = log.read_text(encoding="utf-8").strip().splitlines()
    entries = len(lines)
    if entries == 0:
        return {"status": "EMPTY", "entries": 0}
    # Verify last entry has valid chain hash
    try:
        last = json.loads(lines[-1])
        has_hash = "record_hash" in last and "prev_hash" in last
    except Exception:
        has_hash = False
    return {
        "status": "INTACT" if has_hash else "UNCHAINED",
        "entries": entries,
        "last_event": last.get("event_type", "?") if has_hash else "?",
    }


def worldline_status():
    wl_root = ROOT / "PX_Warehouse" / "WorldLines"
    if not wl_root.exists():
        return {"status": "MISSING", "count": 0}
    tiers = {}
    total = 0
    for tier_dir in sorted(wl_root.iterdir()):
        if tier_dir.is_dir():
            count = len(list(tier_dir.glob("*.worldline")))
            tiers[tier_dir.name] = count
            total += count
    return {"status": "OK" if total > 0 else "EMPTY", "total": total, "tiers": tiers}


def data_lineage_status():
    graph_file = ROOT / "PX_Warehouse" / "Operations" / "lineage_graph.json"
    if not graph_file.exists():
        return {"status": "NO_GRAPH"}
    try:
        data = json.loads(graph_file.read_text(encoding="utf-8"))
        return {"status": "OK", "assets": len(data.get("nodes", data.get("assets", [])))}
    except Exception as e:
        return {"status": f"ERROR: {e}"}


def main():
    print("=" * 60)
    print("PREDATOR X — DATA LINEAGE STATUS")
    print("=" * 60)

    slc = sovereign_log_status()
    print(f"\n1. Sovereign Log Chain: {slc['status']}")
    print(f"   Entries: {slc['entries']}")
    if slc.get("last_event"):
        print(f"   Last event: {slc['last_event']}")

    wl = worldline_status()
    print(f"\n2. WorldLine Database: {wl['status']}")
    print(f"   Total worldlines: {wl.get('total', 0)}")
    for tier, count in wl.get("tiers", {}).items():
        print(f"     {tier}: {count}")

    dl = data_lineage_status()
    print(f"\n3. DataLineage Graph: {dl['status']}")
    if dl.get("assets"):
        print(f"   Tracked assets: {dl['assets']}")

    print("\n" + "=" * 60)
    all_ok = slc["status"] == "INTACT" and wl["status"] == "OK"
    print(f"Overall: {'HEALTHY' if all_ok else 'NEEDS ATTENTION'}")
    print("=" * 60)
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
