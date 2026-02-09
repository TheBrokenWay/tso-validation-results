#!/usr/bin/env python3
"""
Diagnose why FTO failure rate is ~99.7% and why many candidates are "pending evaluation".

Run from repo root: python PX_Executive/diagnose_fto_failure_rate.py
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def main() -> int:
    print("=" * 60)
    print("  FTO FAILURE RATE DIAGNOSIS")
    print("=" * 60)

    # 1. Queue size ("candidates pending evaluation")
    from PX_Warehouse.warehouse_layout import get_queue_path
    queue_path = get_queue_path("prv_24h_queue.json", REPO_ROOT)
    queue_size = 0
    type_n = 0
    type_r = 0
    if queue_path.exists():
        try:
            data = json.loads(queue_path.read_text(encoding="utf-8"))
            candidates = data.get("candidates", []) if isinstance(data, dict) else []
            queue_size = len(candidates)
            for c in candidates:
                if c.get("type") == "N":
                    type_n += 1
                elif c.get("type") == "R":
                    type_r += 1
        except Exception as e:
            print(f"  Queue read error: {e}")
    print(f"\n  1. QUEUE (candidates pending evaluation)")
    print(f"     Path: {queue_path}")
    print(f"     Total candidates: {queue_size}")
    print(f"     Type N (novel): {type_n}  |  Type R (repurposed): {type_r}")

    # 2. Patent index
    from PX_Executive.Patent_Local_Index import index_exists, DEFAULT_INDEX_PATH, get_last_refresh
    index_path = DEFAULT_INDEX_PATH
    exists = index_exists(index_path)
    last_refresh = get_last_refresh(index_path)
    print(f"\n  2. PATENT INDEX (FTO gate)")
    print(f"     Path: {index_path}")
    print(f"     Exists and has patents: {exists}")
    print(f"     Last refresh: {last_refresh}")

    if not exists:
        print(f"\n  >>> ROOT CAUSE: Patent index is missing or empty.")
        print(f"     All FTO checks return LOCAL_INDEX_EMPTY → 100% block.")
        print(f"     Fix: Run patent_data_refresh.py to populate PX_Data/patents/")
        return 0

    # 3. Sample one repurposed candidate through FTO
    if queue_path.exists() and queue_size > 0:
        try:
            data = json.loads(queue_path.read_text(encoding="utf-8"))
            candidates = data.get("candidates", []) if isinstance(data, dict) else []
            rep = next((c for c in candidates if c.get("type") == "R"), None)
            if rep:
                from PX_Executive.PX_Legal_Check import PX_Legal_Check
                legal = PX_Legal_Check(mode="REGULATORY")
                result = legal.check_compound(
                    smiles=rep.get("smiles") or "CCO",
                    iupac_name=rep.get("name") or rep.get("id", ""),
                    compound_id=rep.get("id", ""),
                    novel_invention=False,
                )
                print(f"\n  3. SAMPLE FTO CHECK (first repurposed candidate)")
                print(f"     ID: {rep.get('id')}  Name: {rep.get('name', '')[:40]}")
                print(f"     Status: {result.status}")
                print(f"     Freedom-to-operate: {result.freedom_to_operate}")
                print(f"     Patent hits: {result.patent_hits}")
                print(f"     Blocking patents (recent): {len(result.blocking_patents)}")
                if result.blocking_patents:
                    print(f"     First blocking: {result.blocking_patents[0].get('patent_number')} ({result.blocking_patents[0].get('date')})")
                if result.status == "BLOCKED" and result.patent_hits > 0:
                    print(f"\n  >>> LIKELY CAUSE: Repurposed compounds are searched by NAME + TEXT.")
                    print(f"     Drug names (e.g. CHEMBL ids or names) match patent title/abstract text.")
                    print(f"     Any patent < 20 years → BLOCKED. With 5000+ patents, most names hit something.")
        except Exception as e:
            print(f"\n  3. Sample FTO check error: {e}")
            import traceback
            traceback.print_exc()

    # 4. Why 99.7%?
    print(f"\n  4. WHY ~99.7% FAILURE RATE?")
    print(f"     - Orchestrator runs: IN → PP → E2E → LG (FTO) → OK or FL.")
    print(f"     - Logs show DONE=0 and every item ending in 'FTO blocked' (STAGE=FL).")
    print(f"     - So 100% of processed items are failing at the FTO (stage_legal) gate.")
    print(f"     - '4605 candidates' = number processed in that run (TOT); all failed FTO.")
    print(f"\n  RECOMMENDED FIXES:")
    print(f"     A) RESEARCH mode (default when PRV_LIVE_RESEARCH is not set) only relaxes")
    print(f"        'no search terms'; repurposed still get name/text search → many blocks.")
    print(f"     B) For repurposed: search ONLY by structure (SMILES/InChIKey), not by name/text,")
    print(f"        so only true structural hits block (change novel_invention logic for R).")
    print(f"     C) Or: run with PRV_LIVE_RESEARCH=0 and add FTO bypass for repurposed in RESEARCH")
    print(f"        (e.g. only enforce FTO when PRV_LIVE_RESEARCH=1).")
    print("=" * 60)
    return 0


if __name__ == "__main__":
    sys.exit(main())
