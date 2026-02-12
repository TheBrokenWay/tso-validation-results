#!/usr/bin/env python3
"""
px_finalize.py — Unified Finalization Pipeline

Scans Prv_Dossiers and Novel_Dossiers, runs full finalization checklist
(trial, evidence, Zeus, grade/tier validation), and places into
Finalized_Dossiers/<tier>/.

Usage:
    python PX_Executive/px_finalize.py                # Process all unfinalized
    python PX_Executive/px_finalize.py --dry-run      # Report only
    python PX_Executive/px_finalize.py --limit 100    # Max to process
    python PX_Executive/px_finalize.py --reprocess    # Re-run already finalized
    python PX_Executive/px_finalize.py -v             # Verbose Zeus failure reasons
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def main() -> int:
    parser = argparse.ArgumentParser(description="Unified finalization pipeline")
    parser.add_argument("--dry-run", action="store_true", help="Report only, no writes")
    parser.add_argument("--limit", type=int, default=0, help="Max dossiers to process (0 = all)")
    parser.add_argument("--reprocess", action="store_true", help="Re-process already finalized dossiers")
    parser.add_argument("--verbose", "-v", action="store_true", help="Print Zeus failure details")
    args = parser.parse_args()

    from PX_Warehouse.warehouse_layout import (
        TIERS,
        get_prv_dossier_dir,
        get_finalized_dossier_dir,
    )
    from PX_Warehouse.Finalization_Pipeline import run_finalization, write_finalized_dossier

    wh = REPO_ROOT / "PX_Warehouse"

    # Collect unfinalized dossiers (Novel first for better backfill ordering)
    to_process: list[tuple[Path, str, bool]] = []
    for is_novel, folder in [(True, "Novel_Dossiers"), (False, "Prv_Dossiers")]:
        for tier in TIERS:
            d_dir = wh / folder / tier
            if not d_dir.exists():
                continue
            for f in d_dir.glob("*.json"):
                if f.name.startswith("PATH_TEST"):
                    continue
                item_id = f.stem
                to_process.append((f, item_id, is_novel))

    # Skip already finalized unless --reprocess
    if not args.reprocess:
        finalized_dir = wh / "Finalized_Dossiers"
        seen_finalized: set[str] = set()
        if finalized_dir.exists():
            for t in TIERS:
                td = finalized_dir / t
                if td.exists():
                    for f in td.glob("*.json"):
                        if not f.name.startswith("PATH_TEST"):
                            seen_finalized.add(f.stem)
        to_process = [(p, iid, nov) for p, iid, nov in to_process if iid not in seen_finalized]

    if args.limit > 0:
        to_process = to_process[:args.limit]

    print("=" * 64)
    print("   PREDATOR X — FINALIZATION PIPELINE")
    print("=" * 64)
    print(f"Dossiers: {len(to_process)} | Dry-run: {args.dry_run} | Reprocess: {args.reprocess}")
    print("-" * 64)

    ok = 0
    fail = 0
    zeus_rejected = 0

    for path, item_id, is_novel in to_process:
        try:
            dossier = json.loads(path.read_text(encoding="utf-8"))
        except Exception as e:
            print(f"  Skip {path.name}: read error {e}")
            _log_failure(item_id, str(e), "dossier JSON read")
            fail += 1
            continue

        finalized, tier, err = run_finalization(dossier, item_id, is_novel, REPO_ROOT)
        if err:
            print(f"  {item_id}: {err}")
            _log_failure(item_id, err, "run_finalization")
            fail += 1
            continue

        # Grade/tier consistency check (defense in depth — Finalization_Pipeline also checks)
        _GRADE_RANK = {"DIAMOND_TIER": 3, "GOLD_TIER": 2, "SILVER_TIER": 1, "BRONZE_TIER": 0,
                       "NEEDS_REVIEW": 1, "REJECTED": -1}
        _TIER_RANK = {"Diamond": 3, "Gold": 2, "Silver": 1, "Bronze": 0}
        _RANK_TO_TIER = {3: "Diamond", 2: "Gold", 1: "Silver", 0: "Bronze", -1: "Bronze"}
        grade = ((finalized or {}).get("finalization") or {}).get("discovery_grading", {}).get("grade", "")
        grade_r = _GRADE_RANK.get(grade, 0)
        tier_r = _TIER_RANK.get(tier, 0)
        if grade_r < tier_r:
            tier = _RANK_TO_TIER.get(grade_r, "Bronze")
            if finalized:
                finalized["finalization"]["tier"] = tier

        # Zeus gate: only WRITE is gated
        zeus_authorized = (finalized.get("finalization") or {}).get("zeus_verdict") or {}
        if not zeus_authorized.get("authorized"):
            zeus_rejected += 1
            if not args.dry_run:
                msg = f"  {item_id}: Zeus not authorized (not written)"
                if args.verbose:
                    rationale = zeus_authorized.get("rationale", "")
                    laws = zeus_authorized.get("laws_results") or {}
                    failed_laws = [k for k, v in laws.items() if isinstance(v, dict) and not v.get("passed")]
                    if failed_laws:
                        msg += f" -- {rationale}; failed: {failed_laws}"
                        for k in failed_laws:
                            r = laws.get(k, {})
                            if isinstance(r, dict) and r.get("reason"):
                                msg += f"\n      {k}: {r['reason']}"
                    else:
                        msg += f" -- {rationale}"
                print(msg)
            fail += 1
            continue

        if args.dry_run:
            print(f"  [dry-run] Would finalize {item_id} -> Finalized_Dossiers/{tier}/")
            ok += 1
            continue

        out_path = write_finalized_dossier(finalized, item_id, tier, REPO_ROOT)
        print(f"  {item_id} -> {out_path}")

        # ── Pharma-grade dossier generation (tier-routed) ──
        _generate_pharma_package(finalized, item_id, tier, REPO_ROOT)

        ok += 1

    print("\n" + "=" * 64)
    print("FINALIZATION SUMMARY")
    print("=" * 64)
    print(f"Finalized:     {ok}")
    print(f"Failed/Skip:   {fail}")
    print(f"Zeus rejected: {zeus_rejected}")
    print("=" * 64)

    return 0 if fail == 0 else 1


def _generate_pharma_package(finalized: dict, item_id: str, tier: str, repo_root: Path) -> None:
    """Generate pharma-grade dossier package based on tier routing."""
    tier_upper = tier.upper()
    if tier_upper in ("SILVER", "BRONZE"):
        return  # WorldLine record only

    try:
        from PX_System.foundation.Evidence_Package import generate_pharma_dossier

        # Extract disease_id from finalization metadata
        fin = finalized.get("finalization", {})
        disease_context = fin.get("disease_context", [])
        indication = fin.get("disease_space_anchoring", {})
        disease_id = ""
        if isinstance(disease_context, list) and disease_context:
            disease_id = str(disease_context[0]).strip().lower().replace(" ", "_")
        if not disease_id:
            disease_id = (finalized.get("candidate") or {}).get("indication", "unknown")
            disease_id = disease_id.strip().lower().replace(" ", "_")

        package = generate_pharma_dossier(
            dossier=finalized,
            tier=tier_upper,
            disease_id=disease_id,
            repo_root=str(repo_root),
        )
        if package:
            # Write package alongside finalized dossier
            pkg_dir = repo_root / "PX_Warehouse" / "Finalized_Dossiers" / tier / "packages"
            pkg_dir.mkdir(parents=True, exist_ok=True)
            pkg_path = pkg_dir / f"{item_id}_pharma_package.json"
            pkg_path.write_text(json.dumps(package, indent=2, default=str), encoding="utf-8")
            print(f"    Pharma package ({tier_upper}) -> {pkg_path.name}")
    except Exception as e:
        print(f"    WARN: Pharma package generation failed for {item_id}: {e}", file=sys.stderr)


def _log_failure(item_id: str, error: str, context: str) -> None:
    """Log finalization failure to audit trail."""
    try:
        from PX_System.finalization_log import log_finalization_failure
        log_finalization_failure(
            source_file="px_finalize.py",
            candidate_id=item_id,
            error=error,
            context=context,
        )
    except Exception:
        pass


if __name__ == "__main__":
    sys.exit(main())
