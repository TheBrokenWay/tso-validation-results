#!/usr/bin/env python3
"""
One-shot migration: replace legacy risk_level strings in warehouse JSON with 4-tier taxonomy.
Does not recalculate Diamond; only renames labels for consistency.
Run from repo root: python PX_Warehouse/Operations/scripts/migrate_legacy_labels.py
"""
from pathlib import Path

_SCRIPT_DIR = Path(__file__).resolve().parent
_REPO_ROOT = _SCRIPT_DIR.parents[1]
ROOT = _REPO_ROOT / "PX_Warehouse"


def main() -> int:
    count = 0
    for fpath in ROOT.rglob("*.json"):
        if "Operations" in fpath.parts:
            continue
        if not fpath.is_file():
            continue
        try:
            txt = fpath.read_text(encoding="utf-8")
            if "HIGH_RISK_TRASH" in txt or "SILVER_TIER_RISK" in txt or "GOLD_TIER_SAFE" in txt:
                txt = txt.replace("HIGH_RISK_TRASH", "TOXICITY_FAILURE")
                txt = txt.replace("SILVER_TIER_RISK", "TOXICITY_SILVER")
                txt = txt.replace("GOLD_TIER_SAFE", "TOXICITY_GOLD")
                fpath.write_text(txt, encoding="utf-8")
                count += 1
                print(f"Updated: {fpath.name}")
            # Also migrate 3-tier labels if present
            elif "TOXICITY_PASS" in txt or "TOXICITY_WARNING" in txt:
                txt = txt.replace("TOXICITY_PASS", "TOXICITY_GOLD")
                txt = txt.replace("TOXICITY_WARNING", "TOXICITY_SILVER")
                fpath.write_text(txt, encoding="utf-8")
                count += 1
                print(f"Updated: {fpath.name}")
        except Exception:
            pass
    print(f"MIGRATION COMPLETE: Updated {count} files.")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
