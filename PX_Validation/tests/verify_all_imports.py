"""Verify critical layer imports. Run from repo root with conda env."""
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

MODS = [
    ("PX_Engine.operations.OPE", "run_ope"),
    ("PX_Engine.operations.ADMET", "run_admet"),
    ("PX_Engine.operations.GradingEngine", "GradingEngine"),  # .grade_dossier is method
    ("PX_Warehouse.Finalization_Pipeline", "run_finalization"),
    ("PX_Warehouse.warehouse_layout", "ensure_structure"),
    ("PX_System.foundation.ZeusLaws", "run_zeus_gate"),
    ("PX_Executive.PX_Legal_Check", "PX_Legal_Check"),
    ("PX_Executive.Patent_Local_Index", "search"),
]

def main():
    failed = []
    for m, attr in MODS:
        try:
            mod = __import__(m, fromlist=[attr])
            getattr(mod, attr)
            print("OK", m)
        except Exception as e:
            print("FAIL", m, e)
            failed.append((m, str(e)))
    if failed:
        print("\nFailed:", failed)
        return 1
    print("\nAll critical imports OK.")
    return 0

if __name__ == "__main__":
    sys.exit(main())
