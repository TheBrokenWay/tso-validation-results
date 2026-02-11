import asyncio, json, os, hashlib, sys
from datetime import datetime, timezone
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from PX_Executive.CSA_Pentarchy import CSAPentarchy
from PX_Audit.AAS_Verification import AASVerification

async def pipeline(prop):
    res = {"authorized": False, "error": None, "trace_id": "N/A"}
    try:
        c_r = await CSAPentarchy().evaluate_action(prop)
        a_r = AASVerification().verify_invariants(prop)
        
        res["trace_id"] = f"FDA-GAIP-{hashlib.sha256(str(prop).encode()).hexdigest()[:12].upper()}"
        res["authorized"] = (c_r["status"] == "COHERENT" and a_r["status"] == "SUCCESS" and (not c_r["human_review_required"] or prop.get("human_sign_off_id")))
        res["csa"] = c_r; res["aas"] = a_r
        
        _trace_dir = str(_REPO_ROOT / "PX_Audit")
        os.makedirs(_trace_dir, exist_ok=True)
        with open(os.path.join(_trace_dir, "fda_gaip_trace.log"), "a") as f:
            f.write(json.dumps({"ts": datetime.now(timezone.utc).isoformat(), **res}) + "\n")
        return res
    except Exception as e:
        res["error"] = str(e)
        return res

async def run_master_suite():
    print(">>> [INITIATING MASTER VALIDATION SUITE v1.0.0-GAIP]")
    
    # SMOKE TEST
    p1 = {"task_id": "SMOKE_01", "complexity": 0.1, "energy_delta": 0, "dimensions": 35, "validation_status": "PASSED"}
    res1 = await pipeline(p1)
    assert res1["authorized"] is True
    print("SMOKE TEST: PASSED")

    # PEN TEST (Non-numeric Dimensions Fallback)
    p2 = {"task_id": "PEN_01", "complexity": 0.1, "energy_delta": 0, "dimensions": "INVALID_DIMS", "validation_status": "PASSED"}
    res2 = await pipeline(p2)
    assert res2["authorized"] is False and res2["aas"]["report"]["D_INTEGRITY"] is False
    print("PEN TEST: PASSED (Dimensional Safety Fallback)")

    # CHAOS TEST (Invariant Breach)
    p3 = {"task_id": "CHAOS_01", "energy_delta": 999, "dimensions": 35, "validation_status": "PASSED"}
    res3 = await pipeline(p3)
    assert res3["authorized"] is False
    print("CHAOS TEST: PASSED (U34 Invariant Enforcement)")

    print("\n>>> [MASTER VALIDATION COMPLETE: SYSTEM SEALED]")

if __name__ == "__main__":
    asyncio.run(run_master_suite())
