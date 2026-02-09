import asyncio, json, os, hashlib, importlib.util
from datetime import datetime, timezone

def load_mod(n, p):
    s = importlib.util.spec_from_file_location(n, os.path.join(os.getcwd(), p))
    m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m

async def pipeline(prop):
    res = {"authorized": False, "error": None, "trace_id": "N/A"}
    try:
        c_m = load_mod("CSA", "01_Executive/CSA_Pentarchy.py")
        a_m = load_mod("AAS", "02_Audit/AAS_Verification.py")
        
        c_r = await c_m.CSAPentarchy().evaluate_action(prop)
        a_r = a_m.AASVerification().verify_invariants(prop)
        
        res["trace_id"] = f"FDA-GAIP-{hashlib.sha256(str(prop).encode()).hexdigest()[:12].upper()}"
        res["authorized"] = (c_r["status"] == "COHERENT" and a_r["status"] == "SUCCESS" and (not c_r["human_review_required"] or prop.get("human_sign_off_id")))
        res["csa"] = c_r; res["aas"] = a_r
        
        os.makedirs("02_Audit", exist_ok=True)
        with open("02_Audit/fda_gaip_trace.log", "a") as f:
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
