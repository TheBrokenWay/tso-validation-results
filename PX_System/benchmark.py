import asyncio
import time
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from PX_Executive.CSA_Pentarchy import CSAPentarchy
from PX_Audit.AAS_Verification import AASVerification

async def pipeline(prop):
    # Direct departmental handshake
    c = await CSAPentarchy().evaluate_action(prop)
    a = AASVerification().verify_invariants(prop)
    authorized = (c["status"] == "COHERENT" and a["status"] == "SUCCESS")
    return {"authorized": authorized, "csa": c, "aas": a}

async def run_bench():
    N = 100
    correct = 0
    errors = 0
    lat = []

    # Deterministic Baseline Proposal (Valid Assembly)
    p = {"complexity": 0.1, "energy_delta": 0, "dimensions": 35, "validation_status": "PASSED"}

    print(f">>> [TELEMETRY] Internal Pipeline Benchmark: {N} Cycles...")

    for _ in range(N):
        start = time.perf_counter()
        try:
            r = await pipeline(p)
            if r["authorized"] is True:
                correct += 1
            else:
                errors += 1
        except Exception:
            errors += 1
        lat.append(time.perf_counter() - start)

    avg_lat = sum(lat)/N
    print("\n=== OLYMPUS PIPELINE BENCHMARK (TRUE STATS) ===")
    print(f"Total Cycles:     {N}")
    print(f"Authorized:       {correct}")
    print(f"Errors/Rejects:   {errors}")
    print(f"Accuracy:          {correct/N:.3f}")
    print(f"Avg Latency:       {avg_lat:.6f} sec")
    print(f"Throughput:        {1/avg_lat:.2f} decisions/sec")
    print("================================================")

if __name__ == "__main__":
    asyncio.run(run_bench())
