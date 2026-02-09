import asyncio
import time
import importlib.util
import os
import json
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[1]

def load_mod(n, p):
    path = _REPO_ROOT / p.replace("/", os.sep)
    if not path.exists():
        raise FileNotFoundError(f"Benchmark organ not found: {path}")
    s = importlib.util.spec_from_file_location(n, str(path))
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m

# Load the Sealed Organs (01_Executive = CSA, 02_Audit = AAS)
CSA = load_mod("CSA", "01_Executive/CSA_Pentarchy.py")
AAS = load_mod("AAS", "02_Audit/AAS_Verification.py")

async def pipeline(prop):
    # Direct departmental handshake
    c = await CSA.CSAPentarchy().evaluate_action(prop)
    a = AAS.AASVerification().verify_invariants(prop)
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
