"""
Generate PERFORMANCE_BENCHMARK_REPORT.md

Runs performance benchmarks and compares against baselines.
"""
from __future__ import annotations

import hashlib
import statistics
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

BASELINES = {
    "ope_p50_ms": 5.0,
    "ope_p95_ms": 10.0,
    "ope_p99_ms": 15.0,
    "zeus_p50_ms": 2.0,
    "zeus_p95_ms": 5.0,
    "zeus_p99_ms": 10.0,
    "quint_p99_ms": 5.0,
    "ope_qps": 500,
    "quint_ops": 10000,
}


def _bench_ope(n: int = 500) -> list[float]:
    """Benchmark OPE latency."""
    from PX_Engine.operations.OPE import run_ope
    latencies = []
    for _ in range(n):
        t0 = time.perf_counter()
        run_ope("CCO")
        latencies.append((time.perf_counter() - t0) * 1000)
    return sorted(latencies)


def _bench_zeus(n: int = 500) -> list[float]:
    """Benchmark Zeus gate latency."""
    from PX_System.foundation.ZeusLaws import run_zeus_gate
    dossier = {
        "engines": {
            "ope": {"molecular_weight": 46.07},
            "admet": {"toxicity": {"toxicity_index": 0.005, "risk_level": "TOXICITY_DIAMOND"}},
        },
        "harm_energy": 0.005,
        "trial_outcome_summary": None,
    }
    latencies = []
    for _ in range(n):
        t0 = time.perf_counter()
        run_zeus_gate(dossier)
        latencies.append((time.perf_counter() - t0) * 1000)
    return sorted(latencies)


def _bench_quint(n: int = 2000) -> list[float]:
    """Benchmark QUINT round-trip latency."""
    from PX_System.foundation.quint.kernel import QType
    from PX_System.foundation.quint.converter import ingest, emit, reset_stats
    reset_stats()
    data = {"smiles": "CCO", "molecular_weight": 46.07}
    latencies = []
    for _ in range(n):
        t0 = time.perf_counter()
        frame = ingest(data, qtype=QType.QMOLECULE)
        emit(frame)
        latencies.append((time.perf_counter() - t0) * 1000)
    return sorted(latencies)


def _pN(latencies: list[float], pct: float) -> float:
    """Calculate percentile from sorted latencies."""
    idx = int(len(latencies) * pct / 100)
    idx = min(idx, len(latencies) - 1)
    return latencies[idx]


def _status(actual: float, baseline: float, lower_better: bool = True) -> str:
    """Return PASS/FAIL status."""
    if lower_better:
        return "PASS" if actual <= baseline else "FAIL"
    return "PASS" if actual >= baseline else "FAIL"


def generate(output_dir: Path) -> Path:
    """Generate the performance benchmark report."""
    ts = datetime.now(timezone.utc)
    stamp = ts.strftime("%Y%m%d_%H%M%S")

    # Run benchmarks
    ope_lat = _bench_ope(500)
    zeus_lat = _bench_zeus(500)
    quint_lat = _bench_quint(2000)

    # Compute stats
    ope_p50 = _pN(ope_lat, 50)
    ope_p95 = _pN(ope_lat, 95)
    ope_p99 = _pN(ope_lat, 99)
    ope_avg = statistics.mean(ope_lat)

    zeus_p50 = _pN(zeus_lat, 50)
    zeus_p95 = _pN(zeus_lat, 95)
    zeus_p99 = _pN(zeus_lat, 99)

    quint_p50 = _pN(quint_lat, 50)
    quint_p99 = _pN(quint_lat, 99)

    # Throughput
    ope_elapsed = sum(ope_lat) / 1000
    ope_qps = len(ope_lat) / ope_elapsed if ope_elapsed > 0 else 0

    quint_elapsed = sum(quint_lat) / 1000
    quint_ops = len(quint_lat) / quint_elapsed if quint_elapsed > 0 else 0

    all_pass = (
        ope_p99 <= BASELINES["ope_p99_ms"]
        and zeus_p99 <= BASELINES["zeus_p99_ms"]
        and quint_p99 <= BASELINES["quint_p99_ms"]
        and ope_qps >= BASELINES["ope_qps"]
        and quint_ops >= BASELINES["quint_ops"]
    )

    bench_hash = hashlib.sha256(
        f"{ts.isoformat()}{ope_p99:.4f}{zeus_p99:.4f}{quint_p99:.4f}".encode()
    ).hexdigest()[:16]

    content = f"""# PREDATOR X — PERFORMANCE BENCHMARK REPORT

| Field | Value |
|-------|-------|
| **Report ID** | `PX-BEN-{stamp}-{bench_hash}` |
| **Generated** | {ts.isoformat()} |
| **Overall** | **{"ALL BASELINES MET" if all_pass else "REGRESSION DETECTED"}** |
| **Iterations** | OPE: 500, Zeus: 500, QUINT: 2000 |

---

## 1. OPE Engine Latency

| Metric | Actual | Baseline | Status |
|--------|--------|----------|--------|
| P50 | {ope_p50:.3f} ms | {BASELINES["ope_p50_ms"]} ms | {_status(ope_p50, BASELINES["ope_p50_ms"])} |
| P95 | {ope_p95:.3f} ms | {BASELINES["ope_p95_ms"]} ms | {_status(ope_p95, BASELINES["ope_p95_ms"])} |
| P99 | {ope_p99:.3f} ms | {BASELINES["ope_p99_ms"]} ms | {_status(ope_p99, BASELINES["ope_p99_ms"])} |
| Mean | {ope_avg:.3f} ms | — | — |
| Throughput | {ope_qps:.0f} qps | {BASELINES["ope_qps"]} qps | {_status(ope_qps, BASELINES["ope_qps"], lower_better=False)} |

## 2. Zeus Gate Latency

| Metric | Actual | Baseline | Status |
|--------|--------|----------|--------|
| P50 | {zeus_p50:.3f} ms | {BASELINES["zeus_p50_ms"]} ms | {_status(zeus_p50, BASELINES["zeus_p50_ms"])} |
| P95 | {zeus_p95:.3f} ms | {BASELINES["zeus_p95_ms"]} ms | {_status(zeus_p95, BASELINES["zeus_p95_ms"])} |
| P99 | {zeus_p99:.3f} ms | {BASELINES["zeus_p99_ms"]} ms | {_status(zeus_p99, BASELINES["zeus_p99_ms"])} |

## 3. QUINT Round-Trip Latency

| Metric | Actual | Baseline | Status |
|--------|--------|----------|--------|
| P50 | {quint_p50:.3f} ms | — | — |
| P99 | {quint_p99:.3f} ms | {BASELINES["quint_p99_ms"]} ms | {_status(quint_p99, BASELINES["quint_p99_ms"])} |
| Throughput | {quint_ops:.0f} ops/s | {BASELINES["quint_ops"]} ops/s | {_status(quint_ops, BASELINES["quint_ops"], lower_better=False)} |

## 4. System Profile

| Component | Technology |
|-----------|------------|
| Python | {sys.version.split()[0]} |
| Platform | {sys.platform} |
| Physics | Deterministic (no random in evaluation) |
| Governance | Constitutional (fail-closed) |

---
*Generated by Predator X Enterprise Validation Suite v1.0.0*
"""

    out_path = output_dir / f"PERFORMANCE_BENCHMARK_REPORT_{stamp}.md"
    out_path.write_text(content, encoding="utf-8")
    return out_path
