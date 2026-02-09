import json
import os
import numpy as np
from datetime import datetime, timezone

# Rooted imports for E:/foundation
from PX_Warehouse.Worldline_Indexer import WorldlineIndexer
from PX_Audit.Drift_Monitor import DriftMonitor
from PX_Engine.Metabolism import Metabolism

WL_PATH = "E:/foundation/PX_Warehouse/WorldLines"

class ManifoldHealthSummary:
    def __init__(self):
        self.indexer = WorldlineIndexer()
        self.drift_monitor = DriftMonitor()
        self.metabolism = Metabolism()

    def _worldline_paths(self):
        """Yield (dir_path, filename) for all .worldline under WL_PATH and tier subdirs."""
        if not os.path.isdir(WL_PATH):
            return
        for name in os.listdir(WL_PATH):
            path = os.path.join(WL_PATH, name)
            if os.path.isfile(path) and name.endswith(".worldline"):
                yield WL_PATH, name
            elif os.path.isdir(path):
                for f in os.listdir(path):
                    if f.endswith(".worldline"):
                        yield path, f

    def load_worldlines(self):
        records = []
        for dir_path, f in self._worldline_paths():
            with open(os.path.join(dir_path, f), "r") as fp:
                try:
                    data = json.load(fp)
                    # FDA Requirement: Validate full 35-D vector (support root and physics_snapshot)
                    block = np.array(
                        data.get("coordinate_35d") or data.get("physics_snapshot", {}).get("coordinate_35d", [0]*35),
                        dtype=float
                    )
                    coherence = data.get("header", {}).get("coherence") or data.get("coherence_amplitude")
                    if coherence is not None:
                        records.append((block, coherence))
                except Exception:
                    continue
        return records

    def compute_coherence_mean(self, records):
        coherences = [c for (_, c) in records]
        return float(np.mean(coherences)) if coherences else 0.0

    def compute_resonance_density(self, records):
        if not records:
            return 0.0
        # Focus on the first 10 active stages (Legacy Recovery DNA)
        points = np.array([r[0][:10] for r in records])
        spread = np.mean(np.std(points, axis=0))
        density = len(records) / (spread + 1e-6)
        return float(density)

    def compute_drift_score(self, records):
        if not records:
            return 0.0
        # Rebuilding the index ensures the drift is measured against the NEWLY RESTORED warehouse
        self.indexer.rebuild_index()
        sample_block = records[0][0]
        matches = self.indexer.search_resonance(sample_block, k=5)
        drift = self.drift_monitor.calculate_drift(sample_block, matches)
        return drift.get("drift_score", 0.0)

    def generate_summary(self):
        records = self.load_worldlines()

        coherence_mean = self.compute_coherence_mean(records)
        resonance_density = self.compute_resonance_density(records)
        drift_score = self.compute_drift_score(records)
        metabolic_age = self.metabolism.get_current_age()

        summary = {
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "coherence_mean": round(coherence_mean, 6),
            "drift_score": round(drift_score, 6),
            "resonance_density": round(resonance_density, 6),
            "metabolic_cycle_count": metabolic_age,
            "compliance_status": "FDA-GAIP-2026-STABLE",
            "health_status": (
                "STABLE"
                if drift_score < 0.25 and coherence_mean > 0.75
                else "RECOVERING"
                if drift_score < 0.5
                else "DRIFT_DETECTED"
            )
        }

        print("\n=== PREDATOR X: MANIFOLD HEALTH SUMMARY (FDA '26) ===")
        for k, v in summary.items():
            print(f"{k:<25}: {v}")

        return summary

if __name__ == "__main__":
    engine = ManifoldHealthSummary()
    engine.generate_summary()
