import json
import os
import numpy as np
from PX_Audit.Mural_Network import update_node

_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DEFAULT_DB_PATH = os.path.join(_REPO_ROOT, "PX_Warehouse", "WorldLines")


class ManifoldHealth:
    def __init__(self, db_path=None):
        self.db_path = db_path or DEFAULT_DB_PATH

    def calculate_vitals(self):
        # Default buckets to prevent KeyErrors
        vitals = {"BLUE": [], "RED": [], "TERMINATED": [], "LEGACY": []}
        
        if not os.path.exists(self.db_path):
            return {"status": "EMPTY"}

        for file in os.listdir(self.db_path):
            if file.endswith(".worldline"):
                try:
                    with open(os.path.join(self.db_path, file), "r") as f:
                        data = json.load(f)
                        h = data.get("header", {})
                        # If route is missing or unknown, it goes to LEGACY
                        route = h.get("route", "LEGACY")
                        if route not in vitals:
                            route = "LEGACY"
                            
                        vitals[route].append(h.get("coherence", 0.0))
                except Exception:
                    continue

        all_points = vitals["BLUE"] + vitals["LEGACY"]
        mean_coh = np.mean(all_points) if all_points else 0.0
        
        report = {
            "new_research_count": len(vitals["BLUE"]),
            "legacy_resource_count": len(vitals["LEGACY"]),
            "coherence_mean": round(float(mean_coh), 4),
            "status": "CLEANING_LEGACY" if len(vitals["LEGACY"]) > 0 else "STABLE"
        }
        
        update_node("ManifoldHealth", "1.2.0", f"COH-{report['coherence_mean']}", report["status"])
        return report
