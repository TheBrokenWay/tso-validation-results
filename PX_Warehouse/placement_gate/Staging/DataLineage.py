import json
import os
from datetime import datetime

class DataLineage:
    """
    Automated Data Lineage Graph for asset lifecycle tracking.
    Tracks: Raw → Stratified → Deduplicated → Monetized
    """
    def __init__(self, db_path=r"E:\foundation\PX_Warehouse\Operations\lineage_graph.json"):
        self.db_path = db_path
        self.graph = self._load_graph()

    def record_transition(self, asset_id, source_stage, target_stage, metadata=None):
        """
        Records a transition for an asset.
        """
        if asset_id not in self.graph:
            self.graph[asset_id] = []

        entry = {
            "timestamp": datetime.now().isoformat(),
            "from": source_stage,
            "to": target_stage,
            "metadata": metadata or {}
        }
        self.graph[asset_id].append(entry)
        self._save_graph()

    def get_lineage(self, asset_id):
        return self.graph.get(asset_id, [])

    def _load_graph(self):
        if os.path.exists(self.db_path):
            with open(self.db_path, 'r') as f:
                return json.load(f)
        return {}

    def _save_graph(self):
        os.makedirs(os.path.dirname(self.db_path), exist_ok=True)
        with open(self.db_path, 'w') as f:
            json.dump(self.graph, f, indent=4)

if __name__ == "__main__":
    lineage = DataLineage()
    # Example: Record a discovery to stratification transition
    lineage.record_transition("CHEMBL286791", "RAW_SIMULATION", "GOLD_CORE", {"dominance_score": 0.85})
    print(f"Lineage recorded for CHEMBL286791")
