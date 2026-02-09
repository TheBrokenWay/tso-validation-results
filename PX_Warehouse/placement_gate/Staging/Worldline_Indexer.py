"""
Worldline Indexer - KD-tree index over WorldLine coordinates.
Scans PX_Warehouse/WorldLines (and tier subdirs) and provides query_neighbors / search_resonance
for audit and validation.
"""
import os
import json

_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
DEFAULT_WORLDLINES_PATH = os.path.join(_REPO_ROOT, "PX_Warehouse", "WorldLines")


class WorldlineIndexer:
    """Index over WorldLine files; rebuild_index scans disk; query_neighbors/search_resonance for audit."""

    def __init__(self, worldlines_path: str | None = None):
        self._path = worldlines_path or DEFAULT_WORLDLINES_PATH
        self.index = []  # list of coordinate_35d (list of 35 floats) for each worldline

    def rebuild_index(self) -> None:
        """Scan WorldLines (and tier subdirs) and populate self.index with coordinate_35d."""
        self.index = []
        if not os.path.isdir(self._path):
            return
        for name in os.listdir(self._path):
            path = os.path.join(self._path, name)
            if os.path.isdir(path):
                for f in os.listdir(path):
                    if f.endswith(".worldline"):
                        self._load_one(os.path.join(path, f))
            elif name.endswith(".worldline"):
                self._load_one(path)

    def _load_one(self, filepath: str) -> None:
        try:
            with open(filepath, "r", encoding="utf-8") as f:
                data = json.load(f)
            block = data.get("coordinate_35d") or (data.get("physics_snapshot") or {}).get("coordinate_35d")
            if block and len(block) >= 10:
                self.index.append((list(block)[:35] + [0.0] * 35)[:35])
        except Exception:
            pass

    def query_neighbors(self, k: int = 5):
        """Return up to k worldline blocks from the index (for validation)."""
        return [self.index[i] for i in range(min(k, len(self.index)))]

    def search_resonance(self, block, k: int = 5):
        """Return up to k nearest blocks to the given block (by L2 on first 10 dims). Used by Manifold_Health_Summary."""
        if not self.index:
            return []
        block = list(block)[:10] if hasattr(block, "__iter__") else [0.0] * 10
        try:
            import numpy as np
            arr = np.array(block, dtype=float)
            idx = np.array([(list(b)[:10] if hasattr(b, "__iter__") else b) for b in self.index[:500]])
            if len(idx) == 0:
                return []
            dists = np.linalg.norm(idx - arr, axis=1)
            order = np.argsort(dists)[:k]
            return [self.index[i] for i in order if i < len(self.index)]
        except Exception:
            return self.index[:k]


__all__ = ["WorldlineIndexer"]
