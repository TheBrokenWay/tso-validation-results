import json
import os
from datetime import datetime, timezone
from typing import Dict, Optional


_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

class AASCircuitBreaker:
    CLOSED = "CLOSED"
    OPEN = "OPEN"
    LOCKED = "LOCKED"

    VALID_STATES = {CLOSED, OPEN, LOCKED}

    def __init__(self, state_path: str | None = None):
        self.state_path = state_path or os.path.join(_REPO_ROOT, "PX_Warehouse", "aas_circuit_breaker.json")
        os.makedirs(os.path.dirname(self.state_path), exist_ok=True)
        self.state = self._load_state()

    def _load_state(self) -> Dict:
        if os.path.exists(self.state_path):
            with open(self.state_path, "r", encoding="utf-8") as f:
                data = json.load(f)
                state = data.get("state", self.CLOSED)
                if state in self.VALID_STATES:
                    return data
        return {
            "state": self.CLOSED,
            "reason": "INITIALIZED",
            "updated_at": datetime.now(timezone.utc).isoformat(),
        }

    def _persist(self, state: str, reason: str, meta: Optional[Dict] = None) -> None:
        if state not in self.VALID_STATES:
            raise ValueError(f"Invalid circuit state: {state}")
        payload = {
            "state": state,
            "reason": reason,
            "updated_at": datetime.now(timezone.utc).isoformat(),
        }
        if meta:
            payload["meta"] = meta
        with open(self.state_path, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
        self.state = payload

    def lock(self, reason: str, meta: Optional[Dict] = None) -> None:
        self._persist(self.LOCKED, reason, meta)

    def open(self, reason: str, meta: Optional[Dict] = None) -> None:
        self._persist(self.OPEN, reason, meta)

    def close(self, reason: str, meta: Optional[Dict] = None) -> None:
        self._persist(self.CLOSED, reason, meta)

    def mark_normalized(self, report: Dict) -> None:
        self._persist(self.CLOSED, "NORMALIZATION_COMPLETE", report)

    def is_locked(self) -> bool:
        return self.state.get("state") == self.LOCKED

    def snapshot(self) -> Dict:
        return dict(self.state)
