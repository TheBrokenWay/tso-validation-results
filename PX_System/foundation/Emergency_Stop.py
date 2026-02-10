"""
Emergency Stop Controller
Fail-closed safety switch used across OLYMPUS subsystems.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
from datetime import datetime, timezone
from typing import Any, Dict


try:
    from PX_System.foundation.Sovereign_Log_Chain import append as log_append
    _HAS_LOG_CHAIN = True
except Exception:
    _HAS_LOG_CHAIN = False


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


class EmergencyStop:
    """
    Persistent emergency stop with a local state file.
    """
    STATE_RUNNING = "RUNNING"
    STATE_STOPPED = "STOPPED"

    def __init__(self) -> None:
        data_root = Path(
            os.environ.get("OLYMPUS_CACHE_DIR")
            or os.environ.get("OLYMPUS_DATA_DIR")
            or (Path(__file__).resolve().parents[1] / "temp_cache")
        )
        data_root.mkdir(parents=True, exist_ok=True)
        self._state_path = data_root / "emergency_stop.json"
        self._status = self._load_state()
        self.state = self._status.get("state", self.STATE_RUNNING)

    def _load_state(self) -> Dict[str, Any]:
        if not self._state_path.exists():
            return {
                "state": self.STATE_RUNNING,
                "trigger_reason": None,
                "trigger_severity": None,
                "trigger_time": None,
                "authorized_by": None,
                "source": None,
                "last_checked": _now_iso(),
            }
        try:
            with open(self._state_path, "r", encoding="utf-8") as f:
                data = json.load(f)
            if not isinstance(data, dict):
                raise ValueError("Invalid emergency stop state shape")
            data["last_checked"] = _now_iso()
            return data
        except Exception as e:
            # Fail-closed: if state can't be read, assume STOPPED.
            return {
                "state": self.STATE_STOPPED,
                "trigger_reason": "STATE_READ_ERROR",
                "trigger_severity": "CRITICAL",
                "trigger_time": _now_iso(),
                "authorized_by": "SYSTEM",
                "source": "EMERGENCY_STOP",
                "error": str(e),
                "last_checked": _now_iso(),
            }

    def _save_state(self, data: Dict[str, Any]) -> None:
        tmp_path = self._state_path.with_suffix(".tmp")
        with open(tmp_path, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2)
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp_path, self._state_path)

    def refresh(self) -> Dict[str, Any]:
        self._status = self._load_state()
        self.state = self._status.get("state", self.STATE_RUNNING)
        return self._status

    def is_stopped(self) -> bool:
        status = self.refresh()
        return status.get("state") == self.STATE_STOPPED

    def get_status(self) -> Dict[str, Any]:
        return self.refresh()

    def manual_trigger(self, reason: str, authorized_by: str, severity: str = "CRITICAL", source: str = "MANUAL") -> None:
        status = {
            "state": self.STATE_STOPPED,
            "trigger_reason": str(reason)[:500],
            "trigger_severity": str(severity)[:50],
            "trigger_time": _now_iso(),
            "authorized_by": str(authorized_by)[:100],
            "source": str(source)[:100],
            "last_checked": _now_iso(),
        }
        self._save_state(status)
        self._status = status
        self.state = self.STATE_STOPPED

        if _HAS_LOG_CHAIN:
            try:
                log_append(
                    "EMERGENCY_STOP_TRIGGERED",
                    {
                        "reason": status["trigger_reason"],
                        "authorized_by": status["authorized_by"],
                        "severity": status["trigger_severity"],
                    },
                    {"source": "Emergency_Stop", "function": "manual_trigger"}
                )
            except Exception:
                # Never let audit logging failure re-break emergency stop.
                pass

    def clear(self, authorized_by: str, reason: str = "RESET") -> None:
        status = {
            "state": self.STATE_RUNNING,
            "trigger_reason": str(reason)[:500],
            "trigger_severity": "INFO",
            "trigger_time": _now_iso(),
            "authorized_by": str(authorized_by)[:100],
            "source": "MANUAL",
            "last_checked": _now_iso(),
        }
        self._save_state(status)
        self._status = status
        self.state = self.STATE_RUNNING
