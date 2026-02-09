"""
Intake policy configuration for PRV seed libraries.
Loads PX_Warehouse/Operations/Inputs/intake_policy.json when present.
Used by queue builders and PRV_24H_Orchestrator for queue_file, sources, filters, governance.
"""
from pathlib import Path
import json
from typing import Any, Optional

_DEFAULT_REPO_ROOT = Path(__file__).resolve().parents[2]  # foundation


def get_intake_policy(repo_root: Optional[Path] = None) -> dict[str, Any]:
    """
    Load intake_policy.json from PX_Warehouse/Operations/Inputs.
    Returns the full config dict, or {} if file missing/invalid.
    """
    root = Path(repo_root) if repo_root is not None else _DEFAULT_REPO_ROOT
    path = root / "PX_Warehouse" / "Operations" / "Inputs" / "intake_policy.json"
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return {}


def get_queue_filename(repo_root: Optional[Path] = None) -> str:
    """Return output.queue_file from intake policy, or default 'prv_24h_queue.json'."""
    policy = get_intake_policy(repo_root)
    out = policy.get("output") or {}
    return out.get("queue_file", "prv_24h_queue.json")


__all__ = ["get_intake_policy", "get_queue_filename"]
