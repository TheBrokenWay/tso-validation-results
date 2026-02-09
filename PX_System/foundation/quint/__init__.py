"""
foundation.quint
Architectural bridge for QUINT (Constitutional Core).

Why this exists:
- Canonical QUINT sources live in `01_Executive/quint/` (departmental, governed).
- Python cannot import `01_Executive` as a normal package name.
- Runtime code expects `foundation.quint.*` imports.

Compliance:
- No blanket sys.path injection.
- Explicit, path-based module loading only (importlib).
"""

from __future__ import annotations

import importlib
import importlib.util
import sys
from pathlib import Path
from types import ModuleType


_REPO_ROOT = Path(__file__).resolve().parents[2]
_QUINT_SRC = _REPO_ROOT / "01_Executive" / "quint"


def _load(submodule: str, src_file: Path) -> ModuleType:
    full_name = f"{__name__}.{submodule}"
    if full_name in sys.modules:
        return sys.modules[full_name]
    if not src_file.exists():
        raise ImportError(f"QUINT source missing: {src_file}")
    spec = importlib.util.spec_from_file_location(full_name, src_file)
    if spec is None or spec.loader is None:
        raise ImportError(f"Failed to create spec for {full_name} from {src_file}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[full_name] = mod
    spec.loader.exec_module(mod)
    return mod


# Ensure `foundation.quint.brains` package exists before loading the bridged submodule.
importlib.import_module(f"{__name__}.brains")

# Low-level runtime/compiler must load before brains (brains imports them).
quint_compiler = _load("quint_compiler", _QUINT_SRC / "quint_compiler.py")
quint_runtime = _load("quint_runtime", _QUINT_SRC / "quint_runtime.py")

# Brains subtree must load before consensus (consensus imports it)
_load("brains.domain_brain", _QUINT_SRC / "brains" / "domain_brain.py")

# Preload common QUINT modules so imports like `foundation.quint.converter` resolve.
converter = _load("converter", _QUINT_SRC / "converter.py")
consensus = _load("consensus", _QUINT_SRC / "consensus.py")
executive = _load("executive", _QUINT_SRC / "executive.py")
router = _load("router", _QUINT_SRC / "router.py")
sentinels = _load("sentinels", _QUINT_SRC / "sentinels.py")
quint_registry = _load("quint_registry", _QUINT_SRC / "quint_registry.py")
logging = _load("logging", _QUINT_SRC / "logging.py")

# Orchestrator last (depends on the above in most layouts)
orchestrator = _load("orchestrator", _QUINT_SRC / "orchestrator.py")

__all__ = [
    "converter",
    "consensus",
    "executive",
    "router",
    "sentinels",
    "quint_compiler",
    "quint_runtime",
    "quint_registry",
    "logging",
    "orchestrator",
]

