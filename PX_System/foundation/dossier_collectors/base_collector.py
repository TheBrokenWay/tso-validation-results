"""
Base collector — abstract interface for all dossier section collectors.

Each collector gathers data for one section of the DIAMOND dossier package.
Collectors pull from existing engine results, disease constraint files,
and (future) external APIs.

Constitutional: Python stdlib only. Deterministic. No mock data.
"""

from __future__ import annotations

import sys
from typing import Any, Dict, List, Tuple


class BaseCollector:
    """Base class for all dossier section collectors."""

    section_name: str = "base"

    def __init__(self, compound_data: Dict[str, Any], disease_id: str):
        self.compound = compound_data
        self.disease_id = disease_id
        self.errors: List[str] = []
        self.warnings: List[str] = []

    def collect(self) -> Dict[str, Any]:
        """
        Override in subclass. Returns section data dict.

        Must NOT raise — capture errors in self.errors and return
        partial data with an 'incomplete' flag.
        """
        raise NotImplementedError(f"{self.__class__.__name__}.collect() not implemented")

    def validate(self, data: Dict[str, Any]) -> Tuple[bool, List[str]]:
        """
        Validate collected data meets section requirements.

        Returns:
            (valid, errors) — True if all required fields present and within bounds.
        """
        errors = []
        if not data:
            errors.append(f"{self.section_name}: empty data")
        return len(errors) == 0, errors

    def _get_engine_result(self, engine_name: str) -> Dict[str, Any]:
        """Extract a specific engine result from compound data."""
        stages = self.compound.get("stages", {})
        if engine_name in stages:
            return stages[engine_name]
        # Fallback: engine result at top level
        return self.compound.get(engine_name, {})

    def _get_constraint(self) -> Dict[str, Any]:
        """Load disease constraint data."""
        try:
            from PX_Domain.PRV_Diseases.disease_registry import get_disease_constraint
            frame = get_disease_constraint(self.disease_id)
            if frame is not None:
                return frame.payload
        except (ImportError, Exception) as e:
            self.warnings.append(f"Could not load constraint for {self.disease_id}: {e}")
        return {}

    def _warn(self, msg: str) -> None:
        """Record a warning."""
        self.warnings.append(f"{self.section_name}: {msg}")
        print(f"    WARN [{self.section_name}]: {msg}", file=sys.stderr)

    def _error(self, msg: str) -> None:
        """Record an error."""
        self.errors.append(f"{self.section_name}: {msg}")
        print(f"    ERROR [{self.section_name}]: {msg}", file=sys.stderr)
