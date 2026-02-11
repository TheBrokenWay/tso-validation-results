"""
QUINT Registry — Schema registry for QTypes.

Maintains the canonical field definitions for each QType so that
any component can query what fields a QMOLECULE or QDOSSIER should contain.

The registry is read-only after initialization. Engines register their
output schemas at import time; the registry is frozen before first use.

Constitutional: Python stdlib only.
"""

from __future__ import annotations

from typing import Any, Dict, FrozenSet, List, Optional, Set

from PX_System.foundation.quint.kernel import QType


# ── Schema definition ────────────────────────────────────────────────────────

class FieldDef:
    """Definition of a single field in a QType schema."""
    __slots__ = ("name", "required", "field_type", "description")

    def __init__(
        self,
        name: str,
        required: bool = False,
        field_type: str = "any",
        description: str = "",
    ):
        self.name = name
        self.required = required
        self.field_type = field_type
        self.description = description

    def __repr__(self) -> str:
        req = "required" if self.required else "optional"
        return f"FieldDef({self.name!r}, {req}, {self.field_type})"


class Schema:
    """Schema for a QType — collection of FieldDefs."""
    __slots__ = ("qtype", "fields", "description")

    def __init__(self, qtype: QType, fields: List[FieldDef], description: str = ""):
        self.qtype = qtype
        self.fields = {f.name: f for f in fields}
        self.description = description

    @property
    def required_fields(self) -> FrozenSet[str]:
        return frozenset(n for n, f in self.fields.items() if f.required)

    @property
    def all_fields(self) -> FrozenSet[str]:
        return frozenset(self.fields.keys())

    def check(self, payload: Dict[str, Any]) -> List[str]:
        """Return list of missing required fields."""
        return [n for n in self.required_fields if n not in payload]


# ── Global registry ──────────────────────────────────────────────────────────

_REGISTRY: Dict[QType, Schema] = {}
_FROZEN = False


def register(schema: Schema) -> None:
    """Register a schema. Raises if registry is frozen."""
    if _FROZEN:
        raise RuntimeError("QUINT registry is frozen — cannot register new schemas")
    _REGISTRY[schema.qtype] = schema


def freeze() -> None:
    """Freeze the registry. No more registrations allowed."""
    global _FROZEN
    _FROZEN = True


def is_frozen() -> bool:
    return _FROZEN


def get_schema(qtype: QType) -> Optional[Schema]:
    """Look up the schema for a QType."""
    return _REGISTRY.get(qtype)


def registered_types() -> Set[QType]:
    """Return all registered QTypes."""
    return set(_REGISTRY.keys())


def check_payload(qtype: QType, payload: Dict[str, Any]) -> List[str]:
    """Check a payload against its QType schema. Returns missing required fields."""
    schema = _REGISTRY.get(qtype)
    if schema is None:
        return []  # no schema registered = no requirements
    return schema.check(payload)


# ── Built-in schemas ─────────────────────────────────────────────────────────

def _F(name: str, required: bool = False, ft: str = "any", desc: str = "") -> FieldDef:
    return FieldDef(name, required, ft, desc)


register(Schema(QType.QMOLECULE, [
    _F("smiles", True, "str", "SMILES string"),
    _F("mw", False, "float", "Molecular weight"),
    _F("logp", False, "float", "LogP"),
    _F("hbd", False, "int", "H-bond donors"),
    _F("hba", False, "int", "H-bond acceptors"),
    _F("tpsa", False, "float", "Topological polar surface area"),
    _F("toxicity", False, "float", "Toxicity index"),
    _F("compound_id", False, "str", "Compound identifier"),
], description="Molecular representation"))

register(Schema(QType.QRESULT, [
    _F("engine_id", True, "str", "Engine identifier"),
    _F("status", True, "str", "PASSED/FAILED/REVIEW_REQUIRED"),
    _F("score", False, "float", "Engine score"),
    _F("authorized", False, "bool", "Authorization status"),
    _F("sign_off", False, "dict", "Engine sign-off block"),
], description="Engine computation result"))

register(Schema(QType.QDOSSIER, [
    _F("compound_id", True, "str", "Compound identifier"),
    _F("smiles", False, "str", "SMILES string"),
    _F("tier", False, "str", "DIAMOND/GOLD/SILVER/BRONZE"),
    _F("grade", False, "str", "Grade assignment"),
    _F("zeus_verdict", False, "dict", "Zeus gate verdict"),
    _F("authorization_chain", False, "dict", "12-engine sign-off chain"),
], description="Complete dossier"))

register(Schema(QType.QPHYSICS, [
    _F("dimensions", True, "int", "Physics vector dimensionality"),
    _F("physics_vector", False, "list", "35D state vector"),
    _F("energy", False, "float", "Total energy"),
    _F("global_sum", False, "float", "U34 global sum"),
], description="Physics state vector / worldline snapshot"))

register(Schema(QType.QSIGNAL, [
    _F("signal_type", True, "str", "Signal type identifier"),
    _F("source", True, "str", "Source engine/module"),
    _F("payload", False, "dict", "Signal payload"),
    _F("priority", False, "int", "Signal priority"),
], description="Inter-engine communication"))

register(Schema(QType.QTRIAL, [
    _F("protocol", True, "dict", "Trial protocol"),
    _F("population", False, "dict", "Virtual population"),
    _F("outcome", False, "dict", "Trial outcome summary"),
    _F("arm_count", False, "int", "Number of trial arms"),
], description="Trial simulation data"))

register(Schema(QType.QCONSTRAINT, [
    _F("disease_name", True, "str", "Disease name"),
    _F("pathogen_type", False, "str", "Pathogen type"),
    _F("target_classes", False, "list", "Target classes"),
    _F("ic50_max_um", False, "float", "IC50 threshold"),
    _F("toxicity_threshold", False, "float", "Toxicity threshold"),
], description="Disease constraint / governance rule"))

register(Schema(QType.QRAW, [], description="Untyped passthrough"))


__all__ = [
    "FieldDef",
    "Schema",
    "check_payload",
    "freeze",
    "get_schema",
    "is_frozen",
    "register",
    "registered_types",
]
