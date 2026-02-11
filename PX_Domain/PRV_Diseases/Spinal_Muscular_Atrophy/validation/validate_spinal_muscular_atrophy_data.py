"""
Spinal muscular atrophy data validation.

TODO: Implement disease-specific data validation rules.
Follows the pattern established by Nipah validator.
"""

import json
import os


class SpinalMuscularAtrophyValidationError(Exception):
    pass


def validate_data(data: dict) -> bool:
    """Validate data against Spinal muscular atrophy constraints.

    TODO: Implement disease-specific validation rules.
    """
    if not data:
        raise SpinalMuscularAtrophyValidationError("Empty data")
    return True


def load_constraints() -> dict:
    """Load disease constraint file."""
    config_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "config")
    constraint_path = os.path.join(config_dir, "spinal_muscular_atrophy_constraints.json")
    with open(constraint_path, "r", encoding="utf-8") as f:
        return json.load(f)
