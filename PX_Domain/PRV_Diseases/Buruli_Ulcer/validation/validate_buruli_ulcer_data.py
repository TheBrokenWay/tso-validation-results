"""
Buruli ulcer data validation.

TODO: Implement disease-specific data validation rules.
Follows the pattern established by Nipah validator.
"""

import json
import os


class BuruliUlcerValidationError(Exception):
    pass


def validate_data(data: dict) -> bool:
    """Validate data against Buruli ulcer constraints.

    TODO: Implement disease-specific validation rules.
    """
    if not data:
        raise BuruliUlcerValidationError("Empty data")
    return True


def load_constraints() -> dict:
    """Load disease constraint file."""
    config_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "config")
    constraint_path = os.path.join(config_dir, "buruli_ulcer_constraints.json")
    with open(constraint_path, "r", encoding="utf-8") as f:
        return json.load(f)
