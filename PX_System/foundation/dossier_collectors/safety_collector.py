"""
Safety Profile collector -- Section 04.

Pulls toxicity data, safety margins, and organ toxicity predictions
from ADMET engine results. Classifies into constitutional toxicity tiers.

Constitutional: Python stdlib only. No mock data.
Toxicity tiers (Law L11) are immutable hard limits.
"""

from __future__ import annotations

from typing import Any, Dict

from .base_collector import BaseCollector


class SafetyCollector(BaseCollector):
    """Collects safety profile data from ADMET engine results."""

    section_name = "safety_profile"

    def collect(self) -> Dict[str, Any]:
        try:
            admet = (
                self._get_engine_result("admet")
                or self._get_engine_result("ADMET")
                or {}
            )
            tox_raw = admet.get("toxicity", {})
            tox = tox_raw if isinstance(tox_raw, dict) else {}
            tox_index = tox.get("toxicity_index", admet.get("toxicity_index", 0))

            safety_margin_raw = admet.get("safety_margins", {})
            safety_margin = (
                safety_margin_raw if isinstance(safety_margin_raw, dict) else {}
            )

            # Determine tier classification (Law L11 -- immutable thresholds)
            if tox_index < 0.01:
                tier_class = "TOXICITY_DIAMOND"
            elif tox_index < 0.0200:
                tier_class = "TOXICITY_GOLD"
            elif tox_index < 0.0210:
                tier_class = "TOXICITY_SILVER"
            else:
                tier_class = "TOXICITY_FAILURE"

            return {
                "overall_toxicity_index": tox_index,
                "tier_classification": tier_class,
                "safety_margin_fold": safety_margin.get("therapeutic_index", 0),
                "organ_toxicity": [
                    {
                        "organ": "liver",
                        "prediction_score": tox.get("hepatotoxicity", 0),
                        "threshold": 0.3,
                        "passed": tox.get("hepatotoxicity", 0) < 0.3,
                        "confidence": "MEDIUM",
                        "flags": [],
                    },
                    {
                        "organ": "heart",
                        "prediction_score": tox.get("cardiotoxicity", 0),
                        "threshold": 0.3,
                        "passed": tox.get("cardiotoxicity", 0) < 0.3,
                        "confidence": "MEDIUM",
                        "flags": [],
                    },
                    {
                        "organ": "kidney",
                        "prediction_score": tox.get("nephrotoxicity", 0),
                        "threshold": 0.3,
                        "passed": tox.get("nephrotoxicity", 0) < 0.3,
                        "confidence": "MEDIUM",
                        "flags": [],
                    },
                ],
                "cardiac_safety": {
                    "herg_ic50_uM": tox.get("herg_ic50_uM", 10.0),
                    "herg_threshold_uM": 10.0,
                    "herg_passed": tox.get("herg_ic50_uM", 10.0) >= 10.0,
                    "herg_liability_class": (
                        "LOW"
                        if tox.get("herg_ic50_uM", 10.0) >= 10.0
                        else "HIGH"
                    ),
                },
                "cyp_inhibition": [],
                "genotoxicity": [],
                "safety_confidence": "MEDIUM",
                "recommended_safety_studies": [
                    "hERG patch clamp",
                    "Ames test",
                    "28-day repeat dose toxicity",
                ],
                "monitoring_parameters": [
                    "Liver function tests",
                    "ECG monitoring",
                    "Renal function",
                ],
                "incomplete": not bool(admet),
            }
        except Exception as e:
            self._error(f"Failed to collect safety profile: {e}")
            return {"incomplete": True, "error": str(e)}
