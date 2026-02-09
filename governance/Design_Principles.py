"""
╔══════════════════════════════════════════════════════════════════════════════╗
║ Toxicity_Constitution.py                                                     ║
║ PREDATOR X :: GOVERNANCE & DESIGN PRINCIPLES                                 ║
║ STATUS: IMMUTABLE LAW                                                        ║
╚══════════════════════════════════════════════════════════════════════════════╝

DESIGN PRINCIPLE:
    Predator X evaluates toxicity as a CONTEXTUAL CONSTRAINT, not a 
    disqualifying condition.

CORE DOCTRINE:
    1. A high raw toxicity score (e.g., Hepatotoxicity > 0.2) is NOT a failure.
    2. It is merely a variable to be offset by a wider Therapeutic Window.
    3. Clearance decisions must be based on INTEGRATED DOMINANCE.

    We do not hunt for "safe" molecules. We hunt for "dominant" molecules 
    where the Safety Margin (SM) crushes the risk.

    FORMULA:
    Value = (Efficacy * Exposure * Safety_Margin) / (Risk * Dose_Geometry)

    "The poison is in the dose. The cure is in the geometry."
"""

from typing import Dict, Tuple

class ToxicityPolicy:
    """
    The Governing Logic for all Safety Gates in Predator X.
    This class enforces the 'Smart Offset' doctrine over the 'Zero Harm' fallacy.
    """
    
    # The "Old World" limit (Legacy Pharma). We ignore this if SM is high.
    LEGACY_TOX_LIMIT = 0.021 

    @staticmethod
    def evaluate_dominance(candidate_metrics: Dict) -> Tuple[bool, str]:
        """
        Evaluates clearance based on Integrated Dominance.
        
        Args:
            candidate_metrics: Dict containing:
                - toxicity_index (float): Raw physical toxicity (0.0 - 1.0)
                - safety_margin (float): Ratio of Tox Threshold / Cmax
                - dominance_score (float): Calculated commercial viability
        
        Returns:
            (Authorized [bool], Rationale [str])
        """
        tox = candidate_metrics.get("toxicity_index", 1.0)
        sm = candidate_metrics.get("safety_margin", 1.0)
        
        # 1. THE PERFECT MOLECULE (Low Risk, High Safety)
        if tox < ToxicityPolicy.LEGACY_TOX_LIMIT:
            return True, "AUTHORIZED: CLASSIC_GOLD (Low Toxicity Profile)"

        # 2. THE OFFSET DIAMOND (High Risk, Massive Safety Margin)
        # This is where Predator X beats the industry.
        # We accept toxicity IF the therapeutic geometry allows us to dodge it.
        if sm > 50.0:
            return True, f"AUTHORIZED: OFFSET_DIAMOND (High Tox {tox:.3f} offset by Safety Margin {sm:.1f}x)"

        # 3. THE TRASH (High Risk, No Geometric Escape)
        return False, f"REJECTED: UNMANAGED_RISK (Tox {tox:.3f} with insufficient Safety Margin {sm:.1f}x)"

    @staticmethod
    def get_philosophy() -> str:
        return (
            "Clearance decisions are based on integrated dominance across "
            "efficacy, exposure, safety margins, and therapeutic geometry."
        )

# EXPLICIT WARNING FOR FUTURE ENGINEERS
if __name__ == "__main__":
    print("WARNING: DO NOT REVERT TO BINARY TOXICITY FILTERS.")
    print(ToxicityPolicy.get_philosophy())