"""
Novelty Budget Engine - Deterministic Implementation
Manages computational budget allocation for discovery campaigns
"""
from typing import Dict, Any, Optional


class NoveltyBudgetEngine:
    """
    Novelty budget allocation engine
    
    Provides:
    - Computational budget tracking
    - Novelty scoring for candidate compounds
    - Resource allocation for discovery campaigns
    
    Constitutional Compliance:
        - L51 (Zero Gaps): Fail-closed when budget tracking needed
        - L9 (Axiomatic Grounding): No arbitrary resource allocation
    """
    
    def __init__(self, max_evaluations: int = 10000):
        """
        Initialize novelty budget engine
        
        Args:
            max_evaluations: Maximum compound evaluations per campaign
        """
        self.max_evaluations = max_evaluations
        self.evaluations_used = 0
        self.budget_exhausted = False
    
    def check_budget(self) -> bool:
        """
        Check if computational budget remains
        
        Returns:
            bool: True if budget available, False if exhausted
        """
        return not self.budget_exhausted and self.evaluations_used < self.max_evaluations
    
    def consume_budget(self, amount: int = 1) -> bool:
        """
        Consume computational budget
        
        Args:
            amount: Number of evaluations to consume
        
        Returns:
            bool: True if budget was available and consumed, False otherwise
        """
        if self.evaluations_used + amount > self.max_evaluations:
            self.budget_exhausted = True
            return False
        
        self.evaluations_used += amount
        return True
    
    def calculate_novelty_score(self, compound: Dict[str, Any], reference_set: Optional[list] = None) -> float:
        """
        Calculate novelty score for a compound
        
        Args:
            compound: Compound data dict
            reference_set: Optional reference compounds for novelty comparison
        
        Returns:
            float: Novelty score (0.0 = not novel, 1.0 = highly novel)
        """
        smiles = (compound.get("smiles") or "").strip()
        name = (compound.get("name") or "").strip()

        if reference_set:
            for ref in reference_set:
                if not isinstance(ref, dict):
                    continue
                if smiles and smiles == (ref.get("smiles") or "").strip():
                    return 0.0
                if name and name == (ref.get("name") or "").strip():
                    return 0.0

        length_score = min(1.0, len(smiles) / 200.0) if smiles else 0.1
        alpha_score = min(1.0, sum(1 for c in smiles if c.isalpha()) / 120.0) if smiles else 0.1
        novelty = (0.6 * length_score) + (0.4 * alpha_score)
        return float(f"{max(0.0, min(1.0, novelty)):.8f}")
    
    def get_budget_status(self) -> Dict[str, Any]:
        """
        Get current budget status
        
        Returns:
            dict: Budget allocation status
        """
        return {
            "max_evaluations": self.max_evaluations,
            "evaluations_used": self.evaluations_used,
            "evaluations_remaining": self.max_evaluations - self.evaluations_used,
            "budget_exhausted": self.budget_exhausted,
            "utilization_pct": round((self.evaluations_used / self.max_evaluations) * 100, 1)
        }


__all__ = ["NoveltyBudgetEngine"]
