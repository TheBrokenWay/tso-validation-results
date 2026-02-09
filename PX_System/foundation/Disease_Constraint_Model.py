"""
Disease Constraint Model - Deterministic Implementation
"""
from typing import Dict, Any, List


class DiseaseConstraintModel:
    """
    Disease-specific constraint modeling

    Provides:
    - Target protein validation
    - Disease-specific ADMET thresholds
    - Mechanism of action constraints
    """
    
    def __init__(self, disease_name: str):
        self.disease_name = disease_name
        self.constraints = get_disease_constraints(disease_name)


def create_chagas_dcm() -> Dict[str, Any]:
    """
    Factory method: Create Disease Constraint Model for Chagas Disease (Trypanosoma cruzi)
    
    This model defines the constitutional constraints for 12-hour discovery campaigns
    targeting Chagas Disease, a neglected tropical disease caused by T. cruzi.
    
    Returns:
        dict: Disease-specific constraint model with Lipinski and safety parameters
    
    Constitutional Basis:
        - L51 (Zero Gaps): All thresholds from WHO/DNDi guidelines
        - L10 (Harm is Absolute): Safety margin enforced via selectivity index
        
    References:
        - WHO: Chagas disease (American trypanosomiasis) - Fact Sheet
        - DNDi: Target Product Profile for Chagas Disease
        - Buckner FS et al. (2012) Exp Opin Drug Discov. 7(10): 923-941
    """
    return {
        # Disease Identification
        "disease_id": "TD01",  # Tropical Disease 01
        "disease_name": "Chagas Disease",
        "pathogen": "Trypanosoma cruzi",
        "target_protein": "CYP51",  # Primary drug target (sterol 14α-demethylase)
        
        # Efficacy Constraints (Anti-T. cruzi Activity)
        "ic50_max_um": 10.0,  # Maximum IC50 for T. cruzi inhibition (DNDi TPP)
        "ic50_optimal_um": 1.0,  # Optimal potency target
        
        # Safety Constraints (L10: Harm is Absolute)
        "selectivity_index_min": 10.0,  # SI = CC50/IC50 >= 10 (OSE requirement)
        "cc50_min_um": 50.0,  # Minimum cytotoxicity threshold (mammalian cells)
        
        # Lipinski Constraints (Drug-likeness)
        "molecular_weight_max": 500.0,  # Lipinski Rule of Five
        "logp_min": 0.0,
        "logp_max": 5.0,
        "hbd_max": 5,  # Hydrogen bond donors
        "hba_max": 10,  # Hydrogen bond acceptors
        "lipinski_violations_max": 1,  # Allow 1 violation for tropical disease drugs
        
        # ADMET Constraints
        "oral_bioavailability_min": 30.0,  # F% >= 30% (OPE requirement)
        "half_life_min_hr": 4.0,  # t1/2 >= 4hr for BID dosing (OPE requirement)
        "permeability_flag": True,  # Caco-2 > 10^-6 cm/s preferred
        "metabolic_stability_flag": True,  # Microsomal stability required
        
        # Docking Constraints
        "docking_score_max_kcal": -6.0,  # QVina2 score <= -6.0 kcal/mol (OBE/ODC)
        
        # Legal/Regulatory Constraints
        "prv_eligible": True,  # Priority Review Voucher eligibility (OLE)
        "indication": "TD01",  # Tropical Disease designation
        
        # Discovery Campaign Parameters
        "discovery_mode": "screening",  # screening | optimization | hit-to-lead
        "max_candidates": 10000,  # Maximum compounds to evaluate per campaign
        "level_1_criteria": "PASS",  # Requires all constitutional constraints satisfied
        "level_2_criteria": "CONDITIONAL_PASS",  # Allow partial passes for optimization
        
        # Assumptions & Citations
        "assumptions": [
            "A_DCM_01: CYP51 is primary drug target for T. cruzi (Buckner 2012)",
            "A_DCM_02: In vitro IC50 correlates with in vivo efficacy (WHO guidelines)",
            "A_DCM_03: Selectivity index >= 10 indicates acceptable safety margin",
            "A_DCM_04: Oral bioavailability essential for tropical disease treatment compliance"
        ],
        "citations": [
            "WHO: Chagas disease fact sheet (2023)",
            "DNDi: Target Product Profile for Chagas Disease",
            "Buckner FS et al. (2012) Exp Opin Drug Discov. 7(10): 923-941"
        ],
        
        # Constitutional Seal
        "constitutional_compliance": {
            "L51_zero_gaps": True,
            "L10_harm_absolute": True,
            "L9_axiomatic_grounding": True,
            "L5_transparency": True
        }
    }


def create_nipah_malaysia_dcm() -> Dict[str, Any]:
    """
    Disease Constraint Model for Nipah virus — Malaysia strain (NiV-M).

    CFR ~35–45%; pig-intermediate; no sustained human-to-human.
    Strain-aware: use for intervention modeling, vaccine/therapeutic response for NiV-M.
    """
    return {
        "disease_id": "NIPAH_MALAYSIA",
        "disease_name": "Nipah virus (Malaysia strain)",
        "pathogen": "Nipah virus",
        "strain": "NiV_Malaysia",
        "target_proteins": ["NiV_F", "NiV_G", "NiV_N", "NiV_P", "NiV_L"],
        "host_factors": ["EPHRIN_B2", "EPHRIN_B3"],
        "cfr_range": [0.35, 0.45],
        "primary_vectors": ["pig"],
        "human_to_human": False,
        "ic50_max_um": 10.0,
        "selectivity_index_min": 10.0,
        "cc50_min_um": 50.0,
        "molecular_weight_max": 500.0,
        "logp_min": 0.0,
        "logp_max": 5.0,
        "lipinski_violations_max": 1,
        "admet_score_min": 0.6,
        "constitutional_compliance": {"L51_zero_gaps": True, "L10_harm_absolute": True},
    }


def create_nipah_bangladesh_dcm() -> Dict[str, Any]:
    """
    Disease Constraint Model for Nipah virus — Bangladesh strain (NiV-B).

    CFR ~70–100%; bat and human; human-to-human transmission.
    Strain-aware: use for intervention modeling, vaccine/therapeutic response for NiV-B.
    """
    return {
        "disease_id": "NIPAH_BANGLADESH",
        "disease_name": "Nipah virus (Bangladesh strain)",
        "pathogen": "Nipah virus",
        "strain": "NiV_Bangladesh",
        "target_proteins": ["NiV_F", "NiV_G", "NiV_N", "NiV_P", "NiV_L"],
        "host_factors": ["EPHRIN_B2", "EPHRIN_B3"],
        "cfr_range": [0.70, 1.00],
        "primary_vectors": ["bat", "human"],
        "human_to_human": True,
        "ic50_max_um": 10.0,
        "selectivity_index_min": 10.0,
        "cc50_min_um": 50.0,
        "molecular_weight_max": 500.0,
        "logp_min": 0.0,
        "logp_max": 5.0,
        "lipinski_violations_max": 1,
        "admet_score_min": 0.6,
        "constitutional_compliance": {"L51_zero_gaps": True, "L10_harm_absolute": True},
    }


def get_nipah_dcm(strain: str) -> Dict[str, Any]:
    """Return Nipah DCM for NiV_Malaysia or NiV_Bangladesh. Raises ValueError if strain unknown."""
    if strain == "NiV_Malaysia":
        return create_nipah_malaysia_dcm()
    if strain == "NiV_Bangladesh":
        return create_nipah_bangladesh_dcm()
    raise ValueError(f"Unknown Nipah strain: {strain!r}. Use NiV_Malaysia or NiV_Bangladesh.")


def get_disease_constraints(disease_name: str) -> Dict[str, Any]:
    """
    Get disease-specific constraints

    Args:
        disease_name: Name of disease (e.g. "Chagas Disease", "NiV_Malaysia", "NiV_Bangladesh")

    Returns:
        dict: Constraint configuration
    """
    # Nipah strain-aware DCMs (wire from Nipah_Analysis)
    if disease_name in ("NiV_Malaysia", "NIPAH_MALAYSIA", "Nipah Malaysia"):
        return create_nipah_malaysia_dcm()
    if disease_name in ("NiV_Bangladesh", "NIPAH_BANGLADESH", "Nipah Bangladesh"):
        return create_nipah_bangladesh_dcm()
    # Default generic constraints
    return {
        "disease": disease_name,
        "ic50_max_um": 10.0,
        "selectivity_index_min": 10.0,
        "lipinski_violations_max": 1,
        "admet_score_min": 0.6
    }


__all__ = [
    "DiseaseConstraintModel",
    "get_disease_constraints",
    "create_chagas_dcm",
    "create_nipah_malaysia_dcm",
    "create_nipah_bangladesh_dcm",
    "get_nipah_dcm",
]
