"""
Data Sources Module
Provides data source configuration and utilities for repurposed molecule intake.
Aligns with intake_policy.json: ChEMBL, DrugBank, ClinicalTrials.gov, FDA_OrangeBook, PubChem, NCATS_Repurposing, ZINC.
Fetch scripts support ~24h-scale intake (PubChem 150k, ZINC 200k, CT.gov 100k, NCATS FRDB + optional --extra).
"""


def get_data_sources():
    """
    Returns available data sources for compound research.
    Includes all seed-library sources from intake policy.
    """
    return {
        "chembl": "ChEMBL Database",
        "drugbank": "DrugBank Database",
        "clinical_trials": "ClinicalTrials.gov",
        "fda_orange_book": "FDA Orange Book",
        "pubchem": "PubChem Database",
        "ncats_repurposing": "NCATS Repurposing",
        "zinc": "ZINC15/20 (commercially purchasable)",
    }


def get_enabled_sources_from_policy(repo_root=None):
    """
    Return source names that are enabled in intake_policy.json.
    Uses Intake_Policy.get_intake_policy when available.
    """
    try:
        from PX_System.foundation.Intake_Policy import get_intake_policy
        policy = get_intake_policy(repo_root)
        sources = policy.get("sources") or {}
        return [name for name, cfg in sources.items() if isinstance(cfg, dict) and cfg.get("enabled")]
    except Exception:
        return list(get_data_sources().keys())


__all__ = ["get_data_sources", "get_enabled_sources_from_policy"]
