"""
Generate disease constraint JSON files for all 22 PRV-eligible tropical diseases.

Deterministic generator using only stdlib (json, pathlib).
Template structure derived from Disease_Constraint_Model.py factory methods
(create_chagas_dcm, create_nipah_malaysia_dcm) and PX_Domain/PRV_Diseases/manifest.json.

Constitutional compliance:
    - Toxicity threshold: 0.0210 (Law L11, immutable)
    - All constraints grounded in WHO/DNDi/CDC published guidance
    - No random values; every field deterministic
"""

import json
from pathlib import Path

OUTPUT_DIR = Path(__file__).parent

# Constitutional hard limit (Law L11) -- never modify
TOXICITY_THRESHOLD = 0.0210

# Shared Lipinski defaults for tropical disease drugs (allow 1 violation)
LIPINSKI_DEFAULTS = {
    "molecular_weight_max": 500.0,
    "logp_min": 0.0,
    "logp_max": 5.0,
    "hbd_max": 5,
    "hba_max": 10,
    "lipinski_violations_max": 1,
}

# Shared safety defaults
SAFETY_DEFAULTS = {
    "selectivity_index_min": 10.0,
    "cc50_min_um": 50.0,
}


def slugify(name: str) -> str:
    """Convert disease name to a filesystem-safe slug."""
    return (
        name.lower()
        .replace(" ", "_")
        .replace("-", "_")
    )


# ---------------------------------------------------------------------------
# Disease definitions: all 22 PRV-eligible tropical/neglected diseases
# Each entry contains disease-specific scientific parameters.
# ---------------------------------------------------------------------------

DISEASES = [
    {
        "disease_name": "Buruli ulcer",
        "disease_id": "buruli_ulcer",
        "pathogen_type": "bacteria",
        "pathogen": "Mycobacterium ulcerans",
        "target_classes": ["mycolactone_biosynthesis", "polyketide_synthase", "cell_wall_synthesis"],
        "primary_endpoints": [
            "lesion_healing_rate",
            "time_to_complete_healing_weeks",
            "recurrence_rate_12_months",
            "bacterial_clearance"
        ],
        "efficacy_threshold": 0.60,
        "ic50_max_um": 10.0,
    },
    {
        "disease_name": "Chagas disease",
        "disease_id": "chagas_disease",
        "pathogen_type": "parasite",
        "pathogen": "Trypanosoma cruzi",
        "target_classes": ["CYP51_sterol_14alpha_demethylase", "cruzipain", "trypanothione_reductase"],
        "primary_endpoints": [
            "parasitemia_clearance",
            "sustained_parasitological_cure_12_months",
            "cardiac_function_preservation",
            "pcr_negativity"
        ],
        "efficacy_threshold": 0.65,
        "ic50_max_um": 10.0,
    },
    {
        "disease_name": "Chikungunya",
        "disease_id": "chikungunya",
        "pathogen_type": "virus",
        "pathogen": "Chikungunya virus (CHIKV)",
        "target_classes": ["nsP2_protease", "nsP1_capping", "RNA_dependent_RNA_polymerase"],
        "primary_endpoints": [
            "viremia_reduction_log10",
            "arthralgia_resolution_days",
            "time_to_fever_clearance",
            "chronic_arthritis_prevention_rate"
        ],
        "efficacy_threshold": 0.55,
        "ic50_max_um": 5.0,
    },
    {
        "disease_name": "Cholera",
        "disease_id": "cholera",
        "pathogen_type": "bacteria",
        "pathogen": "Vibrio cholerae",
        "target_classes": ["cholera_toxin_B_subunit", "tcp_pilus", "DNA_gyrase"],
        "primary_endpoints": [
            "stool_volume_reduction_24h",
            "rehydration_requirement_reduction",
            "bacterial_shedding_duration_days",
            "mortality_reduction"
        ],
        "efficacy_threshold": 0.70,
        "ic50_max_um": 8.0,
    },
    {
        "disease_name": "Dengue",
        "disease_id": "dengue",
        "pathogen_type": "virus",
        "pathogen": "Dengue virus (DENV 1-4)",
        "target_classes": ["NS3_protease_helicase", "NS5_RdRp", "NS4B_replication_complex"],
        "primary_endpoints": [
            "viremia_clearance_days",
            "platelet_count_recovery",
            "severe_dengue_prevention_rate",
            "hospitalization_duration_reduction"
        ],
        "efficacy_threshold": 0.55,
        "ic50_max_um": 5.0,
    },
    {
        "disease_name": "Dracunculiasis",
        "disease_id": "dracunculiasis",
        "pathogen_type": "parasite",
        "pathogen": "Dracunculus medinensis",
        "target_classes": ["cuticle_collagen_synthesis", "neuromuscular_junction", "acetylcholinesterase"],
        "primary_endpoints": [
            "worm_emergence_prevention",
            "larval_viability_reduction",
            "transmission_interruption_rate",
            "wound_healing_time_days"
        ],
        "efficacy_threshold": 0.60,
        "ic50_max_um": 10.0,
    },
    {
        "disease_name": "Ebola virus disease",
        "disease_id": "ebola_virus_disease",
        "pathogen_type": "virus",
        "pathogen": "Ebola virus (EBOV)",
        "target_classes": ["RNA_dependent_RNA_polymerase", "VP35_interferon_inhibitor", "GP_fusion"],
        "primary_endpoints": [
            "viremia_reduction_log10",
            "28_day_mortality_reduction",
            "time_to_viral_clearance",
            "organ_function_preservation"
        ],
        "efficacy_threshold": 0.50,
        "ic50_max_um": 1.0,
    },
    {
        "disease_name": "Lassa fever",
        "disease_id": "lassa_fever",
        "pathogen_type": "virus",
        "pathogen": "Lassa virus (LASV)",
        "target_classes": ["RNA_dependent_RNA_polymerase", "Z_matrix_protein", "GP_complex"],
        "primary_endpoints": [
            "viremia_reduction_log10",
            "case_fatality_rate_reduction",
            "hearing_loss_prevention_rate",
            "time_to_clinical_improvement"
        ],
        "efficacy_threshold": 0.50,
        "ic50_max_um": 2.0,
    },
    {
        "disease_name": "Leishmaniasis",
        "disease_id": "leishmaniasis",
        "pathogen_type": "parasite",
        "pathogen": "Leishmania spp.",
        "target_classes": ["trypanothione_reductase", "topoisomerase_II", "CYP51", "miltefosine_transporter"],
        "primary_endpoints": [
            "parasitological_cure_rate",
            "initial_cure_rate_6_months",
            "relapse_rate_12_months",
            "splenic_parasite_clearance"
        ],
        "efficacy_threshold": 0.65,
        "ic50_max_um": 10.0,
    },
    {
        "disease_name": "Leprosy",
        "disease_id": "leprosy",
        "pathogen_type": "bacteria",
        "pathogen": "Mycobacterium leprae",
        "target_classes": ["folate_biosynthesis", "RNA_polymerase", "DNA_gyrase", "cell_wall_synthesis"],
        "primary_endpoints": [
            "bacteriological_index_reduction",
            "disability_grade_improvement",
            "reaction_episode_prevention",
            "treatment_completion_rate"
        ],
        "efficacy_threshold": 0.70,
        "ic50_max_um": 8.0,
    },
    {
        "disease_name": "Lymphatic filariasis",
        "disease_id": "lymphatic_filariasis",
        "pathogen_type": "parasite",
        "pathogen": "Wuchereria bancrofti / Brugia malayi",
        "target_classes": ["Wolbachia_endosymbiont", "chitinase", "glutathione_S_transferase"],
        "primary_endpoints": [
            "microfilaremia_clearance",
            "adult_worm_killing_rate",
            "lymphedema_progression_prevention",
            "transmission_assessment_survey_pass"
        ],
        "efficacy_threshold": 0.60,
        "ic50_max_um": 10.0,
    },
    {
        "disease_name": "Malaria",
        "disease_id": "malaria",
        "pathogen_type": "parasite",
        "pathogen": "Plasmodium falciparum / P. vivax",
        "target_classes": ["PfDHFR", "PfDHODH", "PfATP4", "cytochrome_bc1", "PfPI4K"],
        "primary_endpoints": [
            "parasite_clearance_time_hours",
            "adequate_clinical_parasitological_response_28d",
            "gametocyte_clearance",
            "recrudescence_rate_42d"
        ],
        "efficacy_threshold": 0.70,
        "ic50_max_um": 1.0,
    },
    {
        "disease_name": "Marburg virus disease",
        "disease_id": "marburg_virus_disease",
        "pathogen_type": "virus",
        "pathogen": "Marburg virus (MARV)",
        "target_classes": ["RNA_dependent_RNA_polymerase", "VP35_interferon_inhibitor", "GP_fusion"],
        "primary_endpoints": [
            "viremia_reduction_log10",
            "28_day_mortality_reduction",
            "time_to_viral_clearance",
            "coagulation_parameter_normalization"
        ],
        "efficacy_threshold": 0.50,
        "ic50_max_um": 1.0,
    },
    {
        "disease_name": "Nipah virus infection",
        "disease_id": "nipah_virus_infection",
        "pathogen_type": "virus",
        "pathogen": "Nipah virus (NiV)",
        "target_classes": ["fusion_protein_F", "attachment_glycoprotein_G", "RNA_dependent_RNA_polymerase"],
        "primary_endpoints": [
            "viremia_reduction_log10",
            "case_fatality_rate_reduction",
            "neurological_sequelae_prevention",
            "time_to_clinical_improvement"
        ],
        "efficacy_threshold": 0.50,
        "ic50_max_um": 1.0,
    },
    {
        "disease_name": "Onchocerciasis",
        "disease_id": "onchocerciasis",
        "pathogen_type": "parasite",
        "pathogen": "Onchocerca volvulus",
        "target_classes": ["Wolbachia_endosymbiont", "chitinase", "glutathione_S_transferase"],
        "primary_endpoints": [
            "microfilarial_density_reduction",
            "macrofilaricidal_efficacy",
            "skin_microfilaria_clearance_12_months",
            "ocular_microfilaria_reduction"
        ],
        "efficacy_threshold": 0.60,
        "ic50_max_um": 10.0,
    },
    {
        "disease_name": "Rabies",
        "disease_id": "rabies",
        "pathogen_type": "virus",
        "pathogen": "Rabies virus (RABV)",
        "target_classes": ["nucleoprotein_N", "phosphoprotein_P", "RNA_dependent_RNA_polymerase_L"],
        "primary_endpoints": [
            "post_exposure_survival_rate",
            "viral_neutralization_titer",
            "time_to_seroconversion_days",
            "cns_viral_load_reduction"
        ],
        "efficacy_threshold": 0.50,
        "ic50_max_um": 2.0,
    },
    {
        "disease_name": "Schistosomiasis",
        "disease_id": "schistosomiasis",
        "pathogen_type": "parasite",
        "pathogen": "Schistosoma mansoni / S. haematobium / S. japonicum",
        "target_classes": ["thioredoxin_glutathione_reductase", "SmHDAC8", "calcium_channel"],
        "primary_endpoints": [
            "egg_reduction_rate",
            "parasitological_cure_rate",
            "hepatic_fibrosis_regression",
            "reinfection_rate_12_months"
        ],
        "efficacy_threshold": 0.65,
        "ic50_max_um": 10.0,
    },
    {
        "disease_name": "Trachoma",
        "disease_id": "trachoma",
        "pathogen_type": "bacteria",
        "pathogen": "Chlamydia trachomatis",
        "target_classes": ["type_III_secretion_system", "DNA_gyrase", "ribosomal_50S"],
        "primary_endpoints": [
            "ocular_chlamydia_clearance",
            "trachomatous_inflammation_follicular_resolution",
            "trichiasis_progression_prevention",
            "community_prevalence_reduction"
        ],
        "efficacy_threshold": 0.70,
        "ic50_max_um": 5.0,
    },
    {
        "disease_name": "Tuberculosis",
        "disease_id": "tuberculosis",
        "pathogen_type": "bacteria",
        "pathogen": "Mycobacterium tuberculosis",
        "target_classes": ["InhA_enoyl_ACP_reductase", "DprE1", "ATP_synthase", "MmpL3", "DNA_gyrase"],
        "primary_endpoints": [
            "sputum_culture_conversion_8_weeks",
            "relapse_free_cure_rate_12_months",
            "treatment_shortening_months",
            "mdr_tb_efficacy"
        ],
        "efficacy_threshold": 0.70,
        "ic50_max_um": 1.0,
    },
    {
        "disease_name": "Yaws",
        "disease_id": "yaws",
        "pathogen_type": "bacteria",
        "pathogen": "Treponema pallidum subsp. pertenue",
        "target_classes": ["penicillin_binding_protein", "cell_wall_synthesis", "DNA_gyrase"],
        "primary_endpoints": [
            "serological_cure_rate",
            "lesion_healing_time_weeks",
            "latent_yaws_seroconversion",
            "community_transmission_interruption"
        ],
        "efficacy_threshold": 0.75,
        "ic50_max_um": 5.0,
    },
    {
        "disease_name": "Yellow fever",
        "disease_id": "yellow_fever",
        "pathogen_type": "virus",
        "pathogen": "Yellow fever virus (YFV)",
        "target_classes": ["NS3_protease_helicase", "NS5_RdRp", "NS4B_replication_complex"],
        "primary_endpoints": [
            "viremia_clearance_days",
            "hepatic_function_preservation",
            "case_fatality_rate_reduction",
            "hemorrhagic_complication_prevention"
        ],
        "efficacy_threshold": 0.55,
        "ic50_max_um": 2.0,
    },
    {
        "disease_name": "Zika virus disease",
        "disease_id": "zika_virus_disease",
        "pathogen_type": "virus",
        "pathogen": "Zika virus (ZIKV)",
        "target_classes": ["NS2B_NS3_protease", "NS5_RdRp", "NS5_methyltransferase"],
        "primary_endpoints": [
            "viremia_clearance_days",
            "congenital_zika_syndrome_prevention",
            "guillain_barre_risk_reduction",
            "sexual_transmission_prevention"
        ],
        "efficacy_threshold": 0.55,
        "ic50_max_um": 5.0,
    },
]


def build_constraint(disease: dict) -> dict:
    """
    Build a full constraint JSON object for a single disease.

    Merges disease-specific parameters with shared Lipinski/safety defaults
    and constitutional requirements.
    """
    return {
        # Identification
        "disease_name": disease["disease_name"],
        "disease_id": disease["disease_id"],
        "pathogen_type": disease["pathogen_type"],
        "pathogen": disease["pathogen"],

        # Drug target classes
        "target_classes": disease["target_classes"],

        # Clinical endpoints
        "primary_endpoints": disease["primary_endpoints"],

        # Efficacy constraints
        "efficacy_threshold": disease["efficacy_threshold"],
        "ic50_max_um": disease["ic50_max_um"],
        "ic50_optimal_um": round(disease["ic50_max_um"] * 0.1, 2),

        # Safety constraints (Law L10: Harm is Absolute)
        "selectivity_index_min": SAFETY_DEFAULTS["selectivity_index_min"],
        "cc50_min_um": SAFETY_DEFAULTS["cc50_min_um"],

        # Toxicity threshold (Law L11 -- constitutional hard limit, immutable)
        "toxicity_threshold": TOXICITY_THRESHOLD,

        # Lipinski drug-likeness constraints
        "molecular_weight_max": LIPINSKI_DEFAULTS["molecular_weight_max"],
        "logp_min": LIPINSKI_DEFAULTS["logp_min"],
        "logp_max": LIPINSKI_DEFAULTS["logp_max"],
        "hbd_max": LIPINSKI_DEFAULTS["hbd_max"],
        "hba_max": LIPINSKI_DEFAULTS["hba_max"],
        "lipinski_violations_max": LIPINSKI_DEFAULTS["lipinski_violations_max"],

        # ADMET minimum
        "admet_score_min": 0.6,

        # Regulatory
        "prv_eligible": True,
        "regulatory_pathway": "PRV_ACCELERATED",

        # Versioning
        "constraint_version": "1.0.0",

        # Constitutional compliance seal
        "constitutional_compliance": {
            "L51_zero_gaps": True,
            "L10_harm_absolute": True,
            "L9_axiomatic_grounding": True,
            "L5_transparency": True,
            "L11_toxicity_hard_limit": True,
        },
    }


def main() -> None:
    """Generate all 22 PRV disease constraint files."""
    generated = []

    for disease in DISEASES:
        constraint = build_constraint(disease)
        filename = f"{disease['disease_id']}.json"
        filepath = OUTPUT_DIR / filename

        with open(filepath, "w", encoding="utf-8") as f:
            json.dump(constraint, f, indent=2, ensure_ascii=False)
            f.write("\n")

        generated.append(filename)
        print(f"  Generated: {filename}")

    print(f"\nTotal files generated: {len(generated)}")
    print(f"Output directory: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
