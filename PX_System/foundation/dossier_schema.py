"""
Dossier Schema v1.0 — Pharma-grade dossier structure.

Only DIAMOND tier compounds receive full dossier packages.
GOLD receives summary reports. SILVER/BRONZE receive WorldLines only.

Constitutional: Python stdlib only. Deterministic.
"""

from __future__ import annotations

import json
import uuid
from dataclasses import dataclass, field, asdict
from datetime import datetime
from enum import Enum
from typing import Any, Dict, List, Optional, Tuple

DOSSIER_SCHEMA_VERSION = "1.0.0"


# ---------------------------------------------------------------------------
# Enumerations
# ---------------------------------------------------------------------------

class Tier(Enum):
    DIAMOND = "DIAMOND"
    GOLD = "GOLD"
    SILVER = "SILVER"
    BRONZE = "BRONZE"


class RiskLevel(Enum):
    LOW = "LOW"
    MEDIUM = "MEDIUM"
    HIGH = "HIGH"


class Recommendation(Enum):
    ADVANCE = "ADVANCE_TO_IND"
    OPTIMIZE = "OPTIMIZE_FURTHER"
    HOLD = "HOLD_FOR_DATA"
    REJECT = "REJECT"


def _json_default(obj):
    """Custom JSON serializer that extracts Enum .value instead of str()."""
    if isinstance(obj, Enum):
        return obj.value
    return str(obj)


# ===========================================================================
# SECTION 00 — EXECUTIVE SUMMARY
# ===========================================================================

@dataclass
class MetricScore:
    """Single scored metric with pass/fail against a target."""

    name: str
    value: float
    target: float
    unit: str
    passed: bool
    display_bar_percent: int  # 0-100


@dataclass
class RiskAssessment:
    """Risk factor in one of four categories."""

    category: str  # efficacy / safety / regulatory / commercial
    level: RiskLevel
    factors: List[str]
    mitigations: List[str]


@dataclass
class GovernanceMetadata:
    """Governance provenance attached to every dossier artifact."""

    zeus_gate_status: str
    trace_id: str
    worldline_id: str
    slc_hash: str
    constitutional_compliance: bool
    laws_verified: int


@dataclass
class ExecutiveSummary:
    """Top-level investment-grade summary of the compound."""

    compound_id: str
    compound_name: Optional[str]
    disease_name: str
    disease_id: str
    generation_date: str
    tier: Tier
    investment_thesis: str
    key_metrics: List[MetricScore]
    recommendation: Recommendation
    confidence_percent: int
    cost_to_ind_usd_low: int
    cost_to_ind_usd_high: int
    months_to_ind: int
    risks: List[RiskAssessment]
    governance: GovernanceMetadata


# ===========================================================================
# SECTION 01 — MOLECULE PROFILE
# ===========================================================================

@dataclass
class StructuralAlert:
    """Flagged structural moiety that may affect safety or patentability."""

    alert_name: str
    alert_type: str  # PAINS / Brenk / toxicophore / reactive
    severity: str
    description: str


@dataclass
class SimilarityMatch:
    """Known drug structurally similar to the candidate."""

    drug_name: str
    chembl_id: Optional[str]
    tanimoto_similarity: float
    therapeutic_class: str
    approval_status: str


@dataclass
class PhysicochemicalProperties:
    """Calculated physicochemical descriptors."""

    logp: float
    logd_7_4: Optional[float]
    pka_acidic: Optional[float]
    pka_basic: Optional[float]
    aqueous_solubility_ug_ml: Optional[float]
    psa: float
    hbd: int
    hba: int
    rotatable_bonds: int
    aromatic_rings: int
    heavy_atoms: int
    fraction_sp3: float
    molar_refractivity: Optional[float]


@dataclass
class DruglikenessAssessment:
    """Multi-rule drug-likeness evaluation."""

    lipinski_violations: int
    lipinski_pass: bool
    veber_pass: bool
    ghose_pass: bool
    lead_likeness_pass: bool
    qed_score: float
    cns_mpo_score: Optional[float]


@dataclass
class MoleculeProfile:
    """Complete molecular identity and property profile."""

    smiles: str
    canonical_smiles: str
    inchi: str
    inchi_key: str
    molecular_formula: str
    molecular_weight: float
    exact_mass: float
    stereochemistry: str
    structure_image_base64: Optional[str]
    physicochemical: PhysicochemicalProperties
    druglikeness: DruglikenessAssessment
    structural_alerts: List[StructuralAlert]
    nearest_known_drugs: List[SimilarityMatch]


# ===========================================================================
# SECTION 02 — TARGET VALIDATION
# ===========================================================================

@dataclass
class TargetInteraction:
    """Single molecular interaction between ligand and target residue."""

    interaction_type: str
    residue: str
    distance_angstrom: float
    strength: str


@dataclass
class BindingAnalysis:
    """Computational binding pose and affinity analysis."""

    binding_site_name: str
    binding_mode: str
    docking_score: float
    docking_software: str
    key_interactions: List[TargetInteraction]
    predicted_kd_nM: Optional[float]


@dataclass
class OffTarget:
    """Predicted off-target activity and selectivity concern."""

    target_name: str
    uniprot_id: str
    predicted_activity: float
    selectivity_ratio: float
    concern_level: str


@dataclass
class TargetValidation:
    """Full target identification, druggability, and binding profile."""

    target_name: str
    gene_symbol: str
    uniprot_id: str
    target_class: str
    organism: str
    druggability_score: float
    expression_profile: Dict[str, str]
    genetic_associations: List[str]
    known_inhibitors: List[str]
    clinical_precedent: bool
    binding_analysis: BindingAnalysis
    selectivity_profile: List[OffTarget]
    validation_confidence: str
    rationale: str


# ===========================================================================
# SECTION 03 — EFFICACY DATA
# ===========================================================================

@dataclass
class ActivityPrediction:
    """Predicted biological activity with confidence interval and model info."""

    activity_type: str
    value_nM: float
    confidence_low_nM: float
    confidence_high_nM: float
    model_name: str
    model_version: str
    training_set_size: int
    applicability_domain: str
    similarity_to_training: float


@dataclass
class EfficacyBenchmark:
    """Reference compound for efficacy comparison."""

    compound_name: str
    activity_nM: float
    source: str
    our_compound_comparison: str


@dataclass
class EfficacyData:
    """Complete efficacy evidence package."""

    primary_activity: ActivityPrediction
    secondary_activities: List[ActivityPrediction]
    benchmarks: List[EfficacyBenchmark]
    percentile_rank: int
    dose_projection_mg: float
    dose_projection_basis: str
    therapeutic_index: float
    resistance_liability: str
    resistance_mutations_known: List[str]
    efficacy_confidence: str
    key_uncertainties: List[str]
    recommended_validation_experiments: List[str]


# ===========================================================================
# SECTION 04 — SAFETY PROFILE
# ===========================================================================

@dataclass
class OrganToxicity:
    """Organ-specific toxicity prediction."""

    organ: str
    prediction_score: float
    threshold: float
    passed: bool
    confidence: str
    flags: List[str]


@dataclass
class CardiacSafety:
    """Cardiac liability assessment including hERG and QT risk."""

    herg_ic50_uM: float
    herg_threshold_uM: float
    herg_passed: bool
    herg_liability_class: str
    qt_prolongation_risk: str
    cardiovascular_safety_margin: float


@dataclass
class CYPInhibition:
    """CYP enzyme inhibition assessment for DDI risk."""

    cyp_enzyme: str
    ic50_uM: float
    inhibition_class: str
    ddi_risk: str


@dataclass
class GenotoxicityPrediction:
    """In-silico genotoxicity assay prediction."""

    assay: str
    prediction: str
    confidence: float


@dataclass
class SafetyProfile:
    """Complete safety and toxicology profile."""

    overall_toxicity_index: float
    tier_classification: str
    safety_margin_fold: float
    organ_toxicity: List[OrganToxicity]
    cardiac_safety: CardiacSafety
    cyp_inhibition: List[CYPInhibition]
    cyp_induction_risk: str
    genotoxicity: List[GenotoxicityPrediction]
    max_tolerated_dose_mg_kg: Optional[float]
    dose_limiting_toxicity: Optional[str]
    safety_confidence: str
    recommended_safety_studies: List[str]
    monitoring_parameters: List[str]


# ===========================================================================
# SECTION 05 — PKPD ANALYSIS
# ===========================================================================

@dataclass
class AbsorptionProfile:
    """Oral absorption and permeability assessment."""

    oral_bioavailability_percent: float
    bioavailability_confidence: float
    caco2_permeability_nm_s: float
    pgp_substrate: bool
    pgp_inhibitor: bool
    food_effect: str


@dataclass
class DistributionProfile:
    """Volume of distribution and tissue penetration."""

    vd_L_kg: float
    plasma_protein_binding_percent: float
    blood_brain_barrier: bool
    bbb_permeability: Optional[float]
    tissue_distribution: Dict[str, str]


@dataclass
class MetabolismProfile:
    """Metabolic pathway and stability assessment."""

    primary_cyp: str
    secondary_cyps: List[str]
    hlm_half_life_min: float
    rlm_half_life_min: Optional[float]
    metabolic_stability: str
    active_metabolites: bool
    metabolic_soft_spots: List[str]


@dataclass
class ExcretionProfile:
    """Clearance and elimination route characterization."""

    clearance_mL_min_kg: float
    primary_route: str
    renal_excretion_percent: Optional[float]
    fecal_excretion_percent: Optional[float]


@dataclass
class PKParameters:
    """Predicted human PK parameters with confidence bounds."""

    cmax_ng_mL: float
    cmax_confidence_percent: float
    tmax_hours: float
    tmax_confidence_hours: float
    auc_0_24_ng_h_mL: float
    auc_confidence_percent: float
    half_life_hours: float
    half_life_confidence_percent: float


@dataclass
class DoseProjection:
    """Human dose prediction for Phase 1 starting dose."""

    target_concentration_ng_mL: float
    target_basis: str
    predicted_dose_mg: float
    dosing_frequency: str
    dose_range_phase1_low_mg: float
    dose_range_phase1_high_mg: float
    max_recommended_starting_dose_mg: float
    mrsd_basis: str


@dataclass
class PKPDAnalysis:
    """Integrated pharmacokinetic-pharmacodynamic analysis."""

    absorption: AbsorptionProfile
    distribution: DistributionProfile
    metabolism: MetabolismProfile
    excretion: ExcretionProfile
    pk_parameters_human_predicted: PKParameters
    dose_projection: DoseProjection
    pk_model_type: str
    pk_confidence: str
    pk_uncertainties: List[str]


# ===========================================================================
# SECTION 06 — MANUFACTURING ASSESSMENT
# ===========================================================================

@dataclass
class SynthesisStep:
    """Single step in the proposed synthetic route."""

    step_number: int
    transformation: str
    reagents: List[str]
    conditions: str
    yield_percent: Optional[float]
    concerns: List[str]


@dataclass
class StartingMaterial:
    """Starting material availability and supply chain assessment."""

    name: str
    cas_number: Optional[str]
    commercial_availability: bool
    estimated_cost_per_kg_usd: float
    suppliers: int
    supply_chain_risk: str


@dataclass
class ManufacturingAssessment:
    """Complete CMC feasibility and cost assessment."""

    synthesis_steps_count: int
    longest_linear_sequence: int
    overall_yield_percent: float
    complexity_score: str
    synthesis_route: List[SynthesisStep]
    chiral_centers: int
    resolution_required: bool
    starting_materials: List[StartingMaterial]
    hazardous_reagents: List[str]
    cryogenic_steps: bool
    high_pressure_steps: bool
    chromatography_required: bool
    api_cost_per_kg_usd_low: float
    api_cost_per_kg_usd_high: float
    formulated_cost_per_dose_usd: float
    cogs_confidence: str
    cogs_drivers: List[str]
    recommended_dosage_form: str
    excipient_compatibility: str
    stability_months_25C: int
    special_requirements: List[str]
    cmc_risk: str
    cmc_risks_identified: List[str]
    cmc_mitigations: List[str]


# ===========================================================================
# SECTION 07 — REGULATORY STRATEGY
# ===========================================================================

@dataclass
class DesignationEligibility:
    """Eligibility assessment for a specific regulatory designation."""

    designation: str
    eligible: bool
    confidence: str
    rationale: str


@dataclass
class RegulatoryPrecedent:
    """Historical regulatory precedent relevant to this program."""

    drug_name: str
    approval_year: int
    indication: str
    pathway_used: str
    relevance: str


@dataclass
class RegulatoryTimeline:
    """Single milestone in the projected regulatory timeline."""

    milestone: str
    month_from_start: int
    dependencies: List[str]


@dataclass
class RegulatoryStrategy:
    """Complete regulatory pathway and strategy."""

    designations: List[DesignationEligibility]
    prv_eligible: bool
    prv_qualifying_disease: str
    prv_value_estimate_usd_millions: float
    recommended_pathway: str
    alternative_pathway: str
    animal_rule_applicable: bool
    biomarker_strategy: str
    primary_endpoint_phase3: str
    surrogate_endpoint_available: bool
    precedents: List[RegulatoryPrecedent]
    relevant_fda_guidances: List[str]
    timeline: List[RegulatoryTimeline]
    regulatory_risk: str
    regulatory_hurdles: List[str]
    regulatory_mitigations: List[str]


# ===========================================================================
# SECTION 08 — COMPETITIVE LANDSCAPE
# ===========================================================================

@dataclass
class ApprovedDrug:
    """Approved drug in the competitive landscape."""

    name: str
    company: str
    approval_year: int
    mechanism: str
    limitations: List[str]
    annual_sales_usd_millions: Optional[float]


@dataclass
class ClinicalCandidate:
    """Drug candidate in clinical development."""

    name: str
    company: str
    phase: str
    mechanism: str
    status: str
    results_summary: Optional[str]


@dataclass
class CompetitiveLandscape:
    """Market landscape, pipeline, and commercial opportunity."""

    disease_prevalence: str
    geographic_distribution: List[str]
    current_standard_of_care: str
    unmet_need_level: str
    approved_drugs: List[ApprovedDrug]
    clinical_pipeline: List[ClinicalCandidate]
    recent_failures: List[str]
    our_advantages: List[str]
    our_disadvantages: List[str]
    market_size_peak_sales_usd_millions: float
    prv_value_usd_millions: float
    stockpile_opportunity_usd_millions: float
    total_opportunity_usd_millions: float
    active_acquirers: List[str]
    recent_transactions: List[str]
    valuation_benchmarks: List[str]


# ===========================================================================
# SECTION 09 — PATENT ANALYSIS
# ===========================================================================

@dataclass
class BlockingPatent:
    """Patent that may block freedom-to-operate."""

    patent_number: str
    assignee: str
    title: str
    claims_summary: str
    expiry_date: str
    geographic_coverage: List[str]
    design_around_possible: bool


@dataclass
class PatentabilityAssessment:
    """Patentability evaluation for a specific claim category."""

    category: str
    patentable: str
    prior_art_concerns: List[str]


@dataclass
class PatentAnalysis:
    """Freedom-to-operate and IP strategy."""

    fto_status: str
    blocking_patents_count: int
    blocking_patents: List[BlockingPatent]
    design_around_required: bool
    patent_landscape_families: int
    key_assignees: List[str]
    patentability: List[PatentabilityAssessment]
    recommended_filings: List[str]
    geographic_strategy: List[str]
    orphan_exclusivity_years: int
    pediatric_exclusivity_possible: bool
    ip_risk: str
    ip_risks_identified: List[str]
    ip_mitigations: List[str]


# ===========================================================================
# SECTION 10 — CLINICAL TRIAL DESIGN
# ===========================================================================

@dataclass
class TrialDesign:
    """Single-phase clinical trial design parameters."""

    phase: str
    study_type: str
    population: str
    n_subjects: int
    design: str
    dose_range_mg: str
    primary_objectives: List[str]
    primary_endpoint: str
    secondary_endpoints: List[str]
    duration_weeks: int
    estimated_cost_usd_millions: float


@dataclass
class BiomarkerStrategy:
    """Biomarker plan for pharmacodynamics, efficacy, and safety."""

    pharmacodynamic_biomarkers: List[str]
    efficacy_biomarkers: List[str]
    safety_biomarkers: List[str]
    companion_diagnostic_required: bool


@dataclass
class PatientPopulation:
    """Target patient population and enrollment strategy."""

    inclusion_criteria_key: List[str]
    exclusion_criteria_key: List[str]
    recruitment_considerations: str
    geographic_strategy: List[str]
    enrollment_rate_per_month: Optional[int]


@dataclass
class ClinicalTrialDesign:
    """Integrated Phase 1/2 clinical development plan."""

    phase1_design: TrialDesign
    phase2_design: TrialDesign
    biomarker_strategy: BiomarkerStrategy
    patient_population: PatientPopulation
    outbreak_dependent: bool
    adaptive_design_elements: List[str]
    compassionate_use_pathway: str
    total_cost_to_phase2_usd_millions: float
    total_months_to_phase2_readout: int


# ===========================================================================
# COMPLETE DOSSIER PACKAGE (DIAMOND tier)
# ===========================================================================

@dataclass
class DossierPackage:
    """
    Full DIAMOND-tier pharmaceutical dossier package.

    Contains all 11 sections required for an IND-enabling decision.
    Only compounds that achieve DIAMOND classification receive this
    complete package. Includes built-in validation and serialization.
    """

    schema_version: str
    package_id: str
    generation_timestamp: str
    tier: Tier
    executive_summary: ExecutiveSummary
    molecule_profile: MoleculeProfile
    target_validation: TargetValidation
    efficacy_data: EfficacyData
    safety_profile: SafetyProfile
    pkpd_analysis: PKPDAnalysis
    manufacturing_assessment: ManufacturingAssessment
    regulatory_strategy: RegulatoryStrategy
    competitive_landscape: CompetitiveLandscape
    patent_analysis: PatentAnalysis
    clinical_trial_design: ClinicalTrialDesign
    governance: GovernanceMetadata
    worldline_id: str
    slc_entries: List[str]

    def validate(self) -> Tuple[bool, List[str]]:
        """Validate required fields and structural integrity.

        Returns:
            Tuple of (is_valid, list_of_error_messages).
            An empty error list means the package is valid.
        """
        errors: List[str] = []

        # Schema version must match
        if self.schema_version != DOSSIER_SCHEMA_VERSION:
            errors.append(
                f"Schema version mismatch: got '{self.schema_version}', "
                f"expected '{DOSSIER_SCHEMA_VERSION}'"
            )

        # Package ID must be present
        if not self.package_id or not self.package_id.strip():
            errors.append("package_id is empty")

        # Generation timestamp must be present
        if not self.generation_timestamp or not self.generation_timestamp.strip():
            errors.append("generation_timestamp is empty")

        # Tier must be DIAMOND for a full dossier package
        if self.tier != Tier.DIAMOND:
            errors.append(
                f"DossierPackage requires DIAMOND tier, got {self.tier.value}"
            )

        # Executive summary checks
        if not self.executive_summary.compound_id:
            errors.append("executive_summary.compound_id is empty")
        if not self.executive_summary.disease_name:
            errors.append("executive_summary.disease_name is empty")
        if not self.executive_summary.disease_id:
            errors.append("executive_summary.disease_id is empty")
        if not self.executive_summary.key_metrics:
            errors.append("executive_summary.key_metrics is empty")
        if not (0 <= self.executive_summary.confidence_percent <= 100):
            errors.append(
                f"executive_summary.confidence_percent out of range: "
                f"{self.executive_summary.confidence_percent}"
            )

        # Molecule profile checks
        if not self.molecule_profile.smiles:
            errors.append("molecule_profile.smiles is empty")
        if not self.molecule_profile.canonical_smiles:
            errors.append("molecule_profile.canonical_smiles is empty")
        if not self.molecule_profile.inchi_key:
            errors.append("molecule_profile.inchi_key is empty")
        if self.molecule_profile.molecular_weight <= 0:
            errors.append(
                f"molecule_profile.molecular_weight must be positive: "
                f"{self.molecule_profile.molecular_weight}"
            )

        # Safety profile: toxicity hard limit (Law L11)
        tox = self.safety_profile.overall_toxicity_index
        if tox >= 0.0210:
            errors.append(
                f"safety_profile.overall_toxicity_index exceeds hard limit: "
                f"{tox} >= 0.0210 (Law L11 TOXICITY_FAILURE)"
            )

        # Safety margin must be positive
        if self.safety_profile.safety_margin_fold <= 0:
            errors.append(
                f"safety_profile.safety_margin_fold must be positive: "
                f"{self.safety_profile.safety_margin_fold}"
            )

        # PKPD checks
        if self.pkpd_analysis.absorption.oral_bioavailability_percent < 0:
            errors.append(
                "pkpd_analysis.absorption.oral_bioavailability_percent is negative"
            )
        if self.pkpd_analysis.pk_parameters_human_predicted.half_life_hours <= 0:
            errors.append(
                "pkpd_analysis.pk_parameters_human_predicted.half_life_hours "
                "must be positive"
            )

        # Manufacturing checks
        if self.manufacturing_assessment.synthesis_steps_count <= 0:
            errors.append(
                "manufacturing_assessment.synthesis_steps_count must be positive"
            )

        # Governance checks
        if not self.governance.zeus_gate_status:
            errors.append("governance.zeus_gate_status is empty")
        if not self.governance.constitutional_compliance:
            errors.append("governance.constitutional_compliance is False")
        if not self.governance.trace_id:
            errors.append("governance.trace_id is empty")
        if not self.governance.slc_hash:
            errors.append("governance.slc_hash is empty")

        # WorldLine ID must be present
        if not self.worldline_id or not self.worldline_id.strip():
            errors.append("worldline_id is empty")

        # Target validation checks
        if not self.target_validation.target_name:
            errors.append("target_validation.target_name is empty")
        if not self.target_validation.gene_symbol:
            errors.append("target_validation.gene_symbol is empty")
        if not (0.0 <= self.target_validation.druggability_score <= 1.0):
            errors.append(
                f"target_validation.druggability_score out of range [0,1]: "
                f"{self.target_validation.druggability_score}"
            )

        # Efficacy checks
        if self.efficacy_data.primary_activity.value_nM <= 0:
            errors.append(
                "efficacy_data.primary_activity.value_nM must be positive"
            )
        if self.efficacy_data.therapeutic_index <= 0:
            errors.append(
                "efficacy_data.therapeutic_index must be positive"
            )

        # Patent / FTO status must be present
        if not self.patent_analysis.fto_status:
            errors.append("patent_analysis.fto_status is empty")

        # Regulatory strategy checks
        if not self.regulatory_strategy.recommended_pathway:
            errors.append("regulatory_strategy.recommended_pathway is empty")

        # Clinical trial design checks
        if self.clinical_trial_design.phase1_design.n_subjects <= 0:
            errors.append(
                "clinical_trial_design.phase1_design.n_subjects must be positive"
            )
        if self.clinical_trial_design.phase2_design.n_subjects <= 0:
            errors.append(
                "clinical_trial_design.phase2_design.n_subjects must be positive"
            )

        is_valid = len(errors) == 0
        return (is_valid, errors)

    def to_json(self) -> str:
        """Serialize the complete dossier package to JSON.

        Uses ``default=str`` to handle Enum members, datetime objects,
        and any other non-primitive types encountered during serialization.
        """
        return json.dumps(asdict(self), indent=2, default=_json_default)

    @classmethod
    def from_dict(cls, data: dict) -> DossierPackage:
        """Reconstruct a DossierPackage from a plain dictionary.

        Handles nested dataclass reconstruction and Enum conversion.
        """
        # --- helper converters -------------------------------------------
        def _enum(enum_cls: type, val: Any) -> Any:
            if isinstance(val, enum_cls):
                return val
            if isinstance(val, str):
                # Try value first, then name
                for member in enum_cls:
                    if member.value == val or member.name == val:
                        return member
            return val

        def _metric_score(d: dict) -> MetricScore:
            return MetricScore(**d)

        def _risk_assessment(d: dict) -> RiskAssessment:
            return RiskAssessment(
                category=d["category"],
                level=_enum(RiskLevel, d["level"]),
                factors=d["factors"],
                mitigations=d["mitigations"],
            )

        def _governance(d: dict) -> GovernanceMetadata:
            return GovernanceMetadata(**d)

        def _executive_summary(d: dict) -> ExecutiveSummary:
            return ExecutiveSummary(
                compound_id=d["compound_id"],
                compound_name=d.get("compound_name"),
                disease_name=d["disease_name"],
                disease_id=d["disease_id"],
                generation_date=d["generation_date"],
                tier=_enum(Tier, d["tier"]),
                investment_thesis=d["investment_thesis"],
                key_metrics=[_metric_score(m) for m in d["key_metrics"]],
                recommendation=_enum(Recommendation, d["recommendation"]),
                confidence_percent=d["confidence_percent"],
                cost_to_ind_usd_low=d["cost_to_ind_usd_low"],
                cost_to_ind_usd_high=d["cost_to_ind_usd_high"],
                months_to_ind=d["months_to_ind"],
                risks=[_risk_assessment(r) for r in d["risks"]],
                governance=_governance(d["governance"]),
            )

        def _structural_alert(d: dict) -> StructuralAlert:
            return StructuralAlert(**d)

        def _similarity_match(d: dict) -> SimilarityMatch:
            return SimilarityMatch(**d)

        def _physicochemical(d: dict) -> PhysicochemicalProperties:
            return PhysicochemicalProperties(**d)

        def _druglikeness(d: dict) -> DruglikenessAssessment:
            return DruglikenessAssessment(**d)

        def _molecule_profile(d: dict) -> MoleculeProfile:
            return MoleculeProfile(
                smiles=d["smiles"],
                canonical_smiles=d["canonical_smiles"],
                inchi=d["inchi"],
                inchi_key=d["inchi_key"],
                molecular_formula=d["molecular_formula"],
                molecular_weight=d["molecular_weight"],
                exact_mass=d["exact_mass"],
                stereochemistry=d["stereochemistry"],
                structure_image_base64=d.get("structure_image_base64"),
                physicochemical=_physicochemical(d["physicochemical"]),
                druglikeness=_druglikeness(d["druglikeness"]),
                structural_alerts=[
                    _structural_alert(a) for a in d["structural_alerts"]
                ],
                nearest_known_drugs=[
                    _similarity_match(m) for m in d["nearest_known_drugs"]
                ],
            )

        def _target_interaction(d: dict) -> TargetInteraction:
            return TargetInteraction(**d)

        def _binding_analysis(d: dict) -> BindingAnalysis:
            return BindingAnalysis(
                binding_site_name=d["binding_site_name"],
                binding_mode=d["binding_mode"],
                docking_score=d["docking_score"],
                docking_software=d["docking_software"],
                key_interactions=[
                    _target_interaction(i) for i in d["key_interactions"]
                ],
                predicted_kd_nM=d.get("predicted_kd_nM"),
            )

        def _off_target(d: dict) -> OffTarget:
            return OffTarget(**d)

        def _target_validation(d: dict) -> TargetValidation:
            return TargetValidation(
                target_name=d["target_name"],
                gene_symbol=d["gene_symbol"],
                uniprot_id=d["uniprot_id"],
                target_class=d["target_class"],
                organism=d["organism"],
                druggability_score=d["druggability_score"],
                expression_profile=d["expression_profile"],
                genetic_associations=d["genetic_associations"],
                known_inhibitors=d["known_inhibitors"],
                clinical_precedent=d["clinical_precedent"],
                binding_analysis=_binding_analysis(d["binding_analysis"]),
                selectivity_profile=[
                    _off_target(o) for o in d["selectivity_profile"]
                ],
                validation_confidence=d["validation_confidence"],
                rationale=d["rationale"],
            )

        def _activity_prediction(d: dict) -> ActivityPrediction:
            return ActivityPrediction(**d)

        def _efficacy_benchmark(d: dict) -> EfficacyBenchmark:
            return EfficacyBenchmark(**d)

        def _efficacy_data(d: dict) -> EfficacyData:
            return EfficacyData(
                primary_activity=_activity_prediction(d["primary_activity"]),
                secondary_activities=[
                    _activity_prediction(a)
                    for a in d["secondary_activities"]
                ],
                benchmarks=[
                    _efficacy_benchmark(b) for b in d["benchmarks"]
                ],
                percentile_rank=d["percentile_rank"],
                dose_projection_mg=d["dose_projection_mg"],
                dose_projection_basis=d["dose_projection_basis"],
                therapeutic_index=d["therapeutic_index"],
                resistance_liability=d["resistance_liability"],
                resistance_mutations_known=d["resistance_mutations_known"],
                efficacy_confidence=d["efficacy_confidence"],
                key_uncertainties=d["key_uncertainties"],
                recommended_validation_experiments=d[
                    "recommended_validation_experiments"
                ],
            )

        def _organ_toxicity(d: dict) -> OrganToxicity:
            return OrganToxicity(**d)

        def _cardiac_safety(d: dict) -> CardiacSafety:
            return CardiacSafety(**d)

        def _cyp_inhibition(d: dict) -> CYPInhibition:
            return CYPInhibition(**d)

        def _genotoxicity(d: dict) -> GenotoxicityPrediction:
            return GenotoxicityPrediction(**d)

        def _safety_profile(d: dict) -> SafetyProfile:
            return SafetyProfile(
                overall_toxicity_index=d["overall_toxicity_index"],
                tier_classification=d["tier_classification"],
                safety_margin_fold=d["safety_margin_fold"],
                organ_toxicity=[
                    _organ_toxicity(o) for o in d["organ_toxicity"]
                ],
                cardiac_safety=_cardiac_safety(d["cardiac_safety"]),
                cyp_inhibition=[
                    _cyp_inhibition(c) for c in d["cyp_inhibition"]
                ],
                cyp_induction_risk=d["cyp_induction_risk"],
                genotoxicity=[
                    _genotoxicity(g) for g in d["genotoxicity"]
                ],
                max_tolerated_dose_mg_kg=d.get("max_tolerated_dose_mg_kg"),
                dose_limiting_toxicity=d.get("dose_limiting_toxicity"),
                safety_confidence=d["safety_confidence"],
                recommended_safety_studies=d["recommended_safety_studies"],
                monitoring_parameters=d["monitoring_parameters"],
            )

        def _absorption(d: dict) -> AbsorptionProfile:
            return AbsorptionProfile(**d)

        def _distribution(d: dict) -> DistributionProfile:
            return DistributionProfile(**d)

        def _metabolism(d: dict) -> MetabolismProfile:
            return MetabolismProfile(**d)

        def _excretion(d: dict) -> ExcretionProfile:
            return ExcretionProfile(**d)

        def _pk_parameters(d: dict) -> PKParameters:
            return PKParameters(**d)

        def _dose_projection(d: dict) -> DoseProjection:
            return DoseProjection(**d)

        def _pkpd_analysis(d: dict) -> PKPDAnalysis:
            return PKPDAnalysis(
                absorption=_absorption(d["absorption"]),
                distribution=_distribution(d["distribution"]),
                metabolism=_metabolism(d["metabolism"]),
                excretion=_excretion(d["excretion"]),
                pk_parameters_human_predicted=_pk_parameters(
                    d["pk_parameters_human_predicted"]
                ),
                dose_projection=_dose_projection(d["dose_projection"]),
                pk_model_type=d["pk_model_type"],
                pk_confidence=d["pk_confidence"],
                pk_uncertainties=d["pk_uncertainties"],
            )

        def _synthesis_step(d: dict) -> SynthesisStep:
            return SynthesisStep(**d)

        def _starting_material(d: dict) -> StartingMaterial:
            return StartingMaterial(**d)

        def _manufacturing(d: dict) -> ManufacturingAssessment:
            return ManufacturingAssessment(
                synthesis_steps_count=d["synthesis_steps_count"],
                longest_linear_sequence=d["longest_linear_sequence"],
                overall_yield_percent=d["overall_yield_percent"],
                complexity_score=d["complexity_score"],
                synthesis_route=[
                    _synthesis_step(s) for s in d["synthesis_route"]
                ],
                chiral_centers=d["chiral_centers"],
                resolution_required=d["resolution_required"],
                starting_materials=[
                    _starting_material(m) for m in d["starting_materials"]
                ],
                hazardous_reagents=d["hazardous_reagents"],
                cryogenic_steps=d["cryogenic_steps"],
                high_pressure_steps=d["high_pressure_steps"],
                chromatography_required=d["chromatography_required"],
                api_cost_per_kg_usd_low=d["api_cost_per_kg_usd_low"],
                api_cost_per_kg_usd_high=d["api_cost_per_kg_usd_high"],
                formulated_cost_per_dose_usd=d["formulated_cost_per_dose_usd"],
                cogs_confidence=d["cogs_confidence"],
                cogs_drivers=d["cogs_drivers"],
                recommended_dosage_form=d["recommended_dosage_form"],
                excipient_compatibility=d["excipient_compatibility"],
                stability_months_25C=d["stability_months_25C"],
                special_requirements=d["special_requirements"],
                cmc_risk=d["cmc_risk"],
                cmc_risks_identified=d["cmc_risks_identified"],
                cmc_mitigations=d["cmc_mitigations"],
            )

        def _designation(d: dict) -> DesignationEligibility:
            return DesignationEligibility(**d)

        def _precedent(d: dict) -> RegulatoryPrecedent:
            return RegulatoryPrecedent(**d)

        def _reg_timeline(d: dict) -> RegulatoryTimeline:
            return RegulatoryTimeline(**d)

        def _regulatory(d: dict) -> RegulatoryStrategy:
            return RegulatoryStrategy(
                designations=[
                    _designation(ds) for ds in d["designations"]
                ],
                prv_eligible=d["prv_eligible"],
                prv_qualifying_disease=d["prv_qualifying_disease"],
                prv_value_estimate_usd_millions=d[
                    "prv_value_estimate_usd_millions"
                ],
                recommended_pathway=d["recommended_pathway"],
                alternative_pathway=d["alternative_pathway"],
                animal_rule_applicable=d["animal_rule_applicable"],
                biomarker_strategy=d["biomarker_strategy"],
                primary_endpoint_phase3=d["primary_endpoint_phase3"],
                surrogate_endpoint_available=d["surrogate_endpoint_available"],
                precedents=[_precedent(p) for p in d["precedents"]],
                relevant_fda_guidances=d["relevant_fda_guidances"],
                timeline=[_reg_timeline(t) for t in d["timeline"]],
                regulatory_risk=d["regulatory_risk"],
                regulatory_hurdles=d["regulatory_hurdles"],
                regulatory_mitigations=d["regulatory_mitigations"],
            )

        def _approved_drug(d: dict) -> ApprovedDrug:
            return ApprovedDrug(**d)

        def _clinical_candidate(d: dict) -> ClinicalCandidate:
            return ClinicalCandidate(**d)

        def _competitive(d: dict) -> CompetitiveLandscape:
            return CompetitiveLandscape(
                disease_prevalence=d["disease_prevalence"],
                geographic_distribution=d["geographic_distribution"],
                current_standard_of_care=d["current_standard_of_care"],
                unmet_need_level=d["unmet_need_level"],
                approved_drugs=[
                    _approved_drug(a) for a in d["approved_drugs"]
                ],
                clinical_pipeline=[
                    _clinical_candidate(c) for c in d["clinical_pipeline"]
                ],
                recent_failures=d["recent_failures"],
                our_advantages=d["our_advantages"],
                our_disadvantages=d["our_disadvantages"],
                market_size_peak_sales_usd_millions=d[
                    "market_size_peak_sales_usd_millions"
                ],
                prv_value_usd_millions=d["prv_value_usd_millions"],
                stockpile_opportunity_usd_millions=d[
                    "stockpile_opportunity_usd_millions"
                ],
                total_opportunity_usd_millions=d[
                    "total_opportunity_usd_millions"
                ],
                active_acquirers=d["active_acquirers"],
                recent_transactions=d["recent_transactions"],
                valuation_benchmarks=d["valuation_benchmarks"],
            )

        def _blocking_patent(d: dict) -> BlockingPatent:
            return BlockingPatent(**d)

        def _patentability(d: dict) -> PatentabilityAssessment:
            return PatentabilityAssessment(**d)

        def _patent_analysis(d: dict) -> PatentAnalysis:
            return PatentAnalysis(
                fto_status=d["fto_status"],
                blocking_patents_count=d["blocking_patents_count"],
                blocking_patents=[
                    _blocking_patent(p) for p in d["blocking_patents"]
                ],
                design_around_required=d["design_around_required"],
                patent_landscape_families=d["patent_landscape_families"],
                key_assignees=d["key_assignees"],
                patentability=[
                    _patentability(p) for p in d["patentability"]
                ],
                recommended_filings=d["recommended_filings"],
                geographic_strategy=d["geographic_strategy"],
                orphan_exclusivity_years=d["orphan_exclusivity_years"],
                pediatric_exclusivity_possible=d[
                    "pediatric_exclusivity_possible"
                ],
                ip_risk=d["ip_risk"],
                ip_risks_identified=d["ip_risks_identified"],
                ip_mitigations=d["ip_mitigations"],
            )

        def _trial_design(d: dict) -> TrialDesign:
            return TrialDesign(**d)

        def _biomarker_strategy(d: dict) -> BiomarkerStrategy:
            return BiomarkerStrategy(**d)

        def _patient_population(d: dict) -> PatientPopulation:
            return PatientPopulation(**d)

        def _clinical_trial_design(d: dict) -> ClinicalTrialDesign:
            return ClinicalTrialDesign(
                phase1_design=_trial_design(d["phase1_design"]),
                phase2_design=_trial_design(d["phase2_design"]),
                biomarker_strategy=_biomarker_strategy(
                    d["biomarker_strategy"]
                ),
                patient_population=_patient_population(
                    d["patient_population"]
                ),
                outbreak_dependent=d["outbreak_dependent"],
                adaptive_design_elements=d["adaptive_design_elements"],
                compassionate_use_pathway=d["compassionate_use_pathway"],
                total_cost_to_phase2_usd_millions=d[
                    "total_cost_to_phase2_usd_millions"
                ],
                total_months_to_phase2_readout=d[
                    "total_months_to_phase2_readout"
                ],
            )

        # --- top-level assembly ------------------------------------------
        return cls(
            schema_version=data["schema_version"],
            package_id=data["package_id"],
            generation_timestamp=data["generation_timestamp"],
            tier=_enum(Tier, data["tier"]),
            executive_summary=_executive_summary(data["executive_summary"]),
            molecule_profile=_molecule_profile(data["molecule_profile"]),
            target_validation=_target_validation(data["target_validation"]),
            efficacy_data=_efficacy_data(data["efficacy_data"]),
            safety_profile=_safety_profile(data["safety_profile"]),
            pkpd_analysis=_pkpd_analysis(data["pkpd_analysis"]),
            manufacturing_assessment=_manufacturing(
                data["manufacturing_assessment"]
            ),
            regulatory_strategy=_regulatory(data["regulatory_strategy"]),
            competitive_landscape=_competitive(
                data["competitive_landscape"]
            ),
            patent_analysis=_patent_analysis(data["patent_analysis"]),
            clinical_trial_design=_clinical_trial_design(
                data["clinical_trial_design"]
            ),
            governance=_governance(data["governance"]),
            worldline_id=data["worldline_id"],
            slc_entries=data["slc_entries"],
        )


# ===========================================================================
# TIER-SPECIFIC OUTPUTS
# ===========================================================================

@dataclass
class GoldSummaryReport:
    """
    Summary report for GOLD-tier compounds.

    Contains passing/failing metrics and guidance on what gaps remain
    before the compound could reach DIAMOND classification.
    """

    compound_id: str
    disease_id: str
    tier: Tier
    generation_date: str
    passing_metrics: List[MetricScore]
    failing_metrics: List[MetricScore]
    gaps_to_diamond: List[str]
    optimization_suggestions: List[str]
    governance: GovernanceMetadata
    worldline_id: str

    def to_json(self) -> str:
        """Serialize the Gold summary report to JSON."""
        return json.dumps(asdict(self), indent=2, default=_json_default)


@dataclass
class WorldLineRecord:
    """
    Minimal record for SILVER/BRONZE tier compounds.

    Every compound attempt (pass or fail) generates a WorldLine record
    capturing the 35D physics state at the point of evaluation. This
    builds a complete failure/success map across the search space.
    """

    worldline_id: str
    compound_id: str
    smiles: str
    disease_id: str
    tier: Tier
    generation_timestamp: str
    physics_vector: List[float]
    global_sum: float
    energy_delta: float
    rejection_stage: str
    rejection_reasons: List[str]
    failing_metrics: Dict[str, float]
    trace_id: str
    slc_hash: str

    def to_json(self) -> str:
        """Serialize the WorldLine record to JSON."""
        return json.dumps(asdict(self), indent=2, default=_json_default)


# ===========================================================================
# Public API
# ===========================================================================

__all__ = [
    # Version
    "DOSSIER_SCHEMA_VERSION",
    # Enums
    "Tier",
    "RiskLevel",
    "Recommendation",
    # Section 00 - Executive Summary
    "MetricScore",
    "RiskAssessment",
    "GovernanceMetadata",
    "ExecutiveSummary",
    # Section 01 - Molecule Profile
    "StructuralAlert",
    "SimilarityMatch",
    "PhysicochemicalProperties",
    "DruglikenessAssessment",
    "MoleculeProfile",
    # Section 02 - Target Validation
    "TargetInteraction",
    "BindingAnalysis",
    "OffTarget",
    "TargetValidation",
    # Section 03 - Efficacy Data
    "ActivityPrediction",
    "EfficacyBenchmark",
    "EfficacyData",
    # Section 04 - Safety Profile
    "OrganToxicity",
    "CardiacSafety",
    "CYPInhibition",
    "GenotoxicityPrediction",
    "SafetyProfile",
    # Section 05 - PKPD Analysis
    "AbsorptionProfile",
    "DistributionProfile",
    "MetabolismProfile",
    "ExcretionProfile",
    "PKParameters",
    "DoseProjection",
    "PKPDAnalysis",
    # Section 06 - Manufacturing Assessment
    "SynthesisStep",
    "StartingMaterial",
    "ManufacturingAssessment",
    # Section 07 - Regulatory Strategy
    "DesignationEligibility",
    "RegulatoryPrecedent",
    "RegulatoryTimeline",
    "RegulatoryStrategy",
    # Section 08 - Competitive Landscape
    "ApprovedDrug",
    "ClinicalCandidate",
    "CompetitiveLandscape",
    # Section 09 - Patent Analysis
    "BlockingPatent",
    "PatentabilityAssessment",
    "PatentAnalysis",
    # Section 10 - Clinical Trial Design
    "TrialDesign",
    "BiomarkerStrategy",
    "PatientPopulation",
    "ClinicalTrialDesign",
    # Complete Packages
    "DossierPackage",
    "GoldSummaryReport",
    "WorldLineRecord",
]
