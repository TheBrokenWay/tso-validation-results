"""
Tests for dossier_schema.py — Pharma-grade dossier structure validation.

Covers:
- All dataclass instantiation and validation
- DossierPackage.validate() rules
- Serialization round-trip (to_json / from_dict)
- Tier routing logic
- Constitutional toxicity hard limit enforcement (Law L11)
- GoldSummaryReport and WorldLineRecord
"""

import json
import os
import sys
import unittest
from datetime import datetime, timezone

_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from PX_System.foundation.dossier_schema import (
    DOSSIER_SCHEMA_VERSION,
    Tier,
    RiskLevel,
    Recommendation,
    MetricScore,
    RiskAssessment,
    GovernanceMetadata,
    ExecutiveSummary,
    StructuralAlert,
    SimilarityMatch,
    PhysicochemicalProperties,
    DruglikenessAssessment,
    MoleculeProfile,
    TargetInteraction,
    BindingAnalysis,
    OffTarget,
    TargetValidation,
    ActivityPrediction,
    EfficacyBenchmark,
    EfficacyData,
    OrganToxicity,
    CardiacSafety,
    CYPInhibition,
    GenotoxicityPrediction,
    SafetyProfile,
    AbsorptionProfile,
    DistributionProfile,
    MetabolismProfile,
    ExcretionProfile,
    PKParameters,
    DoseProjection,
    PKPDAnalysis,
    SynthesisStep,
    StartingMaterial,
    ManufacturingAssessment,
    DesignationEligibility,
    RegulatoryPrecedent,
    RegulatoryTimeline,
    RegulatoryStrategy,
    ApprovedDrug,
    ClinicalCandidate,
    CompetitiveLandscape,
    BlockingPatent,
    PatentabilityAssessment,
    PatentAnalysis,
    TrialDesign,
    BiomarkerStrategy,
    PatientPopulation,
    ClinicalTrialDesign,
    DossierPackage,
    GoldSummaryReport,
    WorldLineRecord,
)


# ── Helpers: Build valid test fixtures ─────────────────────────────────

def _governance() -> GovernanceMetadata:
    return GovernanceMetadata(
        zeus_gate_status="PASSED",
        trace_id="PRD-V-TEST0001",
        worldline_id="WL-0001",
        slc_hash="abc123def456",
        constitutional_compliance=True,
        laws_verified=5,
    )


def _metric(name="Tox", val=0.005, target=0.02, passed=True) -> MetricScore:
    return MetricScore(
        name=name, value=val, target=target, unit="index",
        passed=passed, display_bar_percent=90,
    )


def _risk() -> RiskAssessment:
    return RiskAssessment(
        category="safety", level=RiskLevel.LOW,
        factors=["Low toxicity"], mitigations=["GLP tox study"],
    )


def _executive_summary() -> ExecutiveSummary:
    return ExecutiveSummary(
        compound_id="CPD-001", compound_name="TestMol",
        disease_name="Nipah virus infection", disease_id="nipah_virus_infection",
        generation_date=datetime.now(timezone.utc).isoformat(),
        tier=Tier.DIAMOND, investment_thesis="Promising candidate.",
        key_metrics=[_metric()], recommendation=Recommendation.ADVANCE,
        confidence_percent=75,
        cost_to_ind_usd_low=8_000_000, cost_to_ind_usd_high=15_000_000,
        months_to_ind=18, risks=[_risk()], governance=_governance(),
    )


def _physicochemical() -> PhysicochemicalProperties:
    return PhysicochemicalProperties(
        logp=2.5, logd_7_4=2.0, pka_acidic=None, pka_basic=None,
        aqueous_solubility_ug_ml=50.0, psa=80.0, hbd=2, hba=5,
        rotatable_bonds=4, aromatic_rings=2, heavy_atoms=25,
        fraction_sp3=0.4, molar_refractivity=None,
    )


def _druglikeness() -> DruglikenessAssessment:
    return DruglikenessAssessment(
        lipinski_violations=0, lipinski_pass=True, veber_pass=True,
        ghose_pass=True, lead_likeness_pass=True,
        qed_score=0.75, cns_mpo_score=None,
    )


def _molecule_profile() -> MoleculeProfile:
    return MoleculeProfile(
        smiles="CC(=O)Oc1ccccc1C(=O)O", canonical_smiles="CC(=O)Oc1ccccc1C(=O)O",
        inchi="InChI=1S/TEST", inchi_key="TESTINCHIKEY",
        molecular_formula="C9H8O4", molecular_weight=180.16,
        exact_mass=180.042, stereochemistry="achiral",
        structure_image_base64=None, physicochemical=_physicochemical(),
        druglikeness=_druglikeness(), structural_alerts=[], nearest_known_drugs=[],
    )


def _binding_analysis() -> BindingAnalysis:
    return BindingAnalysis(
        binding_site_name="ATP pocket", binding_mode="competitive",
        docking_score=-8.5, docking_software="PX_Dock",
        key_interactions=[], predicted_kd_nM=50.0,
    )


def _target_validation() -> TargetValidation:
    return TargetValidation(
        target_name="NiV-G glycoprotein", gene_symbol="NiV-G",
        uniprot_id="Q9IH62", target_class="viral_surface",
        organism="Nipah henipavirus", druggability_score=0.72,
        expression_profile={"lung": "high"}, genetic_associations=["NiV entry"],
        known_inhibitors=["m102.4"], clinical_precedent=False,
        binding_analysis=_binding_analysis(),
        selectivity_profile=[], validation_confidence="HIGH", rationale="Essential for viral entry.",
    )


def _activity_prediction() -> ActivityPrediction:
    return ActivityPrediction(
        activity_type="IC50", value_nM=100.0,
        confidence_low_nM=50.0, confidence_high_nM=200.0,
        model_name="PX_QSAR", model_version="3.0",
        training_set_size=5000, applicability_domain="INSIDE",
        similarity_to_training=0.85,
    )


def _efficacy_data() -> EfficacyData:
    return EfficacyData(
        primary_activity=_activity_prediction(),
        secondary_activities=[], benchmarks=[],
        percentile_rank=85, dose_projection_mg=100.0,
        dose_projection_basis="IC50 x 10", therapeutic_index=50.0,
        resistance_liability="LOW", resistance_mutations_known=[],
        efficacy_confidence="HIGH", key_uncertainties=[],
        recommended_validation_experiments=["in vitro NiV pseudovirus assay"],
    )


def _cardiac_safety() -> CardiacSafety:
    return CardiacSafety(
        herg_ic50_uM=30.0, herg_threshold_uM=10.0, herg_passed=True,
        herg_liability_class="LOW", qt_prolongation_risk="LOW",
        cardiovascular_safety_margin=300.0,
    )


def _safety_profile(tox_index=0.005) -> SafetyProfile:
    return SafetyProfile(
        overall_toxicity_index=tox_index, tier_classification="TOXICITY_DIAMOND",
        safety_margin_fold=50.0, organ_toxicity=[],
        cardiac_safety=_cardiac_safety(), cyp_inhibition=[],
        cyp_induction_risk="LOW", genotoxicity=[],
        max_tolerated_dose_mg_kg=100.0, dose_limiting_toxicity=None,
        safety_confidence="HIGH", recommended_safety_studies=["Ames test"],
        monitoring_parameters=["LFT"],
    )


def _absorption() -> AbsorptionProfile:
    return AbsorptionProfile(
        oral_bioavailability_percent=65.0, bioavailability_confidence=0.7,
        caco2_permeability_nm_s=25.0, pgp_substrate=False,
        pgp_inhibitor=False, food_effect="MINIMAL",
    )


def _distribution() -> DistributionProfile:
    return DistributionProfile(
        vd_L_kg=1.2, plasma_protein_binding_percent=85.0,
        blood_brain_barrier=False, bbb_permeability=None,
        tissue_distribution={"liver": "high"},
    )


def _metabolism() -> MetabolismProfile:
    return MetabolismProfile(
        primary_cyp="CYP3A4", secondary_cyps=["CYP2D6"],
        hlm_half_life_min=45.0, rlm_half_life_min=30.0,
        metabolic_stability="MODERATE", active_metabolites=False,
        metabolic_soft_spots=["aromatic ring"],
    )


def _excretion() -> ExcretionProfile:
    return ExcretionProfile(
        clearance_mL_min_kg=5.0, primary_route="hepatic",
        renal_excretion_percent=20.0, fecal_excretion_percent=60.0,
    )


def _pk_parameters() -> PKParameters:
    return PKParameters(
        cmax_ng_mL=500.0, cmax_confidence_percent=20.0,
        tmax_hours=2.0, tmax_confidence_hours=0.5,
        auc_0_24_ng_h_mL=4000.0, auc_confidence_percent=25.0,
        half_life_hours=8.0, half_life_confidence_percent=15.0,
    )


def _dose_projection() -> DoseProjection:
    return DoseProjection(
        target_concentration_ng_mL=200.0, target_basis="IC50 coverage",
        predicted_dose_mg=100.0, dosing_frequency="BID",
        dose_range_phase1_low_mg=10.0, dose_range_phase1_high_mg=200.0,
        max_recommended_starting_dose_mg=50.0, mrsd_basis="NOAEL/10",
    )


def _pkpd_analysis() -> PKPDAnalysis:
    return PKPDAnalysis(
        absorption=_absorption(), distribution=_distribution(),
        metabolism=_metabolism(), excretion=_excretion(),
        pk_parameters_human_predicted=_pk_parameters(),
        dose_projection=_dose_projection(),
        pk_model_type="1-compartment", pk_confidence="MEDIUM",
        pk_uncertainties=["Limited allometric data"],
    )


def _manufacturing() -> ManufacturingAssessment:
    return ManufacturingAssessment(
        synthesis_steps_count=5, longest_linear_sequence=4,
        overall_yield_percent=25.0, complexity_score="MODERATE",
        synthesis_route=[], chiral_centers=1, resolution_required=False,
        starting_materials=[], hazardous_reagents=[],
        cryogenic_steps=False, high_pressure_steps=False,
        chromatography_required=True,
        api_cost_per_kg_usd_low=5000, api_cost_per_kg_usd_high=15000,
        formulated_cost_per_dose_usd=2.50, cogs_confidence="LOW",
        cogs_drivers=["Chromatography"], recommended_dosage_form="Oral tablet",
        excipient_compatibility="GOOD", stability_months_25C=24,
        special_requirements=[], cmc_risk="MEDIUM",
        cmc_risks_identified=["Chiral resolution"], cmc_mitigations=["Asymmetric synthesis"],
    )


def _regulatory() -> RegulatoryStrategy:
    return RegulatoryStrategy(
        designations=[], prv_eligible=True,
        prv_qualifying_disease="Nipah virus infection",
        prv_value_estimate_usd_millions=100.0,
        recommended_pathway="Accelerated Approval + Animal Rule",
        alternative_pathway="Standard NDA", animal_rule_applicable=True,
        biomarker_strategy="Viral load reduction",
        primary_endpoint_phase3="Survival at Day 28",
        surrogate_endpoint_available=True, precedents=[],
        relevant_fda_guidances=["Animal Rule Guidance"],
        timeline=[], regulatory_risk="MEDIUM",
        regulatory_hurdles=["Outbreak-dependent trial"],
        regulatory_mitigations=["Ring vaccination study design"],
    )


def _competitive() -> CompetitiveLandscape:
    return CompetitiveLandscape(
        disease_prevalence="Rare, epidemic",
        geographic_distribution=["Southeast Asia", "South Asia"],
        current_standard_of_care="Supportive care only",
        unmet_need_level="CRITICAL", approved_drugs=[], clinical_pipeline=[],
        recent_failures=[], our_advantages=["First-in-class"],
        our_disadvantages=["Outbreak-dependent trials"],
        market_size_peak_sales_usd_millions=50.0,
        prv_value_usd_millions=100.0, stockpile_opportunity_usd_millions=200.0,
        total_opportunity_usd_millions=350.0,
        active_acquirers=["Gilead", "Emergent BioSolutions"],
        recent_transactions=[], valuation_benchmarks=[],
    )


def _patent_analysis() -> PatentAnalysis:
    return PatentAnalysis(
        fto_status="CLEAR", blocking_patents_count=0, blocking_patents=[],
        design_around_required=False, patent_landscape_families=12,
        key_assignees=[], patentability=[], recommended_filings=["Composition of matter"],
        geographic_strategy=["US", "EU", "India"],
        orphan_exclusivity_years=7, pediatric_exclusivity_possible=True,
        ip_risk="LOW", ip_risks_identified=[], ip_mitigations=[],
    )


def _trial_design(phase="Phase 1") -> TrialDesign:
    return TrialDesign(
        phase=phase, study_type="open-label dose escalation",
        population="Healthy volunteers", n_subjects=30,
        design="3+3 dose escalation", dose_range_mg="10-200",
        primary_objectives=["Safety", "PK"],
        primary_endpoint="DLT rate", secondary_endpoints=["PK parameters"],
        duration_weeks=12, estimated_cost_usd_millions=3.0,
    )


def _clinical_trial_design() -> ClinicalTrialDesign:
    return ClinicalTrialDesign(
        phase1_design=_trial_design("Phase 1"),
        phase2_design=_trial_design("Phase 2"),
        biomarker_strategy=BiomarkerStrategy(
            pharmacodynamic_biomarkers=["Viral load"],
            efficacy_biomarkers=["Seroconversion"],
            safety_biomarkers=["ALT", "AST"],
            companion_diagnostic_required=False,
        ),
        patient_population=PatientPopulation(
            inclusion_criteria_key=["Confirmed NiV exposure"],
            exclusion_criteria_key=["Pregnancy"],
            recruitment_considerations="Outbreak-dependent",
            geographic_strategy=["Malaysia", "Bangladesh", "India"],
            enrollment_rate_per_month=5,
        ),
        outbreak_dependent=True,
        adaptive_design_elements=["Sample size re-estimation"],
        compassionate_use_pathway="Emergency IND",
        total_cost_to_phase2_usd_millions=12.0,
        total_months_to_phase2_readout=36,
    )


def _full_dossier_package(**overrides) -> DossierPackage:
    """Build a complete valid DossierPackage for testing."""
    defaults = dict(
        schema_version=DOSSIER_SCHEMA_VERSION,
        package_id="PKG-TEST-001",
        generation_timestamp=datetime.now(timezone.utc).isoformat(),
        tier=Tier.DIAMOND,
        executive_summary=_executive_summary(),
        molecule_profile=_molecule_profile(),
        target_validation=_target_validation(),
        efficacy_data=_efficacy_data(),
        safety_profile=_safety_profile(),
        pkpd_analysis=_pkpd_analysis(),
        manufacturing_assessment=_manufacturing(),
        regulatory_strategy=_regulatory(),
        competitive_landscape=_competitive(),
        patent_analysis=_patent_analysis(),
        clinical_trial_design=_clinical_trial_design(),
        governance=_governance(),
        worldline_id="WL-0001",
        slc_entries=["SLC-001", "SLC-002"],
    )
    defaults.update(overrides)
    return DossierPackage(**defaults)


# ── Test Classes ───────────────────────────────────────────────────────

class TestSchemaVersion(unittest.TestCase):
    def test_schema_version_is_string(self):
        self.assertIsInstance(DOSSIER_SCHEMA_VERSION, str)

    def test_schema_version_format(self):
        parts = DOSSIER_SCHEMA_VERSION.split(".")
        self.assertEqual(len(parts), 3)
        for p in parts:
            self.assertTrue(p.isdigit())


class TestEnums(unittest.TestCase):
    def test_tier_values(self):
        self.assertEqual(Tier.DIAMOND.value, "DIAMOND")
        self.assertEqual(Tier.GOLD.value, "GOLD")
        self.assertEqual(Tier.SILVER.value, "SILVER")
        self.assertEqual(Tier.BRONZE.value, "BRONZE")

    def test_risk_levels(self):
        self.assertEqual(len(RiskLevel), 3)

    def test_recommendation_values(self):
        self.assertEqual(Recommendation.ADVANCE.value, "ADVANCE_TO_IND")
        self.assertEqual(Recommendation.REJECT.value, "REJECT")


class TestDataclassInstantiation(unittest.TestCase):
    """Verify all dataclasses can be instantiated with valid data."""

    def test_metric_score(self):
        m = _metric()
        self.assertEqual(m.name, "Tox")
        self.assertTrue(m.passed)

    def test_governance_metadata(self):
        g = _governance()
        self.assertTrue(g.constitutional_compliance)

    def test_molecule_profile(self):
        mp = _molecule_profile()
        self.assertAlmostEqual(mp.molecular_weight, 180.16)

    def test_target_validation(self):
        tv = _target_validation()
        self.assertGreater(tv.druggability_score, 0)

    def test_efficacy_data(self):
        ed = _efficacy_data()
        self.assertGreater(ed.therapeutic_index, 0)

    def test_safety_profile(self):
        sp = _safety_profile()
        self.assertLess(sp.overall_toxicity_index, 0.0210)

    def test_pkpd_analysis(self):
        pk = _pkpd_analysis()
        self.assertGreater(pk.pk_parameters_human_predicted.half_life_hours, 0)

    def test_manufacturing_assessment(self):
        ma = _manufacturing()
        self.assertGreater(ma.synthesis_steps_count, 0)

    def test_regulatory_strategy(self):
        rs = _regulatory()
        self.assertTrue(rs.prv_eligible)

    def test_competitive_landscape(self):
        cl = _competitive()
        self.assertEqual(cl.unmet_need_level, "CRITICAL")

    def test_patent_analysis(self):
        pa = _patent_analysis()
        self.assertEqual(pa.fto_status, "CLEAR")

    def test_clinical_trial_design(self):
        ct = _clinical_trial_design()
        self.assertTrue(ct.outbreak_dependent)


class TestDossierPackageValidation(unittest.TestCase):
    """Test DossierPackage.validate() enforcement."""

    def test_valid_package_passes(self):
        pkg = _full_dossier_package()
        valid, errors = pkg.validate()
        self.assertTrue(valid, f"Validation failed: {errors}")
        self.assertEqual(len(errors), 0)

    def test_wrong_schema_version_fails(self):
        pkg = _full_dossier_package(schema_version="0.0.0")
        valid, errors = pkg.validate()
        self.assertFalse(valid)
        self.assertTrue(any("Schema version" in e for e in errors))

    def test_empty_package_id_fails(self):
        pkg = _full_dossier_package(package_id="")
        valid, errors = pkg.validate()
        self.assertFalse(valid)
        self.assertTrue(any("package_id" in e for e in errors))

    def test_non_diamond_tier_fails(self):
        pkg = _full_dossier_package(tier=Tier.GOLD)
        valid, errors = pkg.validate()
        self.assertFalse(valid)
        self.assertTrue(any("DIAMOND" in e for e in errors))

    def test_toxicity_hard_limit_enforcement(self):
        """Law L11: tox >= 0.0210 must fail validation."""
        pkg = _full_dossier_package(safety_profile=_safety_profile(tox_index=0.0210))
        valid, errors = pkg.validate()
        self.assertFalse(valid)
        self.assertTrue(any("0.0210" in e for e in errors))

    def test_toxicity_just_below_limit_passes(self):
        """Law L11: tox = 0.0209 should not trigger hard limit."""
        pkg = _full_dossier_package(safety_profile=_safety_profile(tox_index=0.0209))
        valid, errors = pkg.validate()
        # Should not contain toxicity hard-limit error
        tox_errors = [e for e in errors if "0.0210" in e]
        self.assertEqual(len(tox_errors), 0)

    def test_empty_smiles_fails(self):
        mp = _molecule_profile()
        mp.smiles = ""
        pkg = _full_dossier_package(molecule_profile=mp)
        valid, errors = pkg.validate()
        self.assertFalse(valid)
        self.assertTrue(any("smiles" in e for e in errors))

    def test_zero_molecular_weight_fails(self):
        mp = _molecule_profile()
        mp.molecular_weight = 0
        pkg = _full_dossier_package(molecule_profile=mp)
        valid, errors = pkg.validate()
        self.assertFalse(valid)
        self.assertTrue(any("molecular_weight" in e for e in errors))

    def test_empty_worldline_id_fails(self):
        pkg = _full_dossier_package(worldline_id="")
        valid, errors = pkg.validate()
        self.assertFalse(valid)
        self.assertTrue(any("worldline_id" in e for e in errors))

    def test_governance_compliance_false_fails(self):
        gov = _governance()
        gov.constitutional_compliance = False
        pkg = _full_dossier_package(governance=gov)
        valid, errors = pkg.validate()
        self.assertFalse(valid)
        self.assertTrue(any("constitutional_compliance" in e for e in errors))

    def test_negative_safety_margin_fails(self):
        sp = _safety_profile()
        sp.safety_margin_fold = -1.0
        pkg = _full_dossier_package(safety_profile=sp)
        valid, errors = pkg.validate()
        self.assertFalse(valid)
        self.assertTrue(any("safety_margin" in e for e in errors))

    def test_zero_trial_subjects_fails(self):
        ct = _clinical_trial_design()
        ct.phase1_design.n_subjects = 0
        pkg = _full_dossier_package(clinical_trial_design=ct)
        valid, errors = pkg.validate()
        self.assertFalse(valid)
        self.assertTrue(any("n_subjects" in e for e in errors))


class TestSerialization(unittest.TestCase):
    """Test JSON round-trip serialization."""

    def test_to_json_produces_valid_json(self):
        pkg = _full_dossier_package()
        j = pkg.to_json()
        parsed = json.loads(j)
        self.assertIsInstance(parsed, dict)
        self.assertEqual(parsed["schema_version"], DOSSIER_SCHEMA_VERSION)

    def test_to_json_has_all_sections(self):
        pkg = _full_dossier_package()
        parsed = json.loads(pkg.to_json())
        for section in [
            "executive_summary", "molecule_profile", "target_validation",
            "efficacy_data", "safety_profile", "pkpd_analysis",
            "manufacturing_assessment", "regulatory_strategy",
            "competitive_landscape", "patent_analysis", "clinical_trial_design",
        ]:
            self.assertIn(section, parsed, f"Missing section: {section}")

    def test_from_dict_round_trip(self):
        pkg = _full_dossier_package()
        j = pkg.to_json()
        parsed = json.loads(j)
        reconstructed = DossierPackage.from_dict(parsed)
        self.assertEqual(reconstructed.package_id, pkg.package_id)
        self.assertEqual(reconstructed.tier, Tier.DIAMOND)
        self.assertEqual(
            reconstructed.molecule_profile.molecular_weight,
            pkg.molecule_profile.molecular_weight,
        )
        self.assertEqual(
            reconstructed.safety_profile.overall_toxicity_index,
            pkg.safety_profile.overall_toxicity_index,
        )

    def test_from_dict_preserves_nested_structures(self):
        pkg = _full_dossier_package()
        parsed = json.loads(pkg.to_json())
        reconstructed = DossierPackage.from_dict(parsed)
        self.assertEqual(
            reconstructed.pkpd_analysis.absorption.oral_bioavailability_percent,
            65.0,
        )
        self.assertEqual(
            reconstructed.target_validation.binding_analysis.docking_score,
            -8.5,
        )

    def test_reconstructed_validates(self):
        pkg = _full_dossier_package()
        parsed = json.loads(pkg.to_json())
        reconstructed = DossierPackage.from_dict(parsed)
        valid, errors = reconstructed.validate()
        self.assertTrue(valid, f"Reconstructed package validation failed: {errors}")


class TestGoldSummaryReport(unittest.TestCase):
    def test_instantiation(self):
        report = GoldSummaryReport(
            compound_id="CPD-002", disease_id="ebola_virus_disease",
            tier=Tier.GOLD,
            generation_date=datetime.now(timezone.utc).isoformat(),
            passing_metrics=[_metric()],
            failing_metrics=[_metric("Selectivity", 5.0, 10.0, False)],
            gaps_to_diamond=["Improve selectivity > 10x"],
            optimization_suggestions=["Add fluorine at C-3"],
            governance=_governance(), worldline_id="WL-0002",
        )
        self.assertEqual(report.tier, Tier.GOLD)

    def test_to_json(self):
        report = GoldSummaryReport(
            compound_id="CPD-002", disease_id="ebola",
            tier=Tier.GOLD,
            generation_date=datetime.now(timezone.utc).isoformat(),
            passing_metrics=[], failing_metrics=[],
            gaps_to_diamond=[], optimization_suggestions=[],
            governance=_governance(), worldline_id="WL-0002",
        )
        j = report.to_json()
        parsed = json.loads(j)
        self.assertEqual(parsed["tier"], "GOLD")


class TestWorldLineRecord(unittest.TestCase):
    def test_instantiation(self):
        record = WorldLineRecord(
            worldline_id="WL-9999", compound_id="CPD-999",
            smiles="C1=CC=CC=C1", disease_id="malaria",
            tier=Tier.BRONZE,
            generation_timestamp=datetime.now(timezone.utc).isoformat(),
            physics_vector=[0.1] * 35, global_sum=1.0, energy_delta=0.05,
            rejection_stage="OPE", rejection_reasons=["MW > 500"],
            failing_metrics={"molecular_weight": 550.0},
            trace_id="PRD-V-TEST9999", slc_hash="xyz789",
        )
        self.assertEqual(len(record.physics_vector), 35)
        self.assertEqual(record.tier, Tier.BRONZE)

    def test_to_json(self):
        record = WorldLineRecord(
            worldline_id="WL-0003", compound_id="CPD-003",
            smiles="CCO", disease_id="dengue", tier=Tier.SILVER,
            generation_timestamp=datetime.now(timezone.utc).isoformat(),
            physics_vector=[0.0] * 35, global_sum=0.0, energy_delta=0.0,
            rejection_stage="ADMET", rejection_reasons=["Tox borderline"],
            failing_metrics={"toxicity": 0.0205}, trace_id="PRD-V-T0003",
            slc_hash="abc",
        )
        j = record.to_json()
        parsed = json.loads(j)
        self.assertEqual(parsed["tier"], "SILVER")
        self.assertEqual(len(parsed["physics_vector"]), 35)


class TestConstitutionalCompliance(unittest.TestCase):
    """Verify constitutional rules are enforced in the schema."""

    def test_toxicity_boundary_at_0210_fails(self):
        """Law L11: exactly at boundary is FAILURE."""
        pkg = _full_dossier_package(safety_profile=_safety_profile(tox_index=0.0210))
        valid, errors = pkg.validate()
        self.assertFalse(valid)

    def test_toxicity_above_0210_fails(self):
        """Law L11: above boundary is FAILURE."""
        pkg = _full_dossier_package(safety_profile=_safety_profile(tox_index=0.03))
        valid, errors = pkg.validate()
        self.assertFalse(valid)

    def test_toxicity_diamond_tier(self):
        """tox < 0.01 is DIAMOND tier."""
        pkg = _full_dossier_package(safety_profile=_safety_profile(tox_index=0.005))
        valid, errors = pkg.validate()
        self.assertTrue(valid, f"Should pass: {errors}")

    def test_no_mock_data_in_schema(self):
        """Schema module should not contain mock/dummy/placeholder strings."""
        import inspect
        source = inspect.getsource(sys.modules["PX_System.foundation.dossier_schema"])
        for forbidden in ["mock", "dummy", "placeholder", "TODO", "FIXME"]:
            self.assertNotIn(
                forbidden.lower(), source.lower(),
                f"Found '{forbidden}' in dossier_schema.py source"
            )


if __name__ == "__main__":
    unittest.main()
