# ============================================================
#  PREDATOR X — FINALIZATION SPECIFICATION (FULL + COMPLETE)
#  Drop-in governance contract for PX_Warehouse
# ============================================================
#
# GRADING METRICS (required for every finalized dossier):
#   - soc_benchmarking_score   — Standard-of-Care comparison vs best-in-class for disease space
#   - novelty_fingerprint_score — IP readiness; uniqueness vs known SMILES (patent eligibility)
#   - synthetic_accessibility_score — Economic scalability; lab-feasibility (e.g. SA proxy)
#
# This ensures:
#   - Future developers understand the contract
#   - The spec matches the implementation
#   - The system remains self-documenting
#
# ============================================================

# Version tag: recorded on every finalized dossier for re-finalization and governance evolution
FINALIZATION_VERSION = "3.0.0"

FINALIZATION_SPEC = {

    # --------------------------------------------------------
    # 0. Spec version (recorded on each finalized dossier)
    # --------------------------------------------------------
    "finalization_version": FINALIZATION_VERSION,

    # --------------------------------------------------------
    # 1. Definition of a Fully Finalized Dossier
    # --------------------------------------------------------
    "finalized_requires": [
        # Promotion
        "promotion_stage",              # NOV → DISC → REP
        # Mechanistic
        "moa_hypothesis",
        "disease_context",
        "prv_category",
        "prv_rationale",
        "prv_eligibility_scoring",
        # Lineage + ancestry
        "lineage_assignment",
        "ancestral_trace.full_chain",
        "ancestral_trace.first_seen",
        "ancestral_trace.worldline_ids",
        # Physics + engines
        "physics_map_id",
        "physics_map_hash",
        "engines",
        # Trial simulation binding
        "trial_simulation_id",
        "trial_simulation_hash",
        "trial_outcome_summary",
        # Worldline + grading
        "worldline_generated",
        "worldline_path",
        "discovery_grading",
        # Grading metrics (required): Competitive Intelligence, IP Readiness, Economic Scalability
        # See GRADING METRICS block at top of file. Contract: all three must be present.
        "soc_benchmarking_score",
        "novelty_fingerprint_score",
        "synthetic_accessibility_score",
        # Authorization chain (v3.0.0: 12-engine sign-off blocks)
        "authorization_chain",
        "authorization_chain.all_engines_authorized",
        "authorization_chain.authorization_count",
        # Constitutional
        "constitutional_risk_overrides",
        "zeus_verdict",
        "constitutional_seal",
        # Tier + sealing
        "tier",
        "rep_seal",
        # ALCOA + traceability
        "alcoa_metadata",
        "causal_trace_log",
        "fda_compliance",
        # Governance audit
        "finalization_version",
        "spec_compliance_report"
    ],

    # --------------------------------------------------------
    # 1b. Grading metrics contract (self-documenting)
    # --------------------------------------------------------
    "grading_metrics": {
        "description": "Required for every finalized dossier. Ensures future developers understand the contract, spec matches implementation, and system remains self-documenting.",
        "required": [
            "soc_benchmarking_score",
            "novelty_fingerprint_score",
            "synthetic_accessibility_score"
        ],
        "soc_benchmarking_score": "Standard-of-Care comparison vs best-in-class for disease space",
        "novelty_fingerprint_score": "IP readiness; uniqueness vs known SMILES (patent eligibility basis)",
        "synthetic_accessibility_score": "Economic scalability; lab-feasibility (e.g. SA_Score proxy)"
    },

    # --------------------------------------------------------
    # 2. Zeus Constitutional Gate (Full Form)
    # --------------------------------------------------------
    "zeus_gate": {
        "laws_required": [
            "L1_HARM_LAW",
            "U27_STABILITY_LAW",
            "U34_GLOBAL_SUM",
            "L11_DETERMINISTIC_ENGINE"
        ],
        "toxicity_threshold": 0.0210,
        "harm_energy_threshold": 0.0210,
        "verdict_field": "zeus_verdict"
    },

    # --------------------------------------------------------
    # 3. Folder Structure Requirements
    # --------------------------------------------------------
    "warehouse_structure": {
        "roots": [
            "Novel_Dossiers",
            "Prv_Dossiers",
            "WorldLines",
            "Finalized_Dossiers",
            "TrialSimulations",
            "Operations",
            "Calibration_Molecules"
        ],
        "tiers": [
            "Bronze",
            "Silver",
            "Diamond",
            "Gold"
        ]
    },

    # --------------------------------------------------------
    # 4. Automatic Promotion Pipeline
    # --------------------------------------------------------
    "auto_promotion": {
        "enabled": True,
        "watch_folders": [
            "Novel_Dossiers",
            "Prv_Dossiers"
        ],
        "pipeline_stages": [
            "NOV",
            "DISC",
            "REP",
            "FINALIZED"
        ],
        "actions_per_stage": {
            "NOV": [
                "validate_alcoa",
                "validate_engines",
                "assign_item_id"
            ],
            "DISC": [
                "generate_moa",
                "anchor_disease_space",
                "compute_prv_score"
            ],
            "REP": [
                "assign_lineage",
                "extend_ancestral_trace",
                "bind_physics_map",
                "bind_trial_simulation",
                "grade_discovery",
                "apply_constitutional_overrides"
            ],
            "FINALIZED": [
                "run_zeus_gate",
                "compute_tier",
                "seal_rep",
                "write_finalized_dossier",
                "write_worldline",
                "log_operations"
            ]
        }
    },

    # --------------------------------------------------------
    # 5. Warehouse Cleanliness Invariant
    # --------------------------------------------------------
    "cleanliness_invariant": {
        "rule": "NO_NEW_PRODUCTION_IF_UNFINALIZED_EXIST",
        "description": (
            "New dossiers may NOT be produced if any dossier exists in "
            "Novel_Dossiers or Prv_Dossiers that does not have a corresponding "
            "entry in Finalized_Dossiers."
        )
    },

    # --------------------------------------------------------
    # 6. Backfill / Bootstrap Requirements
    # --------------------------------------------------------
    "bootstrap": {
        "script": "PX_Executive/px_finalize.py",
        "behavior": [
            "scan_all_existing_dossiers",
            "skip_if_already_finalized",
            "run_full_finalization",
            "write_finalized_output",
            "assert_all_finalized_before_new_production"
        ]
    },

    # --------------------------------------------------------
    # 7. Finalization Checklist (Hard Gate)
    # --------------------------------------------------------
    "finalization_checklist": {
        "must_pass": [
            "all_required_fields_present",
            "tier_computed",
            "rep_seal_generated",
            "worldline_written",
            "ancestral_trace_complete",
            "physics_map_bound",
            "trial_simulation_bound",
            "grading_metrics_complete"
        ],
        "failure_behavior": "DO_NOT_WRITE_FINALIZED_DOSSIER",
        "note_zeus_gate": "zeus_gate_passed is NOT in must_pass; only the WRITE to Finalized_Dossiers is gated by Zeus. Full dossier is always produced; caller gates write on zeus_verdict.authorized."
    },

    # --------------------------------------------------------
    # 8. Authorization Chain (v3.0.0)
    # --------------------------------------------------------
    "authorization_chain_spec": {
        "description": "Every finalized dossier must contain an authorization_chain proving all 12 mandatory engines executed and signed off.",
        "required_engines": [
            "OPE_V3_DETERMINISTIC",
            "OBE_V3_DETERMINISTIC",
            "OCE_V3_DETERMINISTIC",
            "OLE_V3_DETERMINISTIC",
            "OME_V3_DETERMINISTIC",
            "OSE_V3_DETERMINISTIC",
            "OPE_ADMET_V3_DETERMINISTIC",
            "PKPD_EMAX_V2.0",
            "DOSE_OPTIMIZER_V2.1",
            "TRIAL_ENGINE_V1",
            "GRADING_ENGINE_V1",
        ],
        "sign_off_fields": [
            "engine_id", "version", "timestamp", "execution_time_ms",
            "status", "authorized", "laws_checked", "laws_passed",
            "laws_failed", "reason", "input_hash", "output_hash", "signature",
        ],
        "chain_summary_fields": [
            "authorization_chain",
            "all_engines_authorized",
            "authorization_count",
        ],
    },

    # --------------------------------------------------------
    # 9. Disease Constraint Binding (v3.0.0)
    # --------------------------------------------------------
    "disease_constraint_spec": {
        "description": "Each dossier targeting a PRV-eligible disease should reference a validated disease constraint file from PX_Domain/PRV_Diseases/.",
        "constraint_dir": "PX_Domain/PRV_Diseases/",
        "manifest": "PX_Domain/PRV_Diseases/manifest.json",
        "required_fields": [
            "disease_name", "disease_id", "pathogen_type", "target_classes",
            "ic50_max_um", "toxicity_threshold", "prv_eligible",
            "constitutional_compliance",
        ],
    },

    # --------------------------------------------------------
    # 10. External Data Sources (v3.0.0)
    # --------------------------------------------------------
    "external_data_sources": {
        "clinicaltrials": {
            "client": "PX_Data/clinicaltrials/client.py",
            "refresh_interval_days": 7,
            "data_file": "PX_Data/clinicaltrials/clinicaltrials_interventions.csv",
        },
        "fda": {
            "client": "PX_Data/fda/client.py",
            "refresh_interval_days": 7,
            "data_files": [
                "PX_Data/fda/fda_prv_approvals.jsonl",
                "PX_Data/fda/fda_prv_guidance.jsonl",
            ],
        },
        "patents": {
            "client": "PX_Data/patents/client.py",
            "refresh_interval_days": 7,
            "data_file": "PX_Data/patents/patent_index.db",
        },
    },
}

# ============================================================
#  RECOMMENDATION: FULL PIPELINE EXECUTION FOR ALL DOSSIERS
#  (No single-metric rejection; full multi-factor evaluation)
# ============================================================

FULL_PIPELINE_RECOMMENDATION = {
    "principle": "NO_SINGLE_METRIC_REJECTION",

    "description": (
        "The finalization pipeline must NEVER reject a dossier based on a "
        "single failing metric or early-stage condition. All dossiers—whether "
        "they ultimately pass or fail—must be processed through the ENTIRE "
        "finalization pipeline. This includes: mechanistic inference, disease "
        "anchoring, PRV scoring, lineage extension, ancestral trace, physics "
        "map binding, trial simulation binding, discovery grading, SoC "
        "benchmarking, novelty fingerprint scoring, synthetic accessibility "
        "scoring, compliance reporting, and worldline generation."
    ),

    "rationale": (
        "Predator X is a governed, transparent system. No black-box behavior "
        "is permitted. Every dossier must produce a complete record of ALL "
        "evaluations, regardless of outcome. A dossier that fails Zeus or any "
        "other gate must still contain a full set of metrics and a complete "
        "spec_compliance_report so that failures are explainable, auditable, "
        "and reproducible."
    ),

    "required_behavior": [
        "Do NOT return early from run_finalization() on Zeus failure.",
        "Compute ALL metrics before deciding pass/fail.",
        "Generate a full spec_compliance_report for EVERY dossier.",
        "Record ALL failing checks, not just the first failure.",
        "Apply this logic uniformly to Novel_Dossiers and Prv_Dossiers.",
        "Only the WRITE to Finalized_Dossiers is gated by Zeus.",
        "Failed dossiers must still produce worldlines and compliance reports.",
        "Failed dossiers must still be fully traceable and auditable."
    ],

    "implementation_note": (
        "The Zeus gate runs at the END of the pipeline (IMPLEMENTATION_ROADMAP). "
        "All metrics are computed first; then Zeus is evaluated. Only the caller "
        "gates the WRITE on zeus_verdict.authorized. Dossiers rejected by Zeus "
        "still contain full grading, SoC benchmarking, novelty scoring, "
        "synthetic accessibility scoring, physics map binding, trial binding, "
        "lineage, ancestry, and compliance reporting."
    )
}

# ============================================================
#  IMPLEMENTATION ROADMAP — FULL PIPELINE EXECUTION MODEL
#  (No early exits, no single-metric rejection, full auditability)
# ============================================================

IMPLEMENTATION_ROADMAP = {
    "principle": "FULL_PIPELINE_EVALUATION_FOR_ALL_DOSSIERS",

    "description": (
        "Finalization_Pipeline.py must be restructured so that EVERY dossier "
        "runs through ALL stages of the pipeline, regardless of toxicity, "
        "risk flags, or early failures. No single metric may reject a dossier. "
        "Only the final WRITE to Finalized_Dossiers is gated by Zeus. "
        "All dossiers—pass or fail—must produce a complete, auditable profile."
    ),

    "inference_stage": {
        "required_behavior": [
            "Compute MoA for all dossiers.",
            "Compute disease anchoring for all dossiers.",
            "Compute PRV scoring for all dossiers.",
            "No early exits allowed."
        ]
    },

    "simulation_stage": {
        "required_behavior": [
            "Bind physics maps for all dossiers.",
            "Bind trial outcomes for all dossiers (or mark as virtual).",
            "Simulation must run regardless of toxicity or risk level."
        ]
    },

    "grading_stage": {
        "required_behavior": [
            "Generate SoC Benchmarking score for all dossiers.",
            "Generate Novelty Fingerprint score for all dossiers.",
            "Generate Synthetic Accessibility score for all dossiers.",
            "Generate discovery grading metrics for all dossiers.",
            "No dossier may skip grading due to early failures."
        ]
    },

    "governance_stage": {
        "required_behavior": [
            "Run the Zeus Gate at the END of the pipeline.",
            "Zeus determines folder placement: Finalized_Dossiers vs Operations.",
            "Zeus must NOT stop the pipeline early.",
            "All dossiers must produce a full spec_compliance_report.",
            "All dossiers must include finalization_version."
        ]
    },

    "cleanliness_invariant": {
        "rule": "NO_RAW_DOSSIERS_LEFT_BEHIND",
        "description": (
            "Every dossier in Novel_Dossiers and Prv_Dossiers must be processed "
            "through the full pipeline. No dossier may remain unprocessed or "
            "partially processed. This ensures a complete, auditable record for "
            "the entire inventory."
        )
    },

    "rejection_logic": {
        "required_behavior": [
            "Rejection must be multi-factor, not single-factor.",
            "All failures must be recorded in spec_compliance_report.",
            "Failed dossiers must still generate worldlines.",
            "Failed dossiers must still generate grading metrics."
        ]
    },
}

# --------------------------------------------------------
#  DOSSIER COMPLIANCE CHECKLIST (must have before compliant)
#  Existing finalized dossiers must be reprocessed to add these.
# --------------------------------------------------------
DOSSIER_COMPLIANCE_CHECKLIST = [
    "A. Add the three grading metrics: soc_benchmarking_score, novelty_fingerprint_score, synthetic_accessibility_score",
    "B. Add the version tag: finalization_version: \"1.0.0\"",
    "C. Add the compliance report: spec_compliance_report: {...}",
    "D. Fix discovery grading: resolve dict vs float comparison; populate discovery_grading.metrics with real values",
    "E. Re-run the full finalization pipeline with Zeus at the end so all metrics are computed regardless of pass/fail",
]

# ============================================================
# END OF SPEC
# ============================================================

__all__ = [
    "FINALIZATION_SPEC",
    "FINALIZATION_VERSION",
    "FULL_PIPELINE_RECOMMENDATION",
    "IMPLEMENTATION_ROADMAP",
    "DOSSIER_COMPLIANCE_CHECKLIST",
]
