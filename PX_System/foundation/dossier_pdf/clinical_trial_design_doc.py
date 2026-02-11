"""Clinical Trial Design document generator (Section 10)."""
from __future__ import annotations

from .base_generator import BaseDocumentGenerator


class ClinicalTrialDesignDocument(BaseDocumentGenerator):
    """Generates the 10_CLINICAL_TRIAL_DESIGN report."""

    section_name = "clinical_trial_design"
    section_number = "10"

    def generate(self) -> str:
        data = self.sections_data.get("clinical_trial_design", {})
        phase1 = data.get("phase1", {})
        phase2 = data.get("phase2", {})
        biomarkers = data.get("biomarkers", [])
        cost_timeline = data.get("cost_timeline", {})
        endpoints = data.get("endpoints", {})
        population = data.get("population", {})

        out = self._header("CLINICAL TRIAL DESIGN")

        out += self._section_header("PHASE 1 TRIAL DESIGN")
        out += self._kv("Design", phase1.get("design", "N/A")) + "\n"
        out += self._kv("Dose Escalation", phase1.get("dose_escalation", "N/A")) + "\n"
        out += self._kv("Starting Dose", phase1.get("starting_dose", "N/A")) + "\n"
        out += self._kv("Max Planned Dose", phase1.get("max_dose", "N/A")) + "\n"
        out += self._kv("Cohort Size", phase1.get("cohort_size", "N/A")) + "\n"
        out += self._kv("Number of Cohorts", phase1.get("num_cohorts", "N/A")) + "\n"
        out += self._kv("Total Patients", phase1.get("total_patients", "N/A")) + "\n"
        out += self._kv("Duration", phase1.get("duration", "N/A")) + "\n"
        out += self._kv("Primary Objective", phase1.get("primary_objective", "N/A")) + "\n"
        p1_criteria = phase1.get("inclusion_criteria", [])
        if p1_criteria:
            out += "\n  Key Inclusion Criteria:\n"
            for c in p1_criteria:
                out += f"    - {c}\n"
        p1_exclusion = phase1.get("exclusion_criteria", [])
        if p1_exclusion:
            out += "\n  Key Exclusion Criteria:\n"
            for c in p1_exclusion:
                out += f"    - {c}\n"

        out += self._section_header("PHASE 2 TRIAL DESIGN")
        out += self._kv("Design", phase2.get("design", "N/A")) + "\n"
        out += self._kv("Randomization", phase2.get("randomization", "N/A")) + "\n"
        out += self._kv("Blinding", phase2.get("blinding", "N/A")) + "\n"
        out += self._kv("Arms", phase2.get("arms", "N/A")) + "\n"
        out += self._kv("Patients per Arm", phase2.get("patients_per_arm", "N/A")) + "\n"
        out += self._kv("Total Patients", phase2.get("total_patients", "N/A")) + "\n"
        out += self._kv("Duration", phase2.get("duration", "N/A")) + "\n"
        out += self._kv("Primary Objective", phase2.get("primary_objective", "N/A")) + "\n"
        out += self._kv("Statistical Power", phase2.get("statistical_power", "N/A")) + "\n"
        out += self._kv("Comparator", phase2.get("comparator", "N/A")) + "\n"

        out += self._section_header("STUDY ENDPOINTS")
        primary_ep = endpoints.get("primary", [])
        if primary_ep:
            out += "  Primary Endpoints:\n"
            for ep in primary_ep:
                out += f"    - {ep.get('name', 'N/A')}\n"
                out += f"      Timeframe: {ep.get('timeframe', 'N/A')}\n"
                out += f"      Measure: {ep.get('measure', 'N/A')}\n"
        secondary_ep = endpoints.get("secondary", [])
        if secondary_ep:
            out += "\n  Secondary Endpoints:\n"
            for ep in secondary_ep:
                out += f"    - {ep.get('name', 'N/A')}\n"
                out += f"      Timeframe: {ep.get('timeframe', 'N/A')}\n"
        exploratory_ep = endpoints.get("exploratory", [])
        if exploratory_ep:
            out += "\n  Exploratory Endpoints:\n"
            for ep in exploratory_ep:
                out += f"    - {ep.get('name', 'N/A')}\n"

        out += self._section_header("TARGET POPULATION")
        out += self._kv("Indication", population.get("indication", "N/A")) + "\n"
        out += self._kv("Age Range", population.get("age_range", "N/A")) + "\n"
        out += self._kv("Geographic Regions", population.get("regions", "N/A")) + "\n"
        out += self._kv("Enrollment Strategy", population.get("enrollment_strategy", "N/A")) + "\n"
        out += self._kv("Diversity Requirements", population.get("diversity", "N/A")) + "\n"

        out += self._section_header("BIOMARKER STRATEGY")
        if biomarkers:
            rows = [[b.get("name", "N/A"), b.get("type", "N/A"),
                      b.get("purpose", "N/A"), b.get("validation_status", "N/A")]
                     for b in biomarkers]
            out += self._table(["Biomarker", "Type", "Purpose", "Validation"],
                               rows, [18, 14, 22, 16])
        else:
            out += "  No biomarkers defined.\n"

        out += self._section_header("COST & TIMELINE ESTIMATES")
        out += self._kv("Phase 1 Cost", f"${cost_timeline.get('phase1_cost_usd', 0):,}") + "\n"
        out += self._kv("Phase 1 Duration", cost_timeline.get("phase1_duration", "N/A")) + "\n"
        out += self._kv("Phase 2 Cost", f"${cost_timeline.get('phase2_cost_usd', 0):,}") + "\n"
        out += self._kv("Phase 2 Duration", cost_timeline.get("phase2_duration", "N/A")) + "\n"
        out += self._kv("Total Clinical Cost", f"${cost_timeline.get('total_cost_usd', 0):,}") + "\n"
        out += self._kv("Total Duration", cost_timeline.get("total_duration", "N/A")) + "\n"
        out += self._kv("CRO Strategy", cost_timeline.get("cro_strategy", "N/A")) + "\n"

        out += self._footer()
        return out
