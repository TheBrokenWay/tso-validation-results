# Trial-Integrated Dossier Spec (Final Additions)

Exact spec of what must be added to every finalized dossier to make it a **complete in-silico trial dossier**. The Finalization_Pipeline **runs the full virtual trial** (TrialEngine + DoseOptimizer_v2 + VirtualEfficacyAnalytics) when a dossier has OPE/ADMET but no trial_result/virtual_efficacy, so every finalized dossier is upgraded to the full in-silico trial standard.

---

## 1. Virtual Trial Engine Execution

Dossiers must include (when trial has been run):

- **Cohort simulation:** cohort_size, demographic distribution, variability model, PK/PD variability curves
- **Dose–response simulation:** dose_levels, predicted responder_rate per dose, predicted toxicity per dose, PTA curve
- **Arm-to-arm comparison:** candidate vs SoC, effect size, superiority / non-inferiority metrics
- **Trial outcome summary (actual values):** responder_rate, pta, toxicity_rate, simulated efficacy distribution, uncertainty intervals

---

## 2. Trial Binding Fields (real values when trial run)

| Field | Requirement |
|-------|-------------|
| **trial_simulation_id** | Real ID for each simulation run (not `NOT_RUN`) |
| **trial_simulation_hash** | Hash of the simulation output |
| **trial_outcome_summary** | Full object with: `status: "COMPLETED"`, `responder_rate`, `pta`, `toxicity_rate`, `effect_size`, `uncertainty_bounds`, `cohort_statistics` |

When trial was not run: keep placeholders but use the same structure (status `VIRTUAL_TRIAL_NOT_RUN`, optional empty/null fields).

---

## 3. Grading Engine Integration With Trial Results

- **Trial-aware grading metrics:** trial_responder_rate, trial_pta, trial_variability_cv, trial_effect_size
- **Updated grade:** GOLD / SILVER / BRONZE based on trial outcomes, or NEEDS_REVIEW with trial-based reasoning

---

## 4. Trial-Integrated Tier Assignment

- **Tier recalculation using trial outcomes:** GOLD if trial efficacy + safety meet thresholds; SILVER if moderate; BRONZE if minimal; NEEDS_REVIEW if ambiguous
- When trial_outcome_summary.status is COMPLETED, tier may be derived from trial metrics in addition to discovery physics

---

## 5. Trial-Integrated Worldline

- **worldline_trial_event:** timestamp, simulation_id, trial_hash, outcome_summary (added to worldline or finalization when trial COMPLETED)

---

## 6. Trial-Integrated PRV Eligibility

- **trial_adjusted_prv_score:** efficacy_weight, safety_weight, novelty_weight, SoC comparison weight

---

## 7. Trial-Integrated SoC Benchmarking

- **effect_size_vs_soc:** positive / neutral / negative, magnitude, confidence interval

---

## 8. Trial-Integrated Synthetic Accessibility (optional)

- **trial_adjusted_sa_score:** when trial suggests high dose requirements, SA may be adjusted

---

## 9. Trial-Integrated Constitutional Review (Zeus)

- **Zeus review of trial outcomes:** harm_energy under simulated dosing, toxicity under simulated dosing, stability under variability (when trial data present)

---

## 10. Finalization Block Must Include Trial Results

- **trial_results section:** full cohort summary, dose–response table, PTA curve, toxicity curve, effect size, uncertainty intervals

---

## Clean Checklist (for IDE)

```
ADD TO EVERY FINALIZED DOSSIER:

1. Virtual Trial Engine Execution:
   - cohort simulation, dose–response simulation, PTA curve,
   - responder_rate distribution, toxicity distribution,
   - arm-to-arm comparison vs SoC, uncertainty intervals

2. Populate trial binding fields:
   - trial_simulation_id (real), trial_simulation_hash (real),
   - trial_outcome_summary: responder_rate, pta, toxicity_rate, effect_size,
     uncertainty_bounds, cohort_statistics, status: "COMPLETED" (or NOT_RUN)

3. Integrate trial results into grading:
   - trial_responder_rate, trial_pta, trial_variability_cv, trial_effect_size, updated grade

4. Recompute tier using trial outcomes when status COMPLETED.

5. Add worldline_trial_event: simulation_id, hash, outcome_summary

6. Add trial_adjusted_prv_score (efficacy/safety/novelty/SoC weights).

7. Add effect_size_vs_soc with confidence interval.

8. Add trial_adjusted_sa_score (optional).

9. Add Zeus review of trial outcomes when present.

10. Add trial_results block: cohort summary, dose–response table, PTA curve,
    toxicity curve, effect size, uncertainty intervals.
```
