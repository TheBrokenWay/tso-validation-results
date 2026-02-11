# Warehouse Rejection Analysis — Processed Data That Did Not Pass

## Summary

Analysis of all **processed** (Prv_Dossiers + Novel_Dossiers) dossiers that **did not** reach Finalized_Dossiers shows a single dominant cause: **Zeus gate rejection on toxicity/harm threshold 0.0210**.

- **Counts (current warehouse):** 131 Prv + 221 Novel = 352 total in tiered folders; 21 finalized → **331 unfinalized**.
- **Sample run (80 unfinalized):** 80/80 Zeus rejected; 0 passed; 0 early errors.
- **Failed laws:** L11_DETERMINISTIC_ENGINE (toxicity_index ≥ 0.021), L1_HARM_LAW (harm_energy ≥ 0.021), and when virtual trial ran: TRIAL_TOXICITY, TRIAL_HARM_ENERGY.

## Commonality of Failures

| Cause | Explanation |
|-------|-------------|
| **L11 / L1** | Every rejected dossier has `toxicity_index` and/or `harm_energy` ≥ **0.0210**. Zeus hard-fails single-metric above this. |
| **Toxicity distribution (rejected)** | In the sampled 80: min 0.027, max 0.38, **median 0.23**. Most fall in [0.10, 0.50). |
| **Whole-profile override** | Few if any have TOXICITY_DIAMOND or (tox &gt; 0.02 and safety_margin &gt; 50), so the whole-profile escape rarely applies. |

So the **commonality is**: our **Zeus bar is 2.1%** (0.021), while the **bulk of discovery-stage candidates** have predicted toxicity in the **5–40%** range. They are correctly processed (OPE/ADMET, grading, trial, finalization) but **intentionally not written** to Finalized_Dossiers because they fail the constitutional threshold.

## Are Our Standards Too Strict?

- **For a regulatory / filing bar:** 0.021 is **appropriate** (conservative, safety-first).
- **For discovery-stage funnel:** 0.021 is **much stricter** than:
  - Our **own grading schema**: GOLD toxicity_max = **0.20** (20%), SILVER = 0.35, BRONZE = 0.50.
  - **Typical industry** preclinical filters: often 10–30% predicted risk is considered acceptable for progression; 2.1% is late-stage / high-safety.

So we are **not** mis-calibrated for “finalized = filing-ready” — we are **very strict** for “what gets into Finalized_Dossiers.” That can **impede the process** if the intent is to also use Finalized_Dossiers as the main funnel for discovery (everything else effectively “lost” from the pipeline).

## Alignment With Industry

- **Industry:** Preclinical toxicity filters often allow 10–30% (or higher) for progression; 2.1% is at the tight end.
- **Our design:** Finalized_Dossiers = Zeus-authorized = 2.1% bar. Discovery grading (GradingSchema_Discovery) allows up to 20–50% by tier.
- **Gap:** Candidates that are **Silver/Bronze by discovery grading** (tox &lt; 0.35 or &lt; 0.5) are **all Zeus-rejected** if tox ≥ 0.021, so they never reach Finalized_Dossiers.

## Recalibration Options (Accurate Results, Don’t Impede Process)

1. **Keep 0.021 for “finalized = filing-ready” (no change to governance lock)**  
   - Do **not** loosen 0.021 for any regulatory-facing or filing path.  
   - Add an **audit path** for Zeus-rejected dossiers (e.g. `Operations/Zeus_Rejected` or `Operations/Discovery_Rejected`) so they are not “lost”: same pipeline output, not written to Finalized_Dossiers, but stored for review and potential future recalibration.

2. **Two-tier bar (recommended)**  
   - **Finalized_Dossiers:** unchanged — only Zeus-authorized at **0.021** (current behavior).  
   - **Discovery_Accepted (or similar):** allow a **discovery bar** (e.g. tox &lt; 0.05 or &lt; 0.10) for a **separate** output path (e.g. `Finalized_Dossiers/Discovery_Accepted` or `Operations/Discovery_Accepted`).  
   - Ensures: (a) we don’t impede discovery funnel, (b) we keep a strict, accurate bar for “finalized” and regulatory use, (c) recalibration is explicit and traceable.

3. **Align Zeus with grading for discovery-only path**  
   - For a **discovery-only** pipeline (no filing): allow Zeus “pass” when dossier tier is Silver or better (e.g. tox &lt; 0.10) and use a **relaxed** single-metric threshold (e.g. 0.05) only for that path.  
   - Keep **0.021** for any path that feeds regulatory/filing.

4. **Operational**  
   - Run `PX_Warehouse/analyze_warehouse_rejections.py` periodically (e.g. weekly or after big batches) and keep rejection rate + toxicity/harm distributions in view so recalibration stays **data-driven**.

## Conclusion

- **Why processed data didn’t pass:** Almost entirely **Zeus L11/L1** (toxicity_index / harm_energy ≥ 0.021). No indication of bugs or wrong standards for a “filing-ready” bar.  
- **Standards:** 0.021 is **strict** vs industry and vs our own discovery grading; that is by design for Finalized_Dossiers.  
- **Impedance:** Yes — discovery-stage candidates with tox in 5–40% are all excluded from Finalized_Dossiers.  
- **Recalibration:** Keep 0.021 for finalized/filing; add a **Discovery_Accepted** (or Zeus_Rejected audit) path and/or a **discovery bar** (e.g. 0.05 or 0.10) for a separate output so we get **accurate results** (strict bar where it matters) without **impeding** the discovery process.
