"""
GradingEngine.py
WeightedGradingEngine v2.0 — 100-Point Multi-Parameter Optimization (MPO) Grading

Replaces the binary 5-criterion all-or-nothing system with a weighted scoring
system across 6 categories:
  - Efficacy (30 pts): IC50/EC50, Emax
  - Safety (25 pts): toxicity_index, safety_margin, hERG risk
  - PK (20 pts): bioavailability, half-life, clearance, volume of distribution
  - Drug-likeness (10 pts): Lipinski violations, QED
  - Manufacturability (10 pts): synthetic steps, yield
  - FTO (5 pts): freedom-to-operate status

Tier thresholds:
  DIAMOND_TIER  >= 85%
  GOLD_TIER     70-84%
  SILVER_TIER   50-69%
  BRONZE_TIER   < 50%

Hard limits (immediate REJECTED — cannot be overridden by scoring):
  - toxicity_index >= 0.0210 (Law L11)
  - herg_risk > 0.5
  - zeus_gate_passed == False (when present)
  - constitutional_violations > 0

Updated: 2026-02-11
"""

import json
import sys
import time
from pathlib import Path
from typing import Dict, Any, Tuple, List, Optional
from datetime import datetime, UTC
from PX_System.foundation.sign_off import create_sign_off


# ── Constants ──

TOXICITY_HARD_LIMIT = 0.0210  # Law L11 — immutable
HERG_RISK_LIMIT = 0.5

TIER_THRESHOLDS = {
    "DIAMOND_TIER": 85.0,
    "GOLD_TIER": 70.0,
    "SILVER_TIER": 50.0,
    # Below 50 → BRONZE_TIER
}

MAX_SCORE = 100
CATEGORY_WEIGHTS = {
    "efficacy": 30,
    "safety": 25,
    "pk": 20,
    "druglikeness": 10,
    "manufacturability": 10,
    "fto": 5,
}


class WeightedGradingEngine:
    """100-point weighted MPO grading engine for pharmaceutical dossiers."""

    TOXICITY_LIMIT = TOXICITY_HARD_LIMIT  # backward compat attribute

    def __init__(self, verbose: bool = True):
        self.verbose = verbose
        self.version = "v2.0-weighted"
        self.grading_history: List[Dict[str, Any]] = []

    # ── Main Entry Point ──

    def grade_dossier(self, dossier_path) -> Dict[str, Any]:
        """
        Grade a complete dossier using weighted MPO scoring.

        Args:
            dossier_path: Path to dossier JSON file, or dossier dict in memory.

        Returns:
            Grading result with grade, percentage, category_scores, reasoning, sign_off.
        """
        _t0 = time.monotonic()

        if isinstance(dossier_path, dict):
            dossier = dossier_path
            dossier_ref = "<in-memory>"
        else:
            with open(dossier_path, "r", encoding="utf-8") as f:
                dossier = json.load(f)
            dossier_ref = str(dossier_path)

        # 1. Hard limit checks (immediate REJECTED)
        hard_limits = self._check_hard_limits(dossier)
        if hard_limits["rejected"]:
            grade = "REJECTED"
            percentage = 0.0
            total_score = 0
            category_scores = {}
            reasoning = self._build_reasoning(category_scores, grade, hard_limits)
            metrics = self._extract_raw_inputs(dossier)
        else:
            # 2. Extract raw inputs
            inputs = self._extract_raw_inputs(dossier)

            # 3. Score each category
            category_scores = {
                "efficacy": self._score_efficacy(inputs),
                "safety": self._score_safety(inputs),
                "pk": self._score_pk(inputs),
                "druglikeness": self._score_druglikeness(inputs, dossier),
                "manufacturability": self._score_manufacturability(inputs),
                "fto": self._score_fto(inputs),
            }

            # 4. Sum scores
            total_score = sum(cs["score"] for cs in category_scores.values())
            percentage = round((total_score / MAX_SCORE) * 100, 1)

            # 5. Assign tier
            grade = self._assign_tier(percentage)
            metrics = inputs

            # 6. Build reasoning
            reasoning = self._build_reasoning(category_scores, grade, hard_limits)

        _elapsed_ms = int((time.monotonic() - _t0) * 1000)

        laws_passed = []
        laws_failed = []
        if not hard_limits.get("toxicity_failed"):
            laws_passed.append("L11")
        else:
            laws_failed.append("L11")
        grade_ok = grade != "REJECTED"
        if grade_ok:
            laws_passed.append("GRADE_RULES")
        else:
            laws_failed.append("GRADE_RULES")

        sign_off = create_sign_off(
            engine_id="GRADING_ENGINE_V2",
            version="v2.0-weighted",
            inputs={},
            outputs={"grade": grade, "percentage": percentage},
            laws_checked=["L11", "GRADE_RULES"],
            laws_results={"L11": not hard_limits.get("toxicity_failed", False),
                          "GRADE_RULES": grade_ok},
            execution_time_ms=_elapsed_ms,
        )

        result = {
            "grading_version": "WeightedGradingEngine_v2.0",
            "grade": grade,
            "percentage": percentage,
            "total_score": total_score,
            "max_score": MAX_SCORE,
            "hard_limit_checks": hard_limits,
            "category_scores": category_scores,
            "reasoning": reasoning,
            "metrics": metrics,
            "sign_off": sign_off,
            "timestamp": datetime.now(UTC).isoformat(),
            "grading_engine_version": self.version,
        }

        self.grading_history.append({
            "dossier": dossier_ref,
            "grade": grade,
            "timestamp": result["timestamp"],
        })

        return result

    # ── Hard Limits (Safety Shield) ──

    def _check_hard_limits(self, dossier: Dict[str, Any]) -> Dict[str, Any]:
        """Check hard limits that cause immediate REJECTED. Returns check results."""
        result = {
            "rejected": False,
            "toxicity_failed": False,
            "herg_failed": False,
            "zeus_failed": False,
            "constitutional_failed": False,
            "rejection_reasons": [],
        }

        # Extract toxicity
        tox = self._resolve_toxicity(dossier)
        if tox is not None and tox >= TOXICITY_HARD_LIMIT:
            result["rejected"] = True
            result["toxicity_failed"] = True
            result["rejection_reasons"].append(
                f"Law L11: toxicity {tox:.6f} >= {TOXICITY_HARD_LIMIT}"
            )

        # Extract hERG risk
        herg = self._resolve_herg(dossier)
        if herg is not None and herg > HERG_RISK_LIMIT:
            result["rejected"] = True
            result["herg_failed"] = True
            result["rejection_reasons"].append(
                f"hERG risk {herg:.4f} > {HERG_RISK_LIMIT}"
            )

        # Zeus gate (only if present — zeus runs after grading in px_prv.py)
        zeus = dossier.get("zeus_verdict")
        if isinstance(zeus, dict) and "authorized" in zeus:
            if not zeus["authorized"]:
                result["rejected"] = True
                result["zeus_failed"] = True
                result["rejection_reasons"].append(
                    f"Zeus gate not authorized: {zeus.get('rationale', 'unknown')}"
                )

        # Constitutional violations
        cv = dossier.get("constitutional_violations")
        if isinstance(cv, list) and len(cv) > 0:
            result["rejected"] = True
            result["constitutional_failed"] = True
            result["rejection_reasons"].append(
                f"{len(cv)} constitutional violation(s)"
            )

        return result

    # ── Raw Input Extraction ──

    def _extract_raw_inputs(self, dossier: Dict[str, Any]) -> Dict[str, Any]:
        """Extract all scoring-relevant fields from the dossier."""
        inputs: Dict[str, Any] = {}
        engines = dossier.get("engines") or {}

        # OPE fields
        ope = dossier.get("ope") or engines.get("ope") or \
              (dossier.get("inputs") or {}).get("ope_analysis") or {}

        inputs["ec50"] = self._safe_float(ope.get("ec50"))
        inputs["binding_affinity_nM"] = self._safe_float(ope.get("binding_affinity_nM"))
        inputs["emax"] = self._safe_float(ope.get("emax"))
        # OPE returns emax as 0-1 fraction; scoring expects 0-100 percentage
        if inputs["emax"] is not None and inputs["emax"] <= 1.0:
            inputs["emax"] = inputs["emax"] * 100.0
        inputs["molecular_weight"] = self._safe_float(ope.get("molecular_weight"))
        inputs["logp"] = self._safe_float(ope.get("logp"))
        inputs["hbd"] = self._safe_float(ope.get("hbd"))
        inputs["hba"] = self._safe_float(ope.get("hba"))

        # Binding affinity score (normalized 0-1)
        ba = ope.get("binding_affinity")
        if isinstance(ba, dict) and "score" in ba:
            inputs["binding_affinity_score"] = self._safe_float(ba["score"])
        elif isinstance(ba, (int, float)):
            inputs["binding_affinity_score"] = float(ba) if ba <= 1.0 else ba / 100.0
        else:
            inputs["binding_affinity_score"] = None

        # ADMET
        admet = dossier.get("admet") or engines.get("admet") or \
                (dossier.get("inputs") or {}).get("admet_analysis") or {}

        # Toxicity
        inputs["toxicity_index"] = self._resolve_toxicity(dossier)

        # Safety margin — from ADMET toxicity sub-dict or dose_optimization
        tox_block = admet.get("toxicity") if isinstance(admet.get("toxicity"), dict) else {}
        inputs["safety_margin"] = self._safe_float(tox_block.get("safety_margin"))
        # ADMET safety_margin formula produces tiny fractions (EC50*10/Cmax);
        # selectivity_index is the actual fold-change safety window
        if inputs["safety_margin"] is not None and inputs["safety_margin"] < 1.0:
            si = self._safe_float(tox_block.get("selectivity_index"))
            if si is not None and si >= 1.0:
                inputs["safety_margin"] = si
        if inputs["safety_margin"] is None:
            dose_opt = dossier.get("dose_optimization") or engines.get("dose_optimization") or {}
            best = dose_opt.get("best_regimen") or {}
            inputs["safety_margin"] = self._safe_float(best.get("safety_margin"))

        # hERG
        inputs["herg_risk"] = self._resolve_herg(dossier)

        # Absorption
        absorption = admet.get("absorption") if isinstance(admet.get("absorption"), dict) else {}
        oral_bio = self._safe_float(absorption.get("oral_bioavailability_percent"))
        if oral_bio is not None:
            inputs["bioavailability"] = oral_bio / 100.0
        else:
            inputs["bioavailability"] = None

        # Metabolism
        metabolism = admet.get("metabolism") if isinstance(admet.get("metabolism"), dict) else {}
        inputs["half_life_h"] = self._safe_float(metabolism.get("half_life_h"))
        if inputs["half_life_h"] is None:
            ome = engines.get("ome") or {}
            inputs["half_life_h"] = self._safe_float(ome.get("half_life_h"))

        clearance_L_h = self._safe_float(metabolism.get("clearance_L_per_h"))
        if clearance_L_h is not None:
            # Convert L/h to mL/min/kg: × 1000/60/70
            inputs["clearance_mL_min_kg"] = (clearance_L_h * 1000.0) / (60.0 * 70.0)
        else:
            inputs["clearance_mL_min_kg"] = None

        # Distribution
        distribution = admet.get("distribution") if isinstance(admet.get("distribution"), dict) else {}
        vd = self._safe_float(distribution.get("predicted_vd_L_per_kg"))
        if vd is None:
            vd_L = self._safe_float(distribution.get("vd_L"))
            if vd_L is not None:
                vd = vd_L / 70.0
        inputs["vd_L_per_kg"] = vd

        # QED — from OPE (computed via RDKit) or top-level dossier
        qed_val = self._safe_float(ope.get("qed"))
        if qed_val is None:
            qed_val = self._safe_float(dossier.get("qed"))
        inputs["qed"] = qed_val

        # Lipinski violations — compute from OPE
        inputs["lipinski_violations"] = self._count_lipinski_violations(inputs)

        # Synthesis steps — not computed in current pipeline
        inputs["synthesis_steps"] = self._safe_float(dossier.get("synthesis_steps"))

        # FTO
        ole = engines.get("ole") or {}
        fto = ole.get("freedom_to_operate")
        if isinstance(fto, bool):
            inputs["fto_status"] = "clear" if fto else "blocked"
        elif isinstance(fto, str):
            inputs["fto_status"] = fto.lower()
        else:
            inputs["fto_status"] = None

        # Constraint violations
        cv = dossier.get("constraint_violations")
        inputs["constraint_violations"] = cv if isinstance(cv, list) else []

        return inputs

    # ── Category Scorers ──

    def _score_efficacy(self, inputs: Dict[str, Any]) -> Dict[str, Any]:
        """Score efficacy (30 points max): IC50/EC50 + Emax."""
        max_pts = CATEGORY_WEIGHTS["efficacy"]
        score = 0.0

        # IC50/EC50 sub-score (up to 20 pts out of 30)
        ec50 = inputs.get("ec50")
        ba_nM = inputs.get("binding_affinity_nM")
        potency_nM = ec50 if ec50 is not None else ba_nM

        ic50_score = 0.0
        if potency_nM is not None:
            if potency_nM <= 10:
                ic50_score = 20.0    # Exceptional
            elif potency_nM <= 50:
                ic50_score = 17.0    # Excellent
            elif potency_nM <= 100:
                ic50_score = 14.0    # Good
            elif potency_nM <= 500:
                ic50_score = 10.0    # Moderate
            elif potency_nM <= 1000:
                ic50_score = 6.0     # Weak
            else:
                ic50_score = 2.0     # Very weak
        else:
            # Missing data: partial credit
            ic50_score = 10.0  # 50% of max sub-score

        # Emax sub-score (up to 10 pts out of 30)
        emax = inputs.get("emax")
        emax_score = 0.0
        if emax is not None:
            if emax >= 90:
                emax_score = 10.0
            elif emax >= 70:
                emax_score = 8.0
            elif emax >= 50:
                emax_score = 6.0
            elif emax >= 30:
                emax_score = 4.0
            else:
                emax_score = 2.0
        else:
            emax_score = 5.0  # 50% of max sub-score

        score = ic50_score + emax_score

        return {
            "score": round(score, 1),
            "max": max_pts,
            "percentage": round((score / max_pts) * 100, 1),
            "breakdown": {"ic50": round(ic50_score, 1), "emax": round(emax_score, 1)},
        }

    def _score_safety(self, inputs: Dict[str, Any]) -> Dict[str, Any]:
        """Score safety (25 points max): toxicity_index + safety_margin + hERG."""
        max_pts = CATEGORY_WEIGHTS["safety"]
        score = 0.0

        # Toxicity sub-score (up to 15 pts)
        tox = inputs.get("toxicity_index")
        tox_score = 0.0
        if tox is not None:
            if tox < 0.005:
                tox_score = 15.0     # Exceptional safety
            elif tox < 0.010:
                tox_score = 13.0     # Excellent
            elif tox < 0.015:
                tox_score = 10.0     # Good
            elif tox < 0.020:
                tox_score = 7.0      # Acceptable
            elif tox < TOXICITY_HARD_LIMIT:
                tox_score = 4.0      # Marginal (but below hard limit)
            else:
                tox_score = 0.0      # Hard limit violation (should be caught upstream)
        else:
            tox_score = 5.0  # Missing: pessimistic partial credit

        # Safety margin sub-score (up to 7 pts)
        sm = inputs.get("safety_margin")
        sm_score = 0.0
        if sm is not None:
            if sm >= 100:
                sm_score = 7.0
            elif sm >= 50:
                sm_score = 5.0
            elif sm >= 20:
                sm_score = 3.0
            elif sm >= 10:
                sm_score = 2.0
            else:
                sm_score = 1.0
        else:
            sm_score = 3.0  # Missing: partial credit

        # hERG sub-score (up to 3 pts)
        herg = inputs.get("herg_risk")
        herg_score = 0.0
        if herg is not None:
            if herg <= 0.1:
                herg_score = 3.0
            elif herg <= 0.3:
                herg_score = 2.0
            elif herg <= HERG_RISK_LIMIT:
                herg_score = 1.0
            else:
                herg_score = 0.0  # Hard limit (caught upstream)
        else:
            herg_score = 2.0  # Missing: assume moderate (pass if absent)

        score = tox_score + sm_score + herg_score

        return {
            "score": round(score, 1),
            "max": max_pts,
            "percentage": round((score / max_pts) * 100, 1),
            "breakdown": {
                "toxicity": round(tox_score, 1),
                "safety_margin": round(sm_score, 1),
                "herg": round(herg_score, 1),
            },
        }

    def _score_pk(self, inputs: Dict[str, Any]) -> Dict[str, Any]:
        """Score PK (20 points max): bioavailability + half_life + clearance + Vd."""
        max_pts = CATEGORY_WEIGHTS["pk"]
        score = 0.0

        # Bioavailability sub-score (up to 8 pts)
        bio = inputs.get("bioavailability")
        bio_score = 0.0
        if bio is not None:
            if bio >= 0.80:
                bio_score = 8.0
            elif bio >= 0.60:
                bio_score = 6.0
            elif bio >= 0.40:
                bio_score = 4.0
            elif bio >= 0.20:
                bio_score = 2.0
            else:
                bio_score = 1.0
        else:
            bio_score = 4.0  # Missing: partial credit

        # Half-life sub-score (up to 6 pts) — ideal: 4-12h for oral dosing
        hl = inputs.get("half_life_h")
        hl_score = 0.0
        if hl is not None:
            if 4 <= hl <= 12:
                hl_score = 6.0   # Ideal range
            elif 2 <= hl < 4 or 12 < hl <= 24:
                hl_score = 4.0   # Acceptable
            elif 1 <= hl < 2 or 24 < hl <= 48:
                hl_score = 2.0   # Suboptimal
            else:
                hl_score = 1.0   # Poor (too short or too long)
        else:
            hl_score = 3.0  # Missing: partial credit

        # Clearance sub-score (up to 3 pts) — mL/min/kg, low is better for most drugs
        cl = inputs.get("clearance_mL_min_kg")
        cl_score = 0.0
        if cl is not None:
            if cl <= 5:
                cl_score = 3.0   # Low clearance (good)
            elif cl <= 15:
                cl_score = 2.0   # Moderate
            elif cl <= 30:
                cl_score = 1.0   # High
            else:
                cl_score = 0.5   # Very high (poor)
        else:
            cl_score = 1.5  # Missing: partial credit

        # Volume of distribution sub-score (up to 3 pts)
        vd = inputs.get("vd_L_per_kg")
        vd_score = 0.0
        if vd is not None:
            if 0.3 <= vd <= 3.0:
                vd_score = 3.0   # Ideal
            elif 0.1 <= vd < 0.3 or 3.0 < vd <= 10.0:
                vd_score = 2.0   # Acceptable
            else:
                vd_score = 1.0   # Extreme
        else:
            vd_score = 1.5  # Missing: partial credit

        score = bio_score + hl_score + cl_score + vd_score

        return {
            "score": round(score, 1),
            "max": max_pts,
            "percentage": round((score / max_pts) * 100, 1),
            "breakdown": {
                "bioavailability": round(bio_score, 1),
                "half_life": round(hl_score, 1),
                "clearance": round(cl_score, 1),
                "vd": round(vd_score, 1),
            },
        }

    def _score_druglikeness(self, inputs: Dict[str, Any],
                            dossier: Dict[str, Any]) -> Dict[str, Any]:
        """Score drug-likeness (10 points max): Lipinski + QED - constraint violations."""
        max_pts = CATEGORY_WEIGHTS["druglikeness"]
        score = 0.0

        # Lipinski sub-score (up to 4 pts)
        violations = inputs.get("lipinski_violations", 0)
        if violations == 0:
            lip_score = 4.0
        elif violations == 1:
            lip_score = 3.0
        elif violations == 2:
            lip_score = 1.5
        else:
            lip_score = 0.0

        # QED sub-score (up to 6 pts)
        qed = inputs.get("qed")
        qed_score = 0.0
        if qed is not None:
            if qed >= 0.7:
                qed_score = 6.0
            elif qed >= 0.5:
                qed_score = 4.0
            elif qed >= 0.3:
                qed_score = 2.0
            else:
                qed_score = 1.0
        else:
            qed_score = 3.0  # Missing: 50% partial credit

        score = lip_score + qed_score

        # Constraint violations: -1 point per violation, max -3
        cv = inputs.get("constraint_violations", [])
        cv_penalty = min(len(cv), 3)
        score = max(0.0, score - cv_penalty)

        return {
            "score": round(score, 1),
            "max": max_pts,
            "percentage": round((score / max_pts) * 100, 1),
            "breakdown": {
                "lipinski": round(lip_score, 1),
                "qed": round(qed_score, 1),
                "constraint_penalty": -cv_penalty,
            },
        }

    def _score_manufacturability(self, inputs: Dict[str, Any]) -> Dict[str, Any]:
        """Score manufacturability (10 points max): synthesis steps + yield."""
        max_pts = CATEGORY_WEIGHTS["manufacturability"]
        score = 0.0

        # Synthesis steps sub-score (up to 7 pts)
        steps = inputs.get("synthesis_steps")
        steps_score = 0.0
        if steps is not None:
            if steps <= 3:
                steps_score = 7.0
            elif steps <= 5:
                steps_score = 5.0
            elif steps <= 8:
                steps_score = 3.0
            elif steps <= 12:
                steps_score = 1.5
            else:
                steps_score = 0.5
        else:
            steps_score = 3.5  # Missing: 50% partial credit

        # Yield sub-score (up to 3 pts) — not computed yet
        yield_score = 1.5  # Missing: partial credit

        score = steps_score + yield_score

        return {
            "score": round(score, 1),
            "max": max_pts,
            "percentage": round((score / max_pts) * 100, 1),
            "breakdown": {
                "synthesis_steps": round(steps_score, 1),
                "yield": round(yield_score, 1),
            },
        }

    def _score_fto(self, inputs: Dict[str, Any]) -> Dict[str, Any]:
        """Score FTO (5 points max): freedom-to-operate status."""
        max_pts = CATEGORY_WEIGHTS["fto"]

        fto = inputs.get("fto_status")
        if fto == "clear":
            score = 5.0
        elif fto == "blocked":
            score = 0.0
        elif fto == "partial" or fto == "review":
            score = 2.5
        else:
            score = 2.5  # Missing: partial credit (not penalized for absence)

        return {
            "score": round(score, 1),
            "max": max_pts,
            "percentage": round((score / max_pts) * 100, 1),
            "breakdown": {"fto_status": fto or "unknown"},
        }

    # ── Tier Assignment ──

    def _assign_tier(self, percentage: float) -> str:
        """Assign tier based on percentage score."""
        if percentage >= TIER_THRESHOLDS["DIAMOND_TIER"]:
            return "DIAMOND_TIER"
        elif percentage >= TIER_THRESHOLDS["GOLD_TIER"]:
            return "GOLD_TIER"
        elif percentage >= TIER_THRESHOLDS["SILVER_TIER"]:
            return "SILVER_TIER"
        else:
            return "BRONZE_TIER"

    # ── Reasoning Builder ──

    def _build_reasoning(self, category_scores: Dict[str, Dict],
                         grade: str,
                         hard_limits: Dict[str, Any]) -> Dict[str, Any]:
        """Build audit-trail reasoning with strengths/weaknesses."""
        reasoning: Dict[str, Any] = {
            "grade": grade,
            "strengths": [],
            "weaknesses": [],
            "summary": "",
        }

        if hard_limits.get("rejected"):
            reasoning["summary"] = "REJECTED due to hard limit violation(s): " + \
                "; ".join(hard_limits["rejection_reasons"])
            reasoning["weaknesses"] = hard_limits["rejection_reasons"]
            return reasoning

        for cat_name, cs in category_scores.items():
            pct = cs.get("percentage", 0)
            if pct >= 80:
                reasoning["strengths"].append(
                    f"{cat_name}: {cs['score']}/{cs['max']} ({pct}%)"
                )
            elif pct < 50:
                reasoning["weaknesses"].append(
                    f"{cat_name}: {cs['score']}/{cs['max']} ({pct}%)"
                )

        total = sum(cs["score"] for cs in category_scores.values())
        pct = round((total / MAX_SCORE) * 100, 1)

        if grade == "DIAMOND_TIER":
            reasoning["summary"] = f"Exceptional candidate scoring {pct}% — strong across all categories"
        elif grade == "GOLD_TIER":
            reasoning["summary"] = f"Strong candidate scoring {pct}% — suitable for advancement"
        elif grade == "SILVER_TIER":
            reasoning["summary"] = f"Moderate candidate scoring {pct}% — optimization recommended"
        else:
            reasoning["summary"] = f"Weak candidate scoring {pct}% — significant improvement needed"

        return reasoning

    # ── Statistics ──

    def get_grade_statistics(self) -> Dict[str, Any]:
        """Get grading statistics (alias for get_statistics)."""
        return self.get_statistics()

    def get_statistics(self) -> Dict[str, Any]:
        """Get grading statistics."""
        if not self.grading_history:
            return {}

        grade_counts: Dict[str, int] = {}
        for entry in self.grading_history:
            g = entry["grade"]
            grade_counts[g] = grade_counts.get(g, 0) + 1

        return {
            "total_graded": len(self.grading_history),
            "grade_distribution": grade_counts,
            "engine_version": self.version,
        }

    # ── Helpers ──

    def _resolve_toxicity(self, dossier: Dict[str, Any]) -> Optional[float]:
        """Resolve toxicity_index from multiple dossier schemas."""
        engines = dossier.get("engines") or {}
        admet = dossier.get("admet") or engines.get("admet") or \
                (dossier.get("inputs") or {}).get("admet_analysis") or {}

        # Path 1: admet.toxicity.toxicity_index (nested)
        tox_block = admet.get("toxicity")
        if isinstance(tox_block, dict):
            ti = tox_block.get("toxicity_index")
            if isinstance(ti, (int, float)):
                return float(ti)

        # Path 2: admet.toxicity_index (flat)
        ti = admet.get("toxicity_index")
        if isinstance(ti, (int, float)):
            return float(ti)

        # Path 3: admet.toxicity (direct float = toxicity value)
        if isinstance(tox_block, (int, float)):
            return float(tox_block)

        # Path 4: top-level harm_energy
        he = dossier.get("harm_energy")
        if isinstance(he, (int, float)):
            return float(he)

        return None

    def _resolve_herg(self, dossier: Dict[str, Any]) -> Optional[float]:
        """Resolve hERG risk from dossier."""
        engines = dossier.get("engines") or {}
        admet = dossier.get("admet") or engines.get("admet") or {}
        tox_block = admet.get("toxicity")
        if isinstance(tox_block, dict):
            herg = tox_block.get("herg_risk")
            if isinstance(herg, (int, float)):
                return float(herg)
        return None

    def _count_lipinski_violations(self, inputs: Dict[str, Any]) -> int:
        """Count Lipinski Rule-of-5 violations from OPE data."""
        violations = 0
        mw = inputs.get("molecular_weight")
        logp = inputs.get("logp")
        hbd = inputs.get("hbd")
        hba = inputs.get("hba")

        if mw is not None and mw > 500:
            violations += 1
        if logp is not None and logp > 5:
            violations += 1
        if hbd is not None and hbd > 5:
            violations += 1
        if hba is not None and hba > 10:
            violations += 1

        return violations

    @staticmethod
    def _safe_float(val) -> Optional[float]:
        """Safely convert to float, returning None on failure."""
        if val is None:
            return None
        if isinstance(val, (int, float)):
            return float(val)
        try:
            return float(val)
        except (TypeError, ValueError):
            return None

    # Legacy compatibility
    def extract_metrics(self, dossier: Dict[str, Any]) -> Dict[str, float]:
        """Legacy method — returns raw inputs as flat metrics dict."""
        return self._extract_raw_inputs(dossier)

    def classify_metrics(self, metrics: Dict[str, float]) -> Tuple[str, Dict[str, Any]]:
        """Legacy method — grade from metrics dict. Wraps weighted scoring."""
        # Build a minimal dossier from flat metrics for scoring
        fake_dossier = {}
        if "toxicity" in metrics:
            fake_dossier["harm_energy"] = metrics["toxicity"]
        result = self.grade_dossier(fake_dossier)
        return result["grade"], result["reasoning"]


# Backward compatibility alias
GradingEngine = WeightedGradingEngine
