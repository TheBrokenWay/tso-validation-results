"""
GradingEngine.py
Constitutional 5-Tier Grading System for Predator X v2.0-CORE

Classifies computational dossiers into:
- GOLD_TIER: Exceptional candidates (discovery-stage metrics)
- SILVER_TIER: Strong candidates
- BRONZE_TIER: Moderate candidates
- NEEDS_REVIEW: Ambiguous candidates (2-3 criteria met)
- REJECTED: Failed candidates (<2 criteria)

Updated: 2026-01-26 - Discovery-stage thresholds (10-30% range)
"""

import json
import sys
import time
from pathlib import Path
from typing import Dict, Any, Tuple
from datetime import datetime, UTC
from PX_System.foundation.sign_off import create_sign_off


class GradingEngine:
    """Constitutional grading engine for computational dossiers"""
    TOXICITY_LIMIT = 0.0210
    
    # Discovery-Stage Thresholds (Updated 2026-01-26)
    THRESHOLDS = {
        "GOLD_TIER": {
            "pta_min": 30.0,              # 30% PTA (discovery-stage)
            "responder_rate_min": 30.0,   # 30% responder rate
            "toxicity_max": 0.0210,        # Constitutional hard-lock
            "dose_score_min": 0.8,         # High dose optimization
            "variability_cv_max": 30.0,    # Low variability
        },
        "SILVER_TIER": {
            "pta_min": 20.0,              # 20% PTA
            "responder_rate_min": 20.0,   # 20% responder rate
            "toxicity_max": 0.0210,        # Constitutional hard-lock
            "dose_score_min": 0.6,         # Moderate dose optimization
            "variability_cv_max": 50.0,    # Moderate variability
        },
        "BRONZE_TIER": {
            "pta_min": 10.0,              # 10% PTA
            "responder_rate_min": 10.0,   # 10% responder rate
            "toxicity_max": 0.0210,        # Constitutional hard-lock
            "dose_score_min": 0.4,         # Basic dose optimization
            "variability_cv_max": 70.0,    # Higher variability acceptable
        },
    }
    
    def __init__(self, verbose: bool = True):
        """Initialize grading engine"""
        self.verbose = verbose
        self.version = "v2.0-discovery"
        self.grading_history = []
        self.load_schema()
    
    def load_schema(self):
        """Load grading schema from file if exists, otherwise use defaults"""
        schema_path = Path("PX_Engine/operations/GradingSchema_Discovery.json")
        if schema_path.exists():
            with open(schema_path, 'r') as f:
                schema = json.load(f)
                # Update thresholds from schema
                if "GOLD_TIER" in schema:
                    self.THRESHOLDS["GOLD_TIER"]["pta_min"] = schema["GOLD_TIER"]["pta_min"] * 100
                    self.THRESHOLDS["GOLD_TIER"]["responder_rate_min"] = schema["GOLD_TIER"]["responder_min"] * 100
                    self.THRESHOLDS["GOLD_TIER"]["toxicity_max"] = schema["GOLD_TIER"]["toxicity_max"]
                if "SILVER_TIER" in schema:
                    self.THRESHOLDS["SILVER_TIER"]["pta_min"] = schema["SILVER_TIER"]["pta_min"] * 100
                    self.THRESHOLDS["SILVER_TIER"]["responder_rate_min"] = schema["SILVER_TIER"]["responder_min"] * 100
                    self.THRESHOLDS["SILVER_TIER"]["toxicity_max"] = schema["SILVER_TIER"]["toxicity_max"]
                if "BRONZE_TIER" in schema:
                    self.THRESHOLDS["BRONZE_TIER"]["pta_min"] = schema["BRONZE_TIER"]["pta_min"] * 100
                    self.THRESHOLDS["BRONZE_TIER"]["responder_rate_min"] = schema["BRONZE_TIER"]["responder_min"] * 100
                    self.THRESHOLDS["BRONZE_TIER"]["toxicity_max"] = schema["BRONZE_TIER"]["toxicity_max"]
            if self.verbose:
                print(f"✅ Loaded discovery-stage grading schema from {schema_path}")
    
    def extract_metrics(self, dossier: Dict[str, Any]) -> Dict[str, float]:
        """
        Extract key metrics from Evidence Package v3 dossier
        
        Args:
            dossier: Complete Evidence Package v3 JSON
        
        Returns:
            Dictionary of extracted metrics
        """
        metrics = {
            "pta": 0.0,
            "responder_rate": 0.0,
            "toxicity": 0.5,
            "dose_score": 0.5,
            "variability_cv": 0.0,
            "binding_affinity": 0.5,
        }
        
        # Extract PTA (Probability of Target Attainment)
        try:
            if "virtual_efficacy" in dossier:
                ve = dossier["virtual_efficacy"]
                if "pk_pta" in ve:
                    pk_pta = ve["pk_pta"]
                    if "auc_mg_h_per_L" in pk_pta:
                        metrics["pta"] = pk_pta["auc_mg_h_per_L"].get("pta", 0.0)
        except (KeyError, TypeError, AttributeError) as e:
            print(f"    WARN: PTA extraction failed: {e}", file=sys.stderr)
        
        # Extract Responder Rate (effect_ge_threshold or pd_responders path)
        try:
            if "virtual_efficacy" in dossier:
                ve = dossier["virtual_efficacy"]
                if "responder_rate" in ve:
                    rr = ve["responder_rate"]
                    if "effect_ge_threshold" in rr:
                        metrics["responder_rate"] = rr["effect_ge_threshold"] * 100.0
                elif "pd_responders" in ve:
                    pd = ve["pd_responders"]
                    if "max_effect" in pd and "responder_rate" in pd["max_effect"]:
                        metrics["responder_rate"] = pd["max_effect"]["responder_rate"] * 100.0
        except (KeyError, TypeError, AttributeError) as e:
            print(f"    WARN: responder_rate extraction failed: {e}", file=sys.stderr)
        
        # Extract Toxicity Index (from ADMET: toxicity_index or toxicity sub-dict)
        # Resolve ADMET from multiple dossier schemas:
        #   Evidence Package v2: dossier["engines"]["admet"]
        #   Trial Simulation v3: dossier["inputs"]["admet_analysis"]
        #   Normalized:          dossier["admet"]
        try:
            admet = dossier.get("admet")
            if admet is None:
                admet = (dossier.get("engines") or {}).get("admet")
            if admet is None:
                admet = (dossier.get("inputs") or {}).get("admet_analysis")
            if admet is not None:
                if "toxicity_index" in admet:
                    v = admet["toxicity_index"]
                    metrics["toxicity"] = float(v) if isinstance(v, (int, float)) else 0.5
                elif "toxicity" in admet:
                    tox = admet["toxicity"]
                    if isinstance(tox, (int, float)):
                        metrics["toxicity"] = float(tox)
                    elif isinstance(tox, dict) and "toxicity_index" in tox:
                        metrics["toxicity"] = float(tox["toxicity_index"])
                    else:
                        metrics["toxicity"] = 0.5
        except (KeyError, TypeError, AttributeError) as e:
            print(f"    WARN: toxicity extraction failed (defaults to 0.5 → REJECTED): {e}", file=sys.stderr)

        # Fallback: top-level harm_energy (present in PRV_NOV_ and PRV_REP_ dossiers)
        if metrics["toxicity"] == 0.5 and "harm_energy" in dossier:
            try:
                he = dossier["harm_energy"]
                if isinstance(he, (int, float)):
                    metrics["toxicity"] = float(he)
            except (TypeError, ValueError) as e:
                print(f"    WARN: harm_energy fallback failed: {e}", file=sys.stderr)

        # Extract Dose Optimization Score (safety_margin or optimization_score)
        try:
            if "dose_optimization" in dossier:
                do = dossier["dose_optimization"]
                if "best_regimen" in do:
                    best_regimen = do["best_regimen"]
                    if "safety_margin" in best_regimen:
                        metrics["dose_score"] = min(1.0, best_regimen["safety_margin"] / 2.0)
                    elif "optimization_score" in best_regimen:
                        metrics["dose_score"] = float(best_regimen["optimization_score"])
        except (KeyError, TypeError, AttributeError) as e:
            print(f"    WARN: dose_score extraction failed: {e}", file=sys.stderr)
        
        # Extract Variability CV (pkpd.auc.cv or trial_result.arms[0].exposure_summary.auc_mg_h_per_L.sd)
        try:
            if "pkpd" in dossier:
                pkpd = dossier["pkpd"]
                if "auc" in pkpd and "cv" in pkpd["auc"]:
                    metrics["variability_cv"] = pkpd["auc"]["cv"]
            if metrics["variability_cv"] == 0.0 and "trial_result" in dossier:
                arms = dossier["trial_result"].get("arms", [])
                if arms and "exposure_summary" in arms[0]:
                    es = arms[0]["exposure_summary"]
                    auc = es.get("auc_mg_h_per_L", {})
                    if "sd" in auc and "mean" in auc and auc["mean"]:
                        metrics["variability_cv"] = (auc["sd"] / auc["mean"]) * 100.0
        except (KeyError, TypeError, AttributeError, ZeroDivisionError) as e:
            print(f"    WARN: variability_cv extraction failed: {e}", file=sys.stderr)
        
        # Extract Binding Affinity (binding_affinity_nM or score 0-1 or 0-100)
        # Resolve OPE from multiple dossier schemas:
        #   Evidence Package v2: dossier["engines"]["ope"]
        #   Trial Simulation v3: dossier["inputs"]["ope_analysis"]
        #   Normalized:          dossier["ope"]
        try:
            ope = dossier.get("ope")
            if ope is None:
                ope = (dossier.get("engines") or {}).get("ope")
            if ope is None:
                ope = (dossier.get("inputs") or {}).get("ope_analysis")
            if ope is not None:
                ope = ope or {}
                aff = ope.get("binding_affinity")
                if aff is not None:
                    val = float(aff["score"]) if isinstance(aff, dict) and "score" in aff else float(aff)
                    metrics["binding_affinity"] = min(1.0, val) if val <= 1.0 else min(1.0, val / 100.0)
                elif "binding_affinity_nM" in ope:
                    nM = ope["binding_affinity_nM"]
                    metrics["binding_affinity"] = min(1.0, 1.0 - (float(nM) / 1000.0)) if isinstance(nM, (int, float)) else 0.5
        except (KeyError, TypeError, AttributeError) as e:
            print(f"    WARN: binding_affinity extraction failed: {e}", file=sys.stderr)

        # Coerce all metrics to float so classify_metrics never sees dict
        for k in list(metrics.keys()):
            v = metrics[k]
            if isinstance(v, dict):
                metrics[k] = 0.5
            elif not isinstance(v, (int, float)):
                try:
                    metrics[k] = float(v)
                except (TypeError, ValueError):
                    metrics[k] = 0.5
        
        return metrics
    
    def classify_metrics(self, metrics: Dict[str, float]) -> Tuple[str, Dict[str, Any]]:
        """
        Classify metrics into constitutional grade
        
        Args:
            metrics: Extracted metrics dictionary
        
        Returns:
            Tuple of (grade, reasoning)
        """
        # Law L11 Hard-Lock: reject any dossier at or above toxicity limit
        if metrics["toxicity"] >= self.TOXICITY_LIMIT:
            return "REJECTED", {
                "grade": "REJECTED",
                "reason": f"Law L11 Violation: toxicity {metrics['toxicity']:.6f} >= {self.TOXICITY_LIMIT}",
                "criteria_met_count": 0,
            }
        # Check GOLD_TIER
        gold = self.THRESHOLDS["GOLD_TIER"]
        gold_criteria_met = sum([
            metrics["pta"] >= gold["pta_min"],
            metrics["responder_rate"] >= gold["responder_rate_min"],
            metrics["toxicity"] <= gold["toxicity_max"],
            metrics["dose_score"] >= gold["dose_score_min"],
            metrics["variability_cv"] <= gold["variability_cv_max"] or metrics["variability_cv"] == 0.0,
        ])
        
        if gold_criteria_met >= 5:
            return "GOLD_TIER", {
                "grade": "GOLD_TIER",
                "reason": "Exceptional - meets all 5 GOLD criteria",
                "criteria_met_count": gold_criteria_met,
            }
        
        # Check SILVER_TIER
        silver = self.THRESHOLDS["SILVER_TIER"]
        silver_criteria_met = sum([
            metrics["pta"] >= silver["pta_min"],
            metrics["responder_rate"] >= silver["responder_rate_min"],
            metrics["toxicity"] <= silver["toxicity_max"],
            metrics["dose_score"] >= silver["dose_score_min"],
            metrics["variability_cv"] <= silver["variability_cv_max"] or metrics["variability_cv"] == 0.0,
        ])
        
        if silver_criteria_met >= 5:
            return "SILVER_TIER", {
                "grade": "SILVER_TIER",
                "reason": "Strong - meets all 5 SILVER criteria",
                "criteria_met_count": silver_criteria_met,
            }
        
        # Check BRONZE_TIER
        bronze = self.THRESHOLDS["BRONZE_TIER"]
        bronze_criteria_met = sum([
            metrics["pta"] >= bronze["pta_min"],
            metrics["responder_rate"] >= bronze["responder_rate_min"],
            metrics["toxicity"] <= bronze["toxicity_max"],
            metrics["dose_score"] >= bronze["dose_score_min"],
            metrics["variability_cv"] <= bronze["variability_cv_max"] or metrics["variability_cv"] == 0.0,
        ])
        
        if bronze_criteria_met >= 5:
            return "BRONZE_TIER", {
                "grade": "BRONZE_TIER",
                "reason": "Moderate - meets all 5 BRONZE criteria",
                "criteria_met_count": bronze_criteria_met,
            }
        
        # NEEDS_REVIEW: 2-3 criteria met
        if bronze_criteria_met >= 2:
            return "NEEDS_REVIEW", {
                "grade": "NEEDS_REVIEW",
                "reason": f"Ambiguous - meets {bronze_criteria_met}/5 BRONZE criteria",
                "criteria_met_count": bronze_criteria_met,
            }
        
        # REJECTED: <2 criteria met
        return "REJECTED", {
            "grade": "REJECTED",
            "reason": f"Insufficient - meets only {bronze_criteria_met}/5 criteria",
            "criteria_met_count": bronze_criteria_met,
        }
    
    def grade_dossier(self, dossier_path: str | Path | Dict[str, Any]) -> Dict[str, Any]:
        """
        Grade a complete Evidence Package v3 dossier.

        Args:
            dossier_path: Path to dossier JSON file, or dossier dict in memory.

        Returns:
            Grading result with grade, metrics, reasoning, and metadata.
        """
        _t0 = time.monotonic()
        if isinstance(dossier_path, dict):
            dossier = dossier_path
            dossier_ref = "<in-memory>"
        else:
            with open(dossier_path, "r", encoding="utf-8") as f:
                dossier = json.load(f)
            dossier_ref = str(dossier_path)

        metrics = self.extract_metrics(dossier)
        # Spec 3: trial-aware grading metrics (trial_responder_rate, trial_pta, trial_variability_cv, trial_effect_size)
        engines = dossier.get("engines") or {}
        ve = engines.get("virtual_efficacy") or dossier.get("virtual_efficacy") or {}
        trial_result = dossier.get("trial_result") or engines.get("trial_result") or {}
        if ve or trial_result.get("arms"):
            metrics["trial_responder_rate"] = metrics.get("responder_rate")
            metrics["trial_pta"] = metrics.get("pta")
            metrics["trial_variability_cv"] = metrics.get("variability_cv")
            metrics["trial_effect_size"] = ve.get("effect_size") if isinstance(ve, dict) else None
            metrics["trial_toxicity_rate"] = metrics.get("toxicity")
        else:
            metrics["trial_responder_rate"] = None
            metrics["trial_pta"] = None
            metrics["trial_variability_cv"] = None
            metrics["trial_effect_size"] = None
            metrics["trial_toxicity_rate"] = None
        grade, reasoning = self.classify_metrics(metrics)

        # Olympus constraint enforcement: cap grade if disease constraint violations exist
        constraint_violations = dossier.get("constraint_violations")
        if constraint_violations and isinstance(constraint_violations, list) and len(constraint_violations) > 0:
            if grade in ("GOLD_TIER", "DIAMOND_TIER"):
                grade = "SILVER_TIER"
                reasoning["constraint_cap_applied"] = True
                reasoning["constraint_violations_count"] = len(constraint_violations)

        if metrics.get("trial_pta") is not None or metrics.get("trial_toxicity_rate") is not None:
            reasoning["trial_integrated"] = True
            reasoning["trial_responder_rate"] = metrics.get("trial_responder_rate")
            reasoning["trial_pta"] = metrics.get("trial_pta")
            reasoning["trial_effect_size"] = metrics.get("trial_effect_size")
            reasoning["trial_toxicity_rate"] = metrics.get("trial_toxicity_rate")
        result = {
            "grade": grade,
            "metrics": metrics,
            "reasoning": reasoning,
            "timestamp": datetime.now(UTC).isoformat(),
            "grading_engine_version": self.version,
            "thresholds_used": self.THRESHOLDS,
        }
        _elapsed_ms = int((time.monotonic() - _t0) * 1000)
        result["sign_off"] = create_sign_off(
            engine_id="GRADING_ENGINE_V2",
            version="v2.0-discovery",
            inputs={},
            outputs=result,
            laws_checked=["GRADE_RULES"],
            laws_results={"GRADE_RULES": grade != "REJECTED"},
            execution_time_ms=_elapsed_ms,
        )
        self.grading_history.append({
            "dossier": dossier_ref,
            "grade": grade,
            "timestamp": result["timestamp"],
        })
        return result
    
    def get_grade_statistics(self) -> Dict[str, Any]:
        """Get grading statistics (alias for get_statistics)."""
        return self.get_statistics()

    def get_statistics(self) -> Dict[str, Any]:
        """Get grading statistics"""
        if not self.grading_history:
            return {}
        
        grade_counts = {}
        for entry in self.grading_history:
            grade = entry["grade"]
            grade_counts[grade] = grade_counts.get(grade, 0) + 1
        
        return {
            "total_graded": len(self.grading_history),
            "grade_distribution": grade_counts,
            "engine_version": self.version,
        }
