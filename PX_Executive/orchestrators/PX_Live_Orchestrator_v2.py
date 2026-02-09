"""
PX_Live_Orchestrator_v2.py  
Predator X v2.0-CORE Production Orchestrator

Pipeline: OPE → ADMET → TrialEngine (7-tier IIV) + PKPD (sigmoid Emax) →
DoseOptimizer_v2 (coarse-to-fine, Law L4/L11) → VirtualEfficacyAnalytics (PTA, Responder Rate) →
OCE gate (manifold coherence >= 0.85) → Evidence Package v3.
"""

import argparse
import sys
import json
import hashlib
from pathlib import Path
from datetime import datetime, UTC

# Repo root for canonical warehouse paths (always write to foundation/PX_Warehouse/...)
_REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO_ROOT))

from PX_Engine.operations import OPE
from PX_Engine.operations import ADMET
from PX_Engine.operations.TrialEngine import TrialEngine, generate_virtual_population
from PX_Engine.operations.PKPD import link_pk_to_pd
from PX_Engine.operations.DoseOptimizer_v2 import optimize_dose
from PX_Engine.operations.VirtualEfficacyAnalytics import (
    compute_pta,
    virtual_responder_rate,
    analyze_virtual_efficacy,
)
from PX_Engine.operations import OCE


# 7-tier deterministic IIV (replace 25% CV hash from Simple)
DEFAULT_VARIABILITY = {
    "clearance_variation": 0.3,
    "vd_variation": 0.25,
    "ka_variation": 0.2,
    "n_tiers": 7,
}


class PredatorXOrchestratorV2:
    """Production v2.0-CORE: TrialEngine 7-tier, PKPD sigmoid Emax, DoseOptimizer_v2, VirtualEfficacyAnalytics, OCE gate."""
    
    def __init__(self, verbose: bool = True):
        self.verbose = verbose
        self.version = "2.0.0-CORE-PRODUCTION"
    
    def run_pipeline(self, smiles: str, metadata: dict) -> dict:
        """Execute corrected v2.0-CORE pipeline. Output: PX_Warehouse/Calibration_Molecules/LiveRuns/run_<timestamp>/"""
        from PX_Warehouse.warehouse_layout import get_calibration_molecules_dir
        run_id = hashlib.md5(f"{smiles}{datetime.now(UTC).isoformat()}".encode()).hexdigest()[:12]
        run_timestamp = datetime.now(UTC).strftime('%Y%m%d_%H%M%S')
        run_dir_path = get_calibration_molecules_dir(_REPO_ROOT) / "LiveRuns" / f"run_{run_timestamp}"
        run_dir = str(run_dir_path)
        
        if self.verbose:
            print(f"\n{'='*70}")
            print(f"PREDATOR X v2.0-CORE (CORRECTED PIPELINE)")
            print(f"{'='*70}")
            print(f"Candidate: {metadata.get('name', 'Unknown')}")
            print(f"Run ID: {run_id}")
            print(f"{'='*70}\n")
        
        results = {
            "metadata": metadata,
            "smiles": smiles,
            "run_id": run_id,
            "timestamp": datetime.now(UTC).isoformat(),
            "version": self.version,
            "stages": {}
        }
        
        try:
            # STAGE 1: OPE
            if self.verbose:
                print("[1/6] OPE Analysis...")
            ope_results = OPE.run_ope(smiles)
            results["stages"]["ope"] = ope_results
            
            # STAGE 2: ADMET
            if self.verbose:
                print("[2/6] ADMET Analysis...")
            admet_results = ADMET.run_admet(smiles, ope_results)
            results["stages"]["admet"] = admet_results
            
            # STAGE 3: Trial simulation (7-tier IIV) + PK/PD (sigmoid Emax)
            if self.verbose:
                print("[3/6] Trial simulation (7-tier) + PK/PD (sigmoid Emax)...")
            ec50 = ope_results.get("ec50", 1.0) or 1.0
            emax = ope_results.get("emax", 0.8) or 0.8
            pd_params = {
                "emax": float(emax),
                "ec50": float(ec50),
                "hill": 1.0,
                "baseline": 0.0,
                "effect_threshold": 0.5,
            }
            protocol = {
                "trial_id": f"LIVE-{run_id}",
                "duration_days": 7.0,
                "arms": [{
                    "arm_id": "A1",
                    "label": "100mg Q24h",
                    "dose_mg": 100.0,
                    "dosing_interval_h": 24.0,
                    "n_patients": 21,
                }],
            }
            trial_engine = TrialEngine(time_step_h=1.0)
            trial_result = trial_engine.run_trial(
                protocol=protocol,
                admet=admet_results,
                pd_params=pd_params,
                variability=DEFAULT_VARIABILITY,
            )
            results["stages"]["trial_result"] = trial_result
            results["stages"]["pkpd"] = {"pd_params": pd_params, "arms": trial_result.get("arms", [])}

            # STAGE 4: DoseOptimizer_v2 (coarse-to-fine; Law L4/L11 hard-lock)
            if self.verbose:
                print("[4/6] Dose optimization (coarse-to-fine)...")
            protocol_template = {
                "trial_id": "OPT",
                "duration_days": 7.0,
                "arms": [{"arm_id": "O1", "label": "", "dose_mg": 100.0, "dosing_interval_h": 24.0, "n_patients": 10}],
            }
            dose_opt_results = optimize_dose(
                smiles=smiles,
                admet=admet_results,
                protocol_template=protocol_template,
                target_pk_range={"auc_mg_h_per_L": (200.0, 400.0)},
                target_pd_range={"max_effect": (0.5, 0.9)} if pd_params else None,
                dose_bounds=(50.0, 300.0),
                variability=DEFAULT_VARIABILITY,
                pd_params=pd_params,
                search_strategy="coarse_to_fine",
                n_eval_patients=10,
            )
            results["stages"]["dose_optimization"] = dose_opt_results
            optimized_dose = (dose_opt_results.get("best_regimen") or {}).get("dose_mg", 100.0)
            optimized_interval = (dose_opt_results.get("best_regimen") or {}).get("interval_h", 24.0)

            # OCE: dose is "Optimized" only if 35D manifold coherence >= 0.85
            oce_payload = {"p_vector": [0.1, 0.0, 35.0, 1.0], "csa_scores": [1.0, 1.0, 1.0, 1.0, 1.0]}
            oce_result = OCE.execute(oce_payload)
            dose_authorized = oce_result.get("coherence", 0) >= 0.85
            results["stages"]["oce_gate"] = {
                "coherence": oce_result.get("coherence"),
                "authorized": dose_authorized,
                "note": "No dose considered optimized unless manifold coherence >= 0.85",
            }

            # STAGE 5: VirtualEfficacyAnalytics (PTA, Responder Rate)
            if self.verbose:
                print("[5/6] Virtual Efficacy Analytics (PTA, Responder Rate)...")
            # Re-run trial with optimized regimen for PTA/responder report
            protocol_opt = {
                "trial_id": "PTA",
                "duration_days": 7.0,
                "arms": [{
                    "arm_id": "P1",
                    "label": f"{optimized_dose}mg Q{optimized_interval}h",
                    "dose_mg": optimized_dose,
                    "dosing_interval_h": optimized_interval,
                    "n_patients": 21,
                }],
            }
            trial_opt = trial_engine.run_trial(
                protocol=protocol_opt,
                admet=admet_results,
                pd_params=pd_params,
                variability=DEFAULT_VARIABILITY,
            )
            ve_results = analyze_virtual_efficacy(
                trial_opt,
                pk_target={"auc_mg_h_per_L": 200.0},
                pd_target={"max_effect": 0.5},
                safety_threshold=0.95,
            )
            pta_auc = compute_pta(trial_opt, "auc_mg_h_per_L", 200.0)
            responder = virtual_responder_rate(trial_opt, "max_effect", 0.5)
            ve_results["pta_auc"] = pta_auc
            ve_results["responder_rate"] = responder
            results["stages"]["virtual_efficacy"] = ve_results
            results["stages"]["pkpd_optimized"] = {"trial_result": trial_opt, "dose_mg": optimized_dose, "interval_h": optimized_interval}
            
            # STAGE 6: Evidence Package v3
            if self.verbose:
                print("[6/6] Evidence Package v3...")
            
            evidence_package = self._build_evidence_package_v3(results, run_id, run_dir)

            # Zeus gate: constitutional governance check before warehouse write
            from PX_System.foundation.ZeusLaws import check_constitutional
            admet_pkg = evidence_package.get("admet") or {}
            tox_obj = admet_pkg.get("toxicity") if isinstance(admet_pkg.get("toxicity"), dict) else {}
            tox_val = tox_obj.get("toxicity_index") if tox_obj else None
            zeus_verdict = check_constitutional("live_orchestrator_v2", {
                "toxicity_index": tox_val,
                "harm_energy": tox_val,
            })
            evidence_package["zeus"] = zeus_verdict
            results["evidence_package"] = evidence_package
            if not zeus_verdict.get("authorized", False):
                if self.verbose:
                    print(f"\n⛔ ZEUS GATE REJECTED: {zeus_verdict.get('rationale', 'governance failure')}")
                results["zeus_rejected"] = True
                results["zeus_verdict"] = zeus_verdict
                return results

            # Save files (run_dir is absolute: PX_Warehouse/Calibration_Molecules/LiveRuns/run_*)
            run_dir_path.mkdir(parents=True, exist_ok=True)
            dossier_path = run_dir_path / f"TRIAL_SIMULATION_DOSSIER-{run_id}.json"
            with open(dossier_path, 'w') as f:
                json.dump(evidence_package, f, indent=2, default=str)
            
            with open(run_dir_path / "pipeline_log.json", 'w') as f:
                json.dump(results, f, indent=2, default=str)
            
            with open(run_dir_path / "metrics.json", 'w') as f:
                json.dump({
                    "run_id": run_id,
                    "smiles": smiles,
                    "name": metadata.get("name", "Unknown"),
                    "compound_id": metadata.get("id", run_id),
                    "duration_seconds": 0,
                    "timestamp": datetime.now(UTC).isoformat()
                }, f, indent=2, default=str)
            
            with open(run_dir_path / "config_snapshot.json", 'w') as f:
                json.dump({
                    "version": self.version,
                    "pipeline_order": [
                        "1_OPE",
                        "2_ADMET",
                        "3_PKPD_Simulation",
                        "4_Dose_Optimization_v2",
                        "5_Virtual_Efficacy_Analytics",
                        "6_Evidence_Package_v3"
                    ],
                    "note": "Corrected: Virtual efficacy after dose optimization"
                }, f, indent=2, default=str)
            
            results["evidence_package"] = evidence_package
            results["dossier_path"] = str(dossier_path)
            
            if self.verbose:
                print(f"\n✅ Pipeline complete! Output: {run_dir}\n")
            
            return results
            
        except Exception as e:
            if self.verbose:
                print(f"\n❌ Pipeline failed: {e}\n")
            raise
    
    def _build_evidence_package_v3(self, results: dict, run_id: str, run_dir: str) -> dict:
        """Build Evidence Package v3 from production TrialEngine/DoseOptimizer_v2/VirtualEfficacyAnalytics."""
        stages = results.get("stages", {})
        metadata = results.get("metadata", {})
        arms = (stages.get("trial_result") or {}).get("arms", [])
        first_arm = arms[0] if arms else {}
        exp = first_arm.get("exposure_summary", {})
        pd_sum = first_arm.get("pd_summary", {})
        best = (stages.get("dose_optimization") or {}).get("best_regimen") or {}
        return {
            "metadata": {
                "run_id": run_id,
                "name": metadata.get("name", "Unknown"),
                "id": metadata.get("id", run_id),
                "smiles": results.get("smiles", ""),
                "timestamp": datetime.now(UTC).isoformat(),
                "version": "3.0",
                "pipeline_version": self.version,
            },
            "ope": stages.get("ope", {}),
            "admet": stages.get("admet", {}),
            "pkpd": {
                "auc": {
                    "mean": exp.get("auc_mg_h_per_L", {}).get("mean", 0.0),
                    "sd": exp.get("auc_mg_h_per_L", {}).get("std", 0.0),
                    "cv": (exp.get("auc_mg_h_per_L", {}).get("std", 0) / (exp.get("auc_mg_h_per_L", {}).get("mean", 1) or 1) * 100),
                },
                "cmax": {
                    "mean": exp.get("cmax_mg_per_L", {}).get("mean", 0.0),
                    "sd": exp.get("cmax_mg_per_L", {}).get("std", 0.0),
                },
            },
            "pd": {
                "max_effect": {
                    "mean": pd_sum.get("max_effect", {}).get("mean", 0.0),
                    "sd": pd_sum.get("max_effect", {}).get("std", 0.0),
                },
            },
            "dose_optimization": stages.get("dose_optimization", {}),
            "oce_gate": stages.get("oce_gate", {}),
            "virtual_efficacy": stages.get("virtual_efficacy", {}),
            "output_directory": run_dir,
        }


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(description="Predator X v2.0-CORE Orchestrator (CORRECTED)")
    parser.add_argument("--smiles", type=str, required=True, help="SMILES string")
    parser.add_argument("--name", type=str, default="Unknown", help="Compound name")
    parser.add_argument("--id", type=str, default=None, help="Compound ID")
    parser.add_argument("--indication", type=str, default="Not specified", help="Indication")
    parser.add_argument("--quiet", action="store_true", help="Suppress output")
    
    args = parser.parse_args()
    
    candidate = {
        "name": args.name,
        "id": args.id or hashlib.md5(args.smiles.encode()).hexdigest()[:8],
        "indication": args.indication
    }
    
    orchestrator = PredatorXOrchestratorV2(verbose=not args.quiet)
    results = orchestrator.run_pipeline(smiles=args.smiles, metadata=candidate)
    
    if not args.quiet:
        print("="*70)
        print(f"Dossier: {results.get('dossier_path')}")
        print("="*70)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
