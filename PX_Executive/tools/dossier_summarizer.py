"""
dossier_summarizer.py
Extract Key Metrics from Evidence Package v3 Dossiers

Generates partner-facing summary tables from computational outputs.

Usage:
    python PX_Executive/tools/dossier_summarizer.py <folder>
    python PX_Executive/tools/dossier_summarizer.py --batch BATCH_ID
"""

import json
import sys
from pathlib import Path
from glob import glob
import csv

# Add project root to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))


def extract_metrics(dossier_path: Path) -> dict:
    """
    Extract key metrics from Evidence Package v3 dossier.
    
    Args:
        dossier_path: Path to dossier JSON file
    
    Returns:
        Dict with extracted metrics
    """
    try:
        with open(dossier_path, 'r') as f:
            dossier = json.load(f)
        
        # Extract trial result (first arm)
        trial_result = dossier.get("trial_result", {})
        arms = trial_result.get("arms", [])
        
        if not arms:
            return {"error": "No arms found in trial result"}
        
        arm = arms[0]  # First arm
        
        # Extract PK metrics
        exposure = arm.get("exposure_summary", {})
        auc = exposure.get("auc_mg_h_per_L", {})
        cmax = exposure.get("cmax_mg_per_L", {})
        cmin = exposure.get("cmin_steady_state_mg_per_L", {})
        
        # Extract PD metrics
        pd_summary = arm.get("pd_summary", {})
        max_effect = pd_summary.get("max_effect", {})
        auec = pd_summary.get("auec_h", {})
        
        # Extract dose info
        protocol = dossier.get("protocol", {})
        protocol_arms = protocol.get("arms", [{}])
        dose_mg = protocol_arms[0].get("dose_mg", "N/A")
        interval_h = protocol_arms[0].get("dosing_interval_h", "N/A")
        
        # Extract IIV analysis
        iiv = dossier.get("iiv_analysis", {})
        iiv_enabled = iiv.get("enabled", False) if iiv else False
        
        # Extract adaptive analysis
        adaptive = dossier.get("adaptive_analysis", {})
        adaptive_enabled = adaptive.get("enabled", False) if adaptive else False
        
        # Extract dose optimization (if present)
        dose_opt = dossier.get("dose_optimization", {})
        best_regimen = dose_opt.get("best_regimen", {}) if dose_opt else {}
        
        # Extract virtual efficacy (if present)
        efficacy = dossier.get("virtual_efficacy", {})
        pk_pta = efficacy.get("pk_pta", {}) if efficacy else {}
        pd_responders = efficacy.get("pd_responders", {}) if efficacy else {}
        
        # Build metrics dict
        metrics = {
            "dossier_file": dossier_path.name,
            "trial_id": trial_result.get("trial_id", "N/A"),
            "version": dossier.get("version", "N/A"),
            "timestamp": dossier.get("timestamp_utc", "N/A"),
            
            # Protocol
            "protocol_dose_mg": dose_mg,
            "protocol_interval_h": interval_h,
            "n_patients": arm.get("n_patients", "N/A"),
            
            # PK Metrics
            "auc_mean": auc.get("mean", "N/A"),
            "auc_std": auc.get("std", "N/A"),
            "auc_min": auc.get("min", "N/A"),
            "auc_max": auc.get("max", "N/A"),
            "auc_cv_pct": (auc.get("std", 0) / auc.get("mean", 1) * 100) if auc.get("mean") else "N/A",
            
            "cmax_mean": cmax.get("mean", "N/A"),
            "cmax_std": cmax.get("std", "N/A"),
            
            "cmin_mean": cmin.get("mean", "N/A"),
            "cmin_std": cmin.get("std", "N/A"),
            
            # PD Metrics
            "max_effect_mean": max_effect.get("mean", "N/A"),
            "max_effect_std": max_effect.get("std", "N/A"),
            "max_effect_min": max_effect.get("min", "N/A"),
            "max_effect_max": max_effect.get("max", "N/A"),
            
            "auec_mean": auec.get("mean", "N/A"),
            "auec_std": auec.get("std", "N/A"),
            
            # Features
            "iiv_enabled": iiv_enabled,
            "adaptive_enabled": adaptive_enabled,
            
            # Dose Optimization (if present)
            "optimized_dose_mg": best_regimen.get("dose_mg", "N/A"),
            "optimized_interval_h": best_regimen.get("interval_h", "N/A"),
            "dose_score": best_regimen.get("score", "N/A"),
            
            # Virtual Efficacy (if present)
            "pta_auc": pk_pta.get("auc_mg_h_per_L", {}).get("pta", "N/A") if pk_pta else "N/A",
            "responder_rate": pd_responders.get("max_effect", {}).get("responder_rate", "N/A") if pd_responders else "N/A",
        }
        
        return metrics
    
    except Exception as e:
        return {
            "dossier_file": dossier_path.name,
            "error": str(e)
        }


def summarize_folder(folder_path: str, output_csv: str = None) -> list:
    """
    Summarize all dossiers in a folder.
    
    Args:
        folder_path: Path to folder containing dossiers
        output_csv: Optional CSV output path
    
    Returns:
        List of metric dicts
    """
    folder = Path(folder_path)
    
    if not folder.exists():
        print(f"❌ Folder not found: {folder}")
        return []
    
    # Find all dossier files
    dossier_files = list(folder.glob("TRIAL_SIMULATION_DOSSIER-*.json"))
    
    if not dossier_files:
        print(f"⚠️  No dossier files found in: {folder}")
        return []
    
    print(f"\n{'='*70}")
    print(f"DOSSIER SUMMARIZER")
    print(f"{'='*70}")
    print(f"Folder: {folder}")
    print(f"Dossiers found: {len(dossier_files)}")
    print(f"{'='*70}\n")
    
    summaries = []
    for dossier_file in dossier_files:
        print(f"Processing: {dossier_file.name}")
        metrics = extract_metrics(dossier_file)
        summaries.append(metrics)
        
        if "error" in metrics:
            print(f"  ❌ Error: {metrics['error']}")
        else:
            print(f"  ✅ Extracted {len(metrics)} metrics")
    
    # Save to CSV if requested
    if output_csv and summaries:
        csv_path = Path(output_csv)
        csv_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(csv_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=summaries[0].keys())
            writer.writeheader()
            writer.writerows(summaries)
        
        print(f"\n✅ Saved CSV: {csv_path}")
    
    return summaries


def print_summary_table(summaries: list):
    """Print summary table to console"""
    if not summaries:
        return
    
    print(f"\n{'='*70}")
    print(f"SUMMARY TABLE")
    print(f"{'='*70}\n")
    
    for i, summary in enumerate(summaries, 1):
        if "error" in summary:
            print(f"[{i}] {summary['dossier_file']}: ERROR - {summary['error']}")
            continue
        
        print(f"[{i}] {summary['trial_id']}")
        print(f"    Dossier: {summary['dossier_file']}")
        print(f"    Version: {summary['version']}")
        print(f"    ")
        print(f"    Protocol:")
        print(f"      Dose: {summary['protocol_dose_mg']} mg Q{summary['protocol_interval_h']}h")
        print(f"      Patients: {summary['n_patients']}")
        print(f"    ")
        print(f"    PK Metrics:")
        if summary['auc_mean'] != "N/A":
            print(f"      AUC: {summary['auc_mean']:.1f} ± {summary['auc_std']:.1f} mg·h/L (CV: {summary['auc_cv_pct']:.1f}%)")
            print(f"      Cmax: {summary['cmax_mean']:.2f} ± {summary['cmax_std']:.2f} mg/L")
            print(f"      Cmin: {summary['cmin_mean']:.2f} ± {summary['cmin_std']:.2f} mg/L")
        
        print(f"    ")
        print(f"    PD Metrics:")
        if summary['max_effect_mean'] != "N/A":
            print(f"      Max Effect: {summary['max_effect_mean']:.3f} ± {summary['max_effect_std']:.3f}")
            print(f"      AUEC: {summary['auec_mean']:.1f} ± {summary['auec_std']:.1f} h")
        
        if summary['optimized_dose_mg'] != "N/A":
            print(f"    ")
            print(f"    Dose Optimization:")
            print(f"      Best Regimen: {summary['optimized_dose_mg']:.1f} mg Q{summary['optimized_interval_h']:.0f}h")
            print(f"      Score: {summary['dose_score']:.3f}")
        
        if summary['pta_auc'] != "N/A":
            print(f"    ")
            print(f"    Virtual Efficacy:")
            print(f"      PTA (AUC): {summary['pta_auc']*100:.1f}%")
            if summary['responder_rate'] != "N/A":
                print(f"      Responder Rate: {summary['responder_rate']*100:.1f}%")
        
        print()


def main():
    """Main entry point"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Extract key metrics from Evidence Package v3 dossiers",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Summarize dossiers in a folder
  python PX_Executive/tools/dossier_summarizer.py PX_Warehouse/TrialSimulations/LiveRuns/run_20260126_140402
  
  # Summarize batch run
  python PX_Executive/tools/dossier_summarizer.py --batch BATCH_20260126_123456
  
  # Export to CSV
  python PX_Executive/tools/dossier_summarizer.py PX_Warehouse/TrialSimulations/LiveRuns/run_20260126_140402 --csv summary.csv
        """
    )
    
    parser.add_argument(
        "folder",
        nargs='?',
        default=None,
        help="Folder containing dossiers"
    )
    
    parser.add_argument(
        "--batch",
        type=str,
        default=None,
        help="Batch ID (searches in BatchRuns folder)"
    )
    
    parser.add_argument(
        "--csv",
        type=str,
        default=None,
        help="Export to CSV file"
    )
    
    args = parser.parse_args()
    
    # Determine folder
    if args.batch:
        folder = project_root / "PX_Warehouse" / "TrialSimulations" / "BatchRuns" / args.batch
    elif args.folder:
        folder = Path(args.folder)
    else:
        print("❌ Error: Must specify folder or --batch")
        parser.print_help()
        return 1
    
    # Summarize
    summaries = summarize_folder(str(folder), output_csv=args.csv)
    
    if summaries:
        print_summary_table(summaries)
        return 0
    else:
        return 1


if __name__ == "__main__":
    sys.exit(main())
