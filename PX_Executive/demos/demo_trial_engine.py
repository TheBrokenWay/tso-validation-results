"""
TrialEngine Demo - Virtual clinical trial simulation
Run this to see the TrialEngine in action!
"""

from PX_Engine.operations import TrialEngine, run_ope, run_admet


def demo_trial_engine():
    """Demonstrate TrialEngine capabilities with a 2-arm dose comparison study"""
    
    print("=" * 80)
    print("PREDATOR X - TRIAL ENGINE DEMONSTRATION")
    print("Virtual Clinical Trial Simulation (Exposure-Only)")
    print("=" * 80)
    
    # Example: Aspirin dose comparison study
    smiles = "CC(=O)Oc1ccccc1C(=O)O"
    print(f"\n📊 Molecule: Aspirin (SMILES: {smiles})")
    
    # Step 1: Run ADMET analysis
    print("\n[1] Running ADMET Analysis...")
    ope = run_ope(smiles)
    admet = run_admet(smiles, ope)
    print(f"    Hepatotoxicity Risk: {admet['toxicity_flags']['hepatotoxicity_risk']}")
    
    # Step 2: Define trial protocol
    print("\n[2] Defining Trial Protocol...")
    protocol = {
        "trial_id": "TRIAL-PK-ASPIRIN-001",
        "duration_days": 7.0,  # 7-day study
        "arms": [
            {
                "arm_id": "A1",
                "label": "Low Dose (QD 75 mg)",
                "dose_mg": 75.0,
                "dosing_interval_h": 24.0,
                "n_patients": 20,
            },
            {
                "arm_id": "A2",
                "label": "Standard Dose (QD 100 mg)",
                "dose_mg": 100.0,
                "dosing_interval_h": 24.0,
                "n_patients": 20,
            },
            {
                "arm_id": "A3",
                "label": "High Dose (BID 75 mg)",
                "dose_mg": 75.0,
                "dosing_interval_h": 12.0,
                "n_patients": 20,
            },
        ],
    }
    
    print(f"    Trial ID: {protocol['trial_id']}")
    print(f"    Duration: {protocol['duration_days']} days")
    print(f"    Arms: {len(protocol['arms'])}")
    for arm in protocol['arms']:
        print(f"      • {arm['label']}: {arm['n_patients']} virtual patients")
    
    # Step 3: Run trial simulation
    print("\n[3] Running Trial Simulation...")
    engine = TrialEngine(time_step_h=1.0)
    result = engine.run_trial(protocol, admet)
    
    print(f"    Status: {result['constitutional']['status']}")
    print(f"    Engine: {result['constitutional']['engine']}")
    
    # Step 4: Display results
    print("\n[4] Trial Results - Exposure Summary")
    print("=" * 80)
    
    print(f"\n{'Arm':<25} {'Dose':<15} {'Cmax (mg/L)':<20} {'AUC (mg·h/L)':<20}")
    print("-" * 80)
    
    for arm in result["arms"]:
        label = arm["label"]
        dose = f"{arm['dose_mg']} mg"
        if arm['dosing_interval_h'] == 12.0:
            dose += " BID"
        else:
            dose += " QD"
        
        cmax = arm["exposure_summary"]["cmax_mg_per_L"]
        auc = arm["exposure_summary"]["auc_mg_h_per_L"]
        
        cmax_str = f"{cmax['mean']:.4f} ± {(cmax['max']-cmax['min'])/2:.4f}"
        auc_str = f"{auc['mean']:.2f} ± {(auc['max']-auc['min'])/2:.2f}"
        
        print(f"{label:<25} {dose:<15} {cmax_str:<20} {auc_str:<20}")
    
    # Step 5: Detailed statistics for one arm
    print("\n[5] Detailed Statistics - Arm A2 (Standard Dose)")
    print("=" * 80)
    
    arm = result["arms"][1]  # A2
    
    metrics = [
        ("Cmax (mg/L)", arm["exposure_summary"]["cmax_mg_per_L"]),
        ("AUC (mg·h/L)", arm["exposure_summary"]["auc_mg_h_per_L"]),
        ("Cmin SS (mg/L)", arm["exposure_summary"]["cmin_steady_state_mg_per_L"]),
    ]
    
    for metric_name, stats in metrics:
        print(f"\n{metric_name}:")
        print(f"  Mean:   {stats['mean']:.4f}")
        print(f"  Median: {stats['median']:.4f}")
        print(f"  Min:    {stats['min']:.4f}")
        print(f"  Max:    {stats['max']:.4f}")
        print(f"  Range:  {stats['max'] - stats['min']:.4f}")
    
    # Step 6: Comparison between arms
    print("\n[6] Dose Comparison Analysis")
    print("=" * 80)
    
    a1_auc = result["arms"][0]["exposure_summary"]["auc_mg_h_per_L"]["mean"]
    a2_auc = result["arms"][1]["exposure_summary"]["auc_mg_h_per_L"]["mean"]
    a3_auc = result["arms"][2]["exposure_summary"]["auc_mg_h_per_L"]["mean"]
    
    print(f"\nAUC Comparison (Mean):")
    print(f"  A1 (Low Dose QD):      {a1_auc:.2f} mg·h/L")
    print(f"  A2 (Standard Dose QD): {a2_auc:.2f} mg·h/L  (+{((a2_auc/a1_auc-1)*100):.1f}%)")
    print(f"  A3 (High Dose BID):    {a3_auc:.2f} mg·h/L  (+{((a3_auc/a1_auc-1)*100):.1f}%)")
    
    print("\n" + "=" * 80)
    print("✅ TRIAL ENGINE DEMONSTRATION COMPLETE")
    print("=" * 80)
    
    print("\n🎯 CAPABILITIES DEMONSTRATED:")
    print("   ✅ Virtual population generation (20 patients/arm)")
    print("   ✅ Multi-arm trial simulation (3 arms)")
    print("   ✅ Exposure metrics (Cmax, AUC, Cmin)")
    print("   ✅ Statistical summaries (mean, median, min, max)")
    print("   ✅ Dose comparison analysis")
    print("   ✅ Constitutional compliance (L51/L34)")
    
    print("\n📊 TRIAL CAPABILITIES:")
    print("   • Parallel arm designs")
    print("   • Multiple dosing regimens")
    print("   • Virtual patient populations")
    print("   • Exposure-based endpoints")
    print("   • Deterministic simulations")
    
    print("\n🔬 NEXT STEPS:")
    print("   • Add PK/PD models (efficacy endpoints)")
    print("   • Crossover designs")
    print("   • Adaptive dosing")
    print("   • Inter-individual variability")
    print("=" * 80)


if __name__ == "__main__":
    demo_trial_engine()
