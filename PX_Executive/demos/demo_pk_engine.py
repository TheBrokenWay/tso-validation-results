"""
PK Engine Demo - Quick demonstration of the new PK simulation capability
Run this to see the PK engine in action!
"""

from PX_Laboratory import SimulationEngine
from PX_Engine.operations import run_ope, run_admet


def demo_pk_engine():
    """Demonstrate PK engine capabilities"""
    
    print("=" * 80)
    print("PREDATOR X - PK ENGINE DEMONSTRATION")
    print("=" * 80)
    
    # Example: Aspirin
    smiles = "CC(=O)Oc1ccccc1C(=O)O"
    print(f"\n📊 Molecule: Aspirin (SMILES: {smiles})")
    
    # Step 1: Run ADMET
    print("\n[1] Running ADMET Analysis...")
    ope = run_ope(smiles)
    admet = run_admet(smiles, ope)
    print(f"    Hepatotoxicity Risk: {admet['toxicity_flags']['hepatotoxicity_risk']}")
    
    # Step 2: Create PK simulation engine
    print("\n[2] Creating PK Simulation Engine...")
    engine = SimulationEngine(time_step_h=1.0)
    print(f"    Engine: {engine.metadata}")
    
    # Step 3: Simulate single dose (QD - Once Daily)
    print("\n[3] Simulating Single Dose (100 mg, Once Daily)...")
    patient = {"weight_kg": 70.0}
    
    result_qd = engine.simulate_one_compartment(
        dose_mg=100.0,
        duration_h=24.0,
        dosing_interval_h=24.0,
        patient=patient,
        admet=admet,
    )
    
    print(f"    Model: {result_qd['model']}")
    print(f"    Cmax:  {result_qd['summary']['cmax_mg_per_L']:.4f} mg/L")
    print(f"    Tmax:  {result_qd['summary']['tmax_h']:.1f} hours")
    print(f"    AUC:   {result_qd['summary']['auc_mg_h_per_L']:.2f} mg·h/L")
    print(f"    Cmin:  {result_qd['summary']['cmin_steady_state_mg_per_L']:.4f} mg/L")
    
    # Step 4: Simulate BID regimen (Twice Daily)
    print("\n[4] Simulating BID Regimen (50 mg, Twice Daily, 48h)...")
    
    result_bid = engine.simulate_one_compartment(
        dose_mg=50.0,
        duration_h=48.0,
        dosing_interval_h=12.0,
        patient=patient,
        admet=admet,
    )
    
    print(f"    Cmax:  {result_bid['summary']['cmax_mg_per_L']:.4f} mg/L")
    print(f"    AUC:   {result_bid['summary']['auc_mg_h_per_L']:.2f} mg·h/L")
    print(f"    Cmin:  {result_bid['summary']['cmin_steady_state_mg_per_L']:.4f} mg/L")
    
    # Step 5: Plot concentration-time profile (simple text visualization)
    print("\n[5] Concentration-Time Profile (First 12 hours, BID):")
    print("    Time (h) | Concentration (mg/L)")
    print("    " + "-" * 40)
    
    for i in range(0, min(13, len(result_bid['time_grid_h']))):
        t = result_bid['time_grid_h'][i]
        c = result_bid['concentration_mg_per_L'][i]
        bar = "█" * int(c * 20)
        print(f"    {t:6.1f}   | {c:6.4f} {bar}")
    
    # Step 6: Summary
    print("\n" + "=" * 80)
    print("✅ PK ENGINE DEMONSTRATION COMPLETE")
    print("=" * 80)
    print("\n🎯 CAPABILITIES DEMONSTRATED:")
    print("   ✅ ADMET integration")
    print("   ✅ PK simulation (one-compartment)")
    print("   ✅ Multiple dosing regimens (QD, BID)")
    print("   ✅ PK metrics (Cmax, Tmax, AUC, Cmin)")
    print("   ✅ Concentration-time profiles")
    print("\n📊 NEXT STEPS:")
    print("   • Virtual patient populations")
    print("   • Dose optimization")
    print("   • Clinical trial simulations")
    print("=" * 80)


if __name__ == "__main__":
    demo_pk_engine()
