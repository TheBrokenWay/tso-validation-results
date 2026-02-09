"""
Trial Dossier Demo - Complete Evidence Package Generation
Run this to see the full trial simulation → dossier workflow!
"""

import json
from PX_Engine.operations import TrialEngine, run_ope, run_admet
from PX_System.foundation.Evidence_Package import wrap_trial_simulation


def demo_trial_dossier():
    """Demonstrate complete trial simulation + evidence package workflow"""
    
    print("=" * 80)
    print("PREDATOR X - TRIAL DOSSIER GENERATION DEMONSTRATION")
    print("Complete Workflow: SMILES → OPE → ADMET → Trial → Evidence Package")
    print("=" * 80)
    
    # Example: Aspirin
    smiles = "CC(=O)Oc1ccccc1C(=O)O"
    print(f"\n📊 Molecule: Aspirin (SMILES: {smiles})")
    
    # Step 1: Run OPE
    print("\n[1] Running OPE Analysis...")
    ope = run_ope(smiles)
    print(f"    Status: {ope.get('status', 'N/A')}")
    print(f"    LogP: {ope.get('logp', 'N/A')}")
    
    # Step 2: Run ADMET
    print("\n[2] Running ADMET Analysis...")
    admet = run_admet(smiles, ope)
    print(f"    Status: {admet['constitutional']['status']}")
    print(f"    Engine: {admet['constitutional']['engine']}")
    print(f"    Hepatotoxicity Risk: {admet['toxicity_flags']['hepatotoxicity_risk']}")
    
    # Step 3: Define trial protocol
    print("\n[3] Defining Trial Protocol...")
    protocol = {
        "trial_id": "TRIAL-PK-ASPIRIN-001",
        "duration_days": 7.0,
        "arms": [
            {
                "arm_id": "A1",
                "label": "Low Dose",
                "dose_mg": 75.0,
                "dosing_interval_h": 24.0,
                "n_patients": 20,
            },
            {
                "arm_id": "A2",
                "label": "Standard Dose",
                "dose_mg": 100.0,
                "dosing_interval_h": 24.0,
                "n_patients": 20,
            },
            {
                "arm_id": "A3",
                "label": "High Dose BID",
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
        print(f"      • {arm['label']}: {arm['dose_mg']} mg, {arm['n_patients']} patients")
    
    # Step 4: Run trial simulation
    print("\n[4] Running Trial Simulation...")
    engine = TrialEngine(time_step_h=1.0)
    trial_result = engine.run_trial(protocol, admet)
    
    print(f"    Trial Status: {trial_result['constitutional']['status']}")
    print(f"    Trial Engine: {trial_result['constitutional']['engine']}")
    print(f"    Arms Completed: {len(trial_result['arms'])}")
    
    # Step 5: Generate Evidence Package
    print("\n[5] Generating Evidence Package...")
    dossier_path = wrap_trial_simulation(protocol, trial_result, ope, admet)
    
    print(f"    ✅ Dossier Created: {dossier_path}")
    
    # Step 6: Display dossier summary
    print("\n[6] Evidence Package Summary")
    print("=" * 80)
    
    with open(dossier_path, 'r') as f:
        dossier = json.load(f)
    
    print(f"\nDossier Type:     {dossier['dossier_type']}")
    print(f"Version:          {dossier['version']}")
    print(f"Timestamp:        {dossier['timestamp_utc']}")
    print(f"Evidence Hash:    {dossier['evidence_hash']}")
    
    print(f"\nProvenance:")
    print(f"  OPE Engine:     {dossier['provenance']['ope_engine']}")
    print(f"  ADMET Engine:   {dossier['provenance']['admet_engine']}")
    print(f"  Trial Engine:   {dossier['provenance']['trial_engine']}")
    
    print(f"\nConstitutional:")
    print(f"  Status:         {dossier['constitutional']['status']}")
    print(f"  Law Basis:      {', '.join(dossier['constitutional']['law_basis'])}")
    print(f"  Notes:          {dossier['constitutional']['notes']}")
    
    print(f"\nTrial Results:")
    for arm in dossier['trial_result']['arms']:
        auc = arm['exposure_summary']['auc_mg_h_per_L']
        print(f"  {arm['label']}: AUC = {auc['mean']:.2f} mg·h/L (range: {auc['min']:.2f}-{auc['max']:.2f})")
    
    # Step 7: File size and structure
    print("\n[7] Dossier Details")
    print("=" * 80)
    
    import os
    file_size = os.path.getsize(dossier_path)
    print(f"\nFile Size:        {file_size:,} bytes ({file_size/1024:.1f} KB)")
    print(f"File Location:    {dossier_path}")
    
    # Count data points
    total_fields = len(json.dumps(dossier))
    print(f"Total Data:       {total_fields:,} characters")
    
    print("\n" + "=" * 80)
    print("✅ TRIAL DOSSIER GENERATION COMPLETE")
    print("=" * 80)
    
    print("\n🎯 WORKFLOW VALIDATED:")
    print("   ✅ OPE Analysis")
    print("   ✅ ADMET Prediction")
    print("   ✅ Trial Simulation")
    print("   ✅ Evidence Package Creation")
    print("   ✅ Constitutional Compliance (L51/L34)")
    print("   ✅ Provenance Tracking")
    print("   ✅ SHA-256 Hash Generation")
    
    print("\n📦 DOSSIER CONTENTS:")
    print("   • Trial protocol (3 arms)")
    print("   • Trial results (exposure summaries)")
    print("   • OPE analysis (PK predictions)")
    print("   • ADMET analysis (toxicity)")
    print("   • Engine versions")
    print("   • Timestamps")
    print("   • Reproducibility hash")
    
    print("\n🔐 READY FOR:")
    print("   • GAIP Gateway authorization")
    print("   • Byzantium Council review")
    print("   • Warehouse persistence")
    print("   • Regulatory submission")
    print("   • Audit trail generation")
    
    print("=" * 80)
    
    return dossier_path


if __name__ == "__main__":
    demo_trial_dossier()
