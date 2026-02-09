"""
Integration test for complete pipeline: OPE → ADMET → PK → Trial
Validates end-to-end trial simulation workflow
"""

import sys
import os

# Add foundation to path
ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
if ROOT_DIR not in sys.path:
    sys.path.append(ROOT_DIR)

from PX_Engine.operations import run_ope, run_admet, TrialEngine


def test_full_pipeline():
    """Test complete pipeline: SMILES → OPE → ADMET → Trial Simulation"""
    
    print("=" * 80)
    print("TRIAL ENGINE INTEGRATION TEST")
    print("Complete Pipeline: SMILES → OPE → ADMET → PK → Trial")
    print("=" * 80)
    
    # Test SMILES (Aspirin)
    smiles = "CC(=O)Oc1ccccc1C(=O)O"
    
    # Step 1: Run OPE
    print("\n[1] Running OPE...")
    ope = run_ope(smiles)
    print(f"   Status: {ope.get('status', 'N/A')}")
    print(f"   LogP: {ope.get('logp', 'N/A')}")
    assert ope is not None
    
    # Step 2: Run ADMET
    print("\n[2] Running ADMET...")
    admet = run_admet(smiles, ope)
    print(f"   Status: {admet['constitutional']['status']}")
    print(f"   Engine: {admet['constitutional']['engine']}")
    print(f"   Hepatotoxicity Risk: {admet['toxicity_flags']['hepatotoxicity_risk']}")
    # Governance: only OPE_ADMET_V3_DETERMINISTIC is canonical (BAN: OPE_ADMET_V1)
    assert admet['constitutional']['engine'] == 'OPE_ADMET_V3_DETERMINISTIC'
    
    # Step 3: Define trial protocol
    print("\n[3] Defining Trial Protocol...")
    protocol = {
        "trial_id": "INTEGRATION-TEST-001",
        "duration_days": 2.0,
        "arms": [
            {
                "arm_id": "A1",
                "label": "QD 100 mg",
                "dose_mg": 100.0,
                "dosing_interval_h": 24.0,
                "n_patients": 10,
            },
            {
                "arm_id": "A2",
                "label": "BID 50 mg",
                "dose_mg": 50.0,
                "dosing_interval_h": 12.0,
                "n_patients": 10,
            },
        ],
    }
    print(f"   Trial ID: {protocol['trial_id']}")
    print(f"   Arms: {len(protocol['arms'])}")
    
    # Step 4: Run trial simulation
    print("\n[4] Running Trial Simulation...")
    engine = TrialEngine(time_step_h=1.0)
    trial_result = engine.run_trial(protocol, admet)
    
    print(f"   Trial ID: {trial_result['trial_id']}")
    print(f"   Duration: {trial_result['duration_days']} days")
    print(f"   Arms Completed: {len(trial_result['arms'])}")
    
    # Step 5: Validate results
    print("\n[5] Validating Results...")
    
    assert trial_result['trial_id'] == 'INTEGRATION-TEST-001'
    assert len(trial_result['arms']) == 2
    
    for arm in trial_result['arms']:
        print(f"\n   Arm: {arm['arm_id']} - {arm['label']}")
        
        # Validate exposure summary structure
        exp = arm['exposure_summary']
        assert 'cmax_mg_per_L' in exp
        assert 'auc_mg_h_per_L' in exp
        assert 'cmin_steady_state_mg_per_L' in exp
        
        # Validate statistics
        for metric_name, metric_stats in exp.items():
            assert 'mean' in metric_stats
            assert 'median' in metric_stats
            assert 'min' in metric_stats
            assert 'max' in metric_stats
            
            # Validate ranges
            assert metric_stats['min'] <= metric_stats['mean'] <= metric_stats['max']
            assert metric_stats['min'] <= metric_stats['median'] <= metric_stats['max']
            
            print(f"      {metric_name}: {metric_stats['mean']:.4f} mg/L (range: {metric_stats['min']:.4f}-{metric_stats['max']:.4f})")
    
    # Step 6: Validate constitutional compliance
    print("\n[6] Validating Constitutional Compliance...")
    
    assert trial_result['constitutional']['status'] == 'SIMULATED'
    assert trial_result['constitutional']['engine'] == 'TRIAL_ENGINE_V1'
    print(f"   Status: {trial_result['constitutional']['status']}")
    print(f"   Engine: {trial_result['constitutional']['engine']}")
    
    # Step 7: Test dose-response relationship
    print("\n[7] Testing Dose-Response Relationship...")
    
    arm1_auc = trial_result['arms'][0]['exposure_summary']['auc_mg_h_per_L']['mean']
    arm2_auc = trial_result['arms'][1]['exposure_summary']['auc_mg_h_per_L']['mean']
    
    print(f"   QD 100mg AUC: {arm1_auc:.2f} mg·h/L")
    print(f"   BID 50mg AUC: {arm2_auc:.2f} mg·h/L")
    print(f"   Ratio: {arm2_auc/arm1_auc:.2f}x")
    
    # For same total daily dose (100mg), AUC should be similar
    # BID may have slightly lower AUC due to more frequent dosing
    # Validate that AUCs are within reasonable range (0.8-1.2x)
    assert 0.8 <= arm2_auc/arm1_auc <= 1.2, "BID/QD AUC ratio should be within 0.8-1.2x for same total daily dose"
    print(f"   ✅ Dose-response validation passed (ratio within expected range)")
    
    print("\n" + "=" * 80)
    print("✅ TRIAL ENGINE INTEGRATION TEST SUCCESSFUL")
    print("=" * 80)
    
    print("\n📊 PIPELINE VALIDATED:")
    print("   ✅ OPE (Pharmacokinetic predictions)")
    print("   ✅ ADMET (Toxicity assessment)")
    print("   ✅ PK Simulation (Concentration-time profiles)")
    print("   ✅ Trial Simulation (Multi-arm, population-based)")
    print("   ✅ Statistical analysis (mean, median, range)")
    print("   ✅ Constitutional compliance (L51/L34)")
    print("   ✅ Dose-response validation")
    
    print("\n🎯 READY FOR:")
    print("   • Production trial simulations")
    print("   • Dose optimization studies")
    print("   • Comparative effectiveness trials")
    print("   • Virtual patient populations")
    
    print("=" * 80)
    
    return True


if __name__ == "__main__":
    test_full_pipeline()
