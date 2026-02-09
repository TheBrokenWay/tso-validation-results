"""
End-to-end layer check: governance, Nipah, PX_Engine, PX_System, Evidence_Package.
Run from repo root: python run_e2e_layers.py
"""
import os
import sys

ROOT = os.path.dirname(os.path.abspath(__file__))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

def main():
    errors = []

    # 1. Governance + Constitutional
    try:
        from governance.poison_pill_gate import run_poison_pill_gate, require_poison_pill_gate_before_pipeline
        if not run_poison_pill_gate(timeout_seconds=120):
            errors.append("1. Poison pill gate failed")
        else:
            require_poison_pill_gate_before_pipeline()
            print("1. Governance + Constitutional OK")
    except Exception as e:
        errors.append(f"1. Governance: {e}")

    # 2. Nipah adapter (valid file)
    try:
        from pathlib import Path
        from Nipah_Analysis.adapters.nipah_miner_adapter import NipahMinerAdapter
        valid = Path("Nipah_Analysis/data/raw/valid_malaysia.json")
        if valid.exists():
            adapter = NipahMinerAdapter()
            adapter.mine_candidate(str(valid))
            print("2. Nipah adapter OK")
        else:
            print("2. Nipah valid_malaysia.json missing, skip")
    except Exception as e:
        errors.append(f"2. Nipah: {e}")

    # 3. PX_Engine OPE/ADMET + ZeusLaws
    try:
        from PX_Engine.operations import run_ope, run_admet
        from PX_System.foundation.ZeusLaws import check_constitutional
        ope = run_ope("CCO")
        admet = run_admet("CCO", ope)
        tox = admet.get("toxicity_index") or (admet.get("toxicity") or {}).get("toxicity_index", 0)
        verdict = check_constitutional("e2e", {"toxicity_index": tox, "harm_energy": tox})
        assert "authorized" in verdict
        print("3. OPE/ADMET/ZeusLaws OK")
    except Exception as e:
        errors.append(f"3. PX_Engine: {e}")

    # 4. PX_System foundation + VectorCore + TrialEngine
    try:
        from PX_System.foundation.Disease_Constraint_Model import create_chagas_dcm
        from PX_Engine.Vector_Core import VectorCore
        from PX_Engine.operations.TrialEngine import TrialEngine
        dcm = create_chagas_dcm()
        assert "disease_name" in dcm
        core = VectorCore()
        v = core.execute([0.1, 0.0, 35.0, 1.0])
        assert "authorized" in v
        engine = TrialEngine(time_step_h=1.0)
        protocol = {"trial_id": "E2E", "duration_days": 1.0, "arms": [{"arm_id": "A1", "dose_mg": 100, "dosing_interval_h": 24, "n_patients": 2}]}
        ope = run_ope("CCO")
        admet = run_admet("CCO", ope)
        trial = engine.run_trial(protocol, admet)
        assert isinstance(trial, dict)
        print("4. PX_System + VectorCore + TrialEngine OK")
    except Exception as e:
        errors.append(f"4. PX_System/Engine: {e}")

    # 5. Evidence_Package generate_dossier (ADMET path)
    try:
        from PX_System.foundation.Evidence_Package import generate_dossier
        from PX_Constitution.Virtual_Machine import get_vm_fingerprint
        get_vm_fingerprint()
        candidate = {"name": "E2E", "smiles": "CCO"}
        engine_results = {"ope": {}, "admet": {"toxicity_index": 0.01}}
        dossier = generate_dossier(candidate, engine_results)
        assert "harm_energy" in dossier and dossier["harm_energy"] == 0.01
        assert "dossier_version" in dossier
        print("5. Evidence_Package OK")
    except Exception as e:
        errors.append(f"5. Evidence_Package: {e}")

    # 6. wrap_trial_simulation (orchestrator path)
    try:
        from PX_System.foundation.Evidence_Package import wrap_trial_simulation
        from PX_Engine.operations import run_ope, run_admet
        from PX_Engine.operations.TrialEngine import TrialEngine
        ope = run_ope("CCO")
        admet = run_admet("CCO", ope)
        engine = TrialEngine(time_step_h=1.0)
        protocol = {"trial_id": "E2E", "duration_days": 1.0, "arms": [{"arm_id": "A1", "dose_mg": 100, "dosing_interval_h": 24, "n_patients": 2}]}
        trial = engine.run_trial(protocol, admet)
        path = wrap_trial_simulation(protocol, trial, ope, admet, output_dir=os.path.join(ROOT, "PX_Warehouse", "Calibration_Molecules"))
        assert path and os.path.exists(path)
        print("6. wrap_trial_simulation OK")
    except Exception as e:
        errors.append(f"6. wrap_trial_simulation: {e}")

    # 7. TSO_Validator standalone safety gate
    try:
        tso_script = os.path.join(ROOT, "TSO_Validator", "run_validation.py")
        if os.path.exists(tso_script):
            import subprocess
            r = subprocess.run(
                [sys.executable, tso_script],
                capture_output=True, text=True, timeout=60
            )
            if r.returncode != 0:
                errors.append(f"7. TSO_Validator: exit {r.returncode}: {r.stderr or r.stdout}")
            else:
                print("7. TSO_Validator OK")
        else:
            errors.append("7. TSO_Validator: run_validation.py not found")
    except Exception as e:
        errors.append(f"7. TSO_Validator: {e}")

    if errors:
        print("\nFAILED:", errors)
        sys.exit(1)
    print("\nAll layers OK.")
    return 0

if __name__ == "__main__":
    sys.exit(main() or 0)
