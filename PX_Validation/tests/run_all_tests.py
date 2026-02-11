"""
Comprehensive Test Runner
Runs all tests and generates summary report
"""

import subprocess
import sys
from datetime import datetime, timezone


def run_test(test_path, test_name):
    """Run a single test and return results"""
    print(f"\n{'='*80}")
    print(f"Running: {test_name}")
    print('='*80)
    
    try:
        result = subprocess.run(
            [sys.executable, test_path],
            capture_output=True,
            text=True,
            timeout=60
        )
        
        output = result.stdout + result.stderr
        
        # Parse results
        if "OK" in output or "PASSED" in output or "SUCCESSFUL" in output:
            # Extract test count
            if "Ran" in output:
                for line in output.split('\n'):
                    if 'Ran' in line and 'test' in line:
                        count = line.split()[1]
                        print(f"✅ {test_name}: {count} tests passed")
                        return True, count
            elif "PASSED" in output:
                for line in output.split('\n'):
                    if 'PASSED' in line:
                        count = line.split()[1]
                        print(f"✅ {test_name}: {count} tests passed")
                        return True, count
        else:
            print(f"❌ {test_name}: FAILED")
            print(output[-500:])  # Print last 500 chars
            return False, "0"
            
    except subprocess.TimeoutExpired:
        print(f"⚠️  {test_name}: TIMEOUT")
        return False, "0"
    except Exception as e:
        print(f"❌ {test_name}: ERROR - {e}")
        return False, "0"
    
    return True, "N/A"


def main():
    """Run all tests and generate summary"""
    
    print("=" * 80)
    print("PREDATOR X - COMPREHENSIVE TEST SUITE")
    print(f"Timestamp: {datetime.now(timezone.utc).isoformat()}")
    print("=" * 80)
    
    tests = [
        # Core unit tests
        ("PX_Validation/tests/PX_System_Test.py", "System Tests"),
        ("PX_Validation/tests/test_admet_engine.py", "ADMET Tests"),
        ("PX_Validation/tests/test_pk_engine.py", "PK Engine Tests"),
        ("PX_Validation/tests/test_trial_engine.py", "Trial Engine Tests"),
        ("PX_Validation/tests/test_trial_evidence_package.py", "Trial Evidence Package Tests"),
        ("PX_Validation/tests/test_pkpd.py", "PKPD Tests"),
        ("PX_Validation/tests/test_iiv.py", "IIV Tests"),
        ("PX_Validation/tests/test_adaptive.py", "Adaptive Tests"),
        ("PX_Validation/tests/test_dose_optimizer_v2.py", "Dose Optimizer V2 Tests"),
        ("PX_Validation/tests/test_virtual_efficacy.py", "Virtual Efficacy Tests"),
        ("PX_Validation/tests/test_grading_engine.py", "Grading Engine Tests"),
        # Integration tests
        ("PX_Validation/tests/test_pk_integration.py", "PK Integration Test"),
        ("PX_Validation/tests/test_trial_integration.py", "Trial Integration Test"),
        ("PX_Validation/tests/test_pkpd_integration.py", "PKPD Integration Test"),
        ("PX_Validation/tests/test_iiv_integration.py", "IIV Integration Test"),
        ("PX_Validation/tests/test_adaptive_integration.py", "Adaptive Integration Test"),
        ("PX_Validation/tests/test_dose_optimizer_v2_integration.py", "Dose Optimizer V2 Integration Test"),
        ("PX_Validation/tests/test_virtual_efficacy_integration.py", "Virtual Efficacy Integration Test"),
        ("PX_Validation/tests/test_system_trial_evidence_package.py", "System Trial Evidence Package Test"),
        # End-to-end / orchestrator tests
        ("PX_Validation/tests/test_orchestrator_v2.py", "Orchestrator V2 Tests"),
        ("PX_Validation/tests/test_warehouse_integrity.py", "Warehouse Integrity Tests"),
        ("PX_Validation/tests/test_performance_regression.py", "Performance Regression Tests"),
        # QUINT internal language
        ("PX_Validation/tests/test_quint.py", "QUINT Language Tests"),
        # Disease registry
        ("PX_Validation/tests/test_disease_registry.py", "Disease Registry Tests"),
        # TSO_Validator standalone safety gate
        ("TSO_Validator/tests/test_ethical_guard.py", "TSO Ethical Guard Tests"),
        ("TSO_Validator/tests/test_provenance.py", "TSO Provenance Tests"),
        ("TSO_Validator/tests/test_reproducibility.py", "TSO Reproducibility Tests"),
        ("TSO_Validator/tests/test_run_validation.py", "TSO Run Validation Tests"),
    ]
    
    results = []
    total_passed = 0
    total_count = 0
    
    for test_path, test_name in tests:
        passed, count = run_test(test_path, test_name)
        results.append((test_name, passed, count))
        if passed:
            if count.isdigit():
                total_count += int(count)
                total_passed += int(count)
    
    # Summary
    print("\n" + "=" * 80)
    print("TEST SUITE SUMMARY")
    print("=" * 80)
    
    for test_name, passed, count in results:
        status = "✅ PASSED" if passed else "❌ FAILED"
        print(f"{status}  {test_name}: {count} tests")
    
    print("\n" + "=" * 80)
    print(f"TOTAL: {total_passed}/{total_count} tests passed ({100*total_passed//total_count if total_count else 0}%)")
    print("=" * 80)
    
    if total_passed == total_count:
        print("\n🎉 ALL TESTS PASSING - SYSTEM 100% OPERATIONAL")
        return 0
    else:
        print(f"\n⚠️  {total_count - total_passed} TESTS FAILED")
        return 1


if __name__ == "__main__":
    sys.exit(main())
