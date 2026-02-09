"""
PREDATOR X - COMPREHENSIVE SYSTEM TEST
Tests all PX_ directories and components to identify adjustments needed
"""

import sys
import os
import importlib
import traceback
from pathlib import Path
from datetime import datetime, timezone

# Repo root (no hardcoded drive)
ROOT_DIR = str(Path(__file__).resolve().parents[2])
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)

class SystemTester:
    def __init__(self):
        self.results = {
            "passed": [],
            "failed": [],
            "warnings": []
        }
        self.timestamp = datetime.now(timezone.utc).isoformat()
        
    def log(self, level, module, message, error=None):
        entry = {
            "module": module,
            "message": message,
            "timestamp": self.timestamp
        }
        if error:
            entry["error"] = str(error)
            entry["traceback"] = traceback.format_exc()
        
        self.results[level].append(entry)
        
    def test_import(self, module_name, description):
        """Test if a module can be imported"""
        print(f"\n  Testing: {description}")
        print(f"  Module: {module_name}")
        try:
            module = importlib.import_module(module_name)
            print(f"  ✅ SUCCESS")
            self.log("passed", module_name, description)
            return module
        except Exception as e:
            print(f"  ❌ FAILED: {e}")
            self.log("failed", module_name, description, e)
            return None
    
    def test_class_instantiation(self, module, class_name, description, *args, **kwargs):
        """Test if a class can be instantiated"""
        print(f"\n  Testing: {description}")
        print(f"  Class: {class_name}")
        try:
            if module is None:
                raise Exception(f"Module not available")
            cls = getattr(module, class_name)
            instance = cls(*args, **kwargs)
            print(f"  ✅ SUCCESS")
            self.log("passed", f"{module.__name__}.{class_name}", description)
            return instance
        except Exception as e:
            print(f"  ❌ FAILED: {e}")
            self.log("failed", f"{class_name}", description, e)
            return None
    
    def test_file_exists(self, filepath, description):
        """Test if a file exists"""
        print(f"\n  Testing: {description}")
        print(f"  Path: {filepath}")
        if os.path.exists(filepath):
            print(f"  ✅ EXISTS")
            self.log("passed", filepath, description)
            return True
        else:
            print(f"  ❌ NOT FOUND")
            self.log("failed", filepath, description)
            return False
    
    def test_directory_structure(self, base_dir, expected_subdirs, description):
        """Test if expected subdirectories exist"""
        print(f"\n  Testing: {description}")
        print(f"  Base: {base_dir}")
        missing = []
        for subdir in expected_subdirs:
            path = os.path.join(base_dir, subdir)
            if not os.path.exists(path):
                missing.append(subdir)
        
        if not missing:
            print(f"  ✅ ALL SUBDIRECTORIES PRESENT")
            self.log("passed", base_dir, description)
            return True
        else:
            print(f"  ⚠️  MISSING: {', '.join(missing)}")
            self.log("warnings", base_dir, f"Missing subdirectories: {', '.join(missing)}")
            return False
    
    def print_summary(self):
        """Print test summary"""
        print("\n" + "=" * 80)
        print("TEST SUMMARY")
        print("=" * 80)
        print(f"\n✅ PASSED: {len(self.results['passed'])}")
        print(f"❌ FAILED: {len(self.results['failed'])}")
        print(f"⚠️  WARNINGS: {len(self.results['warnings'])}")
        
        if self.results['failed']:
            print("\n" + "=" * 80)
            print("FAILED TESTS - ADJUSTMENTS NEEDED")
            print("=" * 80)
            for i, failure in enumerate(self.results['failed'], 1):
                print(f"\n{i}. MODULE: {failure['module']}")
                print(f"   TEST: {failure['message']}")
                print(f"   ERROR: {failure['error']}")
                if 'traceback' in failure:
                    print(f"   TRACEBACK:\n{failure['traceback']}")
        
        if self.results['warnings']:
            print("\n" + "=" * 80)
            print("WARNINGS - REVIEW RECOMMENDED")
            print("=" * 80)
            for i, warning in enumerate(self.results['warnings'], 1):
                print(f"\n{i}. MODULE: {warning['module']}")
                print(f"   ISSUE: {warning['message']}")
        
        return len(self.results['failed']) == 0


def run_comprehensive_test():
    """Run comprehensive system test on all PX_ directories"""
    tester = SystemTester()
    
    print("=" * 80)
    print("PREDATOR X - COMPREHENSIVE SYSTEM TEST")
    print("Testing all PX_ directories and components")
    print("=" * 80)
    
    # =========================================================================
    # TEST 1: PX_SYSTEM (Foundation Package)
    # =========================================================================
    print("\n\n[TEST 1] PX_SYSTEM - FOUNDATION PACKAGE")
    print("-" * 80)
    
    # Test foundation core
    foundation_core = tester.test_import("PX_System.foundation.core", "Foundation core module")
    
    # Test ZeusLaws (known issue with old import)
    zeus_laws = tester.test_import("PX_System.foundation.ZeusLaws", "ZeusLaws module")
    
    # Test Sovereign Log Chain
    log_chain = tester.test_import("PX_System.foundation.Sovereign_Log_Chain", "Sovereign Log Chain")
    
    # Test operational engines (moved to PX_Engine/operations)
    print("\n  [PX_Engine Operational Engines]")
    obe = tester.test_import("PX_Engine.operations.OBE", "Operational Blocker Engine")
    oce = tester.test_import("PX_Engine.operations.OCE", "Operational Coherence Engine")
    ole = tester.test_import("PX_Engine.operations.OLE", "Operational Logic Engine")
    ome = tester.test_import("PX_Engine.operations.OME", "Operational Momentum Engine")
    ope = tester.test_import("PX_Engine.operations.OPE", "Operational Physics Engine")
    ose = tester.test_import("PX_Engine.operations.OSE", "Operational Status Engine")
    admet = tester.test_import("PX_Engine.operations.ADMET", "ADMET Engine")
    trial_engine = tester.test_import("PX_Engine.operations.TrialEngine", "Trial Engine")
    
    # Test discovery (moved to PX_Discovery)
    discovery = tester.test_import("PX_Discovery.candidate_discovery_engine", 
                                   "Candidate Discovery Engine")
    
    # Test integrations
    net_policy = tester.test_import("PX_System.foundation.integrations.net_policy", "Network Policy")
    smiles_sec = tester.test_import("PX_System.foundation.integrations.smiles_security", "SMILES Security")
    
    # =========================================================================
    # TEST 2: PX_EXECUTIVE (Governance & Pipelines)
    # =========================================================================
    print("\n\n[TEST 2] PX_EXECUTIVE - GOVERNANCE & PIPELINES")
    print("-" * 80)
    
    gaip_gateway_mod = tester.test_import("PX_Executive.GAIP_Gateway", "GAIP Gateway module")
    gaip_gateway = tester.test_class_instantiation(gaip_gateway_mod, "GAIPGateway", 
                                                   "GAIP Gateway instantiation", mode="RESEARCH")
    
    byzantium_mod = tester.test_import("PX_Executive.Byzantium_Council", "Byzantium Council module")
    byzantium = tester.test_class_instantiation(byzantium_mod, "ByzantiumCouncil", 
                                               "Byzantium Council instantiation")
    
    gold_rush_mod = tester.test_import("PX_Executive.Gold_Rush_Miner", "Gold Rush Miner module")
    
    legal_check_mod = tester.test_import("PX_Executive.PX_Legal_Check", "Legal Check module")
    
    commercial_pipeline_mod = tester.test_import("PX_Executive.Sovereign_Commercial_Pipeline", 
                                                 "Sovereign Commercial Pipeline module")
    
    # =========================================================================
    # TEST 3: PX_ENGINE (Vector Physics & Metabolism)
    # =========================================================================
    print("\n\n[TEST 3] PX_ENGINE - VECTOR PHYSICS & METABOLISM")
    print("-" * 80)
    
    vector_core_mod = tester.test_import("PX_Engine.Vector_Core", "Vector Core module")
    
    metabolism_mod = tester.test_import("PX_Engine.Metabolism", "Metabolism module")
    
    trajectory_mod = tester.test_import("PX_Engine.Trajectory_Predictor", "Trajectory Predictor module")
    
    # =========================================================================
    # TEST 4: PX_LABORATORY (Materialization & Synthesis)
    # =========================================================================
    print("\n\n[TEST 4] PX_LABORATORY - MATERIALIZATION & SYNTHESIS")
    print("-" * 80)
    
    sim_engine_mod = tester.test_import("PX_Laboratory.Simulation_Engine", "Simulation Engine module")
    sim_engine = tester.test_class_instantiation(sim_engine_mod, "SimulationEngine", 
                                                "Simulation Engine instantiation")
    
    manufacturing_mod = tester.test_import("PX_Laboratory.Manufacturing_Manifest", 
                                          "Manufacturing Manifest module")
    
    # =========================================================================
    # TEST 5: PX_WAREHOUSE (Data Persistence & WorldLines)
    # =========================================================================
    print("\n\n[TEST 5] PX_WAREHOUSE - DATA PERSISTENCE & WORLDLINES")
    print("-" * 80)
    
    # Check directory structure (canonical: Prv_Dossiers, Novel_Dossiers, Learning_Material, Operations)
    warehouse_dir = os.path.join(ROOT_DIR, "PX_Warehouse")
    tester.test_directory_structure(warehouse_dir,
                                   ["Prv_Dossiers", "Novel_Dossiers", "Learning_Material", "Operations"],
                                   "Warehouse directory structure")

    # Test WorldLine Database (canonical module at PX_Warehouse)
    wl_db_mod = tester.test_import("PX_Warehouse.WorldLine_Database",
                                   "WorldLine Database module")

    # Check canonical dossier root exists (Prv_Dossiers)
    sample_dossier = os.path.join(ROOT_DIR, "PX_Warehouse", "Prv_Dossiers")
    tester.test_file_exists(sample_dossier, "Prv_Dossiers directory")
    
    # =========================================================================
    # TEST 6: PX_AUDIT (Monitoring & Protocols)
    # =========================================================================
    print("\n\n[TEST 6] PX_AUDIT - MONITORING & PROTOCOLS")
    print("-" * 80)
    
    # Check for key audit protocols
    audit_dir = os.path.join(ROOT_DIR, "PX_Audit")
    key_protocols = [
        "Autonomous_Research_Cycle.py",
        "Drift_Monitor.py",
        "Protocol_Zero.py",
        "Final_System_Seal.py"
    ]
    
    for protocol in key_protocols:
        protocol_path = os.path.join(audit_dir, protocol)
        tester.test_file_exists(protocol_path, f"Audit protocol: {protocol}")
    
    # =========================================================================
    # TEST 7: PX_SECURITY (Immune System)
    # =========================================================================
    print("\n\n[TEST 7] PX_SECURITY - IMMUNE SYSTEM")
    print("-" * 80)
    
    red_surface_mod = tester.test_import("PX_Security.RedSurface", "RedSurface module")
    
    predator_immune_mod = tester.test_import("PX_Security.PredatorImmune_Block", 
                                            "PredatorImmune Block module")
    
    ip_lock_file = os.path.join(ROOT_DIR, "PX_Security", "IP_LOCK.json")
    tester.test_file_exists(ip_lock_file, "IP Lock configuration")
    
    # =========================================================================
    # TEST 8: PX_CONSTITUTION (Virtual Machine)
    # =========================================================================
    print("\n\n[TEST 8] PX_CONSTITUTION - VIRTUAL MACHINE")
    print("-" * 80)
    
    vm_mod = tester.test_import("PX_Constitution.Virtual_Machine", "Virtual Machine module")
    
    block_universe_mod = tester.test_import("PX_Constitution.Block_Universe", "Block Universe module")
    
    # =========================================================================
    # TEST 9: PX_VALIDATION (System Validation)
    # =========================================================================
    print("\n\n[TEST 9] PX_VALIDATION - SYSTEM VALIDATION")
    print("-" * 80)
    
    validation_dir = os.path.join(ROOT_DIR, "PX_Validation")
    # Require only tests/ (manual_tests may be archived)
    tester.test_directory_structure(validation_dir, ["tests"],
                                   "Validation directory structure")
    
    system_inventory_mod = tester.test_import("PX_Validation.system_inventory", 
                                             "System Inventory module")
    
    # =========================================================================
    # TEST 10: PX_DISCOVERY (Autonomous Discovery)
    # =========================================================================
    print("\n\n[TEST 10] PX_DISCOVERY - AUTONOMOUS DISCOVERY")
    print("-" * 80)
    
    # Test discovery directory structure
    discovery_dir = os.path.join(ROOT_DIR, "PX_Discovery")
    tester.test_file_exists(discovery_dir, "PX_Discovery directory")
    
    # Test discovery modules
    discovery_engine = tester.test_import("PX_Discovery.candidate_discovery_engine",
                                         "Candidate Discovery Engine")
    arc = tester.test_import("PX_Discovery", "PX_Discovery package")
    
    # =========================================================================
    # TEST 11: INTEGRATION TEST (PX_Live_Orchestrator_v2 — canonical single-candidate runner)
    # =========================================================================
    print("\n\n[TEST 11] INTEGRATION - PX_LIVE_ORCHESTRATOR_V2 (CANONICAL)")
    print("-" * 80)

    orchestrator_file = os.path.join(ROOT_DIR, "PX_Executive", "orchestrators", "PX_Live_Orchestrator_v2.py")
    tester.test_file_exists(orchestrator_file, "Live Orchestrator v2 script (canonical)")
    try:
        import PX_Executive.orchestrators.PX_Live_Orchestrator_v2 as live_mod
        print(f"  ✅ Live Orchestrator v2 importable")
        tester.log("passed", "PX_Executive.orchestrators.PX_Live_Orchestrator_v2", "Import validation")
    except Exception as e:
        print(f"  ❌ Live Orchestrator v2 import error: {e}")
        tester.log("failed", "PX_Executive.orchestrators.PX_Live_Orchestrator_v2", "Import validation", e)
    
    # =========================================================================
    # PRINT SUMMARY
    # =========================================================================
    success = tester.print_summary()
    
    # Save detailed report
    report_path = os.path.join(ROOT_DIR, "PX_SYSTEM_TEST_REPORT.json")
    import json
    with open(report_path, 'w') as f:
        json.dump(tester.results, f, indent=2)
    print(f"\n📄 Detailed report saved: {report_path}")
    
    return success


if __name__ == "__main__":
    success = run_comprehensive_test()
    sys.exit(0 if success else 1)
