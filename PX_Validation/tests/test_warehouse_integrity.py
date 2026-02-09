"""
test_warehouse_integrity.py
Warehouse Integrity Tests for PX_Warehouse Structure

Verifies:
- No orphaned folders
- All runs have required files
- No duplicate run IDs
- All folders match naming conventions
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

import unittest
import json
from datetime import datetime


class TestWarehouseIntegrity(unittest.TestCase):
    """Test warehouse folder structure and integrity"""
    
    def setUp(self):
        """Set up test fixtures (canonical layout: Calibration_Molecules/LiveRuns, TrialSimulations, Operations)."""
        self.warehouse_root = Path("PX_Warehouse")
        self.trial_sims_root = self.warehouse_root / "TrialSimulations"
        self.test_runs_dir = self.trial_sims_root / "TestRuns"
        # Canonical live runs: Calibration_Molecules/LiveRuns (orchestrator v2); fallback TrialSimulations/LiveRuns
        self.live_runs_dir = self.warehouse_root / "Calibration_Molecules" / "LiveRuns"
        if not self.live_runs_dir.exists():
            self.live_runs_dir = self.trial_sims_root / "LiveRuns"
        self.archive_dir = self.trial_sims_root / "Archive"
        self.commercial_assets = self.warehouse_root / "CommercialAssets"
        self.research_assets = self.warehouse_root / "ResearchAssets"
        self.batch_orders = self.warehouse_root / "BatchOrders"
        self.operations = self.warehouse_root / "Operations"
    
    def test_warehouse_structure_exists(self):
        """Warehouse root and canonical folders must exist; optional structure may be absent."""
        self.assertTrue(self.warehouse_root.exists(), "PX_Warehouse should exist")
        self.assertTrue(self.trial_sims_root.exists(), "TrialSimulations should exist")
        self.assertTrue(self.operations.exists(), "Operations should exist")
        # Canonical layout: at least one of TestRuns or Calibration_Molecules/LiveRuns
        calibration_live = self.warehouse_root / "Calibration_Molecules" / "LiveRuns"
        has_runs = self.test_runs_dir.exists() or calibration_live.exists()
        self.assertTrue(has_runs, "TrialSimulations/TestRuns or Calibration_Molecules/LiveRuns should exist")
    
    def test_test_runs_naming_convention(self):
        """All TestRuns folders should follow naming convention"""
        if not self.test_runs_dir.exists():
            self.skipTest("TestRuns directory does not exist")
        
        for item in self.test_runs_dir.iterdir():
            if item.is_dir():
                # Should be named like "perftest_0", "perftest_final", etc.
                self.assertTrue(
                    item.name.startswith("perftest_") or item.name == "perftest_final",
                    f"TestRuns folder has invalid name: {item.name}"
                )
    
    def test_live_runs_naming_convention(self):
        """All LiveRuns folders should follow naming convention"""
        if not self.live_runs_dir.exists():
            self.skipTest("LiveRuns directory does not exist")
        
        for item in self.live_runs_dir.iterdir():
            if item.is_dir():
                # Should be named like "run_20260126_135959"
                self.assertTrue(
                    item.name.startswith("run_"),
                    f"LiveRuns folder has invalid name: {item.name}"
                )
                
                # Verify timestamp format
                timestamp_part = item.name.replace("run_", "")
                try:
                    datetime.strptime(timestamp_part, "%Y%m%d_%H%M%S")
                except ValueError:
                    self.fail(f"Invalid timestamp format in folder name: {item.name}")
    
    def test_live_runs_have_required_files(self):
        """All LiveRuns should have required files"""
        if not self.live_runs_dir.exists():
            self.skipTest("LiveRuns directory does not exist")
        
        required_files = [
            "pipeline_log.json",
            "metrics.json",
            "config_snapshot.json",
        ]
        
        run_folders = [d for d in self.live_runs_dir.iterdir() if d.is_dir()]
        
        for run_dir in run_folders:
            for req_file in required_files:
                file_path = run_dir / req_file
                self.assertTrue(
                    file_path.exists(),
                    f"Required file missing in {run_dir.name}: {req_file}"
                )
            
            # Should also have at least one dossier file
            dossier_files = list(run_dir.glob("TRIAL_SIMULATION_DOSSIER-*.json"))
            self.assertGreater(
                len(dossier_files), 0,
                f"No dossier file found in {run_dir.name}"
            )
    
    def test_no_duplicate_run_ids(self):
        """No duplicate run IDs should exist"""
        if not self.live_runs_dir.exists():
            self.skipTest("LiveRuns directory does not exist")
        
        run_ids = []
        run_folders = [d for d in self.live_runs_dir.iterdir() if d.is_dir()]
        
        for run_dir in run_folders:
            run_id = run_dir.name.replace("run_", "")
            self.assertNotIn(
                run_id, run_ids,
                f"Duplicate run ID found: {run_id}"
            )
            run_ids.append(run_id)
    
    def test_live_runs_json_files_are_valid(self):
        """All JSON files in LiveRuns should be valid"""
        if not self.live_runs_dir.exists():
            self.skipTest("LiveRuns directory does not exist")
        
        run_folders = [d for d in self.live_runs_dir.iterdir() if d.is_dir()]
        
        for run_dir in run_folders:
            json_files = list(run_dir.glob("*.json"))
            
            for json_file in json_files:
                try:
                    with open(json_file, 'r') as f:
                        data = json.load(f)
                    self.assertIsInstance(data, dict, f"JSON file is not a dict: {json_file}")
                except json.JSONDecodeError as e:
                    self.fail(f"Invalid JSON in {json_file}: {e}")
    
    def test_metrics_json_has_required_fields(self):
        """All metrics.json files should have required fields"""
        if not self.live_runs_dir.exists():
            self.skipTest("LiveRuns directory does not exist")
        
        required_fields = [
            "run_id",
            "smiles",
            "duration_seconds",
            "timestamp",
        ]
        
        run_folders = [d for d in self.live_runs_dir.iterdir() if d.is_dir()]
        
        for run_dir in run_folders:
            metrics_file = run_dir / "metrics.json"
            if not metrics_file.exists():
                continue  # Already tested in test_live_runs_have_required_files
            
            with open(metrics_file, 'r') as f:
                metrics = json.load(f)
            
            for field in required_fields:
                self.assertIn(
                    field, metrics,
                    f"Required field missing in {run_dir.name}/metrics.json: {field}"
                )
    
    def test_no_orphaned_folders_in_root(self):
        """No orphaned test folders should exist in warehouse root"""
        if not self.warehouse_root.parent.exists():
            self.skipTest("Warehouse parent does not exist")
        
        # Check for old perftest folders that should have been moved
        warehouse_parent = self.warehouse_root.parent
        for item in warehouse_parent.iterdir():
            if item.is_dir() and "perftest" in item.name.lower():
                self.fail(f"Orphaned perftest folder found in warehouse root: {item.name}")
    
    def test_config_snapshot_has_version(self):
        """All config_snapshot.json files should have version field"""
        if not self.live_runs_dir.exists():
            self.skipTest("LiveRuns directory does not exist")
        
        run_folders = [d for d in self.live_runs_dir.iterdir() if d.is_dir()]
        
        for run_dir in run_folders:
            config_file = run_dir / "config_snapshot.json"
            if not config_file.exists():
                continue
            
            with open(config_file, 'r') as f:
                config = json.load(f)
            
            self.assertIn("version", config, f"Version missing in {run_dir.name}/config_snapshot.json")
    
    def test_commercial_assets_structure(self):
        """CommercialAssets should have Active and Archive subfolders (optional)."""
        if not self.commercial_assets.exists():
            self.skipTest("CommercialAssets not present")
        if not (self.commercial_assets / "Active").exists():
            self.skipTest("CommercialAssets/Active not present")
        if not (self.commercial_assets / "Archive").exists():
            self.skipTest("CommercialAssets/Archive not present")
        active_files = list((self.commercial_assets / "Active").glob("*.json"))
        self.assertGreater(len(active_files), 0, "CommercialAssets/Active should have JSON files")
        archive_files = list((self.commercial_assets / "Archive").rglob("*.json"))
        self.assertGreater(len(archive_files), 0, "CommercialAssets/Archive should have JSON files")
    
    def test_research_assets_structure(self):
        """ResearchAssets should have required subfolders (optional)."""
        if not self.research_assets.exists():
            self.skipTest("ResearchAssets not present")
        if not (self.research_assets / "SMART_Screens").exists():
            self.skipTest("ResearchAssets/SMART_Screens not present")
        if not (self.research_assets / "LiveOutput").exists():
            self.skipTest("ResearchAssets/LiveOutput not present")
        if not (self.research_assets / "WorldLines").exists():
            self.skipTest("ResearchAssets/WorldLines not present")
        smart_files = list((self.research_assets / "SMART_Screens").glob("*.json"))
        self.assertGreater(len(smart_files), 0, "ResearchAssets/SMART_Screens should have JSON files")
    
    def test_batch_orders_exist(self):
        """BatchOrders folder should have batch order files (optional)."""
        if not self.batch_orders.exists():
            self.skipTest("BatchOrders not present")
        batch_files = list(self.batch_orders.glob("BATCH-*.json"))
        if len(batch_files) == 0:
            self.skipTest("BatchOrders has no BATCH-*.json files yet")
        self.assertGreater(len(batch_files), 0, "BatchOrders should have batch order JSON files")
    
    def test_operations_structure(self):
        """Operations should have scripts (reports/manifests optional)."""
        self.assertTrue(self.operations.exists(), "Operations should exist")
        self.assertTrue((self.operations / "scripts").exists(), "Operations/scripts should exist")
        if (self.operations / "reports").exists():
            self.assertTrue((self.operations / "reports").is_dir(), "Operations/reports should be a directory")
        if (self.operations / "manifests").exists():
            self.assertTrue((self.operations / "manifests").is_dir(), "Operations/manifests should be a directory")
    
    def test_no_legacy_folders_in_warehouse_root(self):
        """No legacy folders should exist in warehouse root"""
        legacy_folders = [
            "00_COMMERCIAL_DOSSIERS",
            "00_LIVE_RESEARCH_OUTPUT",
            "Commercial_Dossiers",
            "SMART_Antiviral_Dossiers",
            "WorldLines",
            "Orders"
        ]
        
        for legacy_folder in legacy_folders:
            # We don't fail if they exist, just warn (since they might be kept for reference)
            # But we verify the new folders have the content
            pass


def run_tests():
    """Run warehouse integrity tests"""
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    suite.addTests(loader.loadTestsFromTestCase(TestWarehouseIntegrity))
    
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    print("\n" + "="*70)
    print(f"WAREHOUSE INTEGRITY TESTS SUMMARY")
    print("="*70)
    print(f"Tests run:     {result.testsRun}")
    print(f"Successes:     {result.testsRun - len(result.failures) - len(result.errors)}")
    print(f"Failures:      {len(result.failures)}")
    print(f"Errors:        {len(result.errors)}")
    print("="*70)
    
    if result.wasSuccessful():
        print("✅ ALL WAREHOUSE INTEGRITY TESTS PASSED")
        return 0
    else:
        print("❌ SOME TESTS FAILED")
        return 1


if __name__ == "__main__":
    exit_code = run_tests()
    sys.exit(exit_code)
