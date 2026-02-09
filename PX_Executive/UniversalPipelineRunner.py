"""
UniversalPipelineRunner.py
Universal Pipeline Runner for Predator X v2.0-CORE

Accepts ANY research asset from consolidated warehouse and:
1. Extracts molecule descriptors (SMILES, structure, metadata)
2. Runs through full v2.0-CORE pipeline
3. Applies constitutional grading
4. Sorts output into PRV_Sorted/<GRADE>/
5. Logs all operations comprehensively
"""

import sys
import json
from pathlib import Path
from datetime import datetime, timezone
from typing import Dict, Any, Optional, List
import traceback
import subprocess

# Repo root for canonical warehouse paths
_REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_REPO_ROOT))

from PX_Engine.operations.GradingEngine import GradingEngine


class UniversalPipelineRunner:
    """
    Universal pipeline runner for any warehouse asset
    
    Capabilities:
    - Auto-detect molecule descriptors from any JSON
    - Run full v2.0-CORE pipeline
    - Apply constitutional grading
    - Sort into tiered folders
    - Comprehensive logging
    """
    
    # Canonical warehouse sources (Prv/Novel tiered + trial outputs; legacy paths for backward read)
    VALID_SOURCES = [
        "Prv_Dossiers",
        "Novel_Dossiers",
        "Learning_Material",
        "Calibration_Molecules",
        "Feeder",
        "TrialSimulations/LiveRuns",
        "TrialSimulations/BatchRuns",
        "CommercialAssets/Gold",
        "CommercialAssets/Silver",
        "CommercialAssets/Diamond",
        "CommercialAssets/Bronze",
        "Commercial_Dossiers",
        "ResearchAssets/SMART_Screens",
        "ResearchAssets/LiveOutput",
        "CommercialAssets/Active",
        "CommercialAssets/Archive",
        "TrialSimulations/Archive/Legacy",
    ]
    
    def __init__(self, warehouse_root="PX_Warehouse", logs_dir="PX_LOGS", verbose=True):
        self.warehouse_root = Path(warehouse_root)
        if not self.warehouse_root.is_absolute():
            self.warehouse_root = (_REPO_ROOT / self.warehouse_root).resolve()
        self.logs_dir = Path(logs_dir)
        if not self.logs_dir.is_absolute():
            self.logs_dir = (_REPO_ROOT / self.logs_dir).resolve()
        self.logs_dir.mkdir(parents=True, exist_ok=True)
        self.verbose = verbose
        
        # Sorted output: Learning_Material/PRV_Sorted/<grade>/ (canonical)
        self.sorted_root = self.warehouse_root / "Learning_Material" / "PRV_Sorted"
        for grade in ["GOLD_TIER", "SILVER_TIER", "BRONZE_TIER", "REJECTED", "NEEDS_REVIEW"]:
            (self.sorted_root / grade).mkdir(parents=True, exist_ok=True)
        
        # Initialize grading engine
        self.grading_engine = GradingEngine(verbose=verbose)
        
        # Stats
        self.stats = {
            "total_processed": 0,
            "successful": 0,
            "failed": 0,
            "needs_review": 0,
            "grades": {
                "GOLD_TIER": 0,
                "SILVER_TIER": 0,
                "BRONZE_TIER": 0,
                "REJECTED": 0,
                "NEEDS_REVIEW": 0,
            }
        }
    
    def extract_molecule_descriptor(self, asset_data: Dict[str, Any], asset_path: Path) -> Optional[Dict[str, str]]:
        """
        Extract molecule descriptor from any asset
        
        Args:
            asset_data: JSON data from asset
            asset_path: Path to asset file
        
        Returns:
            Dictionary with 'smiles', 'name', 'id' or None if not found
        """
        descriptor = {}
        
        # Try direct SMILES field
        if "smiles" in asset_data:
            descriptor["smiles"] = asset_data["smiles"]
        elif "SMILES" in asset_data:
            descriptor["smiles"] = asset_data["SMILES"]
        elif "structure" in asset_data and isinstance(asset_data["structure"], str):
            descriptor["smiles"] = asset_data["structure"]
        
        # Try nested fields (for Evidence Package v3)
        if not descriptor.get("smiles"):
            if "metadata" in asset_data:
                metadata = asset_data["metadata"]
                if "smiles" in metadata:
                    descriptor["smiles"] = metadata["smiles"]
                elif "structure" in metadata:
                    descriptor["smiles"] = metadata["structure"]
        
        # Try candidate fields (for WorldLine assets)
        if not descriptor.get("smiles"):
            if "candidate" in asset_data:
                candidate = asset_data["candidate"]
                if "smiles" in candidate:
                    descriptor["smiles"] = candidate["smiles"]
            elif "prv_candidate" in asset_data:
                prv_candidate = asset_data["prv_candidate"]
                if "smiles" in prv_candidate:
                    descriptor["smiles"] = prv_candidate["smiles"]
            elif "candidate_data" in asset_data:
                candidate_data = asset_data["candidate_data"]
                if "prv_candidate" in candidate_data:
                    prv_candidate = candidate_data["prv_candidate"]
                    if "smiles" in prv_candidate:
                        descriptor["smiles"] = prv_candidate["smiles"]
        
        # Extract name
        descriptor["name"] = (
            asset_data.get("name") or 
            asset_data.get("compound_name") or 
            asset_data.get("candidate_id") or 
            (asset_data.get("prv_candidate", {}).get("common_name") if "prv_candidate" in asset_data else None) or
            (asset_data.get("candidate_data", {}).get("prv_candidate", {}).get("common_name") if "candidate_data" in asset_data else None) or
            (asset_data.get("header", {}).get("worldline_id") if "header" in asset_data else None) or
            asset_path.stem
        )
        
        # Extract ID
        descriptor["id"] = (
            asset_data.get("id") or 
            asset_data.get("compound_id") or 
            asset_data.get("candidate_id") or 
            (asset_data.get("prv_candidate", {}).get("chembl_id") if "prv_candidate" in asset_data else None) or
            (asset_data.get("candidate_data", {}).get("prv_candidate", {}).get("chembl_id") if "candidate_data" in asset_data else None) or
            (asset_data.get("header", {}).get("worldline_id") if "header" in asset_data else None) or
            asset_path.stem
        )
        
        if not descriptor.get("smiles"):
            return None
        
        return descriptor
    
    def run_pipeline(self, smiles: str, name: str, compound_id: str) -> Optional[Dict[str, Any]]:
        """
        Run full v2.0-CORE pipeline via PX_Live_Orchestrator_v2.py
        
        Args:
            smiles: SMILES string
            name: Compound name
            compound_id: Compound ID
        
        Returns:
            Results dictionary with dossier path, or None on failure
        """
        try:
            # Call orchestrator as subprocess
            orchestrator_path = Path("PX_Executive/orchestrators/PX_Live_Orchestrator_v2.py")
            
            cmd = [
                sys.executable,
                str(orchestrator_path),
                "--smiles", smiles,
                "--name", name,
                "--id", compound_id
            ]
            
            if self.verbose:
                print(f"🚀 Running pipeline for: {name}")
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300  # 5 minute timeout
            )
            
            if result.returncode != 0:
                if self.verbose:
                    print(f"❌ Pipeline failed with exit code {result.returncode}")
                    print(f"Error: {result.stderr}")
                return None
            
            # Parse output to find run directory and dossier path
            output_lines = result.stdout.split('\n')
            run_dir = None
            dossier_path = None
            
            for line in output_lines:
                if "Dossier:" in line:
                    # Extract dossier path
                    parts = line.split(':', 1)
                    if len(parts) == 2:
                        dossier_path = parts[1].strip()
                        # Extract run_dir from dossier path
                        dossier_p = Path(dossier_path)
                        run_dir = str(dossier_p.parent)
                        break
            
            if not run_dir or not dossier_path:
                if self.verbose:
                    print("⚠️  Could not parse run directory from output")
                return None
            
            # Verify paths exist
            run_path = Path(run_dir)
            dossier_p = Path(dossier_path)
            
            if not run_path.exists():
                if self.verbose:
                    print(f"⚠️  Run directory does not exist: {run_dir}")
                return None
            
            if not dossier_p.exists():
                if self.verbose:
                    print(f"⚠️  Dossier file does not exist: {dossier_path}")
                return None
            
            return {
                "run_directory": str(run_path),
                "dossier_path": str(dossier_p),
                "success": True,
            }
        
        except subprocess.TimeoutExpired:
            if self.verbose:
                print(f"❌ Pipeline timeout (5 minutes)")
            return None
        
        except Exception as e:
            if self.verbose:
                print(f"❌ Pipeline error: {e}")
                traceback.print_exc()
            return None
    
    def sort_output(self, run_directory: str, dossier_path: str, grade_result: Dict[str, Any]) -> bool:
        """
        Sort pipeline output into PRV_Sorted/<GRADE>/
        
        Args:
            run_directory: Path to run directory
            dossier_path: Path to dossier file
            grade_result: Grading result
        
        Returns:
            True if successful, False otherwise
        """
        try:
            grade = grade_result["grade"]
            run_dir = Path(run_directory)
            
            # Create timestamped folder in sorted tier
            timestamp = datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')
            compound_name = Path(dossier_path).stem.replace("TRIAL_SIMULATION_DOSSIER-", "")
            sorted_dir = self.sorted_root / grade / f"{compound_name}_{timestamp}"
            sorted_dir.mkdir(parents=True, exist_ok=True)
            
            # Copy all files from run directory
            import shutil
            for file in run_dir.glob("*"):
                if file.is_file():
                    shutil.copy2(file, sorted_dir / file.name)
            
            # Save grade metadata
            grade_file = sorted_dir / "grade.json"
            with open(grade_file, 'w') as f:
                json.dump(grade_result, f, indent=2)
            
            if self.verbose:
                print(f"📁 Sorted to: {sorted_dir.relative_to(self.warehouse_root)}")
            
            return True
        
        except Exception as e:
            if self.verbose:
                print(f"❌ Error sorting output: {e}")
                traceback.print_exc()
            return False
    
    def log_result(self, asset_path: Path, success: bool, grade: Optional[str] = None, error: Optional[str] = None):
        """Log processing result"""
        timestamp = datetime.now(timezone.utc).isoformat() + "Z"
        
        if success and grade:
            log_file = self.logs_dir / "pipeline_success.log"
            with open(log_file, 'a') as f:
                f.write(f"{timestamp} | SUCCESS | {asset_path.name} | Grade: {grade}\n")
        
        elif error:
            log_file = self.logs_dir / "pipeline_failures.log"
            with open(log_file, 'a') as f:
                f.write(f"{timestamp} | FAILURE | {asset_path.name} | Error: {error}\n")
        
        if grade == "NEEDS_REVIEW":
            log_file = self.logs_dir / "pipeline_needs_review.log"
            with open(log_file, 'a') as f:
                f.write(f"{timestamp} | NEEDS_REVIEW | {asset_path.name}\n")
    
    def process_asset(self, asset_path: Path) -> bool:
        """
        Process a single asset through the universal pipeline
        
        Args:
            asset_path: Path to asset JSON file
        
        Returns:
            True if successful, False otherwise
        """
        if self.verbose:
            print(f"\n{'='*70}")
            print(f"Processing: {asset_path.name}")
            print(f"{'='*70}")
        
        try:
            # Load asset
            with open(asset_path, 'r') as f:
                asset_data = json.load(f)
            
            # Extract molecule descriptor
            descriptor = self.extract_molecule_descriptor(asset_data, asset_path)
            if not descriptor:
                error = "No valid molecule descriptor found"
                if self.verbose:
                    print(f"❌ {error}")
                self.log_result(asset_path, False, error=error)
                self.stats["failed"] += 1
                return False
            
            smiles = descriptor["smiles"]
            name = descriptor["name"]
            compound_id = descriptor["id"]
            
            if self.verbose:
                print(f"✅ Extracted: {name} ({compound_id})")
                print(f"   SMILES: {smiles[:50]}...")
            
            # Run pipeline
            pipeline_result = self.run_pipeline(smiles, name, compound_id)
            if not pipeline_result:
                error = "Pipeline execution failed"
                if self.verbose:
                    print(f"❌ {error}")
                self.log_result(asset_path, False, error=error)
                self.stats["failed"] += 1
                return False
            
            if self.verbose:
                print(f"✅ Pipeline complete")
            
            # Grade dossier (pass Path, not dict)
            dossier_path = Path(pipeline_result["dossier_path"])
            grade_result = self.grading_engine.grade_dossier(dossier_path)
            grade = grade_result["grade"]
            
            # Sort output
            sort_success = self.sort_output(
                pipeline_result["run_directory"],
                pipeline_result["dossier_path"],
                grade_result
            )
            
            if not sort_success:
                error = "Failed to sort output"
                if self.verbose:
                    print(f"⚠️  {error}")
                self.log_result(asset_path, False, error=error)
                self.stats["failed"] += 1
                return False
            
            # Log success
            self.log_result(asset_path, True, grade=grade)
            self.stats["successful"] += 1
            self.stats["grades"][grade] += 1
            
            if grade == "NEEDS_REVIEW":
                self.stats["needs_review"] += 1
            
            if self.verbose:
                print(f"✅ Processing complete: {grade}")
            
            return True
        
        except Exception as e:
            error = f"Unexpected error: {str(e)}"
            if self.verbose:
                print(f"❌ {error}")
                traceback.print_exc()
            self.log_result(asset_path, False, error=error)
            self.stats["failed"] += 1
            return False
        
        finally:
            self.stats["total_processed"] += 1
    
    def process_batch(self, asset_paths: List[Path]) -> Dict[str, Any]:
        """
        Process multiple assets in batch
        
        Args:
            asset_paths: List of paths to asset files
        
        Returns:
            Statistics dictionary
        """
        print(f"\n{'='*70}")
        print(f"UNIVERSAL PIPELINE RUNNER - BATCH PROCESSING")
        print(f"{'='*70}")
        print(f"Assets to process: {len(asset_paths)}")
        print(f"Started: {datetime.now(timezone.utc).isoformat()}Z")
        print(f"{'='*70}\n")
        
        for i, asset_path in enumerate(asset_paths, 1):
            print(f"\n[{i}/{len(asset_paths)}] Processing: {asset_path.name}")
            self.process_asset(asset_path)
        
        # Generate summary
        print(f"\n{'='*70}")
        print(f"BATCH PROCESSING COMPLETE")
        print(f"{'='*70}")
        print(f"Total Processed:     {self.stats['total_processed']}")
        print(f"Successful:          {self.stats['successful']}")
        print(f"Failed:              {self.stats['failed']}")
        print(f"")
        print(f"Grade Distribution:")
        for grade, count in self.stats["grades"].items():
            print(f"  {grade:20s}: {count}")
        print(f"{'='*70}\n")
        
        # Save summary
        self.save_classification_summary()
        
        return self.stats
    
    def save_classification_summary(self):
        """Save classification summary to audit reports"""
        summary_path = Path("PX_Audit/reports/PIPELINE_CLASSIFICATION_SUMMARY.md")
        
        # Read existing summary or create new
        if summary_path.exists():
            with open(summary_path, 'r') as f:
                existing = f.read()
        else:
            existing = "# PIPELINE CLASSIFICATION SUMMARY\n\n"
        
        # Append new entry
        timestamp = datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M:%S UTC')
        entry = f"""
## Batch Run: {timestamp}

**Statistics:**
- Total Processed: {self.stats['total_processed']}
- Successful: {self.stats['successful']}
- Failed: {self.stats['failed']}

**Grade Distribution:**
- GOLD_TIER: {self.stats['grades']['GOLD_TIER']}
- SILVER_TIER: {self.stats['grades']['SILVER_TIER']}
- BRONZE_TIER: {self.stats['grades']['BRONZE_TIER']}
- NEEDS_REVIEW: {self.stats['grades']['NEEDS_REVIEW']}
- REJECTED: {self.stats['grades']['REJECTED']}

---
"""
        
        with open(summary_path, 'w') as f:
            f.write(existing + entry)


if __name__ == "__main__":
    import argparse
    
    # Constitutional gate: poison pill must run before PX_Engine / PX_Warehouse
    from governance.poison_pill_gate import run_poison_pill_gate, require_poison_pill_gate_before_pipeline
    if not run_poison_pill_gate():
        sys.stderr.write("❌ Poison pill gate failed. Pipeline blocked.\n")
        sys.exit(1)
    require_poison_pill_gate_before_pipeline()

    parser = argparse.ArgumentParser(description="Universal Pipeline Runner for Predator X v2.0-CORE")
    parser.add_argument("assets", nargs="+", help="Path(s) to asset JSON file(s)")
    parser.add_argument("--quiet", action="store_true", help="Suppress verbose output")
    
    args = parser.parse_args()
    
    # Create runner
    runner = UniversalPipelineRunner(verbose=not args.quiet)
    
    # Convert paths to Path objects
    asset_paths = [Path(p) for p in args.assets]
    
    # Validate paths
    valid_paths = []
    for path in asset_paths:
        if not path.exists():
            print(f"❌ File not found: {path}")
        elif path.suffix != ".json":
            print(f"❌ Not a JSON file: {path}")
        else:
            valid_paths.append(path)
    
    if not valid_paths:
        print("❌ No valid assets to process")
        sys.exit(1)
    
    # Process batch
    stats = runner.process_batch(valid_paths)
    
    # Exit with appropriate code
    sys.exit(0 if stats["failed"] == 0 else 1)
