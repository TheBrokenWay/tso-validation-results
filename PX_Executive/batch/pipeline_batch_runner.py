"""
pipeline_batch_runner.py
Batch Runner for Predator X v2.0-CORE

Processes multiple molecules through the complete pipeline.
Generates structured outputs in PX_Warehouse/Calibration_Molecules/BatchRuns/

Usage:
    python PX_Executive/batch/pipeline_batch_runner.py
"""

import subprocess
import json
import sys
from pathlib import Path
from datetime import datetime, UTC

# Add project root to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))


# Molecule library
MOLECULES = [
    {
        "name": "Aticaprant",
        "id": "J&J-ATICAPRANT",
        "smiles": "CC1=CC(=CC(=C1)[C@@H]2CCCN2CC3=CC=C(C=C3)OC4=C(C=C(C=C4)C(=O)N)F)C",
        "indication": "Kappa-opioid receptor antagonist (repurposing candidate)"
    },
    {
        "name": "Aspirin",
        "id": "ASPIRIN-001",
        "smiles": "CC(=O)Oc1ccccc1C(=O)O",
        "indication": "Analgesic/Anti-inflammatory"
    },
    {
        "name": "Ibuprofen",
        "id": "IBUPROFEN-001",
        "smiles": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
        "indication": "NSAID"
    },
    {
        "name": "Caffeine",
        "id": "CAFFEINE-001",
        "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "indication": "CNS stimulant"
    },
    {
        "name": "Paracetamol",
        "id": "PARACETAMOL-001",
        "smiles": "CC(=O)Nc1ccc(cc1)O",
        "indication": "Analgesic/Antipyretic"
    },
    {
        "name": "Benzene",
        "id": "BENZENE-001",
        "smiles": "c1ccccc1",
        "indication": "Test compound (simple aromatic)"
    },
    {
        "name": "Ethanol",
        "id": "ETHANOL-001",
        "smiles": "CCO",
        "indication": "Test compound (minimal structure)"
    },
]


def run_pipeline(molecule: dict, batch_id: str, verbose: bool = False) -> dict:
    """
    Run pipeline for a single molecule.
    
    Args:
        molecule: Dict with name, id, smiles, indication
        batch_id: Batch identifier
        verbose: Show pipeline output
    
    Returns:
        Dict with run results and paths
    """
    timestamp = datetime.now(UTC).strftime("%Y%m%d_%H%M%S")
    
    # Build command
    cmd = [
        sys.executable,
        "PX_Executive/orchestrators/PX_Live_Orchestrator_v2.py",
        "--smiles", molecule["smiles"],
        "--name", molecule["name"],
        "--id", molecule["id"],
        "--indication", molecule.get("indication", "Not specified"),
    ]
    
    if not verbose:
        cmd.append("--quiet")
    
    # Run pipeline
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        cwd=project_root
    )
    
    # Parse output to find run directory
    run_dir = None
    for line in result.stdout.split('\n'):
        if "Run directory:" in line:
            run_dir = line.split("Run directory:")[-1].strip()
            break
    
    return {
        "molecule_name": molecule["name"],
        "molecule_id": molecule["id"],
        "smiles": molecule["smiles"],
        "indication": molecule.get("indication", "Not specified"),
        "timestamp": timestamp,
        "run_directory": run_dir,
        "exit_code": result.returncode,
        "success": result.returncode == 0,
        "stdout": result.stdout,
        "stderr": result.stderr,
    }


def run_batch(molecules: list = None, batch_name: str = None, verbose: bool = False) -> dict:
    """
    Run batch of molecules through pipeline.
    
    Args:
        molecules: List of molecule dicts (default: MOLECULES)
        batch_name: Name for this batch (default: auto-generated)
        verbose: Show detailed output
    
    Returns:
        Dict with batch results
    """
    if molecules is None:
        molecules = MOLECULES
    
    if batch_name is None:
        batch_name = f"BATCH_{datetime.now(UTC).strftime('%Y%m%d_%H%M%S')}"
    
    print(f"\n{'='*70}")
    print(f"PREDATOR X v2.0-CORE BATCH RUNNER")
    print(f"{'='*70}")
    print(f"Batch ID: {batch_name}")
    print(f"Molecules: {len(molecules)}")
    print(f"Started: {datetime.now(UTC).isoformat()}")
    print(f"{'='*70}\n")
    
    results = []
    successful = 0
    failed = 0
    
    for i, molecule in enumerate(molecules, 1):
        print(f"[{i}/{len(molecules)}] Processing: {molecule['name']} ({molecule['id']})")
        
        try:
            result = run_pipeline(molecule, batch_name, verbose=verbose)
            results.append(result)
            
            if result["success"]:
                successful += 1
                print(f"    ✅ SUCCESS - {result['run_directory']}")
            else:
                failed += 1
                print(f"    ❌ FAILED - Exit code: {result['exit_code']}")
                if verbose:
                    print(f"    Error: {result['stderr']}")
        
        except Exception as e:
            failed += 1
            print(f"    ❌ EXCEPTION: {e}")
            results.append({
                "molecule_name": molecule["name"],
                "molecule_id": molecule["id"],
                "success": False,
                "error": str(e)
            })
    
    # Create batch summary
    batch_result = {
        "batch_id": batch_name,
        "timestamp": datetime.now(UTC).isoformat(),
        "total_molecules": len(molecules),
        "successful": successful,
        "failed": failed,
        "success_rate": f"{(successful/len(molecules)*100):.1f}%",
        "results": results,
    }
    
    # Save batch summary
    from PX_Warehouse.warehouse_layout import get_calibration_molecules_dir
    output_dir = get_calibration_molecules_dir(project_root) / "BatchRuns"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    summary_file = output_dir / f"{batch_name}_SUMMARY.json"
    with open(summary_file, 'w') as f:
        json.dump(batch_result, f, indent=2, default=str)
    
    print(f"\n{'='*70}")
    print(f"BATCH COMPLETE")
    print(f"{'='*70}")
    print(f"Total:      {len(molecules)}")
    print(f"Successful: {successful}")
    print(f"Failed:     {failed}")
    print(f"Success Rate: {batch_result['success_rate']}")
    print(f"\nBatch Summary: {summary_file}")
    print(f"{'='*70}\n")
    
    return batch_result


def main():
    """Main entry point"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Predator X v2.0-CORE Batch Runner",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run default molecule set
  python PX_Executive/batch/pipeline_batch_runner.py
  
  # Verbose output
  python PX_Executive/batch/pipeline_batch_runner.py --verbose
  
  # Custom batch name
  python PX_Executive/batch/pipeline_batch_runner.py --batch-name "REPURPOSING-SCREEN-JAN2026"
        """
    )
    
    parser.add_argument(
        "--batch-name",
        type=str,
        default=None,
        help="Custom batch name (default: auto-generated)"
    )
    
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Show detailed pipeline output"
    )
    
    args = parser.parse_args()
    
    # Run batch
    result = run_batch(
        molecules=MOLECULES,
        batch_name=args.batch_name,
        verbose=args.verbose
    )
    
    # Exit with appropriate code
    if result["failed"] > 0:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
