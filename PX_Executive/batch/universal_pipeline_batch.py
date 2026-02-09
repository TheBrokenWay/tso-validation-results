"""
universal_pipeline_batch.py
Batch Processing Script for Universal Pipeline Runner

Automatically discovers and processes assets from consolidated warehouse
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from PX_Executive.UniversalPipelineRunner import UniversalPipelineRunner
import argparse


def discover_assets(source_folder: str, limit: int = None) -> list:
    """
    Discover all JSON assets in a source folder (canonical PX_Warehouse paths).
    Uses repo root so paths resolve correctly regardless of cwd.
    """
    repo_root = Path(__file__).resolve().parents[2]  # batch -> PX_Executive -> foundation
    warehouse = repo_root / "PX_Warehouse"
    source_path = warehouse / source_folder

    if not source_path.exists():
        print(f"❌ Source folder not found: {source_path}")
        return []

    json_files = list(source_path.rglob("*.json"))
    worldline_files = list(source_path.rglob("*.worldline"))
    all_files = json_files + worldline_files

    if limit:
        all_files = all_files[:limit]

    return all_files


def main():
    parser = argparse.ArgumentParser(
        description="Batch process warehouse assets through Universal Pipeline Runner"
    )
    parser.add_argument(
        "--source",
        type=str,
        required=True,
        choices=[
            "Prv_Dossiers",
            "Novel_Dossiers",
            "WorldLines",
            "Calibration_Molecules",
            "Learning_Material",
            "Feeder",
            "ResearchAssets/SMART_Screens",
            "ResearchAssets/LiveOutput",
            "ResearchAssets/WorldLines",
            "CommercialAssets/Active",
            "CommercialAssets/Archive",
            "TrialSimulations/Archive/Legacy",
        ],
        help="Source folder in PX_Warehouse (prefer Prv_Dossiers, Novel_Dossiers, WorldLines)"
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Maximum number of assets to process (default: all)"
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress verbose output"
    )
    
    args = parser.parse_args()
    
    # Discover assets
    print(f"Discovering assets in: {args.source}")
    assets = discover_assets(args.source, args.limit)
    
    if not assets:
        print("❌ No assets found")
        sys.exit(1)
    
    print(f"Found {len(assets)} asset(s)")
    
    # Create runner
    runner = UniversalPipelineRunner(verbose=not args.quiet)
    
    # Process batch
    stats = runner.process_batch(assets)
    
    # Exit with appropriate code
    sys.exit(0 if stats["failed"] == 0 else 1)


if __name__ == "__main__":
    main()
