import os
import shutil
from pathlib import Path

def recover_assets():
    gold_dir = Path(r"E:\foundation\PX_Warehouse\CommercialAssets\Gold")
    learning_dir = Path(r"E:\foundation\PX_Warehouse\Learning_Material")
    
    if not learning_dir.exists():
        print(f"Error: Learning directory not found at {learning_dir}")
        return

    os.makedirs(gold_dir, exist_ok=True)

    recovered_count = 0
    
    # We are looking for the dossiers that were moved
    for dossier_file in learning_dir.glob("TRIAL_SIMULATION_DOSSIER-*.json"):
        try:
            dest_path = gold_dir / dossier_file.name
            shutil.move(str(dossier_file), str(dest_path))
            recovered_count += 1
            if recovered_count % 50 == 0:
                print(f"Recovered {recovered_count} files...")
        except Exception as e:
            print(f"Error recovering {dossier_file.name}: {e}")

    print(f"\nRecovery Complete:")
    print(f"Total Assets Recovered: {recovered_count}")

if __name__ == "__main__":
    recover_assets()
