import json
import os
import numpy as np
from PX_Security.PredatorImmune_Block import PredatorImmuneBlock
from PX_Warehouse.WorldLine_Database import WorldLineDatabase

IMMUNE = PredatorImmuneBlock()
WAREHOUSE = WorldLineDatabase()

def scrub_legacy_warehouse():
    path = WAREHOUSE.path
    print(f"\n>>> [SCRUBBER] Scanning Legacy Resources for Repurposing Potential...")
    
    scrubbed_count = 0
    promoted_count = 0

    for file in os.listdir(path):
        if file.endswith(".worldline"):
            file_path = os.path.join(path, file)
            with open(file_path, "r") as f:
                data = json.load(f)
            
            # Skip files that are already part of the new system
            if data.get("header", {}).get("system_version") == "1.2.0-GAIP":
                continue
                
            scrubbed_count += 1
            
            # Extract coordinates from the old file format
            # (Assuming old format has 'coordinate_35d' or similar)
            coords = np.array(data.get("coordinate_35d", np.zeros(35)))
            p_vec = coords[:4].tolist()
            csa_s = coords[4:9].tolist()
            
            # Re-Evaluate using the NEW Constitutional Ideal
            ctx = {"source_id": "LEGACY-SCRUB", "fingerprint": "PX-ROOT-Sovereign"}
            result = IMMUNE.handle_block_request(f"REPURPOSE-{file}", ctx, p_vec, csa_s)
            
            if result['authorized']:
                promoted_count += 1
                # Promotion: Rewrite the file with the NEW header and BLUE route
                WAREHOUSE.record_materialization(
                    task_id=f"REPURPOSED-{file[:8]}",
                    block=result['block'],
                    coherence=result['coherence'],
                    lab_results=data.get("physical_realization", {"status": "LEGACY_CONVERTED"}),
                    route="BLUE"
                )
                # Optionally remove the old unlabelled file
                # os.remove(file_path)

    print(f"\n>>> [SCRUBBER COMPLETE]")
    print(f"    Scanned:  {scrubbed_count}")
    print(f"    Promoted: {promoted_count} (New 'Blue' Repurposed Candidates)")
    print(f"    Rejected: {scrubbed_count - promoted_count}")

if __name__ == "__main__":
    scrub_legacy_warehouse()
