import os
import time
import glob
from datetime import datetime

_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
WAREHOUSE_ROOT = os.path.join(_REPO_ROOT, "PX_Warehouse")
PATHS = {
    "Novel": os.path.join(WAREHOUSE_ROOT, "Novel_Dossiers"),
    "Repurposed": os.path.join(WAREHOUSE_ROOT, "Prv_Dossiers"),
    "WorldLines": os.path.join(WAREHOUSE_ROOT, "WorldLines"),
    "Vault": os.path.join(WAREHOUSE_ROOT, "Finalized_Dossiers")
}

def count_files(directory, extension="*.json"):
    """Deterministic file count across recursive tiers."""
    if not os.path.exists(directory):
        return 0
    return len(glob.glob(os.path.join(directory, "**", extension), recursive=True))

def run_monitor():
    print(f"--- PREDATOR X IDE MONITOR (FULL VAULT) INITIALIZED ---")
    try:
        while True:
            # 1. RAW INVENTORY
            n_count = count_files(PATHS["Novel"])
            r_count = count_files(PATHS["Repurposed"])
            total_raw = n_count + r_count

            # 2. PHYSICS
            wl_count = count_files(PATHS["WorldLines"], "*.*")

            # 3. FULL TIERED VAULT
            diamond = count_files(os.path.join(PATHS["Vault"], "Diamond"))
            gold    = count_files(os.path.join(PATHS["Vault"], "Gold"))
            silver  = count_files(os.path.join(PATHS["Vault"], "Silver"))
            bronze  = count_files(os.path.join(PATHS["Vault"], "Bronze"))
            total_auth = diamond + gold + silver + bronze

            # 4. METRICS
            gap = max(0, total_raw - total_auth)
            yield_pct = (total_auth / total_raw * 100) if total_raw > 0 else 0

            # IDE FORMATTED OUTPUT
            os.system('cls' if os.name == 'nt' else 'clear')
            print(f"[{datetime.now().strftime('%H:%M:%S')}] PREDATOR X | STATUS: CONVERTING")
            print("-" * 45)
            print(f"RAW INVENTORY:    {total_raw} (N:{n_count} R:{r_count})")
            print(f"WORLDLINES:       {wl_count}")
            print("-" * 45)
            print(f"💎 DIAMOND VAULT:  {diamond}")
            print(f"🥇 GOLD VAULT:     {gold}")
            print(f"🥈 SILVER VAULT:   {silver}")
            print(f"🥉 BRONZE VAULT:   {bronze}")
            print("-" * 45)
            print(f"DIAGNOSTIC GAP:   {gap} pending")
            print(f"CURRENT YIELD:    {yield_pct:.2f}%")
            print("-" * 45)
            print("FTO SYNTAX STATUS: [LITERAL QUOTING ACTIVE]")
            
            time.sleep(10)
    except KeyboardInterrupt:
        print("\nMonitor detached.")

if __name__ == "__main__":
    run_monitor()
