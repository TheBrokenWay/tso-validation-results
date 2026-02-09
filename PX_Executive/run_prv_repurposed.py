"""
PX_Executive/run_prv_repurposed.py
WRAPPER: Launches the 24H Orchestrator in REPURPOSED mode.
"""
import os
import sys
import subprocess
from pathlib import Path

# CONFIGURATION
REPO_ROOT = Path(__file__).resolve().parents[1]
ORCHESTRATOR = REPO_ROOT / "PX_Executive" / "PRV_24H_Orchestrator.py"

def main():
    print("⛏️  STARTING REPURPOSED MINER (Authorized Wrapper)...")
    
    # Environment Variables to Force the Orchestrator into correct mode
    env = os.environ.copy()
    env["PRV_MODE"] = "REPURPOSED"
    env["PRV_QUEUE_FILE"] = "prv_24h_queue.json"  # Looks in Feeder
    env["PYTHONUNBUFFERED"] = "1"

    # Launch the Master Engine
    try:
        cmd = [sys.executable, str(ORCHESTRATOR)]
        subprocess.run(cmd, env=env, check=True)
    except KeyboardInterrupt:
        print("\n🛑 Miner Stopped by User.")
    except Exception as e:
        print(f"\n❌ FATAL ERROR: {e}")

if __name__ == "__main__":
    main()