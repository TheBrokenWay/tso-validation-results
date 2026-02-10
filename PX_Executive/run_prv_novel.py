"""
DEPRECATED: Use px_prv.py --type novel instead.

  python PX_Executive/px_prv.py --type novel
  python PX_Executive/px_prv.py --type novel --limit 100

This file is kept for backward compatibility and still works.
"""
import os
import sys
import subprocess
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
ORCHESTRATOR = REPO_ROOT / "PX_Executive" / "PRV_24H_Orchestrator.py"

def main():
    print("✨ STARTING NOVEL INVENTOR (Authorized Wrapper)...")
    
    # Force NOVEL mode
    env = os.environ.copy()
    env["PRV_MODE"] = "NOVEL"
    env["PRV_QUEUE_FILE"] = "prv_24h_queue.json" 
    env["PYTHONUNBUFFERED"] = "1"

    try:
        cmd = [sys.executable, str(ORCHESTRATOR)]
        subprocess.run(cmd, env=env, check=True)
    except KeyboardInterrupt:
        print("\n🛑 Inventor Stopped.")
    except Exception as e:
        print(f"\n❌ FATAL ERROR: {e}")

if __name__ == "__main__":
    main()
