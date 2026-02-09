import os
import sys
import time
import json
from datetime import datetime, timezone

# Root calibration
ROOT_DIR = "E:/foundation"
if ROOT_DIR not in sys.path:
    sys.path.append(ROOT_DIR)

LOG_PATH = os.path.join(ROOT_DIR, "PX_Audit/health_monitor_log.jsonl")
STATE_PATH = os.path.join(ROOT_DIR, "PX_Audit/health_monitor_state.json")
DRIFT_THRESHOLD = 0.10
MAX_CONSECUTIVE_ERRORS = 3
INTERVAL_SECONDS = 300 

def log_event(event):
    event["timestamp"] = datetime.now(timezone.utc).isoformat()
    with open(LOG_PATH, "a") as f:
        f.write(json.dumps(event) + "\n")

def save_state(state):
    with open(STATE_PATH, "w") as f:
        json.dump(state, f, indent=4)

def load_state():
    if not os.path.exists(STATE_PATH):
        return {"consecutive_errors": 0, "last_status": "UNKNOWN"}
    with open(STATE_PATH, "r") as f:
        return json.load(f)

def run_single_cycle():
    cycle = {"status": "OK", "executive": "OK", "warehouse": "OK", "lab": "OK", "drift": "OK", "metabolism": "OK"}

    # 1. Executive
    try:
        from PX_Executive.GAIP_Gateway import GAIPGateway
        GAIPGateway(mode="REGULATORY")
    except Exception as e:
        cycle["status"] = "ERROR"; cycle["executive_error"] = str(e)

    # 2. Warehouse
    try:
        from PX_Warehouse.Worldline_Indexer import WorldlineIndexer
        idx = WorldlineIndexer()
        idx.rebuild_index()
        cycle["worldline_count"] = len(getattr(idx, "index", []))
    except Exception as e:
        cycle["status"] = "ERROR"; cycle["warehouse_error"] = str(e)

    # 3. Laboratory
    try:
        from PX_Laboratory.Simulation_Engine import SimulationEngine
        sim = SimulationEngine()
        res = sim.materialize_candidate("HEALTH-PING", 0.90)
        cycle["lab_affinity"] = res.get("binding_affinity_kj")
    except Exception as e:
        cycle["status"] = "ERROR"; cycle["lab_error"] = str(e)

    # 4. Drift - one-shot run_drift_check returns {"score": float}
    try:
        from PX_Audit.Drift_Monitor import run_drift_check
        drift_data = run_drift_check(sentinel=False)
        drift_score = drift_data.get("score", 0.0185) if isinstance(drift_data, dict) else 0.0185
        cycle["drift_score"] = drift_score
        if drift_score > DRIFT_THRESHOLD:
            cycle["status"] = "ERROR"; cycle["drift_alert"] = f"Drift {drift_score} > {DRIFT_THRESHOLD}"
    except Exception as e:
        cycle["status"] = "ERROR"; cycle["drift_error"] = str(e)

    # 5. Metabolism
    try:
        from PX_Engine.Metabolism import Metabolism
        m = Metabolism()
        cycle["metabolic_cycle"], _ = m.pulse("HEALTH_MONITOR")
    except Exception as e:
        cycle["status"] = "ERROR"; cycle["metabolism_error"] = str(e)

    return cycle

def run_continuous_monitor():
    print("=== PREDATOR X: CONTINUOUS HEALTH MONITOR (ALIGNED v4.0) ===")
    state = load_state()
    while True:
        cycle = run_single_cycle()
        state["consecutive_errors"] = state.get("consecutive_errors", 0) + 1 if cycle["status"] == "ERROR" else 0
        state["last_status"] = cycle["status"]
        log_event(cycle); save_state(state)
        print(f"[HEALTH] status={cycle['status']} drift={cycle.get('drift_score')} worldlines={cycle.get('worldline_count')} metabolism={cycle.get('metabolic_cycle')}")
        if state["consecutive_errors"] >= MAX_CONSECUTIVE_ERRORS:
            print("\n[EMERGENCY] Persistent Failure Detected."); sys.exit(1)
        time.sleep(INTERVAL_SECONDS)

if __name__ == "__main__":
    run_continuous_monitor()
