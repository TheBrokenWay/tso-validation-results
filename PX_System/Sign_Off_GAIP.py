import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

_REPO_ROOT = str(Path(__file__).resolve().parents[1])
if _REPO_ROOT not in sys.path:
    sys.path.append(_REPO_ROOT)

from PX_System.foundation.Sovereign_Log_Chain import append as slc_append

def apply_human_sign_off(trace_id, auditor_name):
    log_path = os.path.join(_REPO_ROOT, "PX_Audit", "fda_gaip_trace.log")
    if not os.path.exists(log_path):
        print("[ERROR] Audit log not found.")
        return

    # 1. Verify Trace Existence
    traces = []
    target_found = False
    with open(log_path, "r") as f:
        for line in f:
            trace = json.loads(line.strip())
            if trace["trace_id"] == trace_id:
                target_found = True
                # Only sign off if not already authorized
                if not trace["authorized"]:
                    trace["human_signoff"] = {
                        "signed_by": auditor_name,
                        "timestamp": datetime.now(timezone.utc).isoformat(),
                        "decision": "AUTHORIZED"
                    }
                    trace["authorized"] = True
                    trace["decision_rationale"] = f"Human sign-off provided by {auditor_name} per GAIP SOP-001"

                    slc_append("GAIP_HUMAN_SIGN_OFF", {
                        "trace_id": trace_id,
                        "auditor": auditor_name,
                        "decision": "AUTHORIZED",
                        "timestamp": trace["human_signoff"]["timestamp"],
                    })
            traces.append(trace)

    if not target_found:
        print(f"[ERROR] Trace {trace_id} not found in audit log.")
        return

    # 2. Persist Updated Audit Chain (Append-Only Logic)
    with open(log_path, "w") as f:
        for t in traces:
            f.write(json.dumps(t) + "\n")

    print(f">>> [GAIP SIGN-OFF COMPLETE] Trace {trace_id} is now AUTHORIZED.")

if __name__ == "__main__":
    t_id = input("Enter Trace ID for Sign-off: ").strip()
    name = input("Enter Auditor Name: ").strip()
    apply_human_sign_off(t_id, name)
