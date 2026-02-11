from PX_Constitution.Virtual_Machine import get_vm_fingerprint
from PX_System.foundation.Sovereign_Log_Chain import append as slc_append

METADATA = {
    "role": "GAIP_Gateway",
    "version": "1.3.0-GAIP",
    "vm_fingerprint": get_vm_fingerprint()
}

class GAIPGateway:
    def __init__(self, mode="REGULATORY"):
        self.mode = mode

    def evaluate(self, csa_res, aas_res, v_state, prop_meta):
        # 1. Check for Regulatory Mode Hard-Stop
        authorized = False
        rationale = "Awaiting HITL Sign-off"

        # 2. Check for High-Risk HITL Requirement
        needs_human = csa_res.get("human_review_required", True)
        has_signoff = prop_meta.get("human_sign_off_id") is not None

        # Logic Gate: In REGULATORY mode, if it needs a human, it MUST have a signoff
        if self.mode == "REGULATORY":
            if not needs_human:
                authorized = True
                rationale = "Automated Regulatory Approval"
            elif has_signoff:
                authorized = True
                rationale = f"Authorized via HITL ID: {prop_meta.get('human_sign_off_id')}"
        else:
            # RESEARCH mode: still requires CSA pass — no unconditional bypass
            if not needs_human or has_signoff:
                authorized = True
                rationale = "Authorized via RESEARCH mode (CSA pass)"
            else:
                rationale = "RESEARCH mode: CSA requires human review, no sign-off present"

        result = {
            "authorized": authorized,
            "rationale": rationale,
            "mode": self.mode,
            "metadata": METADATA
        }

        slc_append("GAIP_GATEWAY_EVALUATION", {
            "authorized": authorized,
            "rationale": rationale,
            "mode": self.mode,
        })

        return result
