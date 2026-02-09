from PX_Constitution.Virtual_Machine import get_vm_fingerprint

METADATA = {
    "role": "ByzantiumCouncil",
    "version": "1.2.0-GAIP",
    "vm_fingerprint": get_vm_fingerprint()
}

class ByzantiumCouncil:
    def __init__(self):
        self.min_quorum = 4

    def decide(self, vector_res, csa_res, aas_res, gaip_res):
        # Consensus logic: All four pillars must return positive
        votes = [
            vector_res.get("authorized", False),
            csa_res.get("status") == "COHERENT",
            aas_res.get("status") == "SUCCESS",
            gaip_res.get("authorized", False)
        ]
        
        consensus = all(votes)
        
        return {
            "authorized": consensus,
            "quorum_score": f"{sum(votes)}/{self.min_quorum}",
            "pillar_states": {
                "physics": votes[0],
                "ethics": votes[1],
                "invariants": votes[2],
                "regulatory": votes[3]
            },
            "metadata": METADATA
        }
