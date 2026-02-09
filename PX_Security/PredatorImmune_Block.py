import numpy as np
from PX_Constitution.Block_Universe import BlockUniverse

class PredatorImmuneBlock:
    def __init__(self):
        self.manifold = BlockUniverse()
        self.auth_threshold = 0.85

    def handle_block_request(self, task_id, context, p_vec, csa_scores, security_score=1.0):
        """
        p_vec: 4-D physics vector; csa_scores: 5-D CSA; security_score: scalar.
        Builds 35-D block via BlockUniverse.project_proposal for coherence.
        """
        import numpy as np
        p_vec = np.asarray(p_vec, dtype=float).flatten()[:4]
        csa_s = np.asarray(csa_scores, dtype=float).flatten()[:5]
        block = self.manifold.project_proposal(p_vec, csa_s, float(security_score))
        coherence = self.manifold.calculate_coherence(block)
        authorized = coherence >= self.auth_threshold
        classification = "trusted" if "PX-ROOT" in context.get("fingerprint", "") else "unknown"
        route = "BLUE" if authorized and classification == "trusted" else "RED"
        return {
            "authorized": authorized,
            "coherence": round(float(coherence), 6),
            "classification": classification,
            "route": route,
            "block": block.tolist(),
        }
