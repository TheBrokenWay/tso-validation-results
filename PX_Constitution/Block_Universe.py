import numpy as np

class BlockUniverse:
    def __init__(self):
        self.manifold_dims = 35
        self.active_dims = 10  # Only measure Physics, Ethics, and Security
        self.ideal_state = np.zeros(self.manifold_dims)
        self.ideal_state[:4] = [0.1, 0.0, 35.0, 1.0]
        self.ideal_state[4:9] = [1.0, 1.0, 1.0, 1.0, 1.0] # Ethical Ideal
        self.ideal_state[9] = 1.0                         # Security Ideal

    def project_proposal(self, p_vector, csa_scores, security_score):
        block = np.zeros(self.manifold_dims)
        block[:4] = p_vector
        block[4:9] = csa_scores
        block[9] = security_score
        return block

    def calculate_coherence(self, block):
        # Only measure distance in the ACTIVE sub-space (0-9)
        # This prevents the 'Empty' 25 dimensions from adding friction
        distance = np.linalg.norm(block[:self.active_dims] - self.ideal_state[:self.active_dims])
        
        # Quadratic Resonance: 1.0 coherence if distance is 0
        return 1.0 / (1.0 + (distance ** 2))
