"""
PX_Engine Package
Vector physics, metabolism, and operational engines
"""

# Core physics engines
from .Vector_Core import VectorCore
from .Metabolism import Metabolism
from .Trajectory_Predictor import TrajectoryPredictor
from .Genesis_Engine import generate_novel_candidates

__all__ = ['VectorCore', 'Metabolism', 'TrajectoryPredictor', 'generate_novel_candidates']
