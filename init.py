"""
Stability-Distillation Hypothesis (SDH) Simulation Package

This package implements a computational simulation of the stability-distillation
hypothesis for the origin of life, as proposed in arXiv:2403.17072.

Author: Cheng Bi
Email: cb236@foxmail.com
"""

from .simulation import StabilityDistillationSimulation
from .stability import get_stability, calculate_stability_heuristic

# Try to import ViennaRNA-dependent functions if available
try:
    from .stability import calculate_stability_vienna
    VIENNA_AVAILABLE = True
except ImportError:
    VIENNA_AVAILABLE = False

__version__ = "1.0.0"
__all__ = [
    "StabilityDistillationSimulation",
    "get_stability",
    "calculate_stability_heuristic",
    "VIENNA_AVAILABLE"
]
