"""
CrossBuild Assessor Configuration Module

Centralized configuration for scoring, visualization, and constants.
"""

from .constants import VEP_CONSEQUENCE_IMPACT
from .scoring_config import (
    IMPACT_TRANSITION_SCORES,
    CLINICAL_OVERRIDE,
    PRIORITY_CATEGORIES,
    BASE_SCORES
)
from .visualization_config import (
    PLOT_COLORS,
    PLOT_STYLE_CONFIG,
    FIGURE_CONFIG
)

__all__ = [
    'VEP_CONSEQUENCE_IMPACT',
    'IMPACT_TRANSITION_SCORES',
    'CLINICAL_OVERRIDE',
    'PRIORITY_CATEGORIES', 
    'BASE_SCORES',
    'PLOT_COLORS',
    'PLOT_STYLE_CONFIG',
    'FIGURE_CONFIG'
]
