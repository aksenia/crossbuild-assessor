"""
CrossBuild Assessor

A comprehensive genomic variant analysis pipeline for evaluating liftover quality 
and prioritizing clinically relevant discordances between genome builds.

Core Components:
- Database Loader: Load liftover + VEP data into SQLite
- Database Analyzer: Liftover quality control analysis  
- Variant Prioritizer: Clinical discrepancy variant prioritization
"""

__version__ = "1.0.0"
__author__ = "aksenia"

# Import main configuration
from .config import (
    VEP_CONSEQUENCE_IMPACT,
    IMPACT_TRANSITION_SCORES,
    CLINICAL_OVERRIDE,
    PRIORITY_CATEGORIES,
    BASE_SCORES,
    PLOT_COLORS,
    PLOT_STYLE_CONFIG,
    FIGURE_CONFIG
)

__all__ = [
    '__version__',
    '__author__',
    'VEP_CONSEQUENCE_IMPACT',
    'IMPACT_TRANSITION_SCORES', 
    'CLINICAL_OVERRIDE',
    'PRIORITY_CATEGORIES',
    'BASE_SCORES',
    'PLOT_COLORS',
    'PLOT_STYLE_CONFIG',
    'FIGURE_CONFIG'
]
