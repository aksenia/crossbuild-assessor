"""
Visualization configuration

Colors, styles, and layout settings for all plots and visualizations.
Colorblind-friendly palettes optimized for clinical review.
"""

# Colorblind-friendly color palettes
PLOT_COLORS = {
    # Priority categories (4 colors)
    'priority': ['#d62728', '#ff7f0e', '#2ca02c', '#1f77b4'],
    
    # Clinical evidence (7 colors)  
    'clinical': ['#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf', '#ff7f0e', '#2ca02c'],
    
    # Discordance types (extended palette)
    'discordance': ['#d62728', '#ff7f0e', '#2ca02c', '#1f77b4', '#9467bd', 
                   '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22'],
    
    # Clinical transitions
    'transitions': ['#d62728', '#ff7f0e', '#2ca02c', '#1f77b4', '#9467bd', 
                   '#8c564b', '#e377c2']
}

# Plot style configuration
PLOT_STYLE_CONFIG = {
    'style': 'default',  # matplotlib style
    'seaborn_palette': 'husl',  # seaborn palette
    'alpha': 0.8,  # transparency for bars
    'grid_alpha': 0.3,  # grid transparency
    'grid_axis': 'y'  # which axis to show grid
}

# Figure layout configuration
FIGURE_CONFIG = {
    # Main 4-plot figure
    'main_figure': {
        'figsize': (20, 12),
        'rows': 2,
        'cols': 2,
        'dpi': 300,
        'bbox_inches': 'tight'
    },
    
    # Individual plot settings
    'subplot_config': {
        'fontsize_title': 12,
        'fontsize_label': 10,
        'fontsize_tick': 9,
        'fontsize_legend': 8,
        'fontweight_title': 'bold',
        'fontweight_label': 'bold',
        'rotation_xticklabels': 45,
        'ha_xticklabels': 'right'
    },
    
    # Main title
    'main_title': {
        'text': 'Enhanced Variant Prioritization Analysis',
        'fontsize': 16,
        'fontweight': 'bold'
    }
}

# Plot-specific configurations
PLOT_CONFIGS = {
    'priority_categories': {
        'title': 'Variant Priority Categories\n(Clinical Review Levels)',
        'xlabel': '',
        'ylabel': 'Variant Count',
        'show_percentages': True,
        'show_counts': True
    },
    
    'clinical_evidence': {
        'title': 'Priority Categories by Clinical Evidence\n(Intelligent Gap Filling)',
        'xlabel': 'Priority Category',
        'ylabel': 'Variant Count',
        'stacked': True,
        'legend_position': (1.05, 1),
        'legend_location': 'upper left'
    },
    
    'discordance_types': {
        'title': 'Variants by Primary Discordance Type\n(Aggregated Categories)',
        'xlabel': 'Primary Discordance Type',
        'ylabel': 'Variant Count',
        'show_counts': True
    },
    
    'clinical_transitions': {
        'title': 'Clinical Evidence Transitions\n(Dynamic vs Static Evidence)',
        'xlabel': 'Clinical Evidence Transition',
        'ylabel': 'Variant Count',
        'show_counts': True
    }
}

# Clinical evidence category order (for consistent plotting)
CLINICAL_EVIDENCE_ORDER = [
    'Clinical Sig Change',
    'Pathogenic', 
    'Benign',
    'Prediction Change',
    'Pred. Pathogenic',
    'Pred. Benign',
    'High Impact',
    'No Evidence'
]

# Discordance type order (by clinical importance)
DISCORDANCE_TYPE_ORDER = [
    'Same Transcript\nConsequence Change',
    'Clinical Significance\nChange', 
    'Impact Level\nChange',
    'Gene Annotation\nChange',
    'Pathogenicity\nPrediction Change',
    'Unmatched\nConsequences',
    'Position\nMismatch',
    'Genotype\nMismatch',
    'Other/Multiple'
]

# Clinical transition order (by clinical concern level)
CLINICAL_TRANSITION_ORDER = [
    'Pathogenic → Benign',
    'Benign → Pathogenic', 
    'Stable Pathogenic',
    'Prediction Changes',
    'High Impact Only',
    'Stable Benign',
    'Other/Unknown'
]

# Update priority category order
PRIORITY_CATEGORY_ORDER = ['CRITICAL', 'HIGH', 'MODERATE', 'LOW', 'CONCORDANT']  
