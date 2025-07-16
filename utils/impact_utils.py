"""
Impact Utilities

Functions for calculating impact magnitudes, numeric values, and transitions
for clinical significance assessment.
"""


def get_impact_numeric_value(impact):
    """Convert impact to numeric value for magnitude calculation"""
    impact_values = {'HIGH': 4, 'MODERATE': 3, 'LOW': 2, 'MODIFIER': 1}
    return impact_values.get(impact, 0)


def calculate_impact_transition_magnitude(hg19_impact, hg38_impact):
    """Calculate magnitude of impact transition for clinical significance"""
    hg19_val = get_impact_numeric_value(hg19_impact)
    hg38_val = get_impact_numeric_value(hg38_impact)
    magnitude = abs(hg19_val - hg38_val)
    
    # Check if transition involves MODERATE or HIGH (clinically significant)
    involves_moderate_or_high = (
        hg19_impact in ['MODERATE', 'HIGH'] or 
        hg38_impact in ['MODERATE', 'HIGH']
    )
    
    return magnitude, involves_moderate_or_high
