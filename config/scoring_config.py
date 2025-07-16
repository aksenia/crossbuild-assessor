"""
Clinical evidence-driven discrepancy scoring configuration

Scoring weights and rules for variant prioritization based on:
- Clinical evidence (pathogenic/benign)
- Impact magnitude transitions
- Functional consequences
- Annotation discordances
"""

# Magnitude-based impact transition scoring (clinically significant)
IMPACT_TRANSITION_SCORES = {
    # HIGH impact transitions (most clinically significant)
    ('HIGH', 'MODERATE'): 10,    # HIGH ↔ MODERATE
    ('HIGH', 'LOW'): 12,         # HIGH ↔ LOW  
    ('HIGH', 'MODIFIER'): 15,    # HIGH ↔ MODIFIER (most significant)
    
    # MODERATE impact transitions
    ('MODERATE', 'LOW'): 8,      # MODERATE ↔ LOW
    ('MODERATE', 'MODIFIER'): 10, # MODERATE ↔ MODIFIER
    
    # LOW impact transitions (annotation noise - minimal weight)
    ('LOW', 'MODIFIER'): 1       # LOW ↔ MODIFIER (reduced from previous versions)
}

# Clinical evidence override factors
CLINICAL_OVERRIDE = {
    'benign_reduction_factor': 0.1,   # 90% score reduction for benign variants
    'pathogenic_boost_factor': 2.0    # 2x score boost for pathogenic variants
}

# Base scoring weights for different discordance types
BASE_SCORES = {
    # Critical functional changes (highest priority)
    'same_transcript_consequence_changes': 10,  # CRITICAL
    
    # Clinical evidence changes
    'clinical_sig_benign_to_pathogenic': 10,
    'clinical_sig_pathogenic_to_benign': 8,
    'clinical_sig_other_change': 5,
    
    # Pathogenicity prediction changes
    'sift_change': 5,
    'polyphen_change': 5,
    
    # Impact-based gene changes (conditional on clinical significance)
    'gene_changes_high_impact': 8,
    'gene_changes_moderate_impact': 4,
    'gene_changes_low_impact': 2,
    'gene_changes_minimal_weight': 0.1,  # When not clinically relevant
    'gene_changes_mixed_impact': 3,      # Mixed or unknown impact
    
    # Other discordance types
    'unmatched_consequences': 4,                    # INVESTIGATE
    'same_consequence_different_transcripts': 3,    # Moderate priority
    'position_mismatch': 3,
    'genotype_mismatch': 3,
    'position_difference_moderate': 2,              # >10bp difference
    'position_difference_large': 3,                 # >100bp difference
    'ref_alt_swap': 2,                             # BCFtools swap occurred
    'has_clinical_data_bonus': 2,                  # Any clinical data available
    'region_mapping_bonus': 1,                     # REGION mapping status
    'high_impact_bonus': 2                         # HIGH impact in either build
}

# Priority categories in order of clinical importance
PRIORITY_CATEGORIES = ['CRITICAL', 'HIGH', 'MODERATE', 'INVESTIGATE', 'LOW']

# Category determination rules
PRIORITY_RULES = {
    'CRITICAL': {
        'description': 'Same transcript, different consequences',
        'condition': 'same_transcript_consequence_changes > 0'
    },
    'HIGH': {
        'description': 'Clinically significant impact transitions OR clinical significance changes',
        'condition': 'clinically_significant_impact_changes OR major_clinical_significance_changes'
    },
    'MODERATE': {
        'description': 'Gene changes with clinical significance OR prediction changes',
        'condition': 'clinically_relevant_gene_changes OR pathogenicity_prediction_changes'
    },
    'INVESTIGATE': {
        'description': 'Unmatched consequences OR non-significant transitions',
        'condition': 'unmatched_consequences OR minor_impact_changes'
    },
    'LOW': {
        'description': 'Benign variants',
        'condition': 'is_benign_variant'
    }
}
