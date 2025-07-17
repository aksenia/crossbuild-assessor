"""
Clinical evidence-first discrepancy scoring configuration

REDESIGNED PRIORITIZATION FOCUS:
- Clinical significance changes (rare but critical for interpretation)
- Impact level transitions (functionally significant)
- Reduced weight for VEP consequence mismatches (common annotation noise)

Scoring weights prioritize changes that affect clinical decision-making.
"""

# Clinical evidence-first impact transition scoring
IMPACT_TRANSITION_SCORES = {
    # HIGH impact transitions (clinically critical)
    ('HIGH', 'MODERATE'): 15,    # HIGH ↔ MODERATE (critical for clinical interpretation)
    ('HIGH', 'LOW'): 18,         # HIGH ↔ LOW (very critical)
    ('HIGH', 'MODIFIER'): 20,    # HIGH ↔ MODIFIER (extremely critical)
    
    # MODERATE impact transitions (clinically significant)
    ('MODERATE', 'LOW'): 10,     # MODERATE ↔ LOW (significant)
    ('MODERATE', 'MODIFIER'): 12, # MODERATE ↔ MODIFIER (significant)
    
    # LOW impact transitions (annotation noise - minimal weight)
    ('LOW', 'MODIFIER'): 1       # LOW ↔ MODIFIER (mostly annotation differences)
}

# Clinical evidence override factors (unchanged - these work well)
CLINICAL_OVERRIDE = {
    'benign_reduction_factor': 0.1,   # 90% score reduction for benign variants
    'pathogenic_boost_factor': 2.0    # 2x score boost for pathogenic variants
}

# REDESIGNED: Clinical evidence-first scoring weights
BASE_SCORES = {
    # CRITICAL: Clinical interpretation changes (highest priority)
    'clinical_sig_pathogenic_benign_change': 20,    # PATHOGENIC↔BENIGN (critical)
    'clinical_sig_vus_to_pathogenic': 15,          # VUS→PATHOGENIC (critical)
    'clinical_sig_benign_to_pathogenic': 20,       # BENIGN→PATHOGENIC (critical)
    'clinical_sig_pathogenic_to_benign': 18,       # PATHOGENIC→BENIGN (critical)
    
    # HIGH: Functionally significant changes
    'clinical_sig_other_change': 8,                # Other clinical changes (VUS transitions, etc.)
    'same_transcript_consequence_changes': 6,      # DEMOTED: Same transcript changes
    
    # MODERATE: Prediction and annotation changes
    'sift_change': 5,                              # SIFT prediction changes
    'polyphen_change': 5,                          # PolyPhen prediction changes
    'gene_changes_high_impact': 4,                 # Gene changes (high impact context)
    'gene_changes_moderate_impact': 3,             # Gene changes (moderate impact context)
    'gene_changes_low_impact': 2,                  # Gene changes (low impact context)
    'gene_changes_minimal_weight': 0.1,            # Gene changes (not clinically relevant)
    'gene_changes_mixed_impact': 2,                # Gene changes (mixed/unknown impact)
    
    # LOW: Technical and annotation issues (reduced weights)
    'unmatched_consequences': 1,                   # DEMOTED: Unmatched consequences
    'same_consequence_different_transcripts': 0.5, # DEMOTED: Different transcripts
    'position_mismatch': 3,                        # Position issues (liftover problems)
    'genotype_mismatch': 3,                        # Genotype issues (liftover problems)
    'position_difference_moderate': 2,             # >10bp difference
    'position_difference_large': 3,                # >100bp difference
    'ref_alt_swap': 2,                            # BCFtools swap occurred
    
    # Bonuses (unchanged)
    'has_clinical_data_bonus': 2,                 # Any clinical data available
    'region_mapping_bonus': 1,                    # REGION mapping status
    'high_impact_bonus': 2                        # HIGH impact in either build
}

# REDESIGNED: Priority categories with clinical focus
PRIORITY_CATEGORIES = ['CRITICAL', 'HIGH', 'MODERATE', 'LOW', 'INVESTIGATE']

# UPDATED: Category determination rules
PRIORITY_RULES = {
    'CRITICAL': {
        'description': 'Clinical interpretation changes (PATHOGENIC↔BENIGN, VUS→PATHOGENIC) OR high impact transitions',
        'condition': 'major_clinical_significance_changes OR high_impact_transitions'
    },
    'HIGH': {
        'description': 'Functionally significant changes (impact transitions, other clinical changes, same transcript)',
        'condition': 'moderate_impact_transitions OR other_clinical_changes OR same_transcript_changes'
    },
    'MODERATE': {
        'description': 'Prediction changes and clinically relevant gene changes',
        'condition': 'pathogenicity_prediction_changes OR clinically_relevant_gene_changes'
    },
    'LOW': {
        'description': 'Technical issues and annotation differences',
        'condition': 'annotation_differences OR liftover_technical_issues'
    },
    'INVESTIGATE': {
        'description': 'Unclear cases requiring further review',
        'condition': 'insufficient_clinical_context OR mixed_evidence'
    }
}