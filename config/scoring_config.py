"""
Priority transcript-based discrepancy scoring configuration

PRIORITY TRANSCRIPT FOCUS:
- HGVS concordance on priority transcript (clinician #1 priority)
- Clinical significance changes (pathogenic<->benign critical)
- De-prioritized VUS and missing data (less clinically actionable)

Scoring weights prioritize HGVS mismatches and clinical decision-making changes.
"""

# Priority transcript-based scoring thresholds (4-category system)
PRIORITY_THRESHOLDS = {
    'CRITICAL': 100,
    'MODERATE': 60,
    'LOW': 20,
    'CONCORDANT': 0
}

# Priority categories (simplified from 5 to 4)
PRIORITY_CATEGORIES = ['CRITICAL', 'MODERATE', 'LOW', 'CONCORDANT']

# Clinician-priority scoring: HGVS mismatches first, VUS/missing data de-prioritized
BASE_SCORES = {
    # CRITICAL (100+ points) - Clinician TOP priorities
    'priority_hgvsc_mismatch': 100,           # #1 PRIORITY - HGVS mismatches
    'pathogenic_benign_flip': 90,             # P↔B, LP↔B, P↔LB (major clinical impact)
    'priority_hgvsp_mismatch': 50,            # Protein-level mismatches
    
    # MODERATE (60-99 points) - Moderate clinical concern
    'transcript_mismatch': 60,                # No priority transcript in one build
    'pathogenic_vus_change': 40,              # P↔VUS, LP↔VUS (some uncertainty)
    'serious_consequence_difference': 35,      # HIGH→MODERATE impact changes
    
    # LOW (20-59 points) - Minor issues  
    'vus_benign_change': 25,                  # VUS↔B, VUS↔LB (less actionable)
    'minor_clinical_changes': 20,             # Likely↔Definite (same direction)
    'gene_changes': 15,
    'impact_changes': 15,
    'prediction_changes': 10,
    
    # Minimal scoring for missing data
    'missing_clinical_data': 5,               # NONE/missing values (de-prioritized)
    'missing_hgvs_data': 5,                   # No HGVS available for comparison

    # Technical issues 
    'position_mismatch': 20,
    'position_difference_moderate': 10,
    'position_difference_large': 15,
    'genotype_mismatch': 20,
    'ref_alt_swap': 10
}

# Clinical significance change classification
CLINICAL_CHANGE_TYPES = {
    'pathogenic_benign_flip': [
        'PATHOGENIC_TO_BENIGN', 'BENIGN_TO_PATHOGENIC',
        'LIKELY_PATHOGENIC_TO_BENIGN', 'BENIGN_TO_LIKELY_PATHOGENIC',
        'PATHOGENIC_TO_LIKELY_BENIGN', 'LIKELY_BENIGN_TO_PATHOGENIC'
    ],
    'pathogenic_vus_change': [
        'PATHOGENIC_TO_VUS', 'VUS_TO_PATHOGENIC',
        'LIKELY_PATHOGENIC_TO_VUS', 'VUS_TO_LIKELY_PATHOGENIC'
    ],
    'vus_benign_change': [
        'VUS_TO_BENIGN', 'BENIGN_TO_VUS',
        'VUS_TO_LIKELY_BENIGN', 'LIKELY_BENIGN_TO_VUS'
    ],
    'minor_clinical_changes': [
        'PATHOGENIC_TO_LIKELY_PATHOGENIC', 'LIKELY_PATHOGENIC_TO_PATHOGENIC',
        'BENIGN_TO_LIKELY_BENIGN', 'LIKELY_BENIGN_TO_BENIGN'
    ]
}

# Priority category determination rules
PRIORITY_RULES = {
    'CRITICAL': {
        'description': 'HGVS mismatches on priority transcript OR major clinical significance changes (Pathogenic↔Benign)',
        'condition': 'priority_hgvs_mismatch OR pathogenic_benign_clinical_changes'
    },
    'MODERATE': {
        'description': 'Priority transcript unavailable OR moderate clinical changes OR serious functional differences',
        'condition': 'transcript_matching_issues OR pathogenic_vus_changes OR consequence_differences'
    },
    'LOW': {
        'description': 'Minor changes (VUS-Benign, prediction changes, gene symbol differences)',
        'condition': 'minor_clinical_changes OR prediction_changes OR annotation_differences'
    },
    'CONCORDANT': {
        'description': 'Perfect priority transcript match with identical HGVS nomenclature',
        'condition': 'same_priority_transcript AND identical_hgvs AND no_clinical_changes'
    }
}

# Legacy support - maintain backward compatibility where needed
CLINICAL_OVERRIDE = {
    'benign_reduction_factor': 0.1,   # Maintain existing benign variant handling
    'pathogenic_boost_factor': 2.0    # Maintain existing pathogenic variant handling
}