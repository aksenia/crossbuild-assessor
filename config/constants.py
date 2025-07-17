"""
VEP Consequence Impact Mappings and Clinical Significance Normalization

Based on official Ensembl VEP consequence severity hierarchy.
Maps consequence terms to impact levels for clinical prioritization.
"""

# VEP consequence severity mapping based on official Ensembl hierarchy
VEP_CONSEQUENCE_IMPACT = {
    # HIGH IMPACT - Critical functional consequences
    'transcript_ablation': 'HIGH',
    'splice_acceptor_variant': 'HIGH',
    'splice_donor_variant': 'HIGH',
    'stop_gained': 'HIGH',
    'frameshift_variant': 'HIGH',
    'stop_lost': 'HIGH',
    'start_lost': 'HIGH',
    'transcript_amplification': 'HIGH',
    'feature_elongation': 'HIGH',
    'feature_truncation': 'HIGH',
    
    # MODERATE IMPACT - Significant functional consequences
    'inframe_insertion': 'MODERATE',
    'inframe_deletion': 'MODERATE',
    'missense_variant': 'MODERATE',
    'protein_altering_variant': 'MODERATE',
    
    # LOW IMPACT - Minor functional consequences
    'splice_donor_5th_base_variant': 'LOW',
    'splice_region_variant': 'LOW',
    'splice_donor_region_variant': 'LOW',
    'splice_polypyrimidine_tract_variant': 'LOW',
    'incomplete_terminal_codon_variant': 'LOW',
    'start_retained_variant': 'LOW',
    'stop_retained_variant': 'LOW',
    'synonymous_variant': 'LOW',
    
    # MODIFIER IMPACT - Minimal or unknown functional consequences
    'coding_sequence_variant': 'MODIFIER',
    'mature_miRNA_variant': 'MODIFIER',
    '5_prime_UTR_variant': 'MODIFIER',
    '3_prime_UTR_variant': 'MODIFIER',
    'non_coding_transcript_exon_variant': 'MODIFIER',
    'intron_variant': 'MODIFIER',
    'NMD_transcript_variant': 'MODIFIER',
    'non_coding_transcript_variant': 'MODIFIER',
    'coding_transcript_variant': 'MODIFIER',
    'upstream_gene_variant': 'MODIFIER',
    'downstream_gene_variant': 'MODIFIER',
    'TFBS_ablation': 'MODIFIER',
    'TFBS_amplification': 'MODIFIER',
    'TF_binding_site_variant': 'MODIFIER',
    'regulatory_region_ablation': 'MODIFIER',
    'regulatory_region_amplification': 'MODIFIER',
    'regulatory_region_variant': 'MODIFIER',
    'intergenic_variant': 'MODIFIER',
    'sequence_variant': 'MODIFIER'
}

# Impact level numeric values for magnitude calculations
IMPACT_NUMERIC_VALUES = {
    'HIGH': 4,
    'MODERATE': 3, 
    'LOW': 2,
    'MODIFIER': 1
}

# Clinical significance normalization mapping (from 1.5M dataset analysis)
CLINICAL_SIGNIFICANCE_NORMALIZATION = {
    # PATHOGENIC (includes likely_pathogenic)
    'pathogenic': 'PATHOGENIC',
    'likely_pathogenic': 'PATHOGENIC', 
    'pathogenic/likely_pathogenic': 'PATHOGENIC',
    'pathogenic/likely_pathogenic/pathogenic': 'PATHOGENIC',
    'pathogenic/pathogenic': 'PATHOGENIC',
    
    # BENIGN (includes likely_benign)
    'benign': 'BENIGN',
    'likely_benign': 'BENIGN',
    'benign/likely_benign': 'BENIGN',
    
    # VUS (uncertain significance + conflicting + low penetrance)
    'uncertain_significance': 'VUS',
    'conflicting_interpretations_of_pathogenicity': 'VUS',
    'uncertain_significance/uncertain_risk_allele': 'VUS',
    'low_penetrance': 'VUS',
    
    # RISK (risk alleles and factors)
    'established_risk_allele': 'RISK',
    'likely_risk_allele': 'RISK',
    'uncertain_risk_allele': 'RISK',
    'risk_factor': 'RISK',
    
    # DRUG_RESPONSE (pharmacogenomics)
    'drug_response': 'DRUG_RESPONSE',
    
    # PROTECTIVE (protective variants)
    'protective': 'PROTECTIVE',
    
    # OTHER (association studies)
    'association': 'OTHER',
    
    # NONE (missing/no data)
    'not_provided': 'NONE',
    'other': 'NONE',
    'no_classifications_from_unflagged_records': 'NONE',
    '': 'NONE',
    '-': 'NONE',
    'nan': 'NONE'
}

# Normalized clinical significance categories (for consistency)
NORMALIZED_CLINICAL_CATEGORIES = [
    'PATHOGENIC', 'BENIGN', 'VUS', 'RISK', 
    'DRUG_RESPONSE', 'PROTECTIVE', 'OTHER', 'NONE'
]

# Legacy pathogenic and benign terms (for backward compatibility)
PATHOGENIC_TERMS = ['pathogenic', 'likely_pathogenic', 'drug_response']
BENIGN_TERMS = ['benign', 'likely_benign']