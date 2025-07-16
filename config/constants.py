"""
VEP consequence impact mapping

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

# Clinical significance terms for pathogenic/benign classification
PATHOGENIC_TERMS = ['pathogenic', 'likely_pathogenic', 'drug_response']
BENIGN_TERMS = ['benign', 'likely_benign']
