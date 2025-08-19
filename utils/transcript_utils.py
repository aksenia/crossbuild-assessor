"""
Transcript Utilities

Functions for processing transcript IDs, genotype extraction, and other
transcript-related data transformations.
"""

import pandas as pd


def extract_genotype_from_alleles(source_alleles):
    """Extract reference and alternative alleles from source_alleles string"""
    if pd.isna(source_alleles) or source_alleles == '':
        return '', ''
    
    # Handle various formats: "A/G", "A,G", "-/G", etc.
    if '/' in source_alleles:
        parts = source_alleles.split('/')
    elif ',' in source_alleles:
        parts = source_alleles.split(',')
    else:
        return source_alleles, ''
    
    if len(parts) >= 2:
        return parts[0].strip(), parts[1].strip()
    else:
        return parts[0].strip(), ''


def format_consequence_relationship(relationship, hg19_consequences, hg38_consequences):
    """Format consequence relationship with unified display format"""
    if not hg19_consequences and not hg38_consequences:
        return f"{relationship}: no consequence data"
    
    # Parse consequence sets
    hg19_set = set(c.strip() for c in hg19_consequences.split(',') if c.strip()) if hg19_consequences else set()
    hg38_set = set(c.strip() for c in hg38_consequences.split(',') if c.strip()) if hg38_consequences else set()
    
    if relationship == 'matched':
        return ', '.join(sorted(hg19_set))
    elif relationship == 'disjoint_consequences':
        hg19_str = ', '.join(sorted(hg19_set))
        hg38_str = ', '.join(sorted(hg38_set))
        return f"{hg19_str} → {hg38_str}"
    elif relationship == 'partial_overlap_consequences':
        shared = hg19_set & hg38_set
        unique_hg19 = hg19_set - hg38_set
        unique_hg38 = hg38_set - hg19_set
        
        shared_str = f"({', '.join(sorted(shared))})" if shared else ""
        unique_hg19_str = ', '.join(sorted(unique_hg19)) if unique_hg19 else ""
        unique_hg38_str = ', '.join(sorted(unique_hg38)) if unique_hg38 else ""
        
        parts = [p for p in [shared_str, unique_hg19_str] if p]
        left_side = ' '.join(parts)
        return f"{left_side} → {unique_hg38_str}"
    elif relationship == 'hg19_subset_of_hg38':
        shared = hg19_set & hg38_set
        unique_hg38 = hg38_set - hg19_set
        shared_str = f"({', '.join(sorted(shared))})" if shared else ""
        unique_hg38_str = ', '.join(sorted(unique_hg38)) if unique_hg38 else ""
        return f"{shared_str} → {unique_hg38_str}"
    elif relationship == 'hg38_subset_of_hg19':
        shared = hg19_set & hg38_set
        unique_hg19 = hg19_set - hg38_set
        shared_str = f"({', '.join(sorted(shared))})" if shared else ""
        unique_hg19_str = ', '.join(sorted(unique_hg19)) if unique_hg19 else ""
        parts = [p for p in [shared_str, unique_hg19_str] if p]
        left_side = ' '.join(parts)
        return f"{left_side} →"
    else:
        return f"{relationship}: {hg19_consequences} vs {hg38_consequences}"

def analyze_consequence_relationships(hg19_transcripts, hg38_transcripts):
    """
    Analyze consequence set relationships between builds
    
    Args:
        hg19_transcripts: dict of transcript_id -> transcript_data
        hg38_transcripts: dict of transcript_id -> transcript_data
        
    Returns:
        tuple: (consequence_relationship, consequence_change)
    """
    import pandas as pd
    from utils.data_utils import clean_string
    from config.constants import VEP_CONSEQUENCE_IMPACT
    
    all_hg19_consequences = set()
    all_hg38_consequences = set()

    # Build consequence sets with comma-splitting
    for data in hg19_transcripts.values():
        consequence_field = data['consequence']
        if pd.notna(consequence_field) and str(consequence_field).strip():
            individual_consequences = [c.strip() for c in str(consequence_field).split(',')]
            for cons in individual_consequences:
                cleaned_cons = clean_string(cons)
                if cleaned_cons:
                    all_hg19_consequences.add(cleaned_cons)

    for data in hg38_transcripts.values():
        consequence_field = data['consequence']
        if pd.notna(consequence_field) and str(consequence_field).strip():
            individual_consequences = [c.strip() for c in str(consequence_field).split(',')]
            for cons in individual_consequences:
                cleaned_cons = clean_string(cons)
                if cleaned_cons:
                    all_hg38_consequences.add(cleaned_cons)
    
    # Calculate relationships
    shared_consequences = all_hg19_consequences & all_hg38_consequences
    unique_hg19 = all_hg19_consequences - all_hg38_consequences  
    unique_hg38 = all_hg38_consequences - all_hg19_consequences

    # Determine relationship 
    if len(all_hg19_consequences) == 0 and len(all_hg38_consequences) == 0:
        consequence_relationship = 'no_consequences'
    elif len(shared_consequences) > 0 and len(unique_hg19) == 0 and len(unique_hg38) == 0:
        consequence_relationship = 'matched'
    elif len(shared_consequences) == 0:  # Truly disjoint
        consequence_relationship = 'disjoint_consequences'
    elif len(unique_hg19) == 0:  # hg19 subset of hg38
        consequence_relationship = 'hg19_subset_of_hg38'
    elif len(unique_hg38) == 0:  # hg38 subset of hg19 
        consequence_relationship = 'hg38_subset_of_hg19'
    else:  # Partial overlap
        consequence_relationship = 'partial_overlap_consequences'

    # Format consequence change directly
    hg19_str = ', '.join(sorted(all_hg19_consequences)) if all_hg19_consequences else ''
    hg38_str = ', '.join(sorted(all_hg38_consequences)) if all_hg38_consequences else ''
    consequence_change = format_consequence_relationship(consequence_relationship, hg19_str, hg38_str)

    return consequence_relationship, consequence_change