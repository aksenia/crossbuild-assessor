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

def analyze_worst_consequence_transcripts(hg19_transcripts, hg38_transcripts, priority_transcript_crossbuild):
    """
    Analyze worst consequence transcripts for each build
    
    Args:
        hg19_transcripts: dict of transcript_id -> transcript_data
        hg38_transcripts: dict of transcript_id -> transcript_data  
        priority_transcript_crossbuild: priority transcript ID (same for both builds)
        
    Returns:
        dict: worst consequence analysis results
    """
    from config.constants import VEP_CONSEQUENCE_IMPACT
    
    def get_consequence_severity_rank(consequence):
        """Get numeric severity rank for consequence (lower = more severe)"""
        impact = VEP_CONSEQUENCE_IMPACT.get(consequence, 'MODIFIER')
        severity_map = {'HIGH': 1, 'MODERATE': 2, 'LOW': 3, 'MODIFIER': 4}
        return severity_map.get(impact, 4)
    
    def find_worst_consequence_transcripts(transcripts):
        """Find all transcripts with worst consequence in a build"""
        if not transcripts:
            return '', []
        
        # Find the most severe consequence across all transcripts
        worst_severity = float('inf')
        worst_consequence = ''
        
        for tx_id, data in transcripts.items():
            consequence = data.get('consequence', '')
            if consequence and consequence not in ['', '-']:
                # Handle comma-separated consequences
                individual_consequences = [c.strip() for c in consequence.split(',')]
                for cons in individual_consequences:
                    severity = get_consequence_severity_rank(cons)
                    if severity < worst_severity:
                        worst_severity = severity
                        worst_consequence = cons
        
        # Find all transcripts with this worst consequence
        worst_consequence_transcripts = []
        for tx_id, data in transcripts.items():
            consequence = data.get('consequence', '')
            if consequence and consequence not in ['', '-']:
                individual_consequences = [c.strip() for c in consequence.split(',')]
                if worst_consequence in individual_consequences:
                    worst_consequence_transcripts.append(tx_id)
        
        return worst_consequence, worst_consequence_transcripts
    
    # Analyze each build
    hg19_worst_consequence, hg19_worst_tx_list = find_worst_consequence_transcripts(hg19_transcripts)
    hg38_worst_consequence, hg38_worst_tx_list = find_worst_consequence_transcripts(hg38_transcripts)
    
    # Format transcript lists
    hg19_worst_consequence_tx = ';'.join(hg19_worst_tx_list) if hg19_worst_tx_list else ''
    hg38_worst_consequence_tx = ';'.join(hg38_worst_tx_list) if hg38_worst_tx_list else ''
    
   # Check if priority transcript has worst consequence (same transcript for both builds)
    hg19_worst_consequence_tx_is_priority = 'YES' if priority_transcript_crossbuild in hg19_worst_tx_list else 'NO'
    hg38_worst_consequence_tx_is_priority = 'YES' if priority_transcript_crossbuild in hg38_worst_tx_list else 'NO'
    
    # Check if worst consequences differ between builds
    has_worst_consequence_difference = 'YES' if hg19_worst_consequence != hg38_worst_consequence else 'NO'
    
    return {
        'hg19_worst_consequence': hg19_worst_consequence,
        'hg38_worst_consequence': hg38_worst_consequence,
        'hg19_worst_consequence_tx': hg19_worst_consequence_tx,
        'hg38_worst_consequence_tx': hg38_worst_consequence_tx,
        'hg19_worst_consequence_tx_is_priority': hg19_worst_consequence_tx_is_priority,
        'hg38_worst_consequence_tx_is_priority': hg38_worst_consequence_tx_is_priority,
        'has_worst_consequence_difference': has_worst_consequence_difference
    }

def get_priority_transcript_data(transcripts_df, priority_transcript_id, column, fallback_annotations_df=None):
    """Get data from priority transcript, with fallback logic"""
    # If no priority transcript selected, return empty
    if not priority_transcript_id or priority_transcript_id == 'NONE':
        return 'NONE'
    
    # Try to get from priority transcript
    priority_rows = transcripts_df[transcripts_df['feature'] == priority_transcript_id]
    if len(priority_rows) > 0:
        value = priority_rows[column].iloc[0]
        if pd.notna(value) and value not in ['', '-']:
            return value
    
    # If priority transcript doesn't have the data, still return empty
    # (We don't want to fall back to random transcripts - stick to priority)
    return 'NONE'