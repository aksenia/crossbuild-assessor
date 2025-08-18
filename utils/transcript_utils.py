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
