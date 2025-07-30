"""
Data Utilities

Functions for data transformation, type conversion, and other general
data processing tasks.
"""

import pandas as pd


def safe_int_convert(series):
    """Convert float to int, handling NaN values"""
    return series.apply(lambda x: int(x) if pd.notna(x) and x != '' else '')

def clean_string(s):
    """
    Clean strings for robust genomics data comparison
    
    Handles whitespace, case differences, and null variations common in VEP/clinical data
    """
    if pd.isna(s) or s in ['', '-', 'nan', 'None']:
        return ''
    return str(s).strip().lower()

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