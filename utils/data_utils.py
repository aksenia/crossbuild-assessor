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