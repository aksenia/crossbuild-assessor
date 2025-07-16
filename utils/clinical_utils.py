"""
Clinical Utilities

Functions for parsing and analyzing clinical significance, pathogenicity predictions,
and other clinical evidence from variant annotations.
"""

import pandas as pd
import re


def is_pathogenic_clinical_significance(clin_sig):
    """Check if clinical significance indicates pathogenic variant"""
    if pd.isna(clin_sig) or clin_sig in ['', '-', 'nan']:
        return False
    
    clin_sig_lower = str(clin_sig).lower()
    pathogenic_terms = ['pathogenic', 'likely_pathogenic', 'drug_response']
    return any(term in clin_sig_lower for term in pathogenic_terms)


def is_benign_clinical_significance(clin_sig):
    """Check if clinical significance indicates benign variant"""
    if pd.isna(clin_sig) or clin_sig in ['', '-', 'nan']:
        return False
    
    clin_sig_lower = str(clin_sig).lower()
    benign_terms = ['benign', 'likely_benign']
    return any(term in clin_sig_lower for term in benign_terms)


def parse_sift_prediction(sift_str):
    """Parse SIFT prediction and score"""
    if pd.isna(sift_str) or sift_str in ['', '-', 'nan']:
        return None, None
    
    sift_str = str(sift_str).lower()
    
    # Extract prediction
    if 'deleterious' in sift_str:
        prediction = 'deleterious'
    elif 'tolerated' in sift_str:
        prediction = 'tolerated'
    else:
        prediction = None
    
    # Extract score (numbers in parentheses)
    score_match = re.search(r'\(([0-9.]+)\)', sift_str)
    score = float(score_match.group(1)) if score_match else None
    
    return prediction, score


def parse_polyphen_prediction(polyphen_str):
    """Parse PolyPhen prediction and score"""
    if pd.isna(polyphen_str) or polyphen_str in ['', '-', 'nan']:
        return None, None
    
    polyphen_str = str(polyphen_str).lower()
    
    # Extract prediction
    if 'probably_damaging' in polyphen_str:
        prediction = 'probably_damaging'
    elif 'possibly_damaging' in polyphen_str:
        prediction = 'possibly_damaging'
    elif 'benign' in polyphen_str:
        prediction = 'benign'
    else:
        prediction = None
    
    # Extract score (numbers in parentheses)
    score_match = re.search(r'\(([0-9.]+)\)', polyphen_str)
    score = float(score_match.group(1)) if score_match else None
    
    return prediction, score
