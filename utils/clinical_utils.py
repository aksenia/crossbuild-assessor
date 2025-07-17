"""
Clinical Utilities

Functions for parsing and analyzing clinical significance, pathogenicity predictions,
and other clinical evidence from variant annotations.
"""

import pandas as pd
import re
from config.constants import CLINICAL_SIGNIFICANCE_NORMALIZATION, PATHOGENIC_TERMS, BENIGN_TERMS


def normalize_clinical_significance(clin_sig):
    """
    Normalize clinical significance to standardized categories with proper handling of mixed categories
    
    Rules:
    - PATHOGENIC: All components must be pathogenic-related (no conflicting/VUS)
    - BENIGN: All components must be benign-related (no conflicting/pathogenic/VUS)
    - VUS: Any conflicting evidence OR mixed categories OR uncertain significance
    - Others: RISK, DRUG_RESPONSE, PROTECTIVE, OTHER, NONE as before
    
    Args:
        clin_sig: Raw clinical significance string
        
    Returns:
        str: Normalized category (PATHOGENIC, BENIGN, VUS, RISK, DRUG_RESPONSE, PROTECTIVE, OTHER, NONE)
    """
    if pd.isna(clin_sig) or clin_sig in ['', '-', 'nan']:
        return 'NONE'
    
    # Clean and normalize the input
    clin_sig_clean = str(clin_sig).lower().strip()
    
    # Direct mapping for simple cases
    if clin_sig_clean in CLINICAL_SIGNIFICANCE_NORMALIZATION:
        return CLINICAL_SIGNIFICANCE_NORMALIZATION[clin_sig_clean]
    
    # Handle complex cases with multiple categories (comma or slash separated)
    # Split by common separators
    components = []
    for separator in [',', '/', '|', ';']:
        if separator in clin_sig_clean:
            components = [comp.strip() for comp in clin_sig_clean.split(separator)]
            break
    
    if not components:
        components = [clin_sig_clean]
    
    # Categorize each component
    pathogenic_components = []
    benign_components = []
    vus_components = []
    risk_components = []
    drug_components = []
    protective_components = []
    other_components = []
    none_components = []
    
    for comp in components:
        comp = comp.strip()
        if not comp:
            continue
            
        # Check for VUS/conflicting indicators first (these override everything)
        if any(term in comp for term in ['uncertain', 'conflicting', 'vus', 'low_penetrance']):
            vus_components.append(comp)
        # Check for pathogenic terms
        elif any(term in comp for term in ['pathogenic']) and not any(term in comp for term in ['benign']):
            pathogenic_components.append(comp)
        # Check for benign terms  
        elif any(term in comp for term in ['benign']) and not any(term in comp for term in ['pathogenic']):
            benign_components.append(comp)
        # Check for risk terms
        elif any(term in comp for term in ['risk']):
            risk_components.append(comp)
        # Check for drug response
        elif 'drug' in comp:
            drug_components.append(comp)
        # Check for protective
        elif 'protective' in comp:
            protective_components.append(comp)
        # Check for association
        elif 'association' in comp:
            other_components.append(comp)
        # Check for none/missing
        elif comp in ['not_provided', 'other', 'no_classifications_from_unflagged_records', '', '-']:
            none_components.append(comp)
        else:
            # Unknown component - treat as VUS for safety
            vus_components.append(comp)
    
    # Apply decision rules
    
    # If ANY VUS/conflicting components exist → VUS
    if vus_components:
        return 'VUS'
    
    # If mixed pathogenic and benign → VUS (conflicting evidence)
    if pathogenic_components and benign_components:
        return 'VUS'
    
    # Pure pathogenic (no benign, no VUS)
    if pathogenic_components and not benign_components and not vus_components:
        return 'PATHOGENIC'
    
    # Pure benign (no pathogenic, no VUS)  
    if benign_components and not pathogenic_components and not vus_components:
        return 'BENIGN'
    
    # Single category dominance
    if risk_components and not (pathogenic_components or benign_components or vus_components):
        return 'RISK'
    
    if drug_components and not (pathogenic_components or benign_components or vus_components):
        return 'DRUG_RESPONSE'
    
    if protective_components and not (pathogenic_components or benign_components or vus_components):
        return 'PROTECTIVE'
    
    if other_components and not (pathogenic_components or benign_components or vus_components):
        return 'OTHER'
    
    if none_components and not (pathogenic_components or benign_components or vus_components):
        return 'NONE'
    
    # Mixed categories with non-clinical components → VUS for safety
    if (pathogenic_components or benign_components) and (risk_components or drug_components or protective_components or other_components):
        return 'VUS'
    
    # Default fallback
    return 'VUS'


def is_pathogenic_clinical_significance(clin_sig):
    """Check if clinical significance indicates pathogenic variant (legacy function)"""
    if pd.isna(clin_sig) or clin_sig in ['', '-', 'nan']:
        return False
    
    # Use normalized category for consistency
    normalized = normalize_clinical_significance(clin_sig)
    return normalized == 'PATHOGENIC'


def is_benign_clinical_significance(clin_sig):
    """Check if clinical significance indicates benign variant (legacy function)"""
    if pd.isna(clin_sig) or clin_sig in ['', '-', 'nan']:
        return False
    
    # Use normalized category for consistency
    normalized = normalize_clinical_significance(clin_sig)
    return normalized == 'BENIGN'


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