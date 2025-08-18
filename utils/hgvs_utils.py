"""
HGVS Analysis with Exact Version Matching using hgvs library
"""

import pandas as pd
import hgvs.parser
import hgvs.normalizer
from hgvs.exceptions import HGVSError, HGVSParseError, HGVSInvalidVariantError


# Initialize hgvs parser (reuse across calls for performance)
_hgvs_parser = None

def get_hgvs_parser():
    """Get singleton HGVS parser"""
    global _hgvs_parser
    if _hgvs_parser is None:
        _hgvs_parser = hgvs.parser.Parser()
    return _hgvs_parser


def normalize_hgvs_string(hgvs_string):
    """
    Normalize HGVS string using hgvs library
    
    Returns:
        tuple: (normalized_string, success_flag, error_message)
    """
    if pd.isna(hgvs_string) or hgvs_string == '' or hgvs_string == '-':
        return '', False, 'missing_data'
    
    try:
        parser = get_hgvs_parser()
        # Parse the HGVS string
        variant = parser.parse_hgvs_variant(str(hgvs_string).strip())
        
        # Convert back to normalized string representation
        normalized = str(variant)
        return normalized, True, None
        
    except (HGVSParseError, HGVSInvalidVariantError) as e:
        # Fallback to original string for comparison
        return str(hgvs_string).strip(), False, f'parse_error: {str(e)}'
    except Exception as e:
        # Any other error
        return str(hgvs_string).strip(), False, f'unknown_error: {str(e)}'


def compare_hgvsc_strings(hgvsc1, hgvsc2):
    """
    Compare two HGVSc strings using hgvs library with exact transcript version matching
    
    Returns:
        tuple: (match_result, details)
        match_result: 'concordant', 'discordant', 'parse_error', 'missing_data'
    """
    if pd.isna(hgvsc1) or pd.isna(hgvsc2):
        return 'missing_data', 'missing values'
    
    if hgvsc1 == '' or hgvsc2 == '' or hgvsc1 == '-' or hgvsc2 == '-':
        return 'missing_data', 'empty values'
    
    try:
        parser = get_hgvs_parser()
        
        # Parse both HGVS strings
        var1 = parser.parse_hgvs_variant(str(hgvsc1).strip())
        var2 = parser.parse_hgvs_variant(str(hgvsc2).strip())
        
        # Compare transcript WITH VERSION and variant change
        # var1.ac = "NM_181798.2", var1.posedit = "479C>T"
        if var1.ac == var2.ac and str(var1.posedit) == str(var2.posedit):
            return 'concordant', 'exact_match'
        else:
            return 'discordant', f'transcript_or_change_mismatch: {var1.ac}:{var1.posedit} vs {var2.ac}:{var2.posedit}'
            
    except Exception as e:
        # Fallback to string comparison if parsing fails
        if str(hgvsc1).strip() == str(hgvsc2).strip():
            return 'concordant', 'string_match'
        else:
            return 'discordant', f'parse_failed: {str(e)}'


def compare_hgvsp_strings(hgvsp1, hgvsp2):
    """
    Compare two HGVSp strings using hgvs library with exact protein version matching
    
    Returns:
        tuple: (match_result, details)
        match_result: 'concordant', 'discordant', 'parse_error', 'missing_data'
    """
    if pd.isna(hgvsp1) or pd.isna(hgvsp2):
        return 'missing_data', 'missing values'
    
    if hgvsp1 == '' or hgvsp2 == '' or hgvsp1 == '-' or hgvsp2 == '-':
        return 'missing_data', 'empty values'
    
    try:
        parser = get_hgvs_parser()
        
        # Parse both HGVSp strings
        var1 = parser.parse_hgvs_variant(str(hgvsp1).strip())
        var2 = parser.parse_hgvs_variant(str(hgvsp2).strip())
        
        # Compare protein WITH VERSION and variant change
        # var1.ac = "NP_861463.1", var1.posedit = "Ala160Val"
        if var1.ac == var2.ac and str(var1.posedit) == str(var2.posedit):
            return 'concordant', 'exact_match'
        else:
            return 'discordant', f'protein_or_change_mismatch: {var1.ac}:{var1.posedit} vs {var2.ac}:{var2.posedit}'
            
    except Exception as e:
        # Fallback to string comparison if parsing fails
        if str(hgvsp1).strip() == str(hgvsp2).strip():
            return 'concordant', 'string_match'
        else:
            return 'discordant', f'parse_failed: {str(e)}'


def compare_hgvs_strings(hgvs1, hgvs2):
    """
    Universal HGVS comparison function with exact version matching
    Automatically detects HGVSc vs HGVSp and uses appropriate comparison
    
    Returns:
        bool: True if HGVS strings represent the exact same change with same transcript/protein version
    """
    if pd.isna(hgvs1) or pd.isna(hgvs2):
        return False
    
    hgvs1_str = str(hgvs1).strip()
    hgvs2_str = str(hgvs2).strip()
    
    if not hgvs1_str or not hgvs2_str or hgvs1_str in ['', '-'] or hgvs2_str in ['', '-']:
        return False
    
    # Detect HGVSc (coding)
    if any(x in hgvs1_str.lower() for x in ['c.', 'g.', 'n.']) or any(x in hgvs2_str.lower() for x in ['c.', 'g.', 'n.']):
        result, _ = compare_hgvsc_strings(hgvs1, hgvs2)
        return result == 'concordant'
    
    # Detect HGVSp (protein)
    elif any(x in hgvs1_str.lower() for x in ['p.']) or any(x in hgvs2_str.lower() for x in ['p.']):
        result, _ = compare_hgvsp_strings(hgvs1, hgvs2)
        return result == 'concordant'
    
    else:
        # Fallback to exact string comparison
        return hgvs1_str == hgvs2_str


def extract_transcript_from_hgvs(hgvs_string):
    """
    Extract transcript ID from HGVS string with exact version preservation
    
    Returns:
        str: Full transcript ID with version (e.g., "NM_181798.2") or None if extraction fails
    """
    if pd.isna(hgvs_string) or hgvs_string == '' or hgvs_string == '-':
        return None
    
    try:
        parser = get_hgvs_parser()
        variant = parser.parse_hgvs_variant(str(hgvs_string).strip())
        
        if hasattr(variant, 'ac') and variant.ac:
            return variant.ac  # Full transcript ID with version
        
    except Exception:
        # Fallback: try to extract from string pattern
        hgvs_str = str(hgvs_string).strip()
        if ':' in hgvs_str:
            transcript_part = hgvs_str.split(':')[0]
            return transcript_part
    
    return None

def analyze_priority_transcript_hgvs(hg19_transcripts_df, hg38_transcripts_df, transcript_crossbuild_status, priority_transcript_crossbuild):
    """
    Analyze HGVS concordance for the selected priority transcript using HGVS package
    
    Args:
        hg19_transcripts_df: hg19 transcript annotations
        hg38_transcripts_df: hg38 transcript annotations
        transcript_crossbuild_status: Status of transcript matching
        priority_transcript_crossbuild: Selected transcript ID with exact version
        
    Returns:
        dict: HGVS analysis results
    """
    # Skip analysis if no matching transcript
    if transcript_crossbuild_status in ["MANE_hg38_Only", "No_Matching_Transcripts", "No_Transcripts"] or priority_transcript_crossbuild == "NONE":
        return {
            'priority_hgvsc_hg19': '',
            'priority_hgvsc_hg38': '',
            'priority_hgvsp_hg19': '',
            'priority_hgvsp_hg38': '',
            'priority_hgvsc_concordance': 'No_Analysis',
            'priority_hgvsp_concordance': 'No_Analysis'
        }
    
    # Extract HGVS data from the exact priority transcript (with version)
    hg19_hgvsc = ''
    hg19_hgvsp = ''
    hg38_hgvsc = ''
    hg38_hgvsp = ''
    
    # Find HGVS data for priority transcript in hg19 (exact version match)
    for _, row in hg19_transcripts_df.iterrows():
        if row['feature'] == priority_transcript_crossbuild:
            hg19_hgvsc = row.get('hgvsc', '') or ''
            hg19_hgvsp = row.get('hgvsp', '') or ''
            break
    
    # Find HGVS data for priority transcript in hg38 (exact version match)
    for _, row in hg38_transcripts_df.iterrows():
        if row['feature'] == priority_transcript_crossbuild:
            hg38_hgvsc = row.get('hgvsc', '') or ''
            hg38_hgvsp = row.get('hgvsp', '') or ''
            break
    
    # Use HGVS package for exact version comparison
    from utils.hgvs_utils import compare_hgvs_strings
    
    # Compare HGVSc with exact version matching
    if hg19_hgvsc and hg38_hgvsc and hg19_hgvsc not in ['', '-'] and hg38_hgvsc not in ['', '-']:
        hgvsc_match = compare_hgvs_strings(hg19_hgvsc, hg38_hgvsc)
        hgvsc_concordance = "Match" if hgvsc_match else "Mismatch"
    elif (hg19_hgvsc and hg19_hgvsc not in ['', '-']) or (hg38_hgvsc and hg38_hgvsc not in ['', '-']):
        # One has data, other doesn't - this is a mismatch
        hgvsc_concordance = "Mismatch"
    else:
        # Neither has HGVSc data
        hgvsc_concordance = "No_Analysis"
    
    # Compare HGVSp with exact version matching
    if hg19_hgvsp and hg38_hgvsp and hg19_hgvsp not in ['', '-'] and hg38_hgvsp not in ['', '-']:
        hgvsp_match = compare_hgvs_strings(hg19_hgvsp, hg38_hgvsp)
        hgvsp_concordance = "Match" if hgvsp_match else "Mismatch"
    elif (hg19_hgvsp and hg19_hgvsp not in ['', '-']) or (hg38_hgvsp and hg38_hgvsp not in ['', '-']):
        # One has data, other doesn't - this is a mismatch
        hgvsp_concordance = "Mismatch"
    else:
        # Neither has HGVSp data
        hgvsp_concordance = "No_Analysis"
    
    return {
        'priority_hgvsc_hg19': hg19_hgvsc,
        'priority_hgvsc_hg38': hg38_hgvsc,
        'priority_hgvsp_hg19': hg19_hgvsp,
        'priority_hgvsp_hg38': hg38_hgvsp,
        'priority_hgvsc_concordance': hgvsc_concordance,
        'priority_hgvsp_concordance': hgvsp_concordance
    }