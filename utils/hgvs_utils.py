"""
HGVSc Analysis with Canonical Transcript Priority using hgvs library
"""

import pandas as pd
import hgvs.parser
import hgvs.normalizer
from hgvs.exceptions import HGVSError, HGVSParseError, HGVSInvalidVariantError
from utils.transcript_utils import normalize_transcript_id


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
    Compare two HGVSc strings using hgvs library normalization
    
    Returns:
        tuple: (match_result, details)
        match_result: 'concordant', 'discordant', 'parse_error', 'missing_data'
        details: str with comparison details
    """
    if pd.isna(hgvsc1) or pd.isna(hgvsc2):
        return 'missing_data', 'One or both HGVSc values missing'
    
    if hgvsc1 == '' or hgvsc2 == '' or hgvsc1 == '-' or hgvsc2 == '-':
        return 'missing_data', 'One or both HGVSc values empty'
    
    # Normalize both strings
    norm1, success1, error1 = normalize_hgvs_string(hgvsc1)
    norm2, success2, error2 = normalize_hgvs_string(hgvsc2)
    
    # If either failed to parse, fall back to string comparison
    if not success1 or not success2:
        simple_match = str(hgvsc1).strip() == str(hgvsc2).strip()
        if simple_match:
            return 'concordant', f'String match (parse issues: {error1 or error2})'
        else:
            return 'parse_error', f'Parse failed: {error1}, {error2}'
    
    # Compare normalized strings
    if norm1 == norm2:
        return 'concordant', f'Normalized match: {norm1}'
    else:
        return 'discordant', f'Normalized mismatch: {norm1} vs {norm2}'


def extract_transcript_from_hgvs(hgvs_string):
    """
    Extract transcript ID from HGVSc string
    
    Returns:
        tuple: (transcript_id, version_stripped_id)
    """
    if pd.isna(hgvs_string) or hgvs_string == '' or hgvs_string == '-':
        return None, None
    
    try:
        parser = get_hgvs_parser()
        variant = parser.parse_hgvs_variant(str(hgvs_string).strip())
        
        if hasattr(variant, 'ac') and variant.ac:
            full_transcript = variant.ac
            base_transcript = normalize_transcript_id(full_transcript)
            return full_transcript, base_transcript
        
    except Exception:
        # Fallback: try to extract from string pattern
        hgvs_str = str(hgvs_string).strip()
        if ':' in hgvs_str:
            transcript_part = hgvs_str.split(':')[0]
            base_transcript = normalize_transcript_id(transcript_part)
            return transcript_part, base_transcript
    
    return None, None


def find_canonical_transcript(tx_data_dict):
    """Find canonical transcript from transcript data"""
    if not tx_data_dict:
        return None
    
    # Priority: CANONICAL == 1 > longest feature_id > first transcript
    canonical_txs = [tx_id for tx_id, data in tx_data_dict.items() if data.get('is_canonical', 0) == 1]
    
    if canonical_txs:
        return canonical_txs[0]  # Return first canonical
    
    # Fallback: longest transcript ID
    if tx_data_dict:
        return max(tx_data_dict.keys(), key=lambda x: len(tx_data_dict[x].get('feature_id', '')))
    
    return None


def analyze_transcript_hgvsc_matches(hg19_tx_data, hg38_tx_data):
    """
    Compare HGVSc for matched transcript pairs using hgvs library
    Enhanced with canonical + non-empty HGVSc + matched transcript analysis
    
    Returns: dict with detailed HGVSc analysis results
    """
    perfect_matches = []
    mismatches = []
    parse_errors = []
    
    # Find matched transcript pairs (same base ID)
    common_transcripts = set(hg19_tx_data.keys()) & set(hg38_tx_data.keys())
    
    for tx_id in common_transcripts:
        hg19_data = hg19_tx_data[tx_id]
        hg38_data = hg38_tx_data[tx_id]
        
        hg19_hgvsc = hg19_data.get('hgvsc', '')
        hg38_hgvsc = hg38_data.get('hgvsc', '')
        
        match_result, details = compare_hgvsc_strings(hg19_hgvsc, hg38_hgvsc)
        
        if match_result == 'concordant':
            perfect_matches.append((tx_id, hg19_hgvsc, hg38_hgvsc, details))
        elif match_result == 'discordant':
            mismatches.append((tx_id, hg19_hgvsc, hg38_hgvsc, details))
        else:  # parse_error or missing_data
            parse_errors.append((tx_id, hg19_hgvsc, hg38_hgvsc, details))
    
    # Find canonical transcript in each build
    canonical_tx_hg19 = find_canonical_transcript(hg19_tx_data)
    canonical_tx_hg38 = find_canonical_transcript(hg38_tx_data)
    
    # ENHANCED: Only include transcripts WITH HGVSc values in HGVSc_rest
    def get_canonical_and_rest_hgvsc(tx_data, canonical_tx):
        canonical_hgvsc = ''
        rest_hgvsc_list = []
        
        for tx_id, data in tx_data.items():
            hgvsc = data.get('hgvsc', '').strip()
            # Only process if HGVSc exists and is not empty/missing
            if hgvsc and hgvsc != '-' and hgvsc != '':
                # Check if HGVSc already contains transcript ID
                if ':' in hgvsc:
                    # HGVSc already has transcript:consequence format
                    hgvsc_display = hgvsc
                else:
                    # HGVSc is just the consequence part, add transcript
                    hgvsc_display = f"{data['feature_id']}:{hgvsc}"
                
                if tx_id == canonical_tx:
                    canonical_hgvsc = hgvsc_display
                else:
                    rest_hgvsc_list.append(hgvsc_display)
        
        rest_hgvsc = '; '.join(rest_hgvsc_list) if rest_hgvsc_list else ''
        return canonical_hgvsc, rest_hgvsc
    
    # Get canonical and rest HGVSc (only those with actual HGVSc values)
    hg19_canonical_hgvsc, hg19_rest_hgvsc = get_canonical_and_rest_hgvsc(hg19_tx_data, canonical_tx_hg19)
    hg38_canonical_hgvsc, hg38_rest_hgvsc = get_canonical_and_rest_hgvsc(hg38_tx_data, canonical_tx_hg38)
    
    # Get HGVSc from ACTUAL canonical transcripts in each build (for legacy compatibility)
    hg19_canonical_hgvsc_legacy = ''
    hg38_canonical_hgvsc_legacy = ''
    hg19_canonical_hgvsp = ''
    hg38_canonical_hgvsp = ''
    
    if canonical_tx_hg19 and canonical_tx_hg19 in hg19_tx_data:
        hg19_canonical_hgvsc_legacy = hg19_tx_data[canonical_tx_hg19].get('hgvsc', '')
        hg19_canonical_hgvsp = hg19_tx_data[canonical_tx_hg19].get('hgvsp', '')
    
    if canonical_tx_hg38 and canonical_tx_hg38 in hg38_tx_data:
        hg38_canonical_hgvsc_legacy = hg38_tx_data[canonical_tx_hg38].get('hgvsc', '')
        hg38_canonical_hgvsp = hg38_tx_data[canonical_tx_hg38].get('hgvsp', '')
    
    # Only compare HGVSc if the SAME transcript is canonical in BOTH builds
    canonical_hgvsc_match = False
    canonical_match_details = ''
    
    if (canonical_tx_hg19 and canonical_tx_hg38 and 
        canonical_tx_hg19 == canonical_tx_hg38 and
        canonical_tx_hg19 in hg19_tx_data and canonical_tx_hg19 in hg38_tx_data):
        # Same transcript is canonical in both builds - compare its HGVSc
        canon_hg19_hgvsc = hg19_tx_data[canonical_tx_hg19].get('hgvsc', '')
        canon_hg38_hgvsc = hg38_tx_data[canonical_tx_hg19].get('hgvsc', '')
        
        match_result, details = compare_hgvsc_strings(canon_hg19_hgvsc, canon_hg38_hgvsc)
        canonical_hgvsc_match = (match_result == 'concordant')
        canonical_match_details = f"Canonical {canonical_tx_hg19}: {details}"
    else:
        canonical_hgvsc_match = False
        if canonical_tx_hg19 != canonical_tx_hg38:
            canonical_match_details = f"Different canonical: hg19={canonical_tx_hg19}, hg38={canonical_tx_hg38}"
        elif not canonical_tx_hg19 and not canonical_tx_hg38:
            canonical_match_details = "No canonical transcripts found"
        else:
            canonical_match_details = f"Canonical only in one build: hg19={canonical_tx_hg19}, hg38={canonical_tx_hg38}"
    
    # ENHANCED: Matched transcript HGVSc analysis - separate concordant/discordant lists
    matched_concordant_list = []
    matched_discordant_list = []
    
    for tx_id, hg19_hgvsc, hg38_hgvsc, details in perfect_matches:
        # Extract just the HGVSc part (remove transcript prefix if present)
        hgvsc_clean = hg19_hgvsc.split(':')[-1] if ':' in hg19_hgvsc else hg19_hgvsc
        matched_concordant_list.append(f"{tx_id}({hgvsc_clean})")
    
    for tx_id, hg19_hgvsc, hg38_hgvsc, details in mismatches:
        # Extract clean HGVSc parts
        hg19_clean = hg19_hgvsc.split(':')[-1] if ':' in hg19_hgvsc else hg19_hgvsc
        hg38_clean = hg38_hgvsc.split(':')[-1] if ':' in hg38_hgvsc else hg38_hgvsc
        matched_discordant_list.append(f"{tx_id}({hg19_clean}→{hg38_clean})")
    
    # Check high impact matches (for partial match rules)
    high_impact_matches = False
    for tx_id, _, _, _ in perfect_matches:
        hg19_impact = hg19_tx_data.get(tx_id, {}).get('impact', '')
        hg38_impact = hg38_tx_data.get(tx_id, {}).get('impact', '')
        if hg19_impact == 'HIGH' and hg38_impact == 'HIGH':
            high_impact_matches = True
            break
    
    # Create detailed summary
    if not hg19_tx_data and not hg38_tx_data:
        summary = 'no transcripts in either build'
    elif not common_transcripts:
        if hg19_tx_data and hg38_tx_data:
            summary = 'disjoint transcript sets - no shared transcripts'
        elif hg19_tx_data:
            summary = 'transcripts only in hg19'
        else:
            summary = 'transcripts only in hg38'
    else:
        # Existing logic for when there are common transcripts
        summary_parts = []
        if perfect_matches:
            summary_parts.append(f"{len(perfect_matches)} matched")
        if mismatches:
            summary_parts.append(f"{len(mismatches)} mismatched")
        if parse_errors:
            summary_parts.append(f"{len(parse_errors)} parse errors")
        
        summary = ', '.join(summary_parts) if summary_parts else 'common transcripts but no HGVSc data'
    
    return {
        'perfect_matches': perfect_matches,
        'mismatches': mismatches,
        'parse_errors': parse_errors,
        'canonical_hgvsc_match': canonical_hgvsc_match,
        'canonical_match_details': canonical_match_details,
        'high_impact_matches': high_impact_matches,
        'summary': summary,
        
        # HGVSc reporting (only transcripts WITH HGVSc)
        'hg19_canonical_hgvsc': hg19_canonical_hgvsc,
        'hg38_canonical_hgvsc': hg38_canonical_hgvsc,
        'hg19_rest_hgvsc': hg19_rest_hgvsc,
        'hg38_rest_hgvsc': hg38_rest_hgvsc,
        
        # Total transcript counts (includes all transcripts, even without HGVSc)
        'hg19_transcript_count': len(hg19_tx_data),
        'hg38_transcript_count': len(hg38_tx_data),
        
        # Matched transcript analysis (only those with HGVSc)
        'matched_transcript_count': len(perfect_matches) + len(mismatches),
        'matched_hgvsc_concordant': '; '.join(matched_concordant_list) if matched_concordant_list else '',
        'matched_hgvsc_discordant': '; '.join(matched_discordant_list) if matched_discordant_list else '',
        
        # Legacy compatibility fields
        'hg19_canonical_transcript': canonical_tx_hg19 or '',
        'hg38_canonical_transcript': canonical_tx_hg38 or '',
        'hg19_canonical_hgvsp': hg19_canonical_hgvsp,
        'hg38_canonical_hgvsp': hg38_canonical_hgvsp
    }