"""
Clinical Scoring Engine

Calculates priority scores with clinical evidence-first approach:
- Clinical significance changes get highest priority
- Impact transitions are functionally significant  
- VEP consequence mismatches get reduced weight (annotation noise)
"""

import pandas as pd
from config.scoring_config import (
    IMPACT_TRANSITION_SCORES,
    CLINICAL_OVERRIDE,
    BASE_SCORES
)
from utils.impact_utils import calculate_impact_transition_magnitude


class ClinicalScorer:
    """Calculate clinical priority scores with clinical evidence-first approach"""
    
    def __init__(self):
        """Initialize clinical scorer with updated configuration"""
        self.impact_transition_scores = IMPACT_TRANSITION_SCORES
        self.clinical_override = CLINICAL_OVERRIDE
        self.base_scores = BASE_SCORES
    
    def calculate_scores_from_analysis(self, vep_analysis_df):
        """Calculate priority scores with REDESIGNED clinical evidence-first approach"""
        print("Calculating priority scores with REDESIGNED clinical evidence-first approach...")
        
        scored_variants = []
        
        for _, row in vep_analysis_df.iterrows():
            # Calculate base priority score
            priority_score = 0
            
            # CRITICAL: Clinical significance changes (highest priority)
            clin_change = row['clin_sig_change']
            if clin_change == 'BENIGN_TO_PATHOGENIC':
                priority_score += self.base_scores['clinical_sig_benign_to_pathogenic']
            elif clin_change == 'PATHOGENIC_TO_BENIGN':
                priority_score += self.base_scores['clinical_sig_pathogenic_to_benign']
            elif clin_change.startswith('VUS_TO_PATHOGENIC'):
                priority_score += self.base_scores['clinical_sig_vus_to_pathogenic']
            elif 'TO' in str(clin_change) and not clin_change.startswith('STABLE_'):
                # Any other directional clinical change
                priority_score += self.base_scores['clinical_sig_other_change']
            
            # CRITICAL: High impact transitions (functionally significant)
            hg19_impact = row['hg19_impact']
            hg38_impact = row['hg38_impact']
            impact_magnitude, is_clinically_significant = calculate_impact_transition_magnitude(hg19_impact, hg38_impact)
            
            if row['impact_changes'] > 0 and is_clinically_significant:
                # Clinically significant impact transitions  
                transition_key = tuple(sorted([hg19_impact, hg38_impact]))
                if transition_key in self.impact_transition_scores:
                    priority_score += self.impact_transition_scores[transition_key]
                else:
                    # Fallback for unexpected transitions
                    if impact_magnitude >= 2:  # HIGH to LOW/MODIFIER or similar
                        priority_score += 15
                    elif impact_magnitude == 1 and ('HIGH' in [hg19_impact, hg38_impact] or 'MODERATE' in [hg19_impact, hg38_impact]):
                        priority_score += 10
            
            # HIGH: Functionally significant changes (demoted from CRITICAL)
            priority_score += row['same_transcript_consequence_changes'] * self.base_scores['same_transcript_consequence_changes']
            
            # MODERATE: Pathogenicity prediction changes
            if row['sift_change']:
                priority_score += self.base_scores['sift_change']
            if row['polyphen_change']:
                priority_score += self.base_scores['polyphen_change']
            
            # MODERATE: Gene changes - CONDITIONAL on clinical significance
            if row['gene_changes'] > 0:
                if is_clinically_significant or row['clin_sig_change'] or row['sift_change'] or row['polyphen_change']:
                    # Gene change is clinically relevant - use impact-based scoring
                    if hg19_impact == 'HIGH' or hg38_impact == 'HIGH':
                        priority_score += row['gene_changes'] * self.base_scores['gene_changes_high_impact']
                    elif hg19_impact == 'MODERATE' or hg38_impact == 'MODERATE':
                        priority_score += row['gene_changes'] * self.base_scores['gene_changes_moderate_impact']
                    elif hg19_impact == 'LOW' or hg38_impact == 'LOW':
                        priority_score += row['gene_changes'] * self.base_scores['gene_changes_low_impact']
                    else:
                        priority_score += row['gene_changes'] * self.base_scores['gene_changes_mixed_impact']
                else:
                    # Gene change likely annotation synonym - minimal weight
                    priority_score += row['gene_changes'] * self.base_scores['gene_changes_minimal_weight']
            
            # LOW: Annotation differences (DEMOTED - reduced weights)
            priority_score += row['unmatched_consequences'] * self.base_scores['unmatched_consequences']
            priority_score += row['same_consequence_different_transcripts'] * self.base_scores['same_consequence_different_transcripts']
            
            # Technical liftover issues (unchanged)
            if row['pos_match'] == 0:
                priority_score += self.base_scores['position_mismatch']
            if row['pos_difference'] > 10:
                priority_score += self.base_scores['position_difference_moderate']
            if row['pos_difference'] > 100:
                priority_score += self.base_scores['position_difference_large']
            if row['gt_match'] == 0:
                priority_score += self.base_scores['genotype_mismatch']
            
            # BCFtools swap handling
            swap_str = str(row['swap']) if pd.notna(row['swap']) else 'NA'
            if swap_str == '1':
                priority_score += self.base_scores['ref_alt_swap']
            
            # Bonuses (unchanged)
            has_clinical_data = (
                not pd.isna(row['hg19_clin_sig']) and row['hg19_clin_sig'] not in ['', '-'] or
                not pd.isna(row['hg38_clin_sig']) and row['hg38_clin_sig'] not in ['', '-']
            )
            if has_clinical_data:
                priority_score += self.base_scores['has_clinical_data_bonus']
            
            if row['mapping_status'] == 'REGION':
                priority_score += self.base_scores['region_mapping_bonus']
            
            if row['hg19_impact'] == 'HIGH' or row['hg38_impact'] == 'HIGH':
                priority_score += self.base_scores['high_impact_bonus']
            
            # CLINICAL EVIDENCE OVERRIDE (unchanged - works well)
            is_low_modifier_impact = (
                row['hg19_impact'] in ['MODIFIER', 'LOW'] and 
                row['hg38_impact'] in ['MODIFIER', 'LOW']
            )
        
            has_benign_evidence = (
                row['hg19_clin_sig_normalized'] == 'BENIGN' or 
                row['hg38_clin_sig_normalized'] == 'BENIGN' or
                'benign' in str(row['hg19_sift']).lower() or
                'benign' in str(row['hg19_polyphen']).lower() or
                'benign' in str(row['hg38_sift']).lower() or
                'benign' in str(row['hg38_polyphen']).lower()
            )
        
            is_benign_variant = is_low_modifier_impact and has_benign_evidence
        
            has_pathogenic_evidence = (
                row['hg19_clin_sig_normalized'] == 'PATHOGENIC' or 
                row['hg38_clin_sig_normalized'] == 'PATHOGENIC' or
                row['hg19_impact'] == 'HIGH' or row['hg38_impact'] == 'HIGH' or
                'deleterious' in str(row['hg19_sift']).lower() or
                'deleterious' in str(row['hg38_sift']).lower() or
                'probably_damaging' in str(row['hg19_polyphen']).lower() or
                'probably_damaging' in str(row['hg38_polyphen']).lower()
            )
            
            # Apply clinical evidence override
            if is_benign_variant:
                priority_score *= self.clinical_override['benign_reduction_factor']
            elif has_pathogenic_evidence:
                priority_score *= self.clinical_override['pathogenic_boost_factor']
            
            # REDESIGNED: Priority category determination with clinical focus
            has_critical_clinical_change = clin_change in ['BENIGN_TO_PATHOGENIC', 'PATHOGENIC_TO_BENIGN'] or clin_change.startswith('VUS_TO_PATHOGENIC')
            has_high_impact_transition = (row['impact_changes'] > 0 and is_clinically_significant and 
                                        'HIGH' in [hg19_impact, hg38_impact])
            
            if has_critical_clinical_change or has_high_impact_transition:
                priority_category = 'CRITICAL'
            elif (row['impact_changes'] > 0 and is_clinically_significant) or \
                 (row['clin_sig_change'] and 'TO' in str(row['clin_sig_change']) and not row['clin_sig_change'].startswith('STABLE_')) or \
                 row['same_transcript_consequence_changes'] > 0:
                priority_category = 'HIGH'
            elif row['sift_change'] or row['polyphen_change'] or \
                 (row['gene_changes'] > 0 and (is_clinically_significant or row['clin_sig_change'] or row['sift_change'] or row['polyphen_change'])):
                priority_category = 'MODERATE'
            elif is_benign_variant or \
                 (row['unmatched_consequences'] > 0 and not (row['clin_sig_change'] or row['sift_change'] or row['polyphen_change'])):
                priority_category = 'LOW'
            else:
                priority_category = 'INVESTIGATE'
            
            # Create summary of discordance types (updated priorities)
            discordance_summary = []
            if has_critical_clinical_change:
                discordance_summary.append(f"CRITICAL clinical significance change: {row['clin_sig_change']}")
            if has_high_impact_transition:
                discordance_summary.append(f"CRITICAL impact transition: {hg19_impact}→{hg38_impact}")
            if row['same_transcript_consequence_changes'] > 0:
                discordance_summary.append(f"Same transcript consequence changes: {row['same_transcript_consequence_changes']}")
            if row['impact_changes'] > 0 and is_clinically_significant:
                discordance_summary.append(f"Impact level changes: {row['impact_changes']}")
            if row['gene_changes'] > 0:
                discordance_summary.append(f"Gene annotation changes: {row['gene_changes']}")
            if row['sift_change']:
                discordance_summary.append(f"SIFT change: {row['sift_change']}")
            if row['polyphen_change']:
                discordance_summary.append(f"PolyPhen change: {row['polyphen_change']}")
            if row['unmatched_consequences'] > 0:
                discordance_summary.append(f"Unmatched consequences (annotation differences)")
            if row['same_consequence_different_transcripts'] > 0:
                discordance_summary.append(f"Same consequence, different transcripts: {row['same_consequence_different_transcripts']}")
            if row['pos_match'] == 0:
                discordance_summary.append(f"Position mismatch")
            if row['gt_match'] == 0:
                discordance_summary.append(f"Genotype mismatch")
            
            # Only include variants with discordances (score > 0)
            if priority_score == 0:
                continue
            
            # Add calculated fields to the row
            variant_info = row.to_dict()
            variant_info['priority_score'] = priority_score
            variant_info['priority_category'] = priority_category
            variant_info['discordance_summary'] = '; '.join(discordance_summary) if discordance_summary else 'Technical issues only'
            
            scored_variants.append(variant_info)
        
        return pd.DataFrame(scored_variants)