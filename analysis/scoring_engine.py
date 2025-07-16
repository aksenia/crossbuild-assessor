"""
Clinical Scoring Engine

Calculates priority scores from VEP analysis data using clinical evidence-driven
scoring with magnitude-based impact transitions and clinical overrides.
"""

import pandas as pd
from config.scoring_config import (
    IMPACT_TRANSITION_SCORES,
    CLINICAL_OVERRIDE,
    BASE_SCORES
)
from utils.impact_utils import calculate_impact_transition_magnitude


class ClinicalScorer:
    """Calculate clinical priority scores from VEP analysis data"""
    
    def __init__(self):
        """Initialize clinical scorer with configuration"""
        self.impact_transition_scores = IMPACT_TRANSITION_SCORES
        self.clinical_override = CLINICAL_OVERRIDE
        self.base_scores = BASE_SCORES
    
    def calculate_scores_from_analysis(self, vep_analysis_df):
        """Calculate priority scores from cached VEP analysis data with magnitude-based impact scoring"""
        print("Calculating priority scores with magnitude-based impact transitions...")
        
        scored_variants = []
        
        for _, row in vep_analysis_df.iterrows():
            # Calculate base priority score
            priority_score = 0
            
            # Position and genotype discordance scoring
            if row['pos_match'] == 0:
                priority_score += self.base_scores['position_mismatch']
            if row['pos_difference'] > 10:
                priority_score += self.base_scores['position_difference_moderate']
            if row['pos_difference'] > 100:
                priority_score += self.base_scores['position_difference_large']
            if row['gt_match'] == 0:
                priority_score += self.base_scores['genotype_mismatch']
            
            # BCFtools liftover swap values: 1 = swapped, -1 = swap failed, NA = no swap needed
            swap_str = str(row['swap']) if pd.notna(row['swap']) else 'NA'
            if swap_str == '1':  # REF/ALT alleles were swapped during liftover
                priority_score += self.base_scores['ref_alt_swap']
            
            # Core functional changes (highest priority)
            priority_score += row['same_transcript_consequence_changes'] * self.base_scores['same_transcript_consequence_changes']  # CRITICAL
            
            # MAGNITUDE-BASED IMPACT SCORING
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
                    if impact_magnitude == 1:
                        if 'HIGH' in [hg19_impact, hg38_impact] and 'MODERATE' in [hg19_impact, hg38_impact]:
                            priority_score += self.impact_transition_scores[('HIGH', 'MODERATE')]
                        elif 'MODERATE' in [hg19_impact, hg38_impact] and 'LOW' in [hg19_impact, hg38_impact]:
                            priority_score += self.impact_transition_scores[('MODERATE', 'LOW')]
                    elif impact_magnitude == 2:
                        if 'HIGH' in [hg19_impact, hg38_impact]:
                            priority_score += self.impact_transition_scores[('HIGH', 'LOW')]
                        elif 'MODERATE' in [hg19_impact, hg38_impact]:
                            priority_score += self.impact_transition_scores[('MODERATE', 'MODIFIER')]
                    elif impact_magnitude == 3:
                        priority_score += self.impact_transition_scores[('HIGH', 'MODIFIER')]
            elif row['impact_changes'] > 0:
                # Non-clinically significant transitions (annotation noise)
                priority_score += self.impact_transition_scores[('LOW', 'MODIFIER')]
            
            priority_score += row['unmatched_consequences'] * self.base_scores['unmatched_consequences']               # INVESTIGATE
            
            # Clinical significance changes
            if row['clin_sig_change'] == 'BENIGN_TO_PATHOGENIC':
                priority_score += self.base_scores['clinical_sig_benign_to_pathogenic']
            elif row['clin_sig_change'] == 'PATHOGENIC_TO_BENIGN':
                priority_score += self.base_scores['clinical_sig_pathogenic_to_benign']
            elif row['clin_sig_change'] == 'OTHER_CHANGE':
                priority_score += self.base_scores['clinical_sig_other_change']
            
            # Pathogenicity prediction changes
            if row['sift_change']:
                priority_score += self.base_scores['sift_change']
            if row['polyphen_change']:
                priority_score += self.base_scores['polyphen_change']
            
            # Gene changes - CONDITIONAL on clinical significance
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
                        priority_score += row['gene_changes'] * self.base_scores['gene_changes_mixed_impact']  # Mixed or unknown
                else:
                    # Gene change likely annotation synonym - minimal weight
                    priority_score += row['gene_changes'] * self.base_scores['gene_changes_minimal_weight']
            
            # Other transcript-level issues
            priority_score += row['same_consequence_different_transcripts'] * self.base_scores['same_consequence_different_transcripts']  # Moderate priority
            
            # Clinical significance bonus (if any clinical data available)
            has_clinical_data = (
                not pd.isna(row['hg19_clin_sig']) and row['hg19_clin_sig'] not in ['', '-'] or
                not pd.isna(row['hg38_clin_sig']) and row['hg38_clin_sig'] not in ['', '-']
            )
            if has_clinical_data:
                priority_score += self.base_scores['has_clinical_data_bonus']
            
            # Mapping status and impact adjustments
            if row['mapping_status'] == 'REGION':
                priority_score += self.base_scores['region_mapping_bonus']
            
            if row['hg19_impact'] == 'HIGH' or row['hg38_impact'] == 'HIGH':
                priority_score += self.base_scores['high_impact_bonus']
            
            # CLINICAL EVIDENCE OVERRIDE
            # Check if variant is benign based on multiple evidence sources
            is_low_modifier_impact = (
                row['hg19_impact'] in ['MODIFIER', 'LOW'] and 
                row['hg38_impact'] in ['MODIFIER', 'LOW']
            )
            
            has_benign_evidence = (
                row['hg19_is_benign'] or row['hg38_is_benign'] or
                'benign' in str(row['hg19_sift']).lower() or
                'benign' in str(row['hg19_polyphen']).lower() or
                'benign' in str(row['hg38_sift']).lower() or
                'benign' in str(row['hg38_polyphen']).lower()
            )
            
            is_benign_variant = is_low_modifier_impact and has_benign_evidence
            
            # Check if variant has pathogenic evidence
            has_pathogenic_evidence = (
                row['hg19_is_pathogenic'] or row['hg38_is_pathogenic'] or
                row['hg19_impact'] == 'HIGH' or row['hg38_impact'] == 'HIGH' or
                'deleterious' in str(row['hg19_sift']).lower() or
                'deleterious' in str(row['hg38_sift']).lower() or
                'probably_damaging' in str(row['hg19_polyphen']).lower() or
                'probably_damaging' in str(row['hg38_polyphen']).lower()
            )
            
            # Apply clinical evidence override
            if is_benign_variant:
                priority_score *= self.clinical_override['benign_reduction_factor']  # 90% reduction for benign variants
            elif has_pathogenic_evidence:
                priority_score *= self.clinical_override['pathogenic_boost_factor']  # 2x boost for pathogenic variants
            
            # Determine priority category with updated logic
            if row['same_transcript_consequence_changes'] > 0:
                priority_category = 'CRITICAL'
            elif (is_clinically_significant and row['impact_changes'] > 0) or \
                 row['clin_sig_change'] in ['BENIGN_TO_PATHOGENIC', 'PATHOGENIC_TO_BENIGN']:
                priority_category = 'HIGH'
            elif (row['gene_changes'] > 0 and is_clinically_significant) or \
                 row['clin_sig_change'] == 'OTHER_CHANGE' or row['sift_change'] or row['polyphen_change']:
                priority_category = 'MODERATE'
            elif row['unmatched_consequences'] > 0 or \
                 (row['impact_changes'] > 0 and not is_clinically_significant):
                priority_category = 'INVESTIGATE'
            elif is_benign_variant:
                priority_category = 'LOW'
            else:
                priority_category = 'MODERATE'
            
            # Create summary of discordance types for this variant
            discordance_summary = []
            if row['same_transcript_consequence_changes'] > 0:
                discordance_summary.append(f"Same transcript consequence changes: {row['same_transcript_consequence_changes']}")
            if row['gene_changes'] > 0:
                discordance_summary.append(f"Gene annotation changes: {row['gene_changes']}")
            if row['impact_changes'] > 0:
                discordance_summary.append(f"Impact level changes: {row['impact_changes']}")
            if row['same_consequence_different_transcripts'] > 0:
                discordance_summary.append(f"Same consequence, different transcripts: {row['same_consequence_different_transcripts']}")
            if row['unmatched_consequences'] > 0:
                discordance_summary.append(f"Unmatched consequences")
            if row['clin_sig_change']:
                discordance_summary.append(f"Clinical significance change: {row['clin_sig_change']}")
            if row['sift_change']:
                discordance_summary.append(f"SIFT change: {row['sift_change']}")
            if row['polyphen_change']:
                discordance_summary.append(f"PolyPhen change: {row['polyphen_change']}")
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
            variant_info['discordance_summary'] = '; '.join(discordance_summary) if discordance_summary else 'Position/genotype issues only'
            
            scored_variants.append(variant_info)
        
        return pd.DataFrame(scored_variants)
