"""
Clinical Scoring Engine

Calculates priority scores with clinical evidence-first approach:
- Clinical significance changes get highest priority
- Impact transitions are functionally significant  
- VEP consequence mismatches get reduced weight (annotation noise)
"""

import pandas as pd
from config.scoring_config import (
    BASE_SCORES,                    # Updated with new scoring
    PRIORITY_THRESHOLDS,            # New: category thresholds
    CLINICAL_CHANGE_TYPES,          # New: clinical change classifications
    CLINICAL_OVERRIDE               # Kept: legacy benign/pathogenic handling
)
from utils.impact_utils import calculate_impact_transition_magnitude


class ClinicalScorer:
    """Calculate clinical priority scores with clinical evidence-first approach"""
    
    def __init__(self):
        """Initialize clinical scorer with updated configuration"""
        self.base_scores = BASE_SCORES
        self.thresholds = PRIORITY_THRESHOLDS
        self.clinical_types = CLINICAL_CHANGE_TYPES
        self.clinical_override = CLINICAL_OVERRIDE

    def assign_priority_category(self, score, score_details):
        """Assign category based on score thresholds with strict CONCORDANT criteria"""
        
        # STRICT CONCORDANT: Only gene annotation differences allowed
        if score == 0:
            return 'CONCORDANT'
        elif score <= self.base_scores['gene_changes']:
            # Check if ONLY gene changes contributed to the score
            gene_only = all('gene changes' in detail.lower() for detail in score_details if detail)
            if gene_only:
                return 'CONCORDANT'
        
        # Regular thresholds for other categories
        for category in ['CRITICAL', 'MODERATE', 'LOW']:
            if score >= self.thresholds[category]:
                return category
        
        return 'LOW'  # Fallback for edge cases
 
    def calculate_scores_from_analysis(self, vep_analysis_df):
        """Calculate priority scores using priority transcript approach"""
        print("Calculating priority scores using priority transcript approach...")
        
        scored_variants = []
        
        for _, row in vep_analysis_df.iterrows():
            # Initialize scoring
            priority_score = 0
            score_details = []

            # ===== TRANSCRIPT MATCHING ISSUES =====
            transcript_status = row.get('transcript_crossbuild_status', '')
            if transcript_status == 'MANE_hg38_Only':
                priority_score += self.base_scores['transcript_mismatch']
                score_details.append(f"Transcript mismatch ({self.base_scores['transcript_mismatch']})")
            elif transcript_status == 'No_Matching_Transcripts':
                priority_score += self.base_scores['no_matching_transcripts']
                score_details.append(f"No matching transcripts ({self.base_scores['no_matching_transcripts']})")
            elif transcript_status == 'No_Transcripts':
                priority_score += self.base_scores['no_transcripts']
                score_details.append(f"No transcripts ({self.base_scores['no_transcripts']})")
            
            # ===== CRITICAL PRIORITY: HGVS concordance on priority transcript =====
            # Only score HGVS issues if we actually have transcripts to analyze
            if transcript_status not in ['No_Matching_Transcripts', 'No_Transcripts']:
                if row.get('priority_hgvsc_concordance') == 'Mismatch':
                    priority_score += self.base_scores['priority_hgvsc_mismatch']
                    score_details.append(f"HGVSc mismatch ({self.base_scores['priority_hgvsc_mismatch']})")
                elif row.get('priority_hgvsc_concordance') == 'No_Analysis':
                    priority_score += self.base_scores['missing_hgvs_data']
                    score_details.append(f"Missing HGVSc data ({self.base_scores['missing_hgvs_data']})")
                
                if row.get('priority_hgvsp_concordance') == 'Mismatch':
                    priority_score += self.base_scores['priority_hgvsp_mismatch']
                    score_details.append(f"HGVSp mismatch ({self.base_scores['priority_hgvsp_mismatch']})")
                elif row.get('priority_hgvsp_concordance') == 'No_Analysis':
                    priority_score += self.base_scores['missing_hgvs_data']
                    score_details.append(f"Missing HGVSp data ({self.base_scores['missing_hgvs_data']})")
            
            # ===== CLINICAL SIGNIFICANCE CHANGES =====
            clin_change = str(row.get('clin_sig_change', ''))
            if clin_change in self.clinical_types['pathogenic_benign_flip']:
                priority_score += self.base_scores['pathogenic_benign_flip']
                score_details.append(f"Pathogenic↔Benign change ({self.base_scores['pathogenic_benign_flip']})")
            elif clin_change in self.clinical_types['pathogenic_vus_change']:
                priority_score += self.base_scores['pathogenic_vus_change']
                score_details.append(f"Pathogenic↔VUS change ({self.base_scores['pathogenic_vus_change']})")
            elif clin_change in self.clinical_types['vus_benign_change']:
                priority_score += self.base_scores['vus_benign_change']
                score_details.append(f"VUS↔Benign change ({self.base_scores['vus_benign_change']})")
            elif clin_change in self.clinical_types['minor_clinical_changes']:
                priority_score += self.base_scores['minor_clinical_changes']
                score_details.append(f"Minor clinical change ({self.base_scores['minor_clinical_changes']})")
            elif clin_change == 'NONE' or clin_change.startswith('STABLE_NONE'):
                priority_score += self.base_scores['missing_clinical_data']
                score_details.append(f"Missing clinical data ({self.base_scores['missing_clinical_data']})")
            
            # ===== WORST CONSEQUENCE DIFFERENCES =====
            if row.get('has_worst_consequence_difference') == 'YES':
                # Analyze severity of difference using impact changes as proxy
                if row.get('impact_changes', 0) > 0:
                    priority_score += self.base_scores['serious_consequence_difference']
                    score_details.append(f"Serious consequence difference ({self.base_scores['serious_consequence_difference']})")
                else:
                    priority_score += self.base_scores['minor_clinical_changes']
                    score_details.append(f"Minor consequence difference ({self.base_scores['minor_clinical_changes']})")
            
            # ===== BASIC ANNOTATION CHANGES =====
            if row.get('gene_changes', 0) > 0:
                priority_score += self.base_scores['gene_changes']
                score_details.append(f"Gene changes ({self.base_scores['gene_changes']})")
            
            if row.get('impact_changes', 0) > 0:
                priority_score += self.base_scores['impact_changes']
                score_details.append(f"Impact changes ({self.base_scores['impact_changes']})")
            
            # ===== PATHOGENICITY PREDICTION CHANGES =====
            if 'TO' in str(row.get('sift_change', '')):
                priority_score += self.base_scores['prediction_changes']
                score_details.append(f"SIFT change ({self.base_scores['prediction_changes']})")
            
            if 'TO' in str(row.get('polyphen_change', '')):
                priority_score += self.base_scores['prediction_changes']
                score_details.append(f"PolyPhen change ({self.base_scores['prediction_changes']})")
            
            # ===== APPLY CLINICAL EVIDENCE OVERRIDE =====
            # Keep existing benign/pathogenic override logic
            is_low_modifier_impact = (
                row.get('hg19_impact') in ['MODIFIER', 'LOW'] and 
                row.get('hg38_impact') in ['MODIFIER', 'LOW']
            )
        
            has_benign_evidence = (
                row.get('hg19_clin_sig_normalized') == 'BENIGN' or 
                row.get('hg38_clin_sig_normalized') == 'BENIGN' or
                'benign' in str(row.get('hg19_sift', '')).lower() or
                'benign' in str(row.get('hg19_polyphen', '')).lower() or
                'benign' in str(row.get('hg38_sift', '')).lower() or
                'benign' in str(row.get('hg38_polyphen', '')).lower()
            )
        
            is_benign_variant = is_low_modifier_impact and has_benign_evidence
        
            has_pathogenic_evidence = (
                row.get('hg19_clin_sig_normalized') == 'PATHOGENIC' or 
                row.get('hg38_clin_sig_normalized') == 'PATHOGENIC' or
                row.get('hg19_impact') == 'HIGH' or row.get('hg38_impact') == 'HIGH' or
                'deleterious' in str(row.get('hg19_sift', '')).lower() or
                'deleterious' in str(row.get('hg38_sift', '')).lower() or
                'probably_damaging' in str(row.get('hg19_polyphen', '')).lower() or
                'probably_damaging' in str(row.get('hg38_polyphen', '')).lower()
            )

            # ===== TECHNICAL LIFTOVER ISSUES =====
            if row.get('pos_match', 1) == 0:
                priority_score += self.base_scores.get('position_mismatch', 20)
                score_details.append(f"Position mismatch ({self.base_scores.get('position_mismatch', 20)})")
            
            pos_diff = row.get('pos_difference', 0)
            if pos_diff > 100:
                priority_score += self.base_scores.get('position_difference_large', 15)
                score_details.append(f"Large position difference ({self.base_scores.get('position_difference_large', 15)})")
            elif pos_diff > 10:
                priority_score += self.base_scores.get('position_difference_moderate', 10)
                score_details.append(f"Moderate position difference ({self.base_scores.get('position_difference_moderate', 10)})")
            
            if row.get('gt_match', 1) == 0:
                priority_score += self.base_scores.get('genotype_mismatch', 20)
                score_details.append(f"Genotype mismatch ({self.base_scores.get('genotype_mismatch', 20)})")
            
            # BCFtools swap handling
            if str(row.get('swap', 'NA')) == '1':
                priority_score += self.base_scores.get('ref_alt_swap', 10)
                score_details.append(f"Ref/Alt swap ({self.base_scores.get('ref_alt_swap', 10)})")
            
            # Apply clinical evidence override
            if is_benign_variant:
                priority_score *= self.clinical_override['benign_reduction_factor']
                score_details.append("Benign evidence override (×0.1)")
            elif has_pathogenic_evidence:
                priority_score *= self.clinical_override['pathogenic_boost_factor']
                score_details.append("Pathogenic evidence boost (×2.0)")
            
            # ===== ASSIGN PRIORITY CATEGORY BASED ON SCORE THRESHOLDS =====
            priority_category = self.assign_priority_category(priority_score, score_details)

            
            # ===== CREATE VARIANT RECORD =====
            variant_info = row.to_dict()
            variant_info['priority_score'] = int(round(priority_score))
            variant_info['priority_category'] = priority_category
            variant_info['score_breakdown'] = '; '.join(score_details) if score_details else 'No significant issues'
            
            scored_variants.append(variant_info)

        # Create the initial dataframe
        result_df = pd.DataFrame(scored_variants)

        # Add Clinical_Change_Direction column using same categorization as Plot 3
        def categorize_clin_sig_transition_for_plotting(row):
            """Reuse the same categorization logic as Plot 3"""
            hg19 = row.get('hg19_clin_sig_normalized', '')
            hg38 = row.get('hg38_clin_sig_normalized', '')
            
            if hg19 == hg38:
                return f'Stable {hg19}' if hg19 else 'Stable NONE'
            elif hg38 == 'PATHOGENIC':
                return '→ PATHOGENIC'
            elif hg19 == 'PATHOGENIC':
                return 'FROM PATHOGENIC →'
            elif hg38 == 'BENIGN':
                return '→ BENIGN'
            elif hg19 == 'BENIGN':
                return 'FROM BENIGN →'
            else:
                return 'Other Transitions'

        result_df['Clinical_Change_Direction'] = result_df.apply(categorize_clin_sig_transition_for_plotting, axis=1)
    

        return result_df