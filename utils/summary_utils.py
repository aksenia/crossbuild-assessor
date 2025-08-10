"""
Summary Data Calculator

Calculates structured summary data from analysis results for both
liftover QC analysis and variant prioritization.
"""

import pandas as pd
from datetime import datetime


class SummaryDataCalculator:
    """Calculate structured summary data from analysis results"""
    
    def calculate_liftover_summary(self, conn):
        """Calculate liftover QC summary data (for db_analyzer.py)"""
        
        # Basic database statistics
        total_variants = pd.read_sql_query("SELECT COUNT(*) as count FROM comparison", conn).iloc[0]['count']
        
        pos_match_rate = pd.read_sql_query("""
            SELECT AVG(pos_match) as rate FROM comparison
        """, conn).iloc[0]['rate']
        
        gt_match_rate = pd.read_sql_query("""
            SELECT AVG(gt_match) as rate FROM comparison
        """, conn).iloc[0]['rate']
        
        # Calculate proper concordant/discordant variants (non-redundant)
        both_match_query = """
            SELECT COUNT(*) as count 
            FROM comparison 
            WHERE pos_match = 1 AND gt_match = 1
        """
        both_match_result = pd.read_sql_query(both_match_query, conn)
        concordant_variants = both_match_result.iloc[0]['count']

        # Discordant = any mismatch (position OR genotype OR both) - non-redundant
        discordant_variants = total_variants - concordant_variants
        
        # Calculate percentages
        position_match_percentage = round(pos_match_rate * 100, 1) if pos_match_rate else 0
        genotype_match_percentage = round(gt_match_rate * 100, 1) if gt_match_rate else 0
        
        # Calculate proper match categories
        match_categories_query = """
            SELECT 
                COUNT(CASE WHEN pos_match = 1 AND gt_match = 1 THEN 1 END) as both_match,
                COUNT(CASE WHEN pos_match = 1 AND gt_match = 0 THEN 1 END) as position_only,
                COUNT(CASE WHEN pos_match = 0 AND gt_match = 1 THEN 1 END) as genotype_only,
                COUNT(CASE WHEN pos_match = 0 AND gt_match = 0 THEN 1 END) as both_mismatch
            FROM comparison
        """
        
        match_stats = pd.read_sql_query(match_categories_query, conn).iloc[0]
        both_match_count = match_stats['both_match']
        position_only_count = match_stats['position_only']
        genotype_only_count = match_stats['genotype_only']
        both_mismatch_count = match_stats['both_mismatch']
        
        # Mapping status breakdown
        mapping_stats = pd.read_sql_query("""
            SELECT mapping_status, COUNT(*) as count,
                   AVG(pos_match) as pos_match_rate,
                   AVG(gt_match) as gt_match_rate
            FROM comparison 
            GROUP BY mapping_status
            ORDER BY count DESC
        """, conn)
        
        # Position differences analysis
        pos_diff_query = """
        SELECT 
            mapping_status,
            source_chrom,
            ABS(COALESCE(liftover_hg38_pos, 0) - COALESCE(bcftools_hg38_pos, 0)) as pos_diff
        FROM comparison 
        WHERE pos_match = 0 AND liftover_hg38_pos IS NOT NULL AND bcftools_hg38_pos IS NOT NULL
        """
        
        diff_df = pd.read_sql_query(pos_diff_query, conn)
        
        # Flip/swap analysis
        flip_swap_query = """
        SELECT flip, swap, COUNT(*) as count
        FROM comparison
        WHERE pos_match = 0 OR gt_match = 0
        GROUP BY flip, swap
        """
        
        flip_swap_df = pd.read_sql_query(flip_swap_query, conn)
        
        # FIXED: Clean flip/swap category names with descriptive labels
        def clean_flip_swap_category(flip, swap):
            """Convert raw flip/swap values to descriptive category names"""
            flip_clean = str(flip).lower() if pd.notna(flip) else 'none'
            swap_clean = str(swap).lower() if pd.notna(swap) else 'none'
            
            # Handle float values that might be stored as strings
            if swap_clean in ['1.0', '1']:
                swap_clean = '1'
            elif swap_clean in ['-1.0', '-1']:
                swap_clean = '-1'
            elif swap_clean in ['na', 'none', 'nan']:
                swap_clean = 'na'
            
            # Create descriptive category names
            if flip_clean == 'no_flip' and swap_clean == 'na':
                return "No Changes Required"
            elif flip_clean == 'flip' and swap_clean == 'na':
                return "Strand Flip Only"
            elif flip_clean == 'no_flip' and swap_clean == '1':
                return "Allele Swap Successful"
            elif flip_clean == 'flip' and swap_clean == '1':
                return "Strand Flip + Allele Swap Successful"
            elif flip_clean == 'no_flip' and swap_clean == '-1':
                return "Allele Swap Failed"
            elif flip_clean == 'flip' and swap_clean == '-1':
                return "Strand Flip + Allele Swap Failed"
            else:
                return f"Other Configuration ({flip_clean}, {swap_clean})"
        
        flip_swap_cleaned = {}
        for _, row in flip_swap_df.iterrows():
            category = clean_flip_swap_category(row['flip'], row['swap'])
            flip_swap_cleaned[category] = flip_swap_cleaned.get(category, 0) + row['count']
        
        return {
            "metadata": {
                "script": "db_analyzer",
                "timestamp": datetime.now().isoformat(),
                "total_variants": total_variants
            },
            "dataset_overview": {
                "total_variants": total_variants,
                "concordant_variants": concordant_variants,
                "discordant_variants": discordant_variants,
                "position_match_percentage": position_match_percentage,
                "genotype_match_percentage": genotype_match_percentage,
                "position_match_rate": round(pos_match_rate, 3) if pos_match_rate else 0,
                "genotype_match_rate": round(gt_match_rate, 3) if gt_match_rate else 0
            },
            "quality_summary": {
                "position_mismatches": total_variants - int(pos_match_rate * total_variants) if pos_match_rate < 1.0 else 0,
                "genotype_mismatches": total_variants - int(gt_match_rate * total_variants) if gt_match_rate < 1.0 else 0,
                "coordinate_discrepancies": len(diff_df),
                "any_mismatch": discordant_variants 
            },
            "match_categories": {
                "both_match": {
                    "count": both_match_count,
                    "percentage": round(both_match_count / total_variants * 100, 1) if total_variants > 0 else 0
                },
                "position_only": {
                    "count": position_only_count,
                    "percentage": round(position_only_count / total_variants * 100, 1) if total_variants > 0 else 0
                },
                "genotype_only": {
                    "count": genotype_only_count,
                    "percentage": round(genotype_only_count / total_variants * 100, 1) if total_variants > 0 else 0
                },
                "both_mismatch": {
                    "count": both_mismatch_count,
                    "percentage": round(both_mismatch_count / total_variants * 100, 1) if total_variants > 0 else 0
                }
            },
            "mapping_status_breakdown": {
                row['mapping_status']: {
                    "count": row['count'],
                    "pos_match_percentage": round(row['pos_match_rate'] * 100, 1),
                    "gt_match_percentage": round(row['gt_match_rate'] * 100, 1),
                    "pos_match_rate": round(row['pos_match_rate'], 3),
                    "gt_match_rate": round(row['gt_match_rate'], 3)
                }
                for _, row in mapping_stats.iterrows()
            },
            "position_differences": {
                "total_mismatches": len(diff_df),
                "statistics": {
                    "mean": round(diff_df['pos_diff'].mean(), 1) if len(diff_df) > 0 else 0,
                    "max": int(diff_df['pos_diff'].max()) if len(diff_df) > 0 else 0,
                    "min": int(diff_df['pos_diff'].min()) if len(diff_df) > 0 else 0,
                    "std": round(diff_df['pos_diff'].std(), 1) if len(diff_df) > 0 else 0
                },
                "by_chromosome": {
                    chrom: {
                        "count": int(group['pos_diff'].count()),
                        "mean": round(group['pos_diff'].mean(), 1),
                        "max": int(group['pos_diff'].max())
                    }
                    for chrom, group in diff_df.groupby('source_chrom')
                } if len(diff_df) > 0 else {}
            },
            "flip_swap_analysis": flip_swap_cleaned
        }
    
    def calculate_prioritization_summary(self, df_full, df_excel):
        """Calculate prioritization summary data (for variant_prioritizer.py)"""
        
        # Priority distribution
        priority_distribution = df_full['priority_category'].value_counts().to_dict() if len(df_full) > 0 else {}
        
        # Clinical transitions analysis
        clinical_transitions = {}
        if len(df_full) > 0 and 'hg19_clin_sig_normalized' in df_full.columns:
            # Stable annotations
            stable_variants = df_full[df_full['hg19_clin_sig_normalized'] == df_full['hg38_clin_sig_normalized']]
            stable_counts = stable_variants['hg19_clin_sig_normalized'].value_counts().to_dict()
            
            # Directional changes
            changing_variants = df_full[df_full['hg19_clin_sig_normalized'] != df_full['hg38_clin_sig_normalized']]
            directional_changes_list = [] 
            
            if len(changing_variants) > 0:
                transition_counts = changing_variants.groupby(['hg19_clin_sig_normalized', 'hg38_clin_sig_normalized']).size()
                
                # Calculate clinical priority for each transition
                def get_clinical_priority(hg19_cat, hg38_cat):
                    """Determine clinical priority for a transition"""
                    # CRITICAL: Direct clinical interpretation changes
                    if (hg19_cat == 'PATHOGENIC' and hg38_cat == 'BENIGN') or \
                       (hg19_cat == 'BENIGN' and hg38_cat == 'PATHOGENIC'):
                        return 'CRITICAL'
                    elif hg19_cat == 'VUS' and hg38_cat == 'PATHOGENIC':
                        return 'CRITICAL'
                    elif hg19_cat == 'PATHOGENIC' and hg38_cat == 'VUS':
                        return 'CRITICAL'
                    
                    # HIGH: Pathogenic evidence changes
                    elif hg19_cat == 'PATHOGENIC' and hg38_cat in ['RISK', 'OTHER', 'NONE']:
                        return 'HIGH'
                    elif hg19_cat in ['RISK', 'OTHER', 'NONE'] and hg38_cat == 'PATHOGENIC':
                        return 'HIGH'
                    
                    # MODERATE: Benign/VUS transitions, new evidence
                    elif (hg19_cat == 'BENIGN' and hg38_cat == 'VUS') or \
                         (hg19_cat == 'VUS' and hg38_cat == 'BENIGN'):
                        return 'MODERATE'
                    elif hg19_cat == 'NONE' and hg38_cat in ['PATHOGENIC', 'BENIGN']:
                        return 'MODERATE'
                    elif hg19_cat in ['PATHOGENIC', 'BENIGN'] and hg38_cat == 'NONE':
                        return 'MODERATE'
                    
                    # LOW: Minor evidence changes
                    elif (hg19_cat == 'BENIGN' and hg38_cat == 'NONE') or \
                         (hg19_cat == 'NONE' and hg38_cat == 'BENIGN'):
                        return 'LOW'
                    elif hg19_cat == 'NONE' and hg38_cat == 'VUS':
                        return 'LOW'
                    elif hg19_cat == 'VUS' and hg38_cat == 'NONE':
                        return 'LOW'
                    
                    # Default for other transitions
                    else:
                        return 'MODERATE'
                
                # Create ordered list of directional changes
                directional_changes_list = []
                for (hg19_cat, hg38_cat), count in transition_counts.items():
                    clinical_priority = get_clinical_priority(hg19_cat, hg38_cat)
                    directional_changes_list.append({
                        "transition": f"{hg19_cat}→{hg38_cat}",
                        "count": count,
                        "clinical_priority": clinical_priority
                    })
                
                # Sort by priority (CRITICAL first), then by count (descending)
                priority_order = {'CRITICAL': 0, 'HIGH': 1, 'MODERATE': 2, 'LOW': 3}
                directional_changes_list.sort(key=lambda x: (priority_order.get(x['clinical_priority'], 4), -x['count']))
                
                # Convert to dictionary for compatibility
                directional_changes = {
                    item['transition']: {
                        "count": item['count'],
                        "clinical_priority": item['clinical_priority']
                    }
                    for item in directional_changes_list
                }
            else:
                directional_changes = {}
            
            clinical_transitions = {
                "stable_annotations": stable_counts,
                "directional_changes": directional_changes,
                "directional_changes_ordered": directional_changes_list,  # For ordered display
                "total_stable": len(stable_variants),
                "total_changing": len(changing_variants)
            }
        
        # Impact transitions
        impact_transitions = {}
        if len(df_full) > 0:
            stable_impact = df_full[df_full['hg19_impact'] == df_full['hg38_impact']]
            changing_impact = df_full[df_full['hg19_impact'] != df_full['hg38_impact']]
            
            if len(changing_impact) > 0:
                impact_transition_counts = changing_impact.groupby(['hg19_impact', 'hg38_impact']).size()
                impact_changes = {f"{hg19}→{hg38}": count for (hg19, hg38), count in impact_transition_counts.items()}
            else:
                impact_changes = {}
            
            impact_transitions = {
                "stable_count": len(stable_impact),
                "stable_percentage": len(stable_impact) / len(df_full) * 100,
                "transitions": impact_changes
            }
        
        # Clinical coverage analysis 
        clinical_coverage = {}
        if len(df_full) > 0 and 'hg19_clin_sig_normalized' in df_full.columns:
            # Build-specific counts with ordered categories
            hg19_counts = df_full['hg19_clin_sig_normalized'].value_counts().to_dict()
            hg38_counts = df_full['hg38_clin_sig_normalized'].value_counts().to_dict()
            
            # Order categories by clinical interest
            category_order = ['PATHOGENIC', 'BENIGN', 'VUS', 'RISK', 'DRUG_RESPONSE', 'PROTECTIVE', 'OTHER', 'NONE']
            
            # Calculate differences between builds in order
            build_comparison_ordered = []
            for category in category_order:
                hg19_count = hg19_counts.get(category, 0)
                hg38_count = hg38_counts.get(category, 0)
                if hg19_count > 0 or hg38_count > 0:  # Only include categories with data
                    build_comparison_ordered.append({
                        "category": category,
                        "hg19_count": hg19_count,
                        "hg38_count": hg38_count,
                        "difference": hg38_count - hg19_count
                    })
            
            # Also create dictionary for backward compatibility
            build_comparison = {
                item["category"]: {
                    "hg19_count": item["hg19_count"],
                    "hg38_count": item["hg38_count"],
                    "difference": item["difference"]
                }
                for item in build_comparison_ordered
            }
            
            # Overall coverage
            has_hg19_clin = (df_full['hg19_clin_sig_normalized'] != 'NONE').sum()
            has_hg38_clin = (df_full['hg38_clin_sig_normalized'] != 'NONE').sum()
            has_any_clin = ((df_full['hg19_clin_sig_normalized'] != 'NONE') | 
                           (df_full['hg38_clin_sig_normalized'] != 'NONE')).sum()
            
            # Pathogenicity predictions - properly exclude missing values including "-"
            has_sift_hg19 = (df_full['hg19_sift'] != '') & (df_full['hg19_sift'] != '-') & (df_full['hg19_sift'].notna())
            has_sift_hg38 = (df_full['hg38_sift'] != '') & (df_full['hg38_sift'] != '-') & (df_full['hg38_sift'].notna())
            has_polyphen_hg19 = (df_full['hg19_polyphen'] != '') & (df_full['hg19_polyphen'] != '-') & (df_full['hg19_polyphen'].notna())
            has_polyphen_hg38 = (df_full['hg38_polyphen'] != '') & (df_full['hg38_polyphen'] != '-') & (df_full['hg38_polyphen'].notna())
            
            has_sift_any = has_sift_hg19 | has_sift_hg38
            has_polyphen_any = has_polyphen_hg19 | has_polyphen_hg38
            has_predictions = has_sift_any | has_polyphen_any
            
            sift_coverage = has_sift_any.sum()
            polyphen_coverage = has_polyphen_any.sum()
            pred_count = has_predictions.sum()
            
            clinical_coverage = {
                "total_variants": len(df_full),
                "clinical_annotations": {
                    "total_with_annotations": has_any_clin,
                    "percentage_with_annotations": round(has_any_clin / len(df_full) * 100, 1),
                    "hg19_with_annotations": has_hg19_clin,
                    "hg38_with_annotations": has_hg38_clin,
                    "build_comparison": build_comparison,
                    "build_comparison_ordered": build_comparison_ordered,  # For ordered display
                    "hg19_distribution": hg19_counts,
                    "hg38_distribution": hg38_counts
                },
                "pathogenicity_predictions": {
                    "with_predictions": pred_count,
                    "percentage_with_predictions": round(pred_count / len(df_full) * 100, 1),
                    "sift_coverage": sift_coverage,
                    "polyphen_coverage": polyphen_coverage,
                    "sift_percentage": round(sift_coverage / len(df_full) * 100, 1),
                    "polyphen_percentage": round(polyphen_coverage / len(df_full) * 100, 1)
                }
            }
        
        # HGVS analysis 
        hgvs_analysis = {}
        if len(df_full) > 0:
            # CANONICAL HGVSc Match Rate
            if 'CANONICAL_HGVSc_Match' in df_full.columns:
                canonical_match_yes = (df_full['CANONICAL_HGVSc_Match'] == 'YES').sum()
                canonical_match_rate = round(canonical_match_yes / len(df_full) * 100, 1)
            else:
                canonical_match_yes = 0
                canonical_match_rate = 0
            
            # Matched transcripts statistics
            if 'matched_transcript_count' in df_full.columns:
                matched_transcript_counts = df_full['matched_transcript_count'].fillna(0)
                avg_matched_transcripts = round(matched_transcript_counts.mean(), 1)
                total_matched_transcripts = matched_transcript_counts.sum()
            else:
                avg_matched_transcripts = 0
                total_matched_transcripts = 0
            
            # Concordant HGVS Rate (among variants with matched transcripts)
            if 'matched_hgvsc_concordant' in df_full.columns and 'matched_hgvsc_discordant' in df_full.columns:
                # Count concordant and discordant matches
                def count_hgvsc_items(series):
                    """Count items in comma-separated lists, handling empty/NaN values"""
                    if pd.isna(series):
                        return 0
                    items = str(series).strip()
                    if items == '' or items == '-':
                        return 0
                    return len([item.strip() for item in items.split(',') if item.strip()])
                
                concordant_counts = df_full['matched_hgvsc_concordant'].apply(count_hgvsc_items)
                discordant_counts = df_full['matched_hgvsc_discordant'].apply(count_hgvsc_items)
                
                total_concordant = concordant_counts.sum()
                total_discordant = discordant_counts.sum()
                total_matched_hgvsc = total_concordant + total_discordant
                
                concordant_hgvs_rate = round(total_concordant / total_matched_hgvsc * 100, 1) if total_matched_hgvsc > 0 else 0
            else:
                total_concordant = 0
                total_discordant = 0
                total_matched_hgvsc = 0
                concordant_hgvs_rate = 0
            
            hgvs_analysis = {
                "canonical_match": {
                    "variants_with_match": canonical_match_yes,
                    "match_rate_percentage": canonical_match_rate,
                    "total_variants": len(df_full)
                },
                "matched_transcripts": {
                    "average_per_variant": avg_matched_transcripts,
                    "total_matched": int(total_matched_transcripts)
                },
                "hgvsc_concordance": {
                    "total_concordant": int(total_concordant),
                    "total_discordant": int(total_discordant),
                    "total_compared": int(total_matched_hgvsc),
                    "concordant_rate_percentage": concordant_hgvs_rate
                }
            }

        # Functional discordances
        functional_discordances = {}
        if len(df_full) > 0:
            functional_discordances = {
                "same_transcript_changes": (df_full['same_transcript_consequence_changes'] > 0).sum(),
                "gene_changes": (df_full['gene_changes'] > 0).sum(),
                "impact_changes": (df_full['impact_changes'] > 0).sum(),
                "unmatched_consequences": (df_full['unmatched_consequences'] > 0).sum(),
                "clinical_significance_changes": (df_full['clin_sig_change'] != '').sum() if 'clin_sig_change' in df_full.columns else 0,
                "pathogenicity_changes": ((df_full['sift_change'] != '') | (df_full['polyphen_change'] != '')).sum() if 'sift_change' in df_full.columns else 0
            }
        
        # Top variants summary
        top_variants_summary = {}
        if len(df_excel) > 0:
            # Count DISTINCT variants with critical issues
            critical_variants_mask = (
                (df_full['same_transcript_consequence_changes'] > 0) | 
                (df_full['clin_sig_change'].isin(['BENIGN_TO_PATHOGENIC', 'PATHOGENIC_TO_BENIGN', 'VUS_TO_PATHOGENIC']))
            )
            critical_variants_count = critical_variants_mask.sum()
            
            # Count clinical significance changes
            clinical_changes = (df_full['clin_sig_change'] != '').sum() if 'clin_sig_change' in df_full.columns else 0
            
            top_variants_summary = {
                "total_in_excel": len(df_excel),
                "critical_count": len(df_excel[df_excel['Priority_Category'] == 'CRITICAL']) if 'Priority_Category' in df_excel.columns else 0,
                "high_count": len(df_excel[df_excel['Priority_Category'] == 'HIGH']) if 'Priority_Category' in df_excel.columns else 0,
                "moderate_count": len(df_excel[df_excel['Priority_Category'] == 'MODERATE']) if 'Priority_Category' in df_excel.columns else 0,
                "avg_priority_score": round(df_excel['Priority_Score'].mean(), 1) if 'Priority_Score' in df_excel.columns else 0,
                "max_priority_score": round(df_excel['Priority_Score'].max(), 1) if 'Priority_Score' in df_excel.columns else 0,
                "clinical_changes_count": len(df_excel[df_excel['Has_Clinical_Change'] == 'YES']) if 'Has_Clinical_Change' in df_excel.columns else 0,
                "critical_variants_distinct": critical_variants_count
            }
        
        return {
            "metadata": {
                "script": "variant_prioritizer",
                "timestamp": datetime.now().isoformat(),
                "total_analyzed": len(df_full),
                "total_output": len(df_excel)
            },
            "dataset_overview": {
                "total_discordant_variants": len(df_full),
                "variants_in_excel_output": len(df_excel)
            },
            "priority_distribution": priority_distribution,
            "clinical_transitions": clinical_transitions,
            "impact_transitions": impact_transitions,
            "clinical_coverage": clinical_coverage,
            "hgvs_analysis": hgvs_analysis,
            "functional_discordances": functional_discordances,
            "top_variants_summary": top_variants_summary
        }