#!/usr/bin/env python3
"""
Genomic Variant Prioritizer for Clinical Review

DESCRIPTION:
    Creates a ranked, Excel-compatible list of problematic variants for clinical 
    geneticist review. Uses pragmatic variant-level transcript analysis to identify 
    functionally significant discordances between genome builds, focusing on 
    clinical prioritization over exhaustive annotation comparison.

VEP IMPACT-WEIGHTED SCORING SYSTEM:
    Based on official VEP consequence severity hierarchy with clinical evidence override:
    
    CLINICAL EVIDENCE OVERRIDE:
    • 90% score reduction for benign variants (LOW/MODIFIER impact + benign evidence)
    • 2x score boost for pathogenic variants (HIGH impact or pathogenic evidence)
    • Falls back to impact/consequence changes when clinical data unavailable
    
    CORE SCORING:
    • HIGH impact gene changes: +12 points (splice variants, stop_gained, frameshift, etc.)
    • MODERATE impact gene changes: +8 points (missense_variant, inframe indels, etc.)
    • LOW impact gene changes: +4 points (splice_region_variant, synonymous_variant, etc.)
    • MODIFIER-only gene changes: +1 point (intron_variant, UTR variants, etc.)
    • Same transcript consequence changes: +10 points (CRITICAL)
    • Impact level changes: +6 points (IMPORTANT)
    • Clinical significance changes: +10 points (benign ↔ pathogenic)
    • SIFT/PolyPhen changes: +5 points each (deleterious ↔ tolerated)
    • Position/genotype mismatches: +3 points each
    • Reduced gene symbol weight (focus on functional changes over synonyms)

USAGE:
    python variant_prioritizer.py --db-path genomic_analysis.db --output-dir results/
    python variant_prioritizer.py -d data.db -o results/ --force  # Recalculate ignoring cache

CACHING:
    • First run: Always calculates (no cache exists yet)
    • Subsequent runs: Uses cache for faster analysis
    • Use --force flag to recalculate and update cache
    • Cache stores raw VEP analysis only (scores calculated on-the-fly)
"""

import sqlite3
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
from pathlib import Path
import numpy as np
import sys
import re

# Import configuration
from config.constants import VEP_CONSEQUENCE_IMPACT
from config.scoring_config import (
    IMPACT_TRANSITION_SCORES,
    CLINICAL_OVERRIDE,
    BASE_SCORES,
    PRIORITY_CATEGORIES
)
from config.visualization_config import (
    PLOT_COLORS,
    PLOT_STYLE_CONFIG,
    FIGURE_CONFIG
)

# Import utilities  
from utils.clinical_utils import (
    is_pathogenic_clinical_significance,
    is_benign_clinical_significance,
    parse_sift_prediction,
    parse_polyphen_prediction
)
from utils.transcript_utils import (
    normalize_transcript_id,
    extract_genotype_from_alleles
)
from utils.impact_utils import (
    get_impact_numeric_value,
    calculate_impact_transition_magnitude
)
from utils.data_utils import safe_int_convert

# Import visualization
from visualization.plot_generator import PrioritizationPlotter

# Import analysis engine
from analysis.variant_processor import VariantProcessor

# Set visualization style from config
plt.style.use(PLOT_STYLE_CONFIG['style'])
sns.set_palette(PLOT_STYLE_CONFIG['seaborn_palette'])

def connect_database(db_path):
    """Connect to SQLite database and verify structure"""
    if not Path(db_path).exists():
        raise FileNotFoundError(f"Database not found: {db_path}")
    
    conn = sqlite3.connect(db_path)
    
    # Verify required tables
    cursor = conn.cursor()
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
    tables = [row[0] for row in cursor.fetchall()]
    
    required_tables = ['comparison', 'hg19_vep', 'hg38_vep']
    missing_tables = [t for t in required_tables if t not in tables]
    if missing_tables:
        raise ValueError(f"Missing required database tables: {missing_tables}")
    
    return conn

def get_impact_numeric_value(impact):
    """Convert impact to numeric value for magnitude calculation"""
    impact_values = {'HIGH': 4, 'MODERATE': 3, 'LOW': 2, 'MODIFIER': 1}
    return impact_values.get(impact, 0)


def format_for_excel(df):
    """Format dataframe for Excel compatibility with enhanced clinical details and proper genotype extraction"""
    print("Formatting data for Excel...")
    
    # Debug: show what we're working with
    print(f"Input dataframe has {len(df)} rows")
    if len(df) > 0:
        print("First row of data:")
        print(df.iloc[0].to_dict())
    
    # CRITICAL FIX: Reset the index so we have 0-based indexing
    df = df.reset_index(drop=True)
    
    # Create a clean output dataframe
    output_df = pd.DataFrame()
    
    # Basic variant information
    output_df['Rank'] = range(1, len(df) + 1)
    output_df['Priority_Score'] = df['priority_score']
    output_df['Priority_Category'] = df['priority_category']
    output_df['Chromosome'] = df['source_chrom']
    output_df['Position_hg19'] = df['source_pos']
    output_df['Alleles'] = df['source_alleles']
    
    # IMPROVED GENOTYPE EXTRACTION (using existing cached data)
    def extract_genotypes_from_cached_data(row):
        """Extract all 4 genotype components from existing cached data"""
        # Extract hg19 REF/ALT from source_alleles
        source_alleles = row['source_alleles']
        if pd.notna(source_alleles) and '/' in str(source_alleles):
            hg19_parts = str(source_alleles).split('/')
            hg19_ref = hg19_parts[0].strip() if len(hg19_parts) > 0 else ''
            hg19_alt = hg19_parts[1].strip() if len(hg19_parts) > 1 else ''
        elif pd.notna(source_alleles) and ',' in str(source_alleles):
            hg19_parts = str(source_alleles).split(',')
            hg19_ref = hg19_parts[0].strip() if len(hg19_parts) > 0 else ''
            hg19_alt = hg19_parts[1].strip() if len(hg19_parts) > 1 else ''
        else:
            hg19_ref = str(source_alleles) if pd.notna(source_alleles) else ''
            hg19_alt = ''
        
        # Extract hg38 REF/ALT from cached bcftools data
        hg38_ref = row.get('bcftools_hg38_ref', '')
        hg38_alt = row.get('bcftools_hg38_alt', '')
        
        return hg19_ref, hg19_alt, hg38_ref, hg38_alt
    
    # Extract genotypes for all rows
    genotype_data = df.apply(extract_genotypes_from_cached_data, axis=1, result_type='expand')
    genotype_data.columns = ['hg19_ref', 'hg19_alt', 'hg38_ref', 'hg38_alt']
    
    output_df['GT_hg19'] = genotype_data['hg19_ref'] + '/' + genotype_data['hg19_alt']
    output_df['GT_hg38'] = genotype_data['hg38_ref'] + '/' + genotype_data['hg38_alt']
    output_df['Mapping_Status'] = df['mapping_status']
    
    output_df['Position_hg38_CrossMap'] = safe_int_convert(df['liftover_hg38_pos'])
    output_df['Position_hg38_bcftools'] = safe_int_convert(df['bcftools_hg38_pos'])
    output_df['Position_Match'] = df['pos_match'].map({1: 'YES', 0: 'NO'})
    output_df['Position_Difference'] = safe_int_convert(df['pos_difference'])
    
    # Genotype information
    output_df['Genotype_Match'] = df['gt_match'].map({1: 'YES', 0: 'NO'})
    output_df['Strand_Flip'] = df['flip'].fillna('')
    output_df['Ref_Alt_Swap'] = df['swap'].fillna('')
    
    # Enhanced transcript analysis with problematic transcript lists
    output_df['Transcript_Pairs_Analyzed'] = df['transcript_pairs_analyzed']
    output_df['Same_Transcript_Consequence_Changes'] = df['same_transcript_consequence_changes']
    output_df['Same_Consequence_Different_Transcripts'] = df['same_consequence_different_transcripts']
    output_df['Unmatched_Consequences'] = df['unmatched_consequences']
    output_df['Gene_Annotation_Changes'] = df['gene_changes']
    output_df['Impact_Level_Changes'] = df['impact_changes']
    output_df['Problematic_Transcripts_hg19'] = df['problematic_transcripts_hg19']
    output_df['Problematic_Transcripts_hg38'] = df['problematic_transcripts_hg38']
    output_df['Discordance_Summary'] = df['discordance_summary']
    
    # Gene information
    output_df['Gene_hg19'] = df['hg19_gene'].fillna('')
    output_df['Gene_hg38'] = df['hg38_gene'].fillna('')
    output_df['Gene_Match'] = (df['hg19_gene'] == df['hg38_gene']).map({True: 'YES', False: 'NO'})
    
    # VEP consequences (representative - from first transcript if available)
    output_df['Consequence_hg19'] = df['hg19_consequence'].fillna('')
    output_df['Consequence_hg38'] = df['hg38_consequence'].fillna('')
    output_df['Consequence_Match'] = (df['hg19_consequence'] == df['hg38_consequence']).map({True: 'YES', False: 'NO'})
    
    # Impact
    output_df['Impact_hg19'] = df['hg19_impact'].fillna('')
    output_df['Impact_hg38'] = df['hg38_impact'].fillna('')
    output_df['Impact_Match'] = (df['hg19_impact'] == df['hg38_impact']).map({True: 'YES', False: 'NO'})
    
    # Enhanced clinical significance tracking (original + normalized)
    output_df['Clinical_Significance_hg19'] = df['hg19_clin_sig'].fillna('')
    output_df['Clinical_Significance_hg38'] = df['hg38_clin_sig'].fillna('')
    output_df['Clinical_Significance_hg19_Normalized'] = df['hg19_clin_sig_normalized'].fillna('')
    output_df['Clinical_Significance_hg38_Normalized'] = df['hg38_clin_sig_normalized'].fillna('')
    output_df['Clinical_Significance_Change'] = df['clin_sig_change'].fillna('')
    
    # Population frequencies
    output_df['gnomAD_Frequency_hg19'] = df['hg19_gnomad_af'].fillna('')
    output_df['gnomAD_Frequency_hg38'] = df['hg38_gnomad_af'].fillna('')
    
    # Enhanced pathogenicity predictions with change tracking
    output_df['SIFT_hg19'] = df['hg19_sift'].fillna('')
    output_df['SIFT_hg38'] = df['hg38_sift'].fillna('')
    output_df['SIFT_Change'] = df['sift_change'].fillna('')
    output_df['PolyPhen_hg19'] = df['hg19_polyphen'].fillna('')
    output_df['PolyPhen_hg38'] = df['hg38_polyphen'].fillna('')
    output_df['PolyPhen_Change'] = df['polyphen_change'].fillna('')
    
    # Summary flags for quick filtering
    output_df['Has_Position_Issue'] = (df['pos_match'] == 0).map({True: 'YES', False: 'NO'})
    output_df['Has_Genotype_Issue'] = (df['gt_match'] == 0).map({True: 'YES', False: 'NO'})
    output_df['Has_Transcript_Consequence_Issue'] = (df['same_transcript_consequence_changes'] > 0).map({True: 'YES', False: 'NO'})
    output_df['Has_Gene_Issue'] = (df['gene_changes'] > 0).map({True: 'YES', False: 'NO'})
    output_df['Has_Unmatched_Consequences'] = (df['unmatched_consequences'] > 0).map({True: 'YES', False: 'NO'})
    output_df['Has_Clinical_Change'] = (df['clin_sig_change'] != '').map({True: 'YES', False: 'NO'})
    output_df['Has_Pathogenicity_Change'] = ((df['sift_change'] != '') | (df['polyphen_change'] != '')).map({True: 'YES', False: 'NO'})
    
    print(f"Output dataframe created with {len(output_df)} rows")
    print("Sample output:")
    if len(output_df) > 0:
        print(output_df[['Rank', 'Chromosome', 'Position_hg19', 'Gene_hg19', 'Priority_Score', 'Priority_Category']].head().to_string(index=False))
    
    return output_df

def create_summary_statistics(df_full, df_excel, output_dir):
    """Create summary statistics file with clinical evidence-driven analysis details for FULL dataset"""
    print("Creating summary statistics...")
    
    summary_file = output_dir / 'variant_prioritization_summary.txt'
    
    with open(summary_file, 'w') as f:
        f.write("Variant Prioritization Summary - Clinical Evidence-Driven Analysis\n")
        f.write("=" * 70 + "\n\n")
        
        f.write("DATASET OVERVIEW:\n")
        f.write("-" * 20 + "\n")
        f.write(f"Total discordant variants analyzed: {len(df_full):,}\n")
        f.write(f"Variants included in Excel output: {len(df_excel):,}\n\n")
        
        # PRIORITY CATEGORY DISTRIBUTION
        f.write("VARIANT PRIORITIZATION:\n")
        f.write("-" * 25 + "\n")
        if len(df_full) > 0:
            category_counts = df_full['priority_category'].value_counts()
            f.write("Priority category distribution:\n")
            for category, count in category_counts.items():
                pct = count / len(df_full) * 100
                f.write(f"  {category}: {count:,} ({pct:.1f}%)\n")
        f.write("\n")
        
        # CLINICAL SIGNIFICANCE TRANSITION ANALYSIS
        f.write("CLINICAL SIGNIFICANCE TRANSITIONS:\n")
        f.write("-" * 40 + "\n")
        if len(df_full) > 0 and 'hg19_clin_sig_normalized' in df_full.columns:
            # Calculate variants with clinical annotations
            has_hg19_clin = df_full['hg19_clin_sig_normalized'] != 'NONE'
            has_hg38_clin = df_full['hg38_clin_sig_normalized'] != 'NONE'
            has_any_clin = has_hg19_clin | has_hg38_clin
            total_with_clin = has_any_clin.sum()
            
            f.write(f"Total variants with clinical annotations: {total_with_clin:,} ({total_with_clin/len(df_full)*100:.1f}%)\n\n")
            
            # Stable annotations
            stable_variants = df_full[df_full['hg19_clin_sig_normalized'] == df_full['hg38_clin_sig_normalized']]
            f.write(f"Stable annotations: {len(stable_variants):,} variants\n")
            
            stable_counts = stable_variants['hg19_clin_sig_normalized'].value_counts()
            for category in ['PATHOGENIC', 'BENIGN', 'VUS', 'RISK', 'DRUG_RESPONSE', 'PROTECTIVE', 'OTHER', 'NONE']:
                if category in stable_counts:
                    f.write(f"  Stable {category}: {stable_counts[category]:,}\n")
            
            # Directional changes
            changing_variants = df_full[df_full['hg19_clin_sig_normalized'] != df_full['hg38_clin_sig_normalized']]
            f.write(f"\nDirectional changes: {len(changing_variants):,} variants\n")
            
            # Key clinical transitions
            pathogenic_to_benign = len(df_full[(df_full['hg19_clin_sig_normalized'] == 'PATHOGENIC') & 
                                             (df_full['hg38_clin_sig_normalized'] == 'BENIGN')])
            benign_to_pathogenic = len(df_full[(df_full['hg19_clin_sig_normalized'] == 'BENIGN') & 
                                             (df_full['hg38_clin_sig_normalized'] == 'PATHOGENIC')])
            vus_transitions = len(df_full[((df_full['hg19_clin_sig_normalized'] == 'VUS') & 
                                        (df_full['hg38_clin_sig_normalized'].isin(['PATHOGENIC', 'BENIGN']))) |
                                       ((df_full['hg38_clin_sig_normalized'] == 'VUS') & 
                                        (df_full['hg19_clin_sig_normalized'].isin(['PATHOGENIC', 'BENIGN'])))])
            
            if pathogenic_to_benign > 0:
                f.write(f"  Pathogenic → Benign: {pathogenic_to_benign:,} (affects interpretation)\n")
            if benign_to_pathogenic > 0:
                f.write(f"  Benign → Pathogenic: {benign_to_pathogenic:,} (affects interpretation)\n")
            if vus_transitions > 0:
                f.write(f"  VUS transitions: {vus_transitions:,}\n")
        f.write("\n")
        
        # IMPACT LEVEL TRANSITION ANALYSIS
        f.write("IMPACT LEVEL TRANSITIONS:\n")
        f.write("-" * 30 + "\n")
        if len(df_full) > 0:
            # Stable impact levels
            stable_impact = df_full[df_full['hg19_impact'] == df_full['hg38_impact']]
            f.write(f"Variants with stable impact: {len(stable_impact):,} ({len(stable_impact)/len(df_full)*100:.1f}%)\n")
            
            # Clinically significant transitions
            changing_impact = df_full[df_full['hg19_impact'] != df_full['hg38_impact']]
            clinically_significant_changes = changing_impact[
                (changing_impact['hg19_impact'].isin(['HIGH', 'MODERATE'])) |
                (changing_impact['hg38_impact'].isin(['HIGH', 'MODERATE']))
            ]
            f.write(f"Clinically significant transitions: {len(clinically_significant_changes):,} variants\n")
            
            # Specific transitions
            high_moderate = len(df_full[((df_full['hg19_impact'] == 'HIGH') & (df_full['hg38_impact'] == 'MODERATE')) |
                                      ((df_full['hg19_impact'] == 'MODERATE') & (df_full['hg38_impact'] == 'HIGH'))])
            high_low_modifier = len(df_full[((df_full['hg19_impact'] == 'HIGH') & (df_full['hg38_impact'].isin(['LOW', 'MODIFIER']))) |
                                          ((df_full['hg38_impact'] == 'HIGH') & (df_full['hg19_impact'].isin(['LOW', 'MODIFIER'])))])
            low_modifier = len(df_full[((df_full['hg19_impact'] == 'LOW') & (df_full['hg38_impact'] == 'MODIFIER')) |
                                     ((df_full['hg19_impact'] == 'MODIFIER') & (df_full['hg38_impact'] == 'LOW'))])
            
            if high_moderate > 0:
                f.write(f"  HIGH ↔ MODERATE: {high_moderate:,} variants\n")
            if high_low_modifier > 0:
                f.write(f"  HIGH ↔ LOW/MODIFIER: {high_low_modifier:,} variants\n")
            if low_modifier > 0:
                f.write(f"Annotation updates (LOW ↔ MODIFIER): {low_modifier:,} variants\n")
        f.write("\n")

        # DEBUG: Check impact data
        f.write("\nDEBUG - Impact Analysis:\n")
        if len(df_full) > 0:
            # Check if impact data exists
            hg19_impacts = df_full['hg19_impact'].value_counts()
            hg38_impacts = df_full['hg38_impact'].value_counts()
            f.write(f"HG19 impact distribution: {dict(hg19_impacts)}\n")
            f.write(f"HG38 impact distribution: {dict(hg38_impacts)}\n")
            
            # Check impact_changes field
            impact_changes_count = (df_full['impact_changes'] > 0).sum()
            f.write(f"Variants with impact_changes > 0: {impact_changes_count}\n")
            
            # Check specific variants with clinical sig changes
            clin_changing = df_full[df_full['hg19_clin_sig_normalized'] != df_full['hg38_clin_sig_normalized']]
            if len(clin_changing) > 0:
                f.write(f"Clinical significance changing variants impact data:\n")
                for idx, row in clin_changing.head(5).iterrows():
                    f.write(f"  Variant {row['source_chrom']}:{row['source_pos']}: "
                        f"hg19_impact={row['hg19_impact']}, hg38_impact={row['hg38_impact']}, "
                        f"impact_changes={row['impact_changes']}\n")
        
        # CLINICAL DATA COVERAGE ASSESSMENT
        f.write("CLINICAL DATA COVERAGE:\n")
        f.write("-" * 25 + "\n")
        if len(df_full) > 0:
            # Clinical significance data
            if 'hg19_clin_sig_normalized' in df_full.columns:
                has_clin_data = (df_full['hg19_clin_sig_normalized'] != 'NONE') | (df_full['hg38_clin_sig_normalized'] != 'NONE')
                clin_data_count = has_clin_data.sum()
                f.write(f"Variants with clinical significance data: {clin_data_count:,} ({clin_data_count/len(df_full)*100:.1f}%)\n")
            
            # Pathogenicity predictions
            has_sift = (df_full['hg19_sift'] != '') | (df_full['hg38_sift'] != '')
            has_polyphen = (df_full['hg19_polyphen'] != '') | (df_full['hg38_polyphen'] != '')
            has_predictions = has_sift | has_polyphen
            pred_count = has_predictions.sum()
            f.write(f"Variants with pathogenicity predictions: {pred_count:,} ({pred_count/len(df_full)*100:.1f}%)\n")
            
            # Clinical evidence distribution
            if 'hg19_clin_sig_normalized' in df_full.columns:
                f.write("\nClinical evidence distribution:\n")
                all_clin_values = pd.concat([df_full['hg19_clin_sig_normalized'], df_full['hg38_clin_sig_normalized']])
                clin_distribution = all_clin_values.value_counts()
                
                high_confidence = clin_distribution.get('PATHOGENIC', 0) + clin_distribution.get('BENIGN', 0)
                uncertain = clin_distribution.get('VUS', 0)
                risk_drug = clin_distribution.get('RISK', 0) + clin_distribution.get('DRUG_RESPONSE', 0)
                no_data = clin_distribution.get('NONE', 0)
                
                f.write(f"  High confidence (PATHOGENIC/BENIGN): {high_confidence:,}\n")
                f.write(f"  Uncertain (VUS): {uncertain:,}\n")
                f.write(f"  Risk/Drug response: {risk_drug:,}\n")
                f.write(f"  No clinical data: {no_data:,}\n")
        f.write("\n")
        
        # FUNCTIONAL ISSUE BREAKDOWN
        f.write("FUNCTIONAL DISCORDANCES:\n")
        f.write("-" * 25 + "\n")
        if len(df_full) > 0:
            same_transcript_issues = (df_full['same_transcript_consequence_changes'] > 0).sum()
            gene_issues = (df_full['gene_changes'] > 0).sum()
            impact_issues = (df_full['impact_changes'] > 0).sum()
            unmatched_issues = (df_full['unmatched_consequences'] > 0).sum()
            
            f.write(f"Same transcript, different consequences: {same_transcript_issues:,} variants (CRITICAL)\n")
            f.write(f"Gene annotation changes: {gene_issues:,} variants\n")
            f.write(f"Impact level changes: {impact_issues:,} variants\n")
            f.write(f"Unmatched consequences: {unmatched_issues:,} variants\n")
        f.write("\n")
        
        # SPECIFIC TRANSITION MATRICES
        f.write("KEY CLINICAL SIGNIFICANCE TRANSITIONS:\n")
        f.write("-" * 40 + "\n")
        if len(df_full) > 0 and 'hg19_clin_sig_normalized' in df_full.columns:
            transition_matrix = pd.crosstab(df_full['hg19_clin_sig_normalized'], df_full['hg38_clin_sig_normalized'], 
                                          rownames=['hg19'], colnames=['hg38'])
            
            # Show key transitions
            key_categories = ['PATHOGENIC', 'BENIGN', 'VUS', 'NONE']
            f.write("Transition matrix (rows=hg19, columns=hg38):\n")
            f.write(f"{'':12}")
            for cat in key_categories:
                if cat in transition_matrix.columns:
                    f.write(f"{cat:>12}")
            f.write("\n")
            
            for hg19_cat in key_categories:
                if hg19_cat in transition_matrix.index:
                    f.write(f"{hg19_cat:12}")
                    for hg38_cat in key_categories:
                        if hg38_cat in transition_matrix.columns:
                            value = transition_matrix.loc[hg19_cat, hg38_cat] if hg38_cat in transition_matrix.columns else 0
                            f.write(f"{value:>12}")
                    f.write("\n")
        f.write("\n")
        
        # QUALITY RECOMMENDATIONS
        f.write("CLINICAL WORKFLOW RECOMMENDATIONS:\n")
        f.write("-" * 35 + "\n")
        
        f.write("CRITICAL VARIANTS:\n")
        f.write("• Same transcript with different consequences\n")
        f.write("• Potential functional impact on protein products\n")
        f.write("• Manual review and potential validation recommended\n\n")
        
        f.write("HIGH PRIORITY VARIANTS:\n")
        f.write("• Clinical significance changes affecting interpretation\n")
        f.write("• Significant impact level transitions\n")
        f.write("• Priority review within clinical workflow\n\n")
        
        f.write("MODERATE PRIORITY VARIANTS:\n")
        f.write("• Pathogenicity prediction changes\n")
        f.write("• Clinically relevant gene changes\n")
        f.write("• Standard review process\n\n")
        
        f.write("INVESTIGATE VARIANTS:\n")
        f.write("• Unmatched consequences requiring investigation\n")
        f.write("• Non-significant transitions\n")
        f.write("• Secondary review or automated filtering\n\n")
        
        f.write("LOW PRIORITY VARIANTS:\n")
        f.write("• Benign variants with minimal changes\n")
        f.write("• Annotation updates with no clinical impact\n")
        f.write("• Informational only\n\n")
        
        # SCORING METHODOLOGY
        f.write("VARIANT DISCREPANCY SCORING METHODOLOGY:\n")
        f.write("-" * 45 + "\n")
        f.write("Clinical evidence drives prioritization over annotation model differences:\n\n")
        
        f.write("MAGNITUDE-BASED IMPACT SCORING:\n")
        f.write("• HIGH ↔ MODERATE transitions: +10 points (clinically significant)\n")
        f.write("• HIGH ↔ LOW/MODIFIER transitions: +12-15 points (highly significant)\n")
        f.write("• MODERATE ↔ LOW/MODIFIER transitions: +8-10 points (significant)\n")
        f.write("• LOW ↔ MODIFIER transitions: +1 point (annotation noise)\n\n")
        
        f.write("CLINICAL EVIDENCE OVERRIDE:\n")
        f.write("• 90% score reduction for benign variants (LOW/MODIFIER + benign evidence)\n")
        f.write("• 2x score boost for pathogenic variants (HIGH impact or pathogenic evidence)\n")
        f.write("• Falls back to impact/consequence changes when clinical data unavailable\n\n")
        
        f.write("CORE SCORING COMPONENTS:\n")
        f.write("• Same transcript consequence changes: +10 points (CRITICAL)\n")
        f.write("• Clinical significance changes: +5-10 points\n")
        f.write("• SIFT/PolyPhen changes: +5 points each\n")
        f.write("• Gene changes: Impact-weighted (+0.1-8 points)\n")
        f.write("• Position/genotype issues: +2-3 points each\n\n")
        
        f.write("DATA PROCESSING NOTES:\n")
        f.write("• All coordinates use VEP normalization (SNVs: original, Indels: original+1)\n")
        f.write("• Clinical significance normalized to 8 categories (PATHOGENIC, BENIGN, VUS, etc.)\n")
        f.write("• VEP analysis results cached for faster subsequent runs\n")
        f.write("• Priority scores calculated fresh for easy recalibration\n")

    print(f"✓ Enhanced summary statistics saved to: {summary_file}")

def main():
    parser = argparse.ArgumentParser(
        description='Create prioritized variant list with clinical evidence-driven analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
USAGE EXAMPLES:
    python variant_prioritizer.py --db-path genomic_analysis.db --output-dir results/
    python variant_prioritizer.py -d genomic_analysis.db -o results/ --max-variants 5000 --no-plots
    python variant_prioritizer.py --db-path /path/to/data.db --output-dir /path/to/output/ --min-score 5
    python variant_prioritizer.py -d data.db -o results/ --force  # Recalculate ignoring cache

OUTPUT FILES:
    • prioritized_variants.csv - Excel-compatible ranked variant list (top variants only)
    • variant_prioritization_plots.png - Visual analysis plots (full dataset)
    • variant_prioritization_summary.txt - Detailed summary report (full dataset)
    • variant_analysis_cache.pkl - Cached VEP analysis results (no scores)

CACHING BEHAVIOR:
    • First run: Always calculates (no cache exists yet)
    • Subsequent runs: Uses cache for faster analysis
    • Use --force flag to recalculate VEP analysis and update cache
    • Priority scores calculated on-the-fly for easy recalibration

CLINICAL EVIDENCE-DRIVEN SCORING:
    • 90% score reduction for benign variants (LOW/MODIFIER + benign evidence)
    • 2x score boost for pathogenic variants (HIGH impact or pathogenic evidence)
    • Reduced gene symbol weight (focus on functional changes over synonyms)
    • Clinical significance changes: +10 points (benign ↔ pathogenic)
    • SIFT/PolyPhen changes: +5 points each (deleterious ↔ tolerated)
    • Falls back to impact/consequence changes when clinical data unavailable
        """
    )
    
    parser.add_argument('--db-path', '-d', required=True,
                       help='SQLite database file path')
    parser.add_argument('--output-dir', '-o', default='prioritization_output',
                       help='Output directory for CSV file (default: prioritization_output)')
    parser.add_argument('--max-variants', '-m', type=int, default=10000,
                       help='Maximum number of variants to output (default: 10000)')
    parser.add_argument('--min-score', '-s', type=int, default=1,
                       help='Minimum priority score to include (default: 1)')
    parser.add_argument('--no-plots', action='store_true',
                       help='Skip plot generation (faster for large datasets)')
    parser.add_argument('--force', action='store_true',
                       help='Force recalculation of VEP analysis (ignore cache)')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Enable verbose output')
    
    args = parser.parse_args()
    
    db_path = args.db_path
    output_dir = Path(args.output_dir)
    print(f"✓ Using database: {db_path}")
    print(f"✓ Output directory: {output_dir}")
    
    if args.force:
        print("✓ Force recalculation enabled (will ignore cache)")
    
    # Create output directory
    output_dir.mkdir(exist_ok=True, parents=True)
    
    create_plots = not args.no_plots
    
    try:
        print(f"\n=== STARTING VARIANT PRIORITIZATION ===")
        print("Using clinical evidence-driven scoring with benign suppression and pathogenic boost")
        
        # Connect to database
        conn = connect_database(db_path)
        print(f"✓ Connected to database: {db_path}")


        # Calculate variant scores - with caching (VEP analysis only)
        cache_file = output_dir / 'variant_analysis_cache.pkl'
        processor = VariantProcessor()
        df_full = processor.process_all_variants(conn, cache_file, args.force)
                
        if len(df_full) == 0:
            print("No problematic variants found.")
            return
        
        # Sort by priority score (highest first)
        df_full = df_full.sort_values('priority_score', ascending=False)
        
        # Show breakdown of FULL dataset
        print(f"\nFull dataset breakdown (all {len(df_full):,} discordant variants):")
        same_transcript_critical_full = (df_full['same_transcript_consequence_changes'] > 0).sum()
        gene_change_critical_full = (df_full['gene_changes'] > 0).sum()
        unmatched_critical_full = (df_full['unmatched_consequences'] > 0).sum()
        clinical_changes_full = (df_full['clin_sig_change'] != '').sum()
        pathogenicity_changes_full = ((df_full['sift_change'] != '') | (df_full['polyphen_change'] != '')).sum()
        
        print(f"  Same transcript consequence changes: {same_transcript_critical_full:,}")
        print(f"  Gene annotation changes: {gene_change_critical_full:,}")
        print(f"  Unmatched consequences: {unmatched_critical_full:,}")
        print(f"  Clinical significance changes: {clinical_changes_full:,}")
        print(f"  Pathogenicity prediction changes: {pathogenicity_changes_full:,}")
        
        # Show priority category breakdown
        print(f"\nPriority category breakdown:")
        category_counts = df_full['priority_category'].value_counts()
        for category, count in category_counts.items():
            pct = count / len(df_full) * 100
            print(f"  {category}: {count:,} ({pct:.1f}%)")
        
        # Filter by minimum score and maximum variants for EXCEL OUTPUT ONLY
        df_excel = df_full[df_full['priority_score'] >= args.min_score]
        df_excel = df_excel.head(args.max_variants)
        
        print(f"\nSelected {len(df_excel):,} variants for Excel output (score >= {args.min_score}, max {args.max_variants:,})")
        
        
        # Create prioritization plots using FULL dataset
        if create_plots:
            plotter = PrioritizationPlotter(PLOT_COLORS, FIGURE_CONFIG)
            plotter.create_all_plots(df_full, conn, output_dir)
        
        # Format EXCEL subset for Excel
        output_df = format_for_excel(df_excel)
        
        # Save to CSV
        csv_file = output_dir / 'prioritized_variants.csv'
        output_df.to_csv(csv_file, index=False, encoding='utf-8')
        
        # Create summary statistics using BOTH full dataset and Excel subset
        create_summary_statistics(df_full, output_df, output_dir)
        
        conn.close()
        
        print(f"\n=== VARIANT PRIORITIZATION COMPLETED ===")
        print(f"✓ Output files created:")
        print(f"  • {csv_file}")
        if create_plots:
            print(f"  • {output_dir / 'variant_prioritization_plots.png'}")
        print(f"  • {output_dir / 'variant_prioritization_summary.txt'}")
        print(f"\nTop 5 priority variants (Excel output):")
        display_cols = ['Rank', 'Chromosome', 'Position_hg19', 'Gene_hg19', 'Priority_Score', 'Priority_Category']
        print(output_df[display_cols].head().to_string(index=False))
        
        if same_transcript_critical_full > 0:
            print(f"\n⚠️  ATTENTION: {same_transcript_critical_full:,} variants with CRITICAL transcript changes found in full dataset")
            print(f"   → Filter by Priority_Category = 'CRITICAL' for priority review")
        
        if clinical_changes_full > 0:
            print(f"\n🏥 CLINICAL: {clinical_changes_full:,} variants with clinical significance changes")
            print(f"   → Filter by Clinical_Significance_Change column for clinical review")
        
        print(f"\nAll coordinates are VEP-normalized (SNVs: original, Indels: original+1)")
        print(f"\nNOTE: Plots and summaries reflect FULL dataset ({len(df_full):,} variants)")
        print(f"      Excel file contains top {len(df_excel):,} variants only")
        print(f"      Clinical evidence-driven scoring with benign suppression applied")
        
    except Exception as e:
        print(f"Error: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1

if __name__ == "__main__":
    main()
