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

from utils.summary_utils import SummaryDataCalculator
from utils.data_utils import format_consequence_relationship
import json

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

def format_consequence_relationship_for_csv(row):
    return format_consequence_relationship(
        row.get('consequence_relationship', 'unknown'),
        row.get('hg19_consequences_set', ''),
        row.get('hg38_consequences_set', '')
    )

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


def create_clinical_csv_output(df, output_dir, max_variants=10000):
    """Create clinical evidence-focused CSV output for variant prioritization"""
    
    if len(df) == 0:
        print("No variants to output")
        return pd.DataFrame()
    
    print(f"Creating clinical evidence-focused CSV output for top {min(len(df), max_variants):,} variants...")
    
    # Take top variants by priority score
    output_df = df.head(max_variants).copy()
    
    # Add rank column
    output_df['Rank'] = range(1, len(output_df) + 1)
    
    # Clinical evidence-focused column selection and renaming
    column_mapping = {
        'Rank': 'Rank',
        'source_chrom': 'Chromosome',
        'source_pos': 'Position_hg19',
        'bcftools_hg38_pos': 'Position_hg38',
        'source_alleles': 'Alleles',
        'hg19_gene': 'Gene_hg19',
        'hg38_gene': 'Gene_hg38',
        'hg19_clin_sig_normalized': 'Clinical_Significance_hg19',
        'hg38_clin_sig_normalized': 'Clinical_Significance_hg38',
        'hg19_impact': 'Impact_hg19',
        'hg38_impact': 'Impact_hg38',
        'hg19_high_impact_consequences': 'High_Impact_Consequences_hg19',  
        'hg38_high_impact_consequences': 'High_Impact_Consequences_hg38',  
        'hg19_sift': 'SIFT_hg19',
        'hg38_sift': 'SIFT_hg38',
        'hg19_polyphen': 'PolyPhen_hg19',
        'hg38_polyphen': 'PolyPhen_hg38',
        'priority_score': 'Priority_Score',
        'priority_category': 'Priority_Category',
        'discordance_summary': 'Discordance_Summary',
        'transcript_relationship': 'Transcript_Relationship',
        'consequence_relationship': 'Consequence_Relationship',
        'problematic_transcripts_hg19': 'Problematic_Transcripts_hg19',
        'problematic_transcripts_hg38': 'Problematic_Transcripts_hg38',
        'same_transcript_consequence_changes': 'Transcript_Changes',
        'gene_changes': 'Gene_Changes',
        'impact_changes': 'Impact_Changes',
        'clin_sig_change': 'Clinical_Change_Direction',
        'sift_change': 'SIFT_Change',
        'polyphen_change': 'PolyPhen_Change',
        'pos_match': 'Position_Match',
        'gt_match': 'Genotype_Match',
        'mapping_status': 'Mapping_Status'
    }
    
    # Create output dataframe with renamed columns
    output_columns = []
    for old_col, new_col in column_mapping.items():
        if old_col in output_df.columns:
            output_df[new_col] = output_df[old_col]
            output_columns.append(new_col)
    
    # Add clinical change indicator
    if 'Clinical_Change_Direction' in output_columns:
        output_df['Has_Clinical_Change'] = output_df['Clinical_Change_Direction'].apply(
            lambda x: 'YES' if x and x != '' and 'STABLE_' not in str(x) else 'NO'
        )
        output_columns.append('Has_Clinical_Change')
    
    # Add impact change indicator
    if 'Impact_hg19' in output_columns and 'Impact_hg38' in output_columns:
        output_df['Has_Impact_Change'] = (output_df['Impact_hg19'] != output_df['Impact_hg38']).apply(
            lambda x: 'YES' if x else 'NO'
        )
        output_columns.append('Has_Impact_Change')
    
    # Add consequence change indicator based on relationship type
    if 'Consequence_Relationship' in output_columns:
        output_df['Has_Consequence_Change'] = output_df['Consequence_Relationship'].apply(
            lambda x: 'YES' if x in ['disjoint_consequences', 'partial_overlap_consequences', 'hg19_subset_of_hg38', 'hg38_subset_of_hg19'] else 'NO'
        )
        output_columns.append('Has_Consequence_Change')

    # Add computed Consequence_Change column
    if 'consequence_relationship' in df.columns:
        def format_consequence_change(row):
            return format_consequence_relationship_for_csv(row)
        
        output_df['Consequence_Change'] = df.apply(format_consequence_change, axis=1)
        # Insert right after Consequence_Relationship
        if 'Consequence_Relationship' in output_columns:
            insert_index = output_columns.index('Consequence_Relationship') + 1
            output_columns.insert(insert_index, 'Consequence_Change')
        else:
            output_columns.append('Consequence_Change')
    
    # Select only the columns we want in the final output
    final_df = output_df[output_columns].copy()
    
    # Clean up data for clinical compatibility
    final_df = final_df.fillna('')
    
    # Convert numeric columns to appropriate types
    numeric_columns = ['Priority_Score', 'Transcript_Changes', 'Gene_Changes', 'Impact_Changes']
    for col in numeric_columns:
        if col in final_df.columns:
            final_df[col] = pd.to_numeric(final_df[col], errors='coerce').fillna(0)
    
    # Save to CSV
    output_file = output_dir / 'prioritized_variants.csv'
    final_df.to_csv(output_file, index=False)
    
    print(f"✓ Clinical evidence CSV saved to: {output_file}")
    print(f"  - Total variants: {len(final_df):,}")
    print(f"  - Columns: {len(final_df.columns)}")
    
    # Print priority distribution
    if 'Priority_Category' in final_df.columns:
        category_counts = final_df['Priority_Category'].value_counts()
        print("  - Priority distribution:")
        for category, count in category_counts.items():
            print(f"    {category}: {count:,}")
    
    return final_df


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
        
        # CLINICAL SIGNIFICANCE TRANSITION ANALYSIS (UPDATED - removed redundant line)
        f.write("CLINICAL SIGNIFICANCE TRANSITIONS:\n")
        f.write("-" * 40 + "\n")
        if len(df_full) > 0 and 'hg19_clin_sig_normalized' in df_full.columns:
            # Stable annotations
            stable_variants = df_full[df_full['hg19_clin_sig_normalized'] == df_full['hg38_clin_sig_normalized']]
            f.write(f"Stable annotations: {len(stable_variants):,} variants\n")
            
            stable_counts = stable_variants['hg19_clin_sig_normalized'].value_counts()
            for category in ['PATHOGENIC', 'BENIGN', 'VUS', 'RISK', 'DRUG_RESPONSE', 'PROTECTIVE', 'OTHER', 'NONE']:
                if category in stable_counts:
                    f.write(f"  Stable {category}: {stable_counts[category]:,}\n")
            
            # Directional changes (hg19→hg38)
            changing_variants = df_full[df_full['hg19_clin_sig_normalized'] != df_full['hg38_clin_sig_normalized']]
            f.write(f"\nClinical significance transitions (hg19→hg38): {len(changing_variants):,} variants\n")
            
            # Show specific directional transitions
            if len(changing_variants) > 0:
                transition_counts = changing_variants.groupby(['hg19_clin_sig_normalized', 'hg38_clin_sig_normalized']).size()
                
                # Critical transitions first
                critical_transitions = [
                    ('BENIGN', 'PATHOGENIC'),
                    ('PATHOGENIC', 'BENIGN'), 
                    ('VUS', 'PATHOGENIC'),
                    ('PATHOGENIC', 'VUS')
                ]
                
                for hg19_cat, hg38_cat in critical_transitions:
                    if (hg19_cat, hg38_cat) in transition_counts:
                        count = transition_counts[(hg19_cat, hg38_cat)]
                        f.write(f"  {hg19_cat}→{hg38_cat}: {count:,} variants (critical)\n")
                
                # Other transitions
                other_transitions = [(hg19, hg38) for hg19, hg38 in transition_counts.index 
                                   if (hg19, hg38) not in critical_transitions]
                
                for hg19_cat, hg38_cat in other_transitions:
                    count = transition_counts[(hg19_cat, hg38_cat)]
                    f.write(f"  {hg19_cat}→{hg38_cat}: {count:,} variants\n")
        f.write("\n")
        
        # IMPACT LEVEL TRANSITION ANALYSIS
        f.write("IMPACT LEVEL TRANSITIONS:\n")
        f.write("-" * 30 + "\n")
        if len(df_full) > 0:
            # Stable impact levels
            stable_impact = df_full[df_full['hg19_impact'] == df_full['hg38_impact']]
            f.write(f"Variants with stable impact: {len(stable_impact):,} ({len(stable_impact)/len(df_full)*100:.1f}%)\n")
            
            # Impact transitions (hg19→hg38)
            changing_impact = df_full[df_full['hg19_impact'] != df_full['hg38_impact']]
            if len(changing_impact) > 0:
                f.write(f"Impact level transitions (hg19→hg38): {len(changing_impact):,} variants\n")
                
                impact_transition_counts = changing_impact.groupby(['hg19_impact', 'hg38_impact']).size()
                for (hg19_impact, hg38_impact), count in impact_transition_counts.items():
                    f.write(f"  {hg19_impact}→{hg38_impact}: {count:,} variants\n")
            else:
                f.write(f"Impact level transitions (hg19→hg38): 0 variants\n")
        f.write("\n")

        # CLINICAL DATA COVERAGE ASSESSMENT (UPDATED - build-specific table)
        f.write("CLINICAL DATA COVERAGE:\n")
        f.write("-" * 25 + "\n")
        if len(df_full) > 0:
            # Total clinical annotations
            if 'hg19_clin_sig_normalized' in df_full.columns:
                has_hg19_clin = df_full['hg19_clin_sig_normalized'] != 'NONE'
                has_hg38_clin = df_full['hg38_clin_sig_normalized'] != 'NONE'
                has_any_clin = has_hg19_clin | has_hg38_clin
                total_with_clin = has_any_clin.sum()
                f.write(f"Total variants with clinical annotations: {total_with_clin:,} ({total_with_clin/len(df_full)*100:.1f}%)\n\n")
            
            # Pathogenicity predictions
            has_sift = (df_full['hg19_sift'] != '') | (df_full['hg38_sift'] != '')
            has_polyphen = (df_full['hg19_polyphen'] != '') | (df_full['hg38_polyphen'] != '')
            has_predictions = has_sift | has_polyphen
            pred_count = has_predictions.sum()
            f.write(f"Variants with pathogenicity predictions: {pred_count:,} ({pred_count/len(df_full)*100:.1f}%)\n\n")
            
            # Clinical evidence distribution by build (TABLE FORMAT)
            if 'hg19_clin_sig_normalized' in df_full.columns:
                f.write("Clinical evidence distribution by build:\n")
                f.write(f"{'Category':<15} {'hg19':<8} {'hg38':<8}\n")
                f.write("-" * 32 + "\n")
                
                categories = ['PATHOGENIC', 'BENIGN', 'VUS', 'RISK', 'DRUG_RESPONSE', 'PROTECTIVE', 'OTHER', 'NONE']
                hg19_counts = df_full['hg19_clin_sig_normalized'].value_counts()
                hg38_counts = df_full['hg38_clin_sig_normalized'].value_counts()
                
                for category in categories:
                    hg19_count = hg19_counts.get(category, 0)
                    hg38_count = hg38_counts.get(category, 0)
                    f.write(f"{category:<15} {hg19_count:<8} {hg38_count:<8}\n")
                
                f.write("-" * 32 + "\n")
                f.write(f"{'Total':<15} {len(df_full):<8} {len(df_full):<8}\n")
        f.write("\n")
        
        # FUNCTIONAL ISSUE BREAKDOWN
        f.write("FUNCTIONAL DISCORDANCES:\n")
        f.write("-" * 25 + "\n")
        if len(df_full) > 0:
            same_transcript_issues = (df_full['same_transcript_consequence_changes'] > 0).sum()
            gene_issues = (df_full['gene_changes'] > 0).sum()
            impact_issues = (df_full['impact_changes'] > 0).sum()
            unmatched_issues = (df_full['unmatched_consequences'] > 0).sum()
            
            f.write(f"Same transcript, different consequences: {same_transcript_issues:,} variants (HIGH)\n")
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
        f.write("• Clinical interpretation changes (hg19→hg38: PATHOGENIC↔BENIGN, VUS→PATHOGENIC)\n")
        f.write("• High impact transitions affecting protein function\n")
        f.write("• Immediate clinical review and potential validation recommended\n\n")
        
        f.write("HIGH PRIORITY VARIANTS:\n")
        f.write("• Functionally significant impact transitions\n")
        f.write("• Other clinical significance changes\n")
        f.write("• Same transcript consequence changes\n")
        f.write("• Priority review within clinical workflow\n\n")
        
        f.write("MODERATE PRIORITY VARIANTS:\n")
        f.write("• Pathogenicity prediction changes\n")
        f.write("• Clinically relevant gene changes\n")
        f.write("• Standard review process\n\n")
        
        f.write("LOW PRIORITY VARIANTS:\n")
        f.write("• Technical liftover issues\n")
        f.write("• Annotation differences between builds\n")
        f.write("• Secondary review or automated filtering\n\n")
        
        f.write("INVESTIGATE VARIANTS:\n")
        f.write("• Unclear cases requiring investigation\n")
        f.write("• Mixed evidence patterns\n")
        f.write("• Case-by-case review\n\n")
        
        # SCORING METHODOLOGY
        f.write("VARIANT DISCREPANCY SCORING METHODOLOGY:\n")
        f.write("-" * 45 + "\n")
        f.write("CLINICAL EVIDENCE-FIRST PRIORITIZATION:\n")
        f.write("Clinical significance changes drive prioritization over VEP annotation differences.\n")
        f.write("Reduced weight for transcript mismatches (common annotation noise between builds).\n\n")
        
        f.write("PRIORITY CATEGORIES:\n")
        f.write("• CRITICAL: Clinical interpretation changes (hg19→hg38: PATHOGENIC↔BENIGN, VUS→PATHOGENIC)\n")
        f.write("            OR high impact transitions (HIGH↔MODERATE/LOW/MODIFIER)\n")
        f.write("• HIGH: Functionally significant (moderate impact transitions, other clinical changes,\n")
        f.write("        same transcript consequence changes)\n")
        f.write("• MODERATE: Prediction changes (SIFT/PolyPhen) and clinically relevant gene changes\n")
        f.write("• LOW: Technical issues (position/genotype) and annotation differences\n")
        f.write("• INVESTIGATE: Unclear cases requiring further review\n\n")
        
        f.write("CLINICAL EVIDENCE-FIRST SCORING WEIGHTS:\n")
        f.write("• Clinical significance directional changes: +15-20 points (CRITICAL)\n")
        f.write("• High impact transitions (HIGH involved): +15-20 points (CRITICAL)\n")
        f.write("• Moderate impact transitions: +10-12 points (HIGH)\n")
        f.write("• Same transcript consequence changes: +6 points (HIGH - demoted)\n")
        f.write("• SIFT/PolyPhen changes: +5 points each (MODERATE)\n")
        f.write("• Gene changes: Impact-weighted (+0.1-4 points, conditional)\n")
        f.write("• Unmatched consequences: +1 point (LOW - demoted from +4)\n")
        f.write("• Different transcripts: +0.5 points (LOW - demoted from +3)\n")
        f.write("• Position/genotype issues: +2-3 points each (technical)\n\n")
        
        f.write("KEY CHANGES FROM PREVIOUS VERSION:\n")
        f.write("• Clinical significance changes now highest priority (was secondary)\n")
        f.write("• Unmatched consequences demoted from +4 to +1 (annotation noise)\n")
        f.write("• Same transcript changes demoted from +10 to +6 (less critical than clinical)\n")
        f.write("• Different transcripts demoted from +3 to +0.5 (minimal clinical impact)\n")
        f.write("• Impact transitions weighted by clinical significance (HIGH > MODERATE > LOW)\n\n")
        
        f.write("CLINICAL EVIDENCE OVERRIDE:\n")
        f.write("• 90% score reduction for benign variants (LOW/MODIFIER + benign evidence)\n")
        f.write("• 2x score boost for pathogenic variants (HIGH impact or pathogenic evidence)\n")
        f.write("• Falls back to impact/consequence changes when clinical data unavailable\n\n")
        
        f.write("RATIONALE:\n")
        f.write("Clinical significance changes are rare but directly affect patient care decisions.\n")
        f.write("VEP consequence mismatches are common but often represent annotation differences\n")
        f.write("between genome builds rather than true functional changes. This redesign focuses\n")
        f.write("clinical review effort on variants most likely to change clinical interpretation.\n\n")
        
        f.write("DATA PROCESSING NOTES:\n")
        f.write("• All coordinates use VEP normalization (SNVs: original, Indels: original+1)\n")
        f.write("• Clinical significance normalized to 8 categories (PATHOGENIC, BENIGN, VUS, etc.)\n")
        f.write("• VEP analysis results cached for faster subsequent runs\n")
        f.write("• Priority scores calculated fresh for easy recalibration\n")

    print(f"✓ Summary statistics saved to: {summary_file}")


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
    • prioritized_variants.csv - Clinical evidence-focused ranked variant list (top variants only)
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
    • Clinical significance changes: +15-20 points (PATHOGENIC↔BENIGN, VUS→PATHOGENIC)
    • Same transcript consequence changes: +6 points (demoted from CRITICAL)
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
    parser.add_argument('--export-json', action='store_true',
                       help='Export structured results as JSON (optional)')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Enable verbose output')
    
    args = parser.parse_args()
    
    try:
        # Validate inputs
        db_path = Path(args.db_path)
        if not db_path.exists():
            print(f"Error: Database file not found: {db_path}")
            return 1
        
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        if args.verbose:
            print(f"Database: {db_path}")
            print(f"Output directory: {output_dir}")
            print(f"Max variants: {args.max_variants:,}")
            print(f"Min score threshold: {args.min_score}")
            print(f"Force recalculation: {args.force}")
        
        # Connect to database
        print("Connecting to database...")
        conn = sqlite3.connect(db_path)
        
        # Check database contents
        table_check = pd.read_sql_query("""
        SELECT name FROM sqlite_master WHERE type='table'
        """, conn)
        
        required_tables = ['comparison', 'hg19_vep', 'hg38_vep']
        missing_tables = [table for table in required_tables if table not in table_check['name'].values]
        
        if missing_tables:
            print(f"Error: Missing required tables: {missing_tables}")
            print("Available tables:", table_check['name'].tolist())
            conn.close()
            return 1
        
        # Set up caching
        cache_file = output_dir / 'variant_analysis_cache.pkl'
        
        # Initialize variant processor
        processor = VariantProcessor()
        
        print("\n" + "="*80)
        print("CROSSBUILD ASSESSOR - VARIANT PRIORITIZATION")
        print("="*80)
        
        # Process all variants (with caching and scoring)
        result_df = processor.process_all_variants(
            conn, 
            cache_file=cache_file, 
            force_recalculate=args.force
        )
        
        if len(result_df) == 0:
            print("\nNo discordant variants found for prioritization.")
            conn.close()
            return 0
        
        print(f"\n✓ Analysis completed: {len(result_df):,} total discordant variants found")
        
        # Filter by minimum score
        if args.min_score > 0:
            filtered_df = result_df[result_df['priority_score'] >= args.min_score].copy()
            print(f"✓ Filtered by min score ({args.min_score}): {len(filtered_df):,} variants remain")
        else:
            filtered_df = result_df.copy()
        
        # Sort by priority score (descending)
        filtered_df = filtered_df.sort_values('priority_score', ascending=False)
        
        # Create clinical evidence-focused CSV output
        output_df = create_clinical_csv_output(filtered_df, output_dir, args.max_variants)
        
        # Generate plots (unless disabled)
        if not args.no_plots and len(result_df) > 0:
            try:
                print("\nGenerating prioritization visualizations...")
                from visualization.plot_generator import PrioritizationPlotter
                from config.visualization_config import PLOT_COLORS, FIGURE_CONFIG
                
                plotter = PrioritizationPlotter(PLOT_COLORS, FIGURE_CONFIG)
                plotter.create_all_plots(result_df, conn, output_dir)
                
            except ImportError as e:
                print(f"Warning: Could not generate plots: {e}")
                print("Install matplotlib and seaborn for visualization support")
            except Exception as e:
                print(f"Warning: Plot generation failed: {e}")
        
        # Create summary statistics using BOTH full dataset and CSV subset
        create_summary_statistics(result_df, output_df, output_dir)
        
        # Optional JSON export
        if args.export_json:
            print("Exporting structured JSON data...")
            calculator = SummaryDataCalculator()
            priority_data = calculator.calculate_prioritization_summary(result_df, output_df)
            
            json_file = output_dir / 'prioritization_results.json'
            with open(json_file, 'w') as f:
                json.dump(priority_data, f, indent=2, default=str)
            print(f"✓ JSON data exported to: {json_file}")
        
        conn.close()
        
        # Final summary
        print("\n" + "="*80)
        print("PRIORITIZATION COMPLETED SUCCESSFULLY")
        print("="*80)
        print(f"📁 Output directory: {output_dir}")
        print(f"📊 Total variants analyzed: {len(result_df):,}")
        print(f"📋 Variants in CSV output: {len(output_df):,}")
        
        if len(output_df) > 0:
            # Show priority category breakdown
            priority_summary = output_df['Priority_Category'].value_counts()
            print(f"\n🎯 Priority distribution:")
            for category, count in priority_summary.items():
                percentage = count / len(output_df) * 100
                print(f"   {category}: {count:,} ({percentage:.1f}%)")
            
            # Show clinical changes
            if 'Has_Clinical_Change' in output_df.columns:
                clinical_changes = (output_df['Has_Clinical_Change'] == 'YES').sum()
                print(f"🔬 Clinical significance changes: {clinical_changes:,}")
        
        print(f"\n✅ Review prioritized_variants.csv for clinical decision support")
        
    except KeyboardInterrupt:
        print("\n⚠️  Analysis interrupted by user")
        return 1
    except Exception as e:
        print(f"\n❌ Error during analysis: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1


if __name__ == "__main__":
    main()