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
    
    # Enhanced clinical significance tracking
    output_df['Clinical_Significance_hg19'] = df['hg19_clin_sig'].fillna('')
    output_df['Clinical_Significance_hg38'] = df['hg38_clin_sig'].fillna('')
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
        
        # FULL DATASET STATISTICS
        f.write("FULL DATASET ANALYSIS (All Discordant Variants):\n")
        f.write("-" * 50 + "\n")
        
        # Priority category distribution for full dataset
        f.write("Priority Category Distribution (Full Dataset):\n")
        if len(df_full) > 0:
            category_counts = df_full['priority_category'].value_counts()
            for category, count in category_counts.items():
                pct = count / len(df_full) * 100
                f.write(f"  {category}: {count:,} ({pct:.1f}%)\n")
        f.write("\n")
        
        # Score distribution for full dataset
        f.write("Priority Score Distribution (Full Dataset):\n")
        if len(df_full) > 0:
            max_score = df_full['priority_score'].max()
            if max_score <= 50:
                intervals = [(0, 5), (5, 10), (10, 15), (15, 20), (20, 30), (30, 50)]
            elif max_score <= 100:
                intervals = [(0, 10), (10, 20), (20, 30), (30, 50), (50, 75), (75, 100)]
            elif max_score <= 200:
                intervals = [(0, 10), (10, 25), (25, 50), (50, 100), (100, 150), (150, 200)]
            else:
                intervals = [(0, 10), (10, 25), (25, 50), (50, 100), (100, 200), (200, 500), (500, int(max_score) + 1)]
            
            for start, end in intervals:
                count = ((df_full['priority_score'] >= start) & (df_full['priority_score'] < end)).sum()
                if count > 0:
                    pct = count / len(df_full) * 100
                    f.write(f"  Score {start}-{end-1}: {count:,} variants ({pct:.1f}%)\n")
        f.write("\n")
        
        # Clinical evidence breakdown for FULL dataset
        f.write("Clinical Evidence Breakdown (Full Dataset):\n")
        if len(df_full) > 0:
            clinical_evidence_issues = (
                (df_full['clin_sig_change'] != '').sum(),
                (df_full['sift_change'] != '').sum(),
                (df_full['polyphen_change'] != '').sum(),
                (df_full['hg19_is_pathogenic'] | df_full['hg38_is_pathogenic']).sum(),
                (df_full['hg19_is_benign'] | df_full['hg38_is_benign']).sum()
            )
            
            f.write(f"  Clinical significance changes: {clinical_evidence_issues[0]:,} variants\n")
            f.write(f"  SIFT prediction changes: {clinical_evidence_issues[1]:,} variants\n")
            f.write(f"  PolyPhen prediction changes: {clinical_evidence_issues[2]:,} variants\n")
            f.write(f"  Variants with pathogenic evidence: {clinical_evidence_issues[3]:,} variants\n")
            f.write(f"  Variants with benign evidence: {clinical_evidence_issues[4]:,} variants\n")
        f.write("\n")
        
        # Functional issue breakdown for FULL dataset
        f.write("Functional Issue Breakdown (Full Dataset):\n")
        if len(df_full) > 0:
            same_transcript_issues = (df_full['same_transcript_consequence_changes'] > 0).sum()
            gene_issues = (df_full['gene_changes'] > 0).sum()
            impact_issues = (df_full['impact_changes'] > 0).sum()
            different_transcript_issues = (df_full['same_consequence_different_transcripts'] > 0).sum()
            unmatched_issues = (df_full['unmatched_consequences'] > 0).sum()
            
            f.write(f"  Same transcript, different consequence: {same_transcript_issues:,} variants (CRITICAL)\n")
            f.write(f"  Gene annotation changes: {gene_issues:,} variants (reduced weight)\n")
            f.write(f"  Impact level changes: {impact_issues:,} variants (HIGH)\n")
            f.write(f"  Unmatched consequences: {unmatched_issues:,} variants (INVESTIGATE)\n")
        f.write("\n")
        f.write("\n")
        
        # Clinical review recommendations
        f.write("CLINICAL REVIEW RECOMMENDATIONS - CLINICAL EVIDENCE-DRIVEN STRATEGY\n")
        f.write("-" * 60 + "\n")
        
        f.write("PRIORITY 1 - CRITICAL:\n")
        f.write("1. Same transcript, different consequences (Score contribution: +10)\n")
        f.write("   - Same RefSeq/Ensembl ID, but different functional impact\n")
        f.write("   - May indicate annotation errors or genome build issues\n")
        f.write("   - Requires clinical review and potential wet lab validation\n\n")
        
        f.write("PRIORITY 2 - HIGH:\n")
        f.write("2. Clinical significance changes (Score contribution: +10 benign↔pathogenic, +8 pathogenic↔benign, +5 other)\n")
        f.write("   - Changes in ClinVar clinical significance between builds\n")
        f.write("   - Critical for clinical interpretation and reporting\n")
        f.write("3. High impact functional changes (Score contribution: VEP IMPACT-WEIGHTED)\n")
        f.write("   - HIGH impact variants get priority regardless of gene symbol differences\n\n")
        
        f.write("VEP IMPACT-BASED SCORING SYSTEM (with Clinical Override):\n")
        f.write("-" * 55 + "\n")
        f.write("CLINICAL EVIDENCE OVERRIDE:\n")
        f.write("• 90% score reduction for benign variants (LOW/MODIFIER + benign evidence)\n")
        f.write("• 2x score boost for pathogenic variants (HIGH impact or pathogenic evidence)\n")
        f.write("• Falls back to impact/consequence changes when clinical data unavailable\n\n")
        
        f.write("CORE SCORING (follows official VEP consequence severity hierarchy):\n\n")
        f.write("HIGH IMPACT (reduced gene weight: +8 points per gene change):\n")
        f.write("• transcript_ablation, splice_acceptor_variant, splice_donor_variant\n")
        f.write("• stop_gained, frameshift_variant, stop_lost, start_lost\n")
        f.write("• transcript_amplification, feature_elongation, feature_truncation\n\n")
        f.write("MODERATE IMPACT (reduced gene weight: +4 points per gene change):\n")
        f.write("• inframe_insertion, inframe_deletion, missense_variant\n")
        f.write("• protein_altering_variant\n\n")
        f.write("LOW IMPACT (reduced gene weight: +2 points per gene change):\n")
        f.write("• splice_donor_5th_base_variant, splice_region_variant\n")
        f.write("• splice_donor_region_variant, splice_polypyrimidine_tract_variant\n")
        f.write("• incomplete_terminal_codon_variant, start_retained_variant\n")
        f.write("• stop_retained_variant, synonymous_variant\n\n")
        f.write("MODIFIER IMPACT (minimal gene weight: +0.5 points per gene change):\n")
        f.write("• coding_sequence_variant, mature_miRNA_variant\n")
        f.write("• 5_prime_UTR_variant, 3_prime_UTR_variant\n")
        f.write("• non_coding_transcript_exon_variant, intron_variant\n")
        f.write("• upstream_gene_variant, downstream_gene_variant\n")
        f.write("• regulatory_region_variant, intergenic_variant, etc.\n\n")
        
        f.write("PATHOGENICITY PREDICTION CHANGES (+5 points each):\n")
        f.write("• SIFT: deleterious ↔ tolerated transitions\n")
        f.write("• PolyPhen: probably_damaging ↔ benign transitions\n\n")
        
        f.write("KEY IMPROVEMENTS:\n")
        f.write("• Reduced gene symbol weight (many differences are synonyms, not functional changes)\n")
        f.write("• Clinical evidence drives prioritization over annotation model differences\n")
        f.write("• Benign variants suppressed regardless of annotation changes\n")
        f.write("• Pathogenic variants boosted for clinical attention\n\n")
        
        f.write("ENHANCED OUTPUT COLUMNS:\n")
        f.write("• GT_hg19, GT_hg38: Actual genotypes from source alleles\n")
        f.write("• Problematic_Transcripts_hg19/hg38: Specific transcript IDs with issues\n")
        f.write("• Clinical_Significance_Change: Tracks benign↔pathogenic transitions\n")
        f.write("• SIFT_Change, PolyPhen_Change: Tracks pathogenicity prediction changes\n")
        f.write("• Priority_Category: Single column replacing multiple boolean flags\n\n")
        
        f.write("CACHING INFORMATION:\n")
        f.write("- VEP analysis results cached for faster subsequent runs\n")
        f.write("- Priority scores calculated on-the-fly for easy recalibration\n")
        f.write("- Use --force flag to recalculate VEP analysis and update cache\n")
        f.write("- Cache location: [output_dir]/variant_analysis_cache.pkl\n")

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
