#!/usr/bin/env python3
"""
Genomic Variant Database Analyzer

DESCRIPTION:
    Performs comprehensive analysis of genomic variant liftover tool performance
    comparing CrossMap and bcftools between genome builds. Focuses on liftover
    quality control with professional visualizations and detailed reports.

EXPECTED DATABASE STRUCTURE:
    SQLite database created by db_loader.py with three tables:
    
    1. comparison: Liftover comparison data
       - source_chrom, source_pos, source_alleles
       - liftover_hg38_pos, bcftools_hg38_pos
       - pos_match, gt_match, flip, swap
       - mapping_status

ANALYSIS MODULES:

1. LIFTOVER ANALYSIS:
   - Position and genotype match analysis
   - Strand flip and reference/alt swap patterns
   - Mapping status performance comparison
   - Statistical breakdowns with percentage calculations

2. POSITION DIFFERENCES ANALYSIS:
   - Coordinate discrepancies between CrossMap and bcftools
   - Distribution of position differences by chromosome
   - Genotype match correlation with position accuracy
   - Chromosome-level proportion analysis

OUTPUT:
    - Liftover analysis visualizations (PNG)
    - Position differences analysis plots (PNG)
    - Comprehensive analysis summary report (TXT)
    - Detailed statistics on liftover tool performance

USAGE:
    python db_analyzer.py --db-path genomic_analysis.db --output-dir results/
    python db_analyzer.py -d genomic_analysis.db -o results/ --no-plots

COORDINATE SYSTEM:
    All coordinates in the database use VEP normalization:
    - SNVs: original coordinates
    - Indels: original coordinates + 1

BCFTOOLS SWAP VALUES:
    The 'swap' column contains bcftools liftover swap status:
    - "NA": No action needed; alleles already matched reference genome
    - "1": REF and ALT alleles were swapped to match reference genome
    - "-1": REF and ALT alleles could not be swapped (ambiguous/invalid)

NOTE: VEP consequence analysis has been moved to variant_prioritizer.py
      to eliminate redundancy and improve memory efficiency.
"""

import sqlite3
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
from pathlib import Path
import numpy as np
import sys

from utils.summary_utils import SummaryDataCalculator
import json

# Set visualization style
plt.style.use('default')
sns.set_palette("husl")

def connect_database(db_path):
    """Connect to SQLite database and verify structure"""
    if not Path(db_path).exists():
        raise FileNotFoundError(f"Database not found: {db_path}")
    
    conn = sqlite3.connect(db_path)
    
    # Verify required tables exist
    cursor = conn.cursor()
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
    tables = [row[0] for row in cursor.fetchall()]
    
    required_tables = ['comparison']
    missing_tables = [t for t in required_tables if t not in tables]
    if missing_tables:
        raise ValueError(f"Missing required database tables: {missing_tables}")
    
    return conn

def analyze_cross_variables(conn, output_dir):
    """Analyze relationships between mapping_status, pos_match, and gt_match with enhanced plots"""
    print("\n=== LIFTOVER TOOL PERFORMANCE ANALYSIS ===")
    
    # Get all comparison data
    detailed_query = """
    SELECT 
        mapping_status,
        pos_match,
        gt_match,
        flip,
        swap
    FROM comparison
    """
    
    detailed_df = pd.read_sql_query(detailed_query, conn)
    
    print(f"Total variants analyzed: {len(detailed_df):,}")
    
    # Define colorblind-friendly palette
    colors_categorical = ['#1f77b4', '#ff7f0e', '#2ca02c', '#9467bd']  # Blue, Orange, Green, Purple
    colors_concordance = ['#2ca02c', '#1f77b4']  # Green for concordant, Blue for discordant
    colors_flip_swap = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']  # Blue, Orange, Green, Red, Purple
    
    # Calculate match categories
    detailed_df['match_category'] = detailed_df.apply(lambda row: 
        'Both Match' if row['pos_match'] == 1 and row['gt_match'] == 1
        else 'Position Only' if row['pos_match'] == 1 and row['gt_match'] == 0
        else 'Genotype Only' if row['pos_match'] == 0 and row['gt_match'] == 1
        else 'Both Mismatch', axis=1
    )
    
    match_counts = detailed_df['match_category'].value_counts()
    category_order = ['Both Match', 'Position Only', 'Genotype Only', 'Both Mismatch']
    match_counts = match_counts.reindex([cat for cat in category_order if cat in match_counts.index])
    
    # Determine if pie chart or bar chart is better
    dominant_category_pct = match_counts.max() / len(detailed_df) * 100
    use_pie_chart = dominant_category_pct < 85  # Switch to bar if one category dominates
    
    # Create comprehensive visualization
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))
    fig.suptitle('Liftover Tool Performance Analysis', fontsize=16, fontweight='bold')
    
    # PLOT 1: Match categories - smart chart selection
    print("1. Creating match categories visualization...")
    
    if use_pie_chart:
        # Use pie chart for balanced data
        wedges, texts, autotexts = axes[0].pie(
            match_counts.values, 
            labels=match_counts.index,
            autopct=lambda pct: f'{pct:.1f}%\n({int(pct/100*len(detailed_df)):,})',
            colors=colors_categorical[:len(match_counts)],
            startangle=90,
            pctdistance=0.85  # Move percentages closer to center
        )
        
        # Improve text formatting with smart positioning
        for i, (autotext, text) in enumerate(zip(autotexts, texts)):
            autotext.set_color('white')
            autotext.set_fontweight('bold')
            autotext.set_fontsize(9)
            
            # Adjust label position for small segments
            if match_counts.iloc[i] / len(detailed_df) < 0.05:  # Small segments
                text.set_fontsize(8)
                # Move small labels further out
                angle = (wedges[i].theta1 + wedges[i].theta2) / 2
                x, y = np.cos(np.radians(angle)), np.sin(np.radians(angle))
                text.set_position((x*1.3, y*1.3))
        
        axes[0].set_title('Variant Match Categories\n(Position & Genotype)', fontsize=12, fontweight='bold')
    else:
        # Use bar chart for heavily skewed data
        bars = axes[0].bar(range(len(match_counts)), match_counts.values, 
                          color=colors_categorical[:len(match_counts)], alpha=0.8)
        
        # Add percentage labels on bars
        for i, (bar, count) in enumerate(zip(bars, match_counts.values)):
            pct = count / len(detailed_df) * 100
            axes[0].text(bar.get_x() + bar.get_width()/2, bar.get_height() + count*0.01,
                        f'{pct:.1f}%\n({count:,})', ha='center', va='bottom', 
                        fontweight='bold', fontsize=10)
        
        axes[0].set_xticks(range(len(match_counts)))
        axes[0].set_xticklabels(match_counts.index, rotation=45, ha='right')
        axes[0].set_ylabel('Variant Count', fontweight='bold')
        axes[0].set_title('Variant Match Categories\n(Position & Genotype)', fontsize=12, fontweight='bold')
        axes[0].grid(True, alpha=0.3, axis='y')
    
    # PLOT 2: Concordance by mapping status (normalized percentages)
    print("2. Creating concordance by mapping status...")
    
    # Define concordance
    detailed_df['concordance'] = detailed_df.apply(lambda row:
        'Concordant' if row['pos_match'] == 1 and row['gt_match'] == 1
        else 'Discordant', axis=1
    )
    
    # Create contingency table and normalize within each mapping status
    concordance_mapping = pd.crosstab(detailed_df['mapping_status'], detailed_df['concordance'])
    concordance_pct = concordance_mapping.div(concordance_mapping.sum(axis=1), axis=0) * 100
    
    # Create percentage-based stacked bar plot
    concordance_pct.plot(
        kind='bar',
        ax=axes[1],
        color=colors_concordance,
        width=0.7,
        stacked=True
    )
    
    axes[1].set_title('Concordance by Mapping Status\n(Normalized Within Status)', fontsize=12, fontweight='bold')
    axes[1].set_xlabel('Mapping Status', fontweight='bold')
    axes[1].set_ylabel('Percentage', fontweight='bold')
    axes[1].legend(title='Concordance', title_fontsize=10, fontsize=10)
    axes[1].tick_params(axis='x', rotation=45)
    axes[1].grid(True, alpha=0.3, axis='y')
    axes[1].set_ylim(0, 100)
    
    # Add percentage labels on bars
    for container in axes[1].containers:
        axes[1].bar_label(container, fmt='%.1f%%', label_type='center', fontweight='bold', color='white')
    
    # PLOT 3: Mismatch details - Flip/Swap analysis (normalized percentages)
    print("3. Creating mismatch details (flip/swap analysis)...")
    
    # Filter to only mismatched variants - avoid SettingWithCopyWarning
    mismatch_mask = (detailed_df['pos_match'] == 0) | (detailed_df['gt_match'] == 0)
    mismatch_df = detailed_df[mismatch_mask].copy()  # Explicit copy to avoid warning
    
    print(f"Mismatched variants: {len(mismatch_df):,}")
    
    if len(mismatch_df) > 0:
        # Create flip/swap categories with CORRECTED LOGIC for float values
        def categorize_flip_swap(row):
            flip_status = row['flip'] if pd.notna(row['flip']) else 'no_flip'
            
            # Handle SWAP values that might be stored as floats (1.0, -1.0) or strings
            swap_val = row['swap']
            if pd.notna(swap_val):
                try:
                    # Convert to float first, then check value
                    swap_float = float(swap_val)
                    if swap_float == 1.0:
                        swap_status = '1'
                    elif swap_float == -1.0:
                        swap_status = '-1'
                    else:
                        swap_status = 'NA'
                except (ValueError, TypeError):
                    # If conversion fails, treat as string
                    swap_str = str(swap_val).strip()
                    if swap_str in ['1', '1.0']:
                        swap_status = '1'
                    elif swap_str in ['-1', '-1.0']:
                        swap_status = '-1'
                    else:
                        swap_status = 'NA'
            else:
                swap_status = 'NA'
            
            # CORRECTED BCFtools swap values interpretation:
            # "NA" = No action needed; alleles already matched reference genome
            # "1" = REF and ALT alleles were swapped to match reference genome  
            # "-1" = REF and ALT alleles could not be swapped (ambiguous/invalid)
            
            if flip_status == 'flip' and swap_status == '1':
                return 'FLIP+SWAP'
            elif flip_status == 'flip':
                return 'FLIP'
            elif swap_status == '1':  # REF/ALT alleles were swapped
                return 'SWAP'
            elif swap_status == '-1':  # Swap attempted but failed (ambiguous)
                return 'SWAP_FAILED'
            else:
                return 'NONE'
        
        mismatch_df['flip_swap_category'] = mismatch_df.apply(categorize_flip_swap, axis=1)
        
        # Debug: Check for SWAP values and categorization results
        swap_debug = mismatch_df['swap'].value_counts()
        category_debug = mismatch_df['flip_swap_category'].value_counts()
        print(f"Debug - SWAP value counts: {swap_debug.to_dict()}")
        print(f"Debug - Category counts: {category_debug.to_dict()}")
        
        # Create contingency table and normalize within each mapping status
        flip_swap_mapping = pd.crosstab(mismatch_df['mapping_status'], mismatch_df['flip_swap_category'])
        flip_swap_pct = flip_swap_mapping.div(flip_swap_mapping.sum(axis=1), axis=0) * 100
        
        # Order categories logically
        flip_swap_order = ['NONE', 'FLIP', 'SWAP', 'SWAP_FAILED', 'FLIP+SWAP']
        flip_swap_pct = flip_swap_pct.reindex(
            columns=[cat for cat in flip_swap_order if cat in flip_swap_pct.columns],
            fill_value=0
        )
        
        # Create percentage-based stacked bar plot
        flip_swap_pct.plot(
            kind='bar',
            ax=axes[2],
            color=colors_flip_swap[:len(flip_swap_pct.columns)],
            width=0.7,
            stacked=True
        )
        
        axes[2].set_title('Mismatch Details by Mapping Status\n(Normalized Within Status)', fontsize=12, fontweight='bold')
        axes[2].set_xlabel('Mapping Status', fontweight='bold')
        axes[2].set_ylabel('Percentage', fontweight='bold')
        axes[2].legend(title='Flip/Swap Status', title_fontsize=10, fontsize=10)
        axes[2].tick_params(axis='x', rotation=45)
        axes[2].grid(True, alpha=0.3, axis='y')
        axes[2].set_ylim(0, 100)
        
        # Add percentage labels on bars
        for container in axes[2].containers:
            axes[2].bar_label(container, fmt='%.1f%%', label_type='center', fontweight='bold', color='white')
            
    else:
        axes[2].text(0.5, 0.5, 'No mismatched variants found', 
                    ha='center', va='center', fontsize=12, transform=axes[2].transAxes)
        axes[2].set_title('Mismatch Details by Mapping Status\n(Flip/Swap Analysis)', fontsize=12, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_dir / 'liftover_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Print summary statistics
    print("\n=== SUMMARY STATISTICS ===")
    print(f"Total variants: {len(detailed_df):,}")
    print("\nMatch category breakdown:")
    for category, count in match_counts.items():
        pct = count / len(detailed_df) * 100
        print(f"  {category}: {count:,} ({pct:.1f}%)")
    
    print(f"\nConcordance breakdown:")
    concordance_summary = detailed_df['concordance'].value_counts()
    for status, count in concordance_summary.items():
        pct = count / len(detailed_df) * 100
        print(f"  {status}: {count:,} ({pct:.1f}%)")
    
    if len(mismatch_df) > 0:
        print(f"\nFlip/Swap analysis (mismatched variants only):")
        print(f"BCFtools SWAP value meanings:")
        print(f"  • NA: No action needed; alleles already matched reference")
        print(f"  • 1: REF and ALT alleles were swapped to match reference")
        print(f"  • -1: REF and ALT alleles could not be swapped (ambiguous/invalid)")
        print(f"")
        flip_swap_summary = mismatch_df['flip_swap_category'].value_counts()
        for category, count in flip_swap_summary.items():
            pct = count / len(mismatch_df) * 100
            print(f"  {category}: {count:,} ({pct:.1f}%)")
    
    print("✓ Liftover analysis completed")

def analyze_position_differences(conn, output_dir):
    """Analyze position differences between liftover tools with enhanced visualizations"""
    print("\n=== POSITION DIFFERENCES ANALYSIS ===")
    
    query = """
    SELECT 
        mapping_status,
        source_chrom,
        liftover_hg38_pos,
        bcftools_hg38_pos,
        gt_match,
        ABS(COALESCE(liftover_hg38_pos, 0) - COALESCE(bcftools_hg38_pos, 0)) as pos_diff
    FROM comparison 
    WHERE pos_match = 0 AND liftover_hg38_pos IS NOT NULL AND bcftools_hg38_pos IS NOT NULL
    """
    
    diff_df = pd.read_sql_query(query, conn)
    
    if len(diff_df) == 0:
        print("No position differences found between liftover tools")
        return
    
    print(f"Position differences found in {len(diff_df):,} variants")
    print(f"\nDifference statistics:")
    print(diff_df['pos_diff'].describe())
    
    print(f"\nDifferences by mapping status:")
    mapping_stats = diff_df.groupby('mapping_status')['pos_diff'].describe()
    print(mapping_stats)
    
    print(f"\nDifferences by chromosome (top 10):")
    chr_stats = diff_df.groupby('source_chrom')['pos_diff'].agg(['count', 'mean', 'max']).round(1)
    print(chr_stats.head(10))
    
    # Create visualization with enhanced plots
    fig, axes = plt.subplots(2, 2, figsize=(18, 12))
    fig.suptitle('Position Differences Between Liftover Tools', fontsize=16, fontweight='bold')
    
    # Define consistent colors
    colors_concordance = ['#2ca02c', '#1f77b4']  # Green for match, Blue for mismatch (colorblind friendly)
    
    # Plot 1: Distribution of differences (log scale)
    axes[0,0].hist(diff_df['pos_diff'], bins=50, alpha=0.7, edgecolor='black', color='#ff7f0e')
    axes[0,0].set_xlabel('Position Difference (bp)', fontweight='bold')
    axes[0,0].set_ylabel('Count', fontweight='bold')
    axes[0,0].set_title('Distribution of Position Differences', fontsize=12, fontweight='bold')
    axes[0,0].set_yscale('log')
    axes[0,0].grid(True, alpha=0.3)
    
    # Plot 2: Box plot by mapping status, colored by genotype match
    if diff_df['mapping_status'].nunique() > 1:
        # Create grouped data for coloring by gt_match
        for i, (status, group) in enumerate(diff_df.groupby('mapping_status')):
            gt_match_data = [group[group['gt_match'] == 1]['pos_diff'].values,
                           group[group['gt_match'] == 0]['pos_diff'].values]
            
            positions = [i*3 + 1, i*3 + 2]  # Group positions
            bp = axes[0,1].boxplot(gt_match_data, positions=positions, 
                                 patch_artist=True, widths=0.6)
            
            # Color by genotype match status
            for patch, color in zip(bp['boxes'], colors_concordance):
                patch.set_facecolor(color)
                patch.set_alpha(0.7)
        
        axes[0,1].set_title('Position Differences by Mapping Status & Genotype Match', fontsize=12, fontweight='bold')
        axes[0,1].set_xlabel('Mapping Status', fontweight='bold')
        axes[0,1].set_ylabel('Position Difference (bp)', fontweight='bold')
        axes[0,1].set_yscale('log')
        
        # Custom x-axis labels
        mapping_statuses = diff_df['mapping_status'].unique()
        axes[0,1].set_xticks([i*3 + 1.5 for i in range(len(mapping_statuses))])
        axes[0,1].set_xticklabels(mapping_statuses)
        
        # Add legend
        legend_elements = [plt.Rectangle((0,0),1,1, facecolor=colors_concordance[0], alpha=0.7, label='GT Match'),
                          plt.Rectangle((0,0),1,1, facecolor=colors_concordance[1], alpha=0.7, label='GT Mismatch')]
        axes[0,1].legend(handles=legend_elements)
    else:
        axes[0,1].text(0.5, 0.5, 'Single mapping\nstatus only', ha='center', va='center', fontsize=12)
        axes[0,1].set_title('Position Differences by Mapping Status & Genotype Match', fontsize=12, fontweight='bold')
    
    # Plot 3: Chromosome proportions (normalized within each chromosome)
    chr_counts = diff_df['source_chrom'].value_counts().head(15)
    
    # Get genotype match proportions for each chromosome (normalized within chromosome)
    chr_gt_data = []
    chr_labels = []
    for chrom in chr_counts.index:
        chrom_data = diff_df[diff_df['source_chrom'] == chrom]
        total_chrom = len(chrom_data)
        gt_match_pct = (chrom_data['gt_match'] == 1).sum() / total_chrom * 100
        gt_mismatch_pct = (chrom_data['gt_match'] == 0).sum() / total_chrom * 100
        chr_gt_data.append([gt_match_pct, gt_mismatch_pct])
        chr_labels.append(f"{chrom}\n(n={total_chrom})")
    
    # Create percentage-based stacked bar chart
    chr_gt_array = np.array(chr_gt_data)
    x_pos = np.arange(len(chr_labels))
    
    axes[1,0].bar(x_pos, chr_gt_array[:, 0], color=colors_concordance[0], label='GT Match', alpha=0.8)
    axes[1,0].bar(x_pos, chr_gt_array[:, 1], bottom=chr_gt_array[:, 0], 
                 color=colors_concordance[1], label='GT Mismatch', alpha=0.8)
    
    axes[1,0].set_title('Position Mismatches by Chromosome\n(GT Match Proportions)', fontsize=12, fontweight='bold')
    axes[1,0].set_xlabel('Chromosome', fontweight='bold')
    axes[1,0].set_ylabel('Percentage', fontweight='bold')
    axes[1,0].set_xticks(x_pos)
    axes[1,0].set_xticklabels(chr_labels, rotation=45, fontsize=9)
    axes[1,0].legend()
    axes[1,0].grid(True, alpha=0.3, axis='y')
    axes[1,0].set_ylim(0, 100)
    
    # Plot 4: GWAS-style genome-wide position plot
    sample_size = min(5000, len(diff_df))
    sample_df = diff_df.sample(sample_size)
    
    # Create chromosome-aware x-axis positions
    chromosome_order = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', 
                       '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']
    
    # Calculate cumulative positions for each chromosome
    chr_positions = {}
    cumulative_pos = 0
    chr_centers = {}
    chr_colors = ['#1f77b4', '#ff7f0e']  # Alternating blue and orange
    
    for i, chrom in enumerate(chromosome_order):
        if chrom in sample_df['source_chrom'].values:
            chr_data = sample_df[sample_df['source_chrom'] == chrom]
            if len(chr_data) > 0:
                max_pos = chr_data['liftover_hg38_pos'].max()
                chr_positions[chrom] = cumulative_pos
                chr_centers[chrom] = cumulative_pos + max_pos / 2
                cumulative_pos += max_pos + 1000000  # Add 1Mb gap between chromosomes
    
    # Plot points with alternating chromosome colors
    for i, chrom in enumerate(chromosome_order):
        if chrom in chr_positions:
            chr_data = sample_df[sample_df['source_chrom'] == chrom]
            if len(chr_data) > 0:
                x_positions = chr_data['liftover_hg38_pos'] + chr_positions[chrom]
                
                # Color by genotype match status
                gt_match_mask = chr_data['gt_match'] == 1
                color = chr_colors[i % 2]
                
                axes[1,1].scatter(x_positions[gt_match_mask], chr_data[gt_match_mask]['pos_diff'],
                                alpha=0.6, s=15, c=color, label='GT Match' if i == 0 else "", marker='o')
                axes[1,1].scatter(x_positions[~gt_match_mask], chr_data[~gt_match_mask]['pos_diff'],
                                alpha=0.6, s=15, c=color, label='GT Mismatch' if i == 0 else "", marker='^')
    
    # Customize x-axis with chromosome labels
    chr_ticks = [chr_centers[chrom] for chrom in chromosome_order if chrom in chr_centers]
    chr_tick_labels = [chrom for chrom in chromosome_order if chrom in chr_centers]
    
    axes[1,1].set_xticks(chr_ticks)
    axes[1,1].set_xticklabels(chr_tick_labels, fontsize=9)
    axes[1,1].set_xlabel('Chromosome', fontweight='bold')
    axes[1,1].set_ylabel('Position Difference (bp)', fontweight='bold')
    axes[1,1].set_title(f'Genome-wide Position Differences\n(Sample: {sample_size:,}, GT Match: ○ vs △)', 
                       fontsize=12, fontweight='bold')
    axes[1,1].grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(output_dir / 'position_differences_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print("✓ Position differences analysis completed")

def generate_comprehensive_report(conn, output_dir):
    """Generate comprehensive analysis summary report"""
    print("\n=== GENERATING COMPREHENSIVE SUMMARY REPORT ===")
    
    # Basic database statistics
    total_variants = pd.read_sql_query("SELECT COUNT(*) as count FROM comparison", conn).iloc[0]['count']
    
    pos_match_rate = pd.read_sql_query("""
        SELECT AVG(pos_match) as rate FROM comparison
    """, conn).iloc[0]['rate']
    
    gt_match_rate = pd.read_sql_query("""
        SELECT AVG(gt_match) as rate FROM comparison
    """, conn).iloc[0]['rate']
    
    print(f"Total variants analyzed: {total_variants:,}")
    print(f"Position match rate: {pos_match_rate:.1%}")
    print(f"Genotype match rate: {gt_match_rate:.1%}")
    
    # Detailed mapping status breakdown
    mapping_stats = pd.read_sql_query("""
        SELECT mapping_status, COUNT(*) as count,
               AVG(pos_match) as pos_match_rate,
               AVG(gt_match) as gt_match_rate
        FROM comparison 
        GROUP BY mapping_status
        ORDER BY count DESC
    """, conn)
    
    print(f"\nBreakdown by mapping status:")
    print(mapping_stats.round(3))
    
    # Get detailed statistics for summary file
    detailed_query = """
    SELECT 
        mapping_status,
        pos_match,
        gt_match,
        flip,
        swap
    FROM comparison
    """
    
    detailed_df = pd.read_sql_query(detailed_query, conn)
    
    # Calculate match categories
    detailed_df['match_category'] = detailed_df.apply(lambda row: 
        'Both Match' if row['pos_match'] == 1 and row['gt_match'] == 1
        else 'Position Only' if row['pos_match'] == 1 and row['gt_match'] == 0
        else 'Genotype Only' if row['pos_match'] == 0 and row['gt_match'] == 1
        else 'Both Mismatch', axis=1
    )
    
    match_counts = detailed_df['match_category'].value_counts()
    category_order = ['Both Match', 'Position Only', 'Genotype Only', 'Both Mismatch']
    match_counts = match_counts.reindex([cat for cat in category_order if cat in match_counts.index])
    
    # Calculate concordance
    detailed_df['concordance'] = detailed_df.apply(lambda row:
        'Concordant' if row['pos_match'] == 1 and row['gt_match'] == 1
        else 'Discordant', axis=1
    )
    concordance_summary = detailed_df['concordance'].value_counts()
    
    # Flip/Swap analysis for mismatched variants
    mismatch_mask = (detailed_df['pos_match'] == 0) | (detailed_df['gt_match'] == 0)
    mismatch_df = detailed_df[mismatch_mask].copy()
    
    flip_swap_summary = None
    if len(mismatch_df) > 0:
        def categorize_flip_swap(row):
            flip_status = row['flip'] if pd.notna(row['flip']) else 'no_flip'
            
            # Handle SWAP values that might be stored as floats (1.0, -1.0) or strings
            swap_val = row['swap']
            if pd.notna(swap_val):
                try:
                    # Convert to float first, then check value
                    swap_float = float(swap_val)
                    if swap_float == 1.0:
                        swap_status = '1'
                    elif swap_float == -1.0:
                        swap_status = '-1'
                    else:
                        swap_status = 'NA'
                except (ValueError, TypeError):
                    # If conversion fails, treat as string
                    swap_str = str(swap_val).strip()
                    if swap_str in ['1', '1.0']:
                        swap_status = '1'
                    elif swap_str in ['-1', '-1.0']:
                        swap_status = '-1'
                    else:
                        swap_status = 'NA'
            else:
                swap_status = 'NA'
            
            # CORRECTED BCFtools swap values interpretation:
            # "NA" = No action needed; alleles already matched reference genome
            # "1" = REF and ALT alleles were swapped to match reference genome  
            # "-1" = REF and ALT alleles could not be swapped (ambiguous/invalid)
            
            if flip_status == 'flip' and swap_status == '1':
                return 'FLIP+SWAP'
            elif flip_status == 'flip':
                return 'FLIP'
            elif swap_status == '1':  # REF/ALT alleles were swapped
                return 'SWAP'
            elif swap_status == '-1':  # Swap attempted but failed (ambiguous)
                return 'SWAP_FAILED'
            else:
                return 'NONE'
        
        mismatch_df['flip_swap_category'] = mismatch_df.apply(categorize_flip_swap, axis=1)
        flip_swap_summary = mismatch_df['flip_swap_category'].value_counts()
    
    # Save comprehensive summary to file
    summary_file = output_dir / 'liftover_analysis_summary.txt'
    with open(summary_file, 'w') as f:
        f.write("GENOMIC VARIANT LIFTOVER TOOL ANALYSIS SUMMARY\n")
        f.write("=" * 60 + "\n\n")
        
        # Database overview
        f.write("DATABASE OVERVIEW\n")
        f.write("-" * 20 + "\n")
        f.write(f"Total variants analyzed: {total_variants:,}\n")
        f.write(f"Position match rate: {pos_match_rate:.1%}\n")
        f.write(f"Genotype match rate: {gt_match_rate:.1%}\n")
        f.write(f"Analysis performed on VEP-normalized coordinates\n\n")
        
        # Mapping status breakdown
        f.write("LIFTOVER TOOL PERFORMANCE\n")
        f.write("-" * 30 + "\n")
        f.write("Breakdown by mapping status:\n")
        f.write(mapping_stats.to_string(index=False))
        f.write("\n\n")
        
        # Match category breakdown
        f.write("VARIANT MATCH CATEGORIES\n")
        f.write("-" * 25 + "\n")
        for category, count in match_counts.items():
            pct = count / len(detailed_df) * 100
            f.write(f"{category}: {count:,} ({pct:.1f}%)\n")
        f.write("\n")
        
        # Concordance breakdown
        f.write("CONCORDANCE ANALYSIS\n")
        f.write("-" * 20 + "\n")
        for status, count in concordance_summary.items():
            pct = count / len(detailed_df) * 100
            f.write(f"{status}: {count:,} ({pct:.1f}%)\n")
        f.write("\n")
        
        # Flip/Swap analysis
        if flip_swap_summary is not None and len(flip_swap_summary) > 0:
            f.write("FLIP/SWAP ANALYSIS (Mismatched Variants Only)\n")
            f.write("-" * 45 + "\n")
            f.write(f"Total mismatched variants: {len(mismatch_df):,}\n\n")
            f.write("BCFtools SWAP value meanings:\n")
            f.write("• NA: No action needed; alleles already matched reference genome\n")
            f.write("• 1: REF and ALT alleles were swapped to match reference genome\n")
            f.write("• -1: REF and ALT alleles could not be swapped (ambiguous/invalid)\n\n")
            f.write("Flip/Swap category breakdown:\n")
            f.write("• FLIP: Strand flip occurred\n")
            f.write("• SWAP: REF/ALT alleles were swapped during liftover\n")
            f.write("• SWAP_FAILED: Swap attempted but failed (ambiguous alleles)\n")
            f.write("• FLIP+SWAP: Both strand flip and allele swap occurred\n")
            f.write("• NONE: No flip or swap operations needed\n\n")
            for category, count in flip_swap_summary.items():
                pct = count / len(mismatch_df) * 100
                f.write(f"{category}: {count:,} ({pct:.1f}%)\n")
            f.write("\n")
        
        # Generated files
        f.write("GENERATED OUTPUT FILES\n")
        f.write("-" * 25 + "\n")
        f.write("• liftover_analysis.png - Liftover tool performance visualization\n")
        f.write("• position_differences_analysis.png - Coordinate discrepancy analysis\n")
        f.write("• liftover_analysis_summary.txt - This report\n\n")
        
        f.write("ANALYSIS FOCUS\n")
        f.write("-" * 15 + "\n")
        f.write("This analysis focuses on liftover tool quality control comparing CrossMap\n")
        f.write("and bcftools performance. For detailed VEP consequence analysis and variant\n")
        f.write("prioritization, use variant_prioritizer.py which provides comprehensive\n")
        f.write("clinical review outputs with priority scoring.\n\n")
        
        f.write("COORDINATE SYSTEM NOTES\n")
        f.write("-" * 25 + "\n")
        f.write("• All coordinates use VEP normalization\n")
        f.write("• SNVs: original input coordinates\n")
        f.write("• Indels: original input coordinates + 1\n")
        f.write("• Position differences calculated on normalized coordinates\n\n")
        
        f.write("BCFTOOLS SWAP VALUE INTERPRETATION\n")
        f.write("-" * 35 + "\n")
        f.write("The 'swap' column in the comparison data contains bcftools liftover status:\n")
        f.write("• NA: No action needed; alleles already matched reference genome\n")
        f.write("• 1: REF and ALT alleles were swapped to match reference genome\n")
        f.write("• -1: REF and ALT alleles could not be swapped (ambiguous/invalid)\n\n")
        f.write("This information is crucial for understanding allele orientation changes\n")
        f.write("during genome build liftover and potential impact on variant interpretation.\n")
    
    print(f"✓ Comprehensive summary report saved to: {summary_file}")
    return summary_file

def main():
    parser = argparse.ArgumentParser(
        description='Analyze genomic variant liftover tool performance',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
USAGE EXAMPLES:
    python db_analyzer.py --db-path genomic_analysis.db --output-dir results/
    python db_analyzer.py -d genomic_analysis.db -o results/ --no-plots
    python db_analyzer.py --db-path /path/to/data.db --output-dir /path/to/output/

BCFTOOLS SWAP VALUES:
    • NA: No action needed; alleles already matched reference genome
    • 1: REF and ALT alleles were swapped to match reference genome  
    • -1: REF and ALT alleles could not be swapped (ambiguous/invalid)

NOTE: VEP consequence analysis has been moved to variant_prioritizer.py
      for better memory efficiency and to eliminate redundancy.
        """
    )
    
    parser.add_argument('--db-path', '-d', required=True,
                       help='SQLite database file path')
    parser.add_argument('--output-dir', '-o', default='analysis_output',
                       help='Output directory for results (default: analysis_output)')
    parser.add_argument('--no-plots', action='store_true',
                       help='Skip plot generation (faster for large datasets)')
    parser.add_argument('--export-json', action='store_true',
                       help='Export structured results as JSON (optional)') 
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Enable verbose output')
    
    args = parser.parse_args()
    
    db_path = args.db_path
    output_dir = Path(args.output_dir)
    print(f"✓ Using database: {db_path}")
    print(f"✓ Output directory: {output_dir}")
    
    # Create output directory
    output_dir.mkdir(exist_ok=True, parents=True)
    
    create_plots = not args.no_plots
    
    try:
        print(f"\n=== STARTING LIFTOVER ANALYSIS ===")
        print("Note: VEP consequence analysis moved to variant_prioritizer.py for better memory efficiency")
        
        # Connect to database
        conn = connect_database(db_path)
        print(f"✓ Connected to database: {db_path}")
        
        # Run liftover analysis modules
        if create_plots:
            print("Running liftover analysis with visualizations...")
            analyze_cross_variables(conn, output_dir)
            analyze_position_differences(conn, output_dir)
        else:
            print("Running analysis without plots...")
        
        # Generate comprehensive summary
        generate_comprehensive_report(conn, output_dir)

        # Optional JSON export
        if args.export_json:
            print("Exporting structured JSON data...")
            calculator = SummaryDataCalculator()
            liftover_data = calculator.calculate_liftover_summary(conn)
            
            json_file = output_dir / 'liftover_analysis.json'
            with open(json_file, 'w') as f:
                json.dump(liftover_data, f, indent=2, default=str)
            print(f"✓ JSON data exported to: {json_file}")
        
        conn.close()
        
        print(f"\n=== LIFTOVER ANALYSIS COMPLETED ===")
        print(f"✓ Results saved in: {output_dir}")
        print(f"\nGenerated files:")
        if create_plots:
            print(f"  • liftover_analysis.png")
            print(f"  • position_differences_analysis.png") 
        print(f"  • liftover_analysis_summary.txt")
        
        print(f"\nNext steps:")
        print(f"  → Run variant_prioritizer.py for VEP consequence analysis and clinical prioritization")
        print(f"  → Use liftover_analysis.png for quality control assessment")
        print(f"  → Review position_differences_analysis.png for coordinate discrepancies")
        
    except Exception as e:
        print(f"Error during analysis: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
