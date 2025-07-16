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

# Set visualization style
plt.style.use('default')
sns.set_palette("husl")

# VEP consequence severity mapping based on official Ensembl hierarchy
VEP_CONSEQUENCE_IMPACT = {
    # HIGH IMPACT
    'transcript_ablation': 'HIGH',
    'splice_acceptor_variant': 'HIGH',
    'splice_donor_variant': 'HIGH',
    'stop_gained': 'HIGH',
    'frameshift_variant': 'HIGH',
    'stop_lost': 'HIGH',
    'start_lost': 'HIGH',
    'transcript_amplification': 'HIGH',
    'feature_elongation': 'HIGH',
    'feature_truncation': 'HIGH',
    
    # MODERATE IMPACT
    'inframe_insertion': 'MODERATE',
    'inframe_deletion': 'MODERATE',
    'missense_variant': 'MODERATE',
    'protein_altering_variant': 'MODERATE',
    
    # LOW IMPACT
    'splice_donor_5th_base_variant': 'LOW',
    'splice_region_variant': 'LOW',
    'splice_donor_region_variant': 'LOW',
    'splice_polypyrimidine_tract_variant': 'LOW',
    'incomplete_terminal_codon_variant': 'LOW',
    'start_retained_variant': 'LOW',
    'stop_retained_variant': 'LOW',
    'synonymous_variant': 'LOW',
    
    # MODIFIER IMPACT
    'coding_sequence_variant': 'MODIFIER',
    'mature_miRNA_variant': 'MODIFIER',
    '5_prime_UTR_variant': 'MODIFIER',
    '3_prime_UTR_variant': 'MODIFIER',
    'non_coding_transcript_exon_variant': 'MODIFIER',
    'intron_variant': 'MODIFIER',
    'NMD_transcript_variant': 'MODIFIER',
    'non_coding_transcript_variant': 'MODIFIER',
    'coding_transcript_variant': 'MODIFIER',
    'upstream_gene_variant': 'MODIFIER',
    'downstream_gene_variant': 'MODIFIER',
    'TFBS_ablation': 'MODIFIER',
    'TFBS_amplification': 'MODIFIER',
    'TF_binding_site_variant': 'MODIFIER',
    'regulatory_region_ablation': 'MODIFIER',
    'regulatory_region_amplification': 'MODIFIER',
    'regulatory_region_variant': 'MODIFIER',
    'intergenic_variant': 'MODIFIER',
    'sequence_variant': 'MODIFIER'
}

def get_impact_numeric_value(impact):
    """Convert impact to numeric value for magnitude calculation"""
    impact_values = {'HIGH': 4, 'MODERATE': 3, 'LOW': 2, 'MODIFIER': 1}
    return impact_values.get(impact, 0)

def calculate_impact_transition_magnitude(hg19_impact, hg38_impact):
    """Calculate magnitude of impact transition for clinical significance"""
    hg19_val = get_impact_numeric_value(hg19_impact)
    hg38_val = get_impact_numeric_value(hg38_impact)
    magnitude = abs(hg19_val - hg38_val)
    
    # Check if transition involves MODERATE or HIGH (clinically significant)
    involves_moderate_or_high = (
        hg19_impact in ['MODERATE', 'HIGH'] or 
        hg38_impact in ['MODERATE', 'HIGH']
    )
    
    return magnitude, involves_moderate_or_high

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

def extract_genotype_from_alleles(source_alleles):
    """Extract reference and alternative alleles from source_alleles string"""
    if pd.isna(source_alleles) or source_alleles == '':
        return '', ''
    
    # Handle various formats: "A/G", "A,G", "-/G", etc.
    if '/' in source_alleles:
        parts = source_alleles.split('/')
    elif ',' in source_alleles:
        parts = source_alleles.split(',')
    else:
        return source_alleles, ''
    
    if len(parts) >= 2:
        return parts[0].strip(), parts[1].strip()
    else:
        return parts[0].strip(), ''

def is_pathogenic_clinical_significance(clin_sig):
    """Check if clinical significance indicates pathogenic variant"""
    if pd.isna(clin_sig) or clin_sig in ['', '-', 'nan']:
        return False
    
    clin_sig_lower = str(clin_sig).lower()
    pathogenic_terms = ['pathogenic', 'likely_pathogenic', 'drug_response']
    return any(term in clin_sig_lower for term in pathogenic_terms)

def is_benign_clinical_significance(clin_sig):
    """Check if clinical significance indicates benign variant"""
    if pd.isna(clin_sig) or clin_sig in ['', '-', 'nan']:
        return False
    
    clin_sig_lower = str(clin_sig).lower()
    benign_terms = ['benign', 'likely_benign']
    return any(term in clin_sig_lower for term in benign_terms)

def parse_sift_prediction(sift_str):
    """Parse SIFT prediction and score"""
    if pd.isna(sift_str) or sift_str in ['', '-', 'nan']:
        return None, None
    
    sift_str = str(sift_str).lower()
    
    # Extract prediction
    if 'deleterious' in sift_str:
        prediction = 'deleterious'
    elif 'tolerated' in sift_str:
        prediction = 'tolerated'
    else:
        prediction = None
    
    # Extract score (numbers in parentheses)
    score_match = re.search(r'\(([0-9.]+)\)', sift_str)
    score = float(score_match.group(1)) if score_match else None
    
    return prediction, score

def parse_polyphen_prediction(polyphen_str):
    """Parse PolyPhen prediction and score"""
    if pd.isna(polyphen_str) or polyphen_str in ['', '-', 'nan']:
        return None, None
    
    polyphen_str = str(polyphen_str).lower()
    
    # Extract prediction
    if 'probably_damaging' in polyphen_str:
        prediction = 'probably_damaging'
    elif 'possibly_damaging' in polyphen_str:
        prediction = 'possibly_damaging'
    elif 'benign' in polyphen_str:
        prediction = 'benign'
    else:
        prediction = None
    
    # Extract score (numbers in parentheses)
    score_match = re.search(r'\(([0-9.]+)\)', polyphen_str)
    score = float(score_match.group(1)) if score_match else None
    
    return prediction, score

def calculate_variant_scores(conn, cache_file=None, force_recalculate=False):
    """Calculate priority scores using clinical evidence-driven analysis with memory optimization and caching"""
    print("Calculating variant priority scores with clinical evidence-driven analysis...")
    
    # Check for cached results (raw VEP analysis only, no scores)
    if cache_file and Path(cache_file).exists() and not force_recalculate:
        print(f"Loading cached VEP analysis from: {cache_file}")
        try:
            cached_vep_analysis = pd.read_pickle(cache_file)
            print(f"Loaded {len(cached_vep_analysis):,} cached VEP analyses")
            
            # Calculate scores on-the-fly (not cached)
            result_df = calculate_scores_from_vep_analysis(cached_vep_analysis)
            return result_df
        except Exception as e:
            print(f"Warning: Could not load cache file ({e}), recalculating...")
    
    print("Calculating fresh VEP analysis...")
    
    def normalize_transcript_id(transcript_id):
        """Strip version numbers for transcript matching (ENST123.4 → ENST123)"""
        if pd.isna(transcript_id) or transcript_id == '' or transcript_id == '-':
            return None
        return str(transcript_id).split('.')[0]
    
    # STEP 1: Get all unique variants (both concordant and discordant) - MEMORY EFFICIENT
    print("Step 1: Identifying all variants for VEP analysis...")
    variant_query = """
    SELECT DISTINCT 
        c.source_chrom,
        c.source_pos,
        c.bcftools_hg38_chrom,
        c.bcftools_hg38_pos,
        c.mapping_status,
        c.pos_match,
        c.gt_match,
        c.flip,
        c.swap,
        c.liftover_hg38_pos,
        c.source_alleles,
        ABS(COALESCE(c.liftover_hg38_pos, 0) - COALESCE(c.bcftools_hg38_pos, 0)) as pos_difference
    FROM comparison c
    WHERE c.bcftools_hg38_chrom IS NOT NULL 
      AND c.bcftools_hg38_pos IS NOT NULL
    """
    
    variants_df = pd.read_sql_query(variant_query, conn)
    print(f"Found {len(variants_df):,} unique variants for VEP analysis")
    
    # STEP 2: Process variants in chunks - MEMORY SAFE
    chunk_size = 10000  # Process 10K variants at a time
    all_vep_analyses = []
    total_chunks = (len(variants_df) + chunk_size - 1) // chunk_size
    
    print(f"Processing {len(variants_df):,} variants in {total_chunks} chunks of {chunk_size:,}...")
    
    for chunk_idx in range(total_chunks):
        start_idx = chunk_idx * chunk_size
        end_idx = min(start_idx + chunk_size, len(variants_df))
        chunk_variants = variants_df.iloc[start_idx:end_idx]
        
        print(f"  Processing chunk {chunk_idx + 1}/{total_chunks} ({len(chunk_variants):,} variants)...")
        
        # Process this chunk of variants
        chunk_analyses = process_variant_chunk(conn, chunk_variants, normalize_transcript_id)
        all_vep_analyses.extend(chunk_analyses)
        
        # Memory cleanup
        del chunk_variants
        
    print(f"Completed VEP analysis of {len(all_vep_analyses):,} variants")
    
    if len(all_vep_analyses) == 0:
        print("No variants found for analysis")
        return pd.DataFrame()
    
    vep_analysis_df = pd.DataFrame(all_vep_analyses)
    
    # Save VEP analysis to cache (without scores)
    if cache_file:
        try:
            vep_analysis_df.to_pickle(cache_file)
            print(f"✓ Cached VEP analysis saved to: {cache_file}")
        except Exception as e:
            print(f"Warning: Could not save cache file ({e})")
    
    # Calculate scores on-the-fly
    result_df = calculate_scores_from_vep_analysis(vep_analysis_df)
    return result_df

def process_variant_chunk(conn, chunk_variants, normalize_transcript_id):
    """Process a chunk of variants with comprehensive VEP analysis"""
    
    vep_analyses = []
    
    for _, variant_row in chunk_variants.iterrows():
        chrom = variant_row['source_chrom']
        pos = variant_row['source_pos']
        hg38_chrom = variant_row['bcftools_hg38_chrom']
        hg38_pos = variant_row['bcftools_hg38_pos']
        
        # Get VEP annotations for this specific variant - TARGETED QUERIES
        hg19_query = """
        SELECT feature_type, consequence, impact, symbol, feature, sift, polyphen, gnomadg_af, clin_sig
        FROM hg19_vep 
        WHERE chr = ? AND pos = ?
        """
        
        hg38_query = """
        SELECT feature_type, consequence, impact, symbol, feature, sift, polyphen, gnomadg_af, clin_sig
        FROM hg38_vep 
        WHERE chr = ? AND pos = ?
        """
        
        hg19_annotations = pd.read_sql_query(hg19_query, conn, params=[chrom, pos])
        hg38_annotations = pd.read_sql_query(hg38_query, conn, params=[hg38_chrom, hg38_pos])
        
        # Apply comprehensive VEP analysis to this variant
        vep_analysis = analyze_single_variant_vep(
            variant_row, hg19_annotations, hg38_annotations, normalize_transcript_id
        )
        
        if vep_analysis:
            vep_analyses.append(vep_analysis)
    
    return vep_analyses

def analyze_single_variant_vep(variant_row, hg19_annotations, hg38_annotations, normalize_transcript_id):
    """Apply comprehensive VEP analysis to a single variant (cached part)"""
    
    chrom = variant_row['source_chrom']
    pos = variant_row['source_pos']
    
    # Extract genotypes
    ref_allele, alt_allele = extract_genotype_from_alleles(variant_row['source_alleles'])
    
    # Focus on transcript annotations only for transcript analysis
    hg19_transcripts_df = hg19_annotations[hg19_annotations['feature_type'] == 'Transcript']
    hg38_transcripts_df = hg38_annotations[hg38_annotations['feature_type'] == 'Transcript']
    
    # Initialize analysis variables
    same_transcript_consequence_changes = 0
    same_consequence_different_transcripts = 0
    unmatched_consequences = 0
    gene_changes = 0
    impact_changes = 0
    problematic_transcripts_hg19 = []
    problematic_transcripts_hg38 = []
    
    if len(hg19_transcripts_df) > 0 or len(hg38_transcripts_df) > 0:
        # Organize transcript data by build
        hg19_transcripts = {}
        hg38_transcripts = {}
        
        for _, row in hg19_transcripts_df.iterrows():
            hg19_base = normalize_transcript_id(row['feature'])
            if hg19_base:
                hg19_transcripts[hg19_base] = {
                    'consequence': row['consequence'],
                    'impact': row['impact'],
                    'symbol': row['symbol'],
                    'feature_id': row['feature']
                }
        
        for _, row in hg38_transcripts_df.iterrows():
            hg38_base = normalize_transcript_id(row['feature'])
            if hg38_base:
                hg38_transcripts[hg38_base] = {
                    'consequence': row['consequence'],
                    'impact': row['impact'],
                    'symbol': row['symbol'],
                    'feature_id': row['feature']
                }
        
        # Strategy 1: Match by transcript ID (highest priority)
        matched_transcripts = []
        unmatched_hg19 = dict(hg19_transcripts)
        unmatched_hg38 = dict(hg38_transcripts)
        
        for tx_id in set(hg19_transcripts.keys()) & set(hg38_transcripts.keys()):
            hg19_data = hg19_transcripts[tx_id]
            hg38_data = hg38_transcripts[tx_id]
            
            consequence_match = hg19_data['consequence'] == hg38_data['consequence']
            impact_match = hg19_data['impact'] == hg38_data['impact']
            gene_match = hg19_data['symbol'] == hg38_data['symbol']
            
            matched_transcripts.append({
                'transcript': tx_id,
                'consequence_match': consequence_match,
                'impact_match': impact_match,
                'gene_match': gene_match,
                'hg19_consequence': hg19_data['consequence'],
                'hg38_consequence': hg38_data['consequence'],
                'hg19_impact': hg19_data['impact'],
                'hg38_impact': hg38_data['impact']
            })
            
            # FIXED LOGIC: Only count actual differences
            if not consequence_match:
                same_transcript_consequence_changes += 1
                problematic_transcripts_hg19.append(f"{hg19_data['feature_id']}({hg19_data['consequence']})")
                problematic_transcripts_hg38.append(f"{hg38_data['feature_id']}({hg38_data['consequence']})")
            
            if not impact_match:
                impact_changes += 1
            if not gene_match:
                gene_changes += 1
            
            # Remove from unmatched
            unmatched_hg19.pop(tx_id, None)
            unmatched_hg38.pop(tx_id, None)
        
        # Strategy 2: Match by consequence type (moderate priority)
        # Group remaining transcripts by consequence
        hg19_by_consequence = {}
        for tx_id, data in unmatched_hg19.items():
            cons = data['consequence']
            if cons not in hg19_by_consequence:
                hg19_by_consequence[cons] = []
            hg19_by_consequence[cons].append(tx_id)
        
        hg38_by_consequence = {}
        for tx_id, data in unmatched_hg38.items():
            cons = data['consequence']
            if cons not in hg38_by_consequence:
                hg38_by_consequence[cons] = []
            hg38_by_consequence[cons].append(tx_id)
        
        # Count matched consequences
        matched_consequences = set(hg19_by_consequence.keys()) & set(hg38_by_consequence.keys())
        same_consequence_different_transcripts = len(matched_consequences)
        
        # Remove matched consequences from remaining
        remaining_hg19 = dict(unmatched_hg19)
        remaining_hg38 = dict(unmatched_hg38)
        
        for consequence in matched_consequences:
            for tx_id in hg19_by_consequence[consequence]:
                remaining_hg19.pop(tx_id, None)
            for tx_id in hg38_by_consequence[consequence]:
                remaining_hg38.pop(tx_id, None)
        
        # Strategy 3: Count unmatched consequences
        unmatched_hg19_consequences = set(data['consequence'] for data in remaining_hg19.values())
        unmatched_hg38_consequences = set(data['consequence'] for data in remaining_hg38.values())
        
        if len(unmatched_hg19_consequences) > 0 or len(unmatched_hg38_consequences) > 0:
            unmatched_consequences = 1
            
            # Add unmatched transcripts to problematic lists
            for tx_id, data in remaining_hg19.items():
                problematic_transcripts_hg19.append(f"{data['feature_id']}({data['consequence']})")
            for tx_id, data in remaining_hg38.items():
                problematic_transcripts_hg38.append(f"{data['feature_id']}({data['consequence']})")
    
    # Get representative VEP information
    hg19_gene = hg19_transcripts_df['symbol'].iloc[0] if len(hg19_transcripts_df) > 0 else (hg19_annotations['symbol'].iloc[0] if len(hg19_annotations) > 0 else '')
    hg38_gene = hg38_transcripts_df['symbol'].iloc[0] if len(hg38_transcripts_df) > 0 else (hg38_annotations['symbol'].iloc[0] if len(hg38_annotations) > 0 else '')
    hg19_consequence = hg19_transcripts_df['consequence'].iloc[0] if len(hg19_transcripts_df) > 0 else (hg19_annotations['consequence'].iloc[0] if len(hg19_annotations) > 0 else '')
    hg38_consequence = hg38_transcripts_df['consequence'].iloc[0] if len(hg38_transcripts_df) > 0 else (hg38_annotations['consequence'].iloc[0] if len(hg38_annotations) > 0 else '')
    hg19_impact = hg19_annotations['impact'].iloc[0] if len(hg19_annotations) > 0 else ''
    hg38_impact = hg38_annotations['impact'].iloc[0] if len(hg38_annotations) > 0 else ''
    
    # Clinical significance and pathogenicity predictions
    hg19_clin_sig = hg19_annotations['clin_sig'].iloc[0] if len(hg19_annotations) > 0 else ''
    hg38_clin_sig = hg38_annotations['clin_sig'].iloc[0] if len(hg38_annotations) > 0 else ''
    hg19_sift = hg19_annotations['sift'].iloc[0] if len(hg19_annotations) > 0 else ''
    hg38_sift = hg38_annotations['sift'].iloc[0] if len(hg38_annotations) > 0 else ''
    hg19_polyphen = hg19_annotations['polyphen'].iloc[0] if len(hg19_annotations) > 0 else ''
    hg38_polyphen = hg38_annotations['polyphen'].iloc[0] if len(hg38_annotations) > 0 else ''
    
    # Parse pathogenicity predictions
    hg19_sift_pred, hg19_sift_score = parse_sift_prediction(hg19_sift)
    hg38_sift_pred, hg38_sift_score = parse_sift_prediction(hg38_sift)
    hg19_polyphen_pred, hg19_polyphen_score = parse_polyphen_prediction(hg19_polyphen)
    hg38_polyphen_pred, hg38_polyphen_score = parse_polyphen_prediction(hg38_polyphen)
    
    # Check for clinical significance changes
    hg19_is_pathogenic = is_pathogenic_clinical_significance(hg19_clin_sig)
    hg38_is_pathogenic = is_pathogenic_clinical_significance(hg38_clin_sig)
    hg19_is_benign = is_benign_clinical_significance(hg19_clin_sig)
    hg38_is_benign = is_benign_clinical_significance(hg38_clin_sig)
    
    clin_sig_change = ''
    if hg19_is_benign and hg38_is_pathogenic:
        clin_sig_change = 'BENIGN_TO_PATHOGENIC'
    elif hg19_is_pathogenic and hg38_is_benign:
        clin_sig_change = 'PATHOGENIC_TO_BENIGN'
    elif hg19_is_benign != hg38_is_benign or hg19_is_pathogenic != hg38_is_pathogenic:
        clin_sig_change = 'OTHER_CHANGE'
    
    # Check for SIFT/PolyPhen changes
    sift_change = ''
    if hg19_sift_pred and hg38_sift_pred and hg19_sift_pred != hg38_sift_pred:
        if (hg19_sift_pred == 'tolerated' and hg38_sift_pred == 'deleterious') or \
           (hg19_sift_pred == 'deleterious' and hg38_sift_pred == 'tolerated'):
            sift_change = f'{hg19_sift_pred.upper()}_TO_{hg38_sift_pred.upper()}'
    
    polyphen_change = ''
    if hg19_polyphen_pred and hg38_polyphen_pred and hg19_polyphen_pred != hg38_polyphen_pred:
        if (hg19_polyphen_pred == 'benign' and hg38_polyphen_pred in ['possibly_damaging', 'probably_damaging']) or \
           (hg19_polyphen_pred in ['possibly_damaging', 'probably_damaging'] and hg38_polyphen_pred == 'benign'):
            polyphen_change = f'{hg19_polyphen_pred.upper()}_TO_{hg38_polyphen_pred.upper()}'
    
    # Store variant VEP analysis (cached part - no scores)
    vep_analysis = {
        'mapping_status': variant_row['mapping_status'],
        'source_chrom': chrom,
        'source_pos': pos,
        'source_alleles': variant_row['source_alleles'],
        'gt_hg19': ref_allele,
        'gt_hg38': alt_allele,
        'bcftools_hg38_ref': variant_row.get('bcftools_hg38_ref', ''),
        'bcftools_hg38_alt': variant_row.get('bcftools_hg38_alt', ''),
        'pos_match': variant_row['pos_match'],
        'gt_match': variant_row['gt_match'],
        'flip': variant_row['flip'],
        'swap': variant_row['swap'],
        'liftover_hg38_pos': variant_row['liftover_hg38_pos'],
        'bcftools_hg38_pos': variant_row['bcftools_hg38_pos'],
        'pos_difference': variant_row['pos_difference'],
        
        # VEP analysis results
        'transcript_pairs_analyzed': len(hg19_transcripts_df) + len(hg38_transcripts_df),
        'same_transcript_consequence_changes': same_transcript_consequence_changes,
        'same_consequence_different_transcripts': same_consequence_different_transcripts,
        'unmatched_consequences': unmatched_consequences,
        'gene_changes': gene_changes,
        'impact_changes': impact_changes,
        'problematic_transcripts_hg19': '; '.join(problematic_transcripts_hg19) if problematic_transcripts_hg19 else '',
        'problematic_transcripts_hg38': '; '.join(problematic_transcripts_hg38) if problematic_transcripts_hg38 else '',
        
        # Representative VEP information
        'hg19_gene': hg19_gene,
        'hg38_gene': hg38_gene,
        'hg19_consequence': hg19_consequence,
        'hg38_consequence': hg38_consequence,
        'hg19_impact': hg19_impact or '',
        'hg38_impact': hg38_impact or '',
        'hg19_clin_sig': hg19_clin_sig or '',
        'hg38_clin_sig': hg38_clin_sig or '',
        'clin_sig_change': clin_sig_change,
        'hg19_gnomad_af': hg19_annotations['gnomadg_af'].iloc[0] if len(hg19_annotations) > 0 else None,
        'hg38_gnomad_af': hg38_annotations['gnomadg_af'].iloc[0] if len(hg38_annotations) > 0 else None,
        'hg19_sift': hg19_sift or '',
        'hg38_sift': hg38_sift or '',
        'sift_change': sift_change,
        'hg19_polyphen': hg19_polyphen or '',
        'hg38_polyphen': hg38_polyphen or '',
        'polyphen_change': polyphen_change,
        
        # Pathogenicity flags
        'hg19_is_pathogenic': hg19_is_pathogenic,
        'hg38_is_pathogenic': hg38_is_pathogenic,
        'hg19_is_benign': hg19_is_benign,
        'hg38_is_benign': hg38_is_benign
    }
    
    return vep_analysis

def get_impact_numeric_value(impact):
    """Convert impact to numeric value for magnitude calculation"""
    impact_values = {'HIGH': 4, 'MODERATE': 3, 'LOW': 2, 'MODIFIER': 1}
    return impact_values.get(impact, 0)

def calculate_impact_transition_magnitude(hg19_impact, hg38_impact):
    """Calculate magnitude of impact transition for clinical significance"""
    hg19_val = get_impact_numeric_value(hg19_impact)
    hg38_val = get_impact_numeric_value(hg38_impact)
    magnitude = abs(hg19_val - hg38_val)
    
    # Check if transition involves MODERATE or HIGH (clinically significant)
    involves_moderate_or_high = (
        hg19_impact in ['MODERATE', 'HIGH'] or 
        hg38_impact in ['MODERATE', 'HIGH']
    )
    
    return magnitude, involves_moderate_or_high

def calculate_scores_from_vep_analysis(vep_analysis_df):
    """Calculate priority scores from cached VEP analysis data with magnitude-based impact scoring"""
    print("Calculating priority scores with magnitude-based impact transitions...")
    
    scored_variants = []
    
    for _, row in vep_analysis_df.iterrows():
        # Calculate base priority score
        priority_score = 0
        
        # Position and genotype discordance scoring
        if row['pos_match'] == 0:
            priority_score += 3
        if row['pos_difference'] > 10:
            priority_score += 2
        if row['pos_difference'] > 100:
            priority_score += 3
        if row['gt_match'] == 0:
            priority_score += 3
        
        # BCFtools liftover swap values: 1 = swapped, -1 = swap failed, NA = no swap needed
        swap_str = str(row['swap']) if pd.notna(row['swap']) else 'NA'
        if swap_str == '1':  # REF/ALT alleles were swapped during liftover
            priority_score += 2
        
        # Core functional changes (highest priority)
        priority_score += row['same_transcript_consequence_changes'] * 10  # CRITICAL
        
        # MAGNITUDE-BASED IMPACT SCORING (NEW)
        hg19_impact = row['hg19_impact']
        hg38_impact = row['hg38_impact']
        impact_magnitude, is_clinically_significant = calculate_impact_transition_magnitude(hg19_impact, hg38_impact)
        
        if row['impact_changes'] > 0 and is_clinically_significant:
            # Clinically significant impact transitions
            if impact_magnitude == 1:
                if 'HIGH' in [hg19_impact, hg38_impact] and 'MODERATE' in [hg19_impact, hg38_impact]:
                    priority_score += 10  # HIGH ↔ MODERATE
                elif 'MODERATE' in [hg19_impact, hg38_impact] and 'LOW' in [hg19_impact, hg38_impact]:
                    priority_score += 8   # MODERATE ↔ LOW
            elif impact_magnitude == 2:
                if 'HIGH' in [hg19_impact, hg38_impact]:
                    priority_score += 12  # HIGH ↔ LOW
                elif 'MODERATE' in [hg19_impact, hg38_impact]:
                    priority_score += 10  # MODERATE ↔ MODIFIER
            elif impact_magnitude == 3:
                priority_score += 15      # HIGH ↔ MODIFIER (most significant)
        elif row['impact_changes'] > 0:
            # Non-clinically significant transitions (annotation noise)
            priority_score += 1           # LOW ↔ MODIFIER (minimal weight)
        
        priority_score += row['unmatched_consequences'] * 4               # INVESTIGATE
        
        # Clinical significance changes
        if row['clin_sig_change'] == 'BENIGN_TO_PATHOGENIC':
            priority_score += 10
        elif row['clin_sig_change'] == 'PATHOGENIC_TO_BENIGN':
            priority_score += 8
        elif row['clin_sig_change'] == 'OTHER_CHANGE':
            priority_score += 5
        
        # Pathogenicity prediction changes
        if row['sift_change']:
            priority_score += 5
        if row['polyphen_change']:
            priority_score += 5
        
        # Gene changes - CONDITIONAL on clinical significance (NEW LOGIC)
        if row['gene_changes'] > 0:
            if is_clinically_significant or row['clin_sig_change'] or row['sift_change'] or row['polyphen_change']:
                # Gene change is clinically relevant - use impact-based scoring
                if hg19_impact == 'HIGH' or hg38_impact == 'HIGH':
                    priority_score += row['gene_changes'] * 8
                elif hg19_impact == 'MODERATE' or hg38_impact == 'MODERATE':
                    priority_score += row['gene_changes'] * 4
                elif hg19_impact == 'LOW' or hg38_impact == 'LOW':
                    priority_score += row['gene_changes'] * 2
                else:
                    priority_score += row['gene_changes'] * 3  # Mixed or unknown
            else:
                # Gene change likely annotation synonym - minimal weight
                priority_score += row['gene_changes'] * 0.1
        
        # Other transcript-level issues
        priority_score += row['same_consequence_different_transcripts'] * 3  # Moderate priority
        
        # Clinical significance bonus (if any clinical data available)
        has_clinical_data = (
            not pd.isna(row['hg19_clin_sig']) and row['hg19_clin_sig'] not in ['', '-'] or
            not pd.isna(row['hg38_clin_sig']) and row['hg38_clin_sig'] not in ['', '-']
        )
        if has_clinical_data:
            priority_score += 2
        
        # Mapping status and impact adjustments
        if row['mapping_status'] == 'REGION':
            priority_score += 1
        
        if row['hg19_impact'] == 'HIGH' or row['hg38_impact'] == 'HIGH':
            priority_score += 2
        
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
            priority_score *= 0.1  # 90% reduction for benign variants
        elif has_pathogenic_evidence:
            priority_score *= 2.0  # 2x boost for pathogenic variants
        
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

def create_prioritization_plots(df, conn, output_dir):
    """Create enhanced visualization plots for variant prioritization analysis"""
    print("Creating enhanced prioritization visualizations...")

    # Define consistent colorblind-friendly palette
    colors_priority = ['#d62728', '#ff7f0e', '#2ca02c', '#1f77b4', '#9467bd']
    colors_clinical = ['#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf', '#ff7f0e', '#2ca02c']

    # Create enhanced visualization with 4 plots
    fig, axes = plt.subplots(2, 2, figsize=(20, 12))
    fig.suptitle('Enhanced Variant Prioritization Analysis', fontsize=16, fontweight='bold')

    # PLOT 1: Priority categories bar chart (keep as is)
    print("1. Creating priority categories bar chart...")

    # Calculate total variants including concordant ones
    total_query = """SELECT COUNT(*) as total FROM comparison"""
    total_result = pd.read_sql_query(total_query, conn)
    total_variants = total_result.iloc[0]['total']

    category_counts = df['priority_category'].value_counts()

    # Add actual concordant variants count
    discordant_count = len(df)
    concordant_count = max(0, total_variants - discordant_count)
    if concordant_count > 0:
        category_counts['CONCORDANT'] = concordant_count

    # Order categories by priority
    category_order = ['CRITICAL', 'HIGH', 'MODERATE', 'INVESTIGATE', 'LOW', 'CONCORDANT']
    category_counts = category_counts.reindex([cat for cat in category_order if cat in category_counts.index])

    # Create bar chart with percentages
    bars = axes[0,0].bar(range(len(category_counts)), category_counts.values,
                      color=colors_priority[:len(category_counts)], alpha=0.8)

    # Add percentage labels on bars
    for i, (bar, count) in enumerate(zip(bars, category_counts.values)):
        pct = count / category_counts.sum() * 100
        axes[0,0].text(bar.get_x() + bar.get_width()/2, bar.get_height() + count*0.01,
                    f'{pct:.1f}%\n({count:,})', ha='center', va='bottom',
                    fontweight='bold', fontsize=9)

    axes[0,0].set_xticks(range(len(category_counts)))
    axes[0,0].set_xticklabels(category_counts.index, rotation=45, ha='right', fontsize=10)
    axes[0,0].set_ylabel('Variant Count', fontweight='bold')
    axes[0,0].set_title('Variant Priority Categories\n(Clinical Review Levels)', fontsize=12, fontweight='bold')
    axes[0,0].grid(True, alpha=0.3, axis='y')

    # PLOT 2: Priority Categories × Clinical Evidence (NEW)
    print("2. Creating priority categories × clinical evidence...")

    def get_clinical_evidence_color(row):
        """Smart gap filling for clinical evidence"""
        # Priority 1: Clinical significance
        if row['clin_sig_change']:
            return 'Clinical Sig Change'
        elif row['hg19_is_pathogenic'] or row['hg38_is_pathogenic']:
            return 'Pathogenic'
        elif row['hg19_is_benign'] or row['hg38_is_benign']:
            return 'Benign'

        # Priority 2: SIFT + PolyPhen (combined)
        elif row['sift_change'] or row['polyphen_change']:
            return 'Prediction Change'
        elif ('deleterious' in str(row['hg19_sift']).lower() or 'deleterious' in str(row['hg38_sift']).lower() or
              'probably_damaging' in str(row['hg19_polyphen']).lower() or 'probably_damaging' in str(row['hg38_polyphen']).lower()):
            return 'Pred. Pathogenic'
        elif ('tolerated' in str(row['hg19_sift']).lower() or 'tolerated' in str(row['hg38_sift']).lower() or
              'benign' in str(row['hg19_polyphen']).lower() or 'benign' in str(row['hg38_polyphen']).lower()):
            return 'Pred. Benign'

        # Priority 3: Impact level fallback
        elif row['hg19_impact'] == 'HIGH' or row['hg38_impact'] == 'HIGH':
            return 'High Impact'
        else:
            return 'No Evidence'

    if len(df) > 0:
        df['clinical_evidence'] = df.apply(get_clinical_evidence_color, axis=1)

        # Create contingency table
        evidence_priority = pd.crosstab(df['priority_category'], df['clinical_evidence'])
        evidence_priority = evidence_priority.reindex(
            index=[cat for cat in category_order[:-1] if cat in evidence_priority.index],  # Exclude CONCORDANT
            fill_value=0
        )

        # Create stacked bar chart
        evidence_priority.plot(kind='bar', ax=axes[0,1], stacked=True,
                              color=colors_clinical[:len(evidence_priority.columns)], alpha=0.8)

        axes[0,1].set_title('Priority Categories by Clinical Evidence\n(Intelligent Gap Filling)',
                           fontsize=12, fontweight='bold')
        axes[0,1].set_xlabel('Priority Category', fontweight='bold')
        axes[0,1].set_ylabel('Variant Count', fontweight='bold')
        axes[0,1].legend(title='Clinical Evidence', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
        axes[0,1].tick_params(axis='x', rotation=45)
        axes[0,1].grid(True, alpha=0.3, axis='y')
    else:
        axes[0,1].text(0.5, 0.5, 'No variants\nfound', ha='center', va='center',
                      fontsize=12, transform=axes[0,1].transAxes)
        axes[0,1].set_title('Priority Categories by Clinical Evidence\n(Intelligent Gap Filling)',
                           fontsize=12, fontweight='bold')

    # PLOT 3: Primary Discordance Types (aggregated, no numbers)
    print("3. Creating primary discordance types visualization...")

    def categorize_discordance_primary(summary):
        """Categorize discordance without numbers for aggregation"""
        if pd.isna(summary) or summary == '':
            return 'Other'

        summary_lower = summary.lower()

        if 'same transcript consequence' in summary_lower:
            return 'Same Transcript\nConsequence Change'
        elif 'gene annotation' in summary_lower:
            return 'Gene Annotation\nChange'
        elif 'impact level' in summary_lower:
            return 'Impact Level\nChange'
        elif 'clinical significance' in summary_lower:
            return 'Clinical Significance\nChange'
        elif 'sift change' in summary_lower or 'polyphen change' in summary_lower:
            return 'Pathogenicity\nPrediction Change'
        elif 'unmatched consequences' in summary_lower:
            return 'Unmatched\nConsequences'
        elif 'position mismatch' in summary_lower:
            return 'Position\nMismatch'
        elif 'genotype mismatch' in summary_lower:
            return 'Genotype\nMismatch'
        else:
            return 'Other/Multiple'

    if len(df) > 0:
        df['discordance_primary'] = df['discordance_summary'].apply(categorize_discordance_primary)
        discordance_counts = df['discordance_primary'].value_counts()

        # Order categories by clinical importance
        discordance_order = [
            'Same Transcript\nConsequence Change',
            'Clinical Significance\nChange',
            'Impact Level\nChange',
            'Gene Annotation\nChange',
            'Pathogenicity\nPrediction Change',
            'Unmatched\nConsequences',
            'Position\nMismatch',
            'Genotype\nMismatch',
            'Other/Multiple'
        ]

        discordance_counts = discordance_counts.reindex([cat for cat in discordance_order if cat in discordance_counts.index])

        bars = axes[1,0].bar(range(len(discordance_counts)), discordance_counts.values,
                            color=colors_priority[:len(discordance_counts)], alpha=0.8)

        axes[1,0].set_xlabel('Primary Discordance Type', fontweight='bold')
        axes[1,0].set_ylabel('Variant Count', fontweight='bold')
        axes[1,0].set_title('Variants by Primary Discordance Type\n(Aggregated Categories)',
                           fontsize=12, fontweight='bold')
        axes[1,0].set_xticks(range(len(discordance_counts)))
        axes[1,0].set_xticklabels(discordance_counts.index, rotation=45, ha='right', fontsize=9)
        axes[1,0].grid(True, alpha=0.3, axis='y')

        # Add count labels on bars
        for i, (bar, count) in enumerate(zip(bars, discordance_counts.values)):
            axes[1,0].text(bar.get_x() + bar.get_width()/2, bar.get_height() + count*0.01,
                          str(count), ha='center', va='bottom', fontweight='bold', fontsize=10)
    else:
        axes[1,0].text(0.5, 0.5, 'No discordant variants\nfound', ha='center', va='center',
                      fontsize=12, transform=axes[1,0].transAxes)
        axes[1,0].set_title('Variants by Primary Discordance Type\n(Aggregated Categories)',
                           fontsize=12, fontweight='bold')

    # PLOT 4: Clinical Evidence Transitions (NEW CATEGORIES)
    print("4. Creating clinical evidence transitions...")

    def categorize_clinical_transitions(row):
        """Categorize variants by clinical evidence transitions"""
        # Check for clinical significance transitions
        if row['clin_sig_change'] == 'PATHOGENIC_TO_BENIGN':
            return 'Pathogenic → Benign'
        elif row['clin_sig_change'] == 'BENIGN_TO_PATHOGENIC':
            return 'Benign → Pathogenic'
        elif row['hg19_is_pathogenic'] and row['hg38_is_pathogenic']:
            return 'Stable Pathogenic'
        elif row['hg19_is_benign'] and row['hg38_is_benign']:
            return 'Stable Benign'
        elif row['sift_change'] or row['polyphen_change']:
            return 'Prediction Changes'
        elif row['hg19_impact'] == 'HIGH' or row['hg38_impact'] == 'HIGH':
            return 'High Impact Only'
        else:
            return 'Other/Unknown'

    if len(df) > 0:
        df['clinical_transitions'] = df.apply(categorize_clinical_transitions, axis=1)
        transition_counts = df['clinical_transitions'].value_counts()

        # Order categories by clinical concern level
        transition_order = [
            'Pathogenic → Benign',
            'Benign → Pathogenic',
            'Stable Pathogenic',
            'Prediction Changes',
            'High Impact Only',
            'Stable Benign',
            'Other/Unknown'
        ]

        transition_counts = transition_counts.reindex([cat for cat in transition_order if cat in transition_counts.index])

        bars = axes[1,1].bar(range(len(transition_counts)), transition_counts.values,
                            color=colors_clinical[:len(transition_counts)], alpha=0.8)

        axes[1,1].set_xlabel('Clinical Evidence Transition', fontweight='bold')
        axes[1,1].set_ylabel('Variant Count', fontweight='bold')
        axes[1,1].set_title('Clinical Evidence Transitions\n(Dynamic vs Static Evidence)',
                           fontsize=12, fontweight='bold')
        axes[1,1].set_xticks(range(len(transition_counts)))
        axes[1,1].set_xticklabels(transition_counts.index, rotation=45, ha='right', fontsize=9)
        axes[1,1].grid(True, alpha=0.3, axis='y')

        # Add count labels on bars
        for i, (bar, count) in enumerate(zip(bars, transition_counts.values)):
            axes[1,1].text(bar.get_x() + bar.get_width()/2, bar.get_height() + count*0.01,
                          str(count), ha='center', va='bottom', fontweight='bold', fontsize=10)
    else:
        axes[1,1].text(0.5, 0.5, 'No variants\nfound', ha='center', va='center',
                      fontsize=12, transform=axes[1,1].transAxes)
        axes[1,1].set_title('Clinical Evidence Transitions\n(Dynamic vs Static Evidence)',
                           fontsize=12, fontweight='bold')

    plt.tight_layout()
    plt.savefig(output_dir / 'variant_prioritization_plots.png', dpi=300, bbox_inches='tight')
    plt.show()

    # Print enhanced summary statistics
    print("\n=== ENHANCED VISUALIZATION SUMMARY ===")
    print(f"Total variants in database: {total_variants:,}")
    print(f"Discordant variants analyzed: {len(df):,}")
    print(f"Concordant variants: {concordant_count:,}")

    if len(df) > 0:
        print("\nPriority category breakdown:")
        for category, count in category_counts.items():
            pct = count / category_counts.sum() * 100
            print(f"  {category}: {count:,} ({pct:.1f}%)")

        if 'clinical_transitions' in df.columns:
            print(f"\nClinical evidence transitions:")
            for transition, count in transition_counts.items():
                pct = count / len(df) * 100
                print(f"  {transition}: {count:,} ({pct:.1f}%)")

        if 'discordance_primary' in df.columns:
            print(f"\nPrimary discordance types:")
            for disc_type, count in discordance_counts.items():
                pct = count / len(df) * 100
                print(f"  {disc_type.replace(chr(10), ' ')}: {count:,} ({pct:.1f}%)")

    print("✓ Enhanced prioritization visualizations completed")

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
    
    # Liftover information - Convert floats to integers properly
    def safe_int_convert(series):
        """Convert float to int, handling NaN values"""
        return series.apply(lambda x: int(x) if pd.notna(x) and x != '' else '')
    
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
        df_full = calculate_variant_scores(conn, cache_file, args.force)
        
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
            create_prioritization_plots(df_full, conn, output_dir)
        
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
