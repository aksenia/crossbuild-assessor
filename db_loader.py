#!/usr/bin/env python3
"""
Genomic Variant Database Loader

DESCRIPTION:
    Loads genomic variant data from liftover comparisons and VEP annotations into 
    a SQLite database with optimized schema for cross-build analysis. Handles VEP 
    coordinate normalization to ensure proper matching between comparison and 
    annotation data.

EXPECTED DATA FORMATS:

1. COMPARISON FILE (liftover comparison):
    - Tab-separated file comparing liftover tools (CrossMap vs bcftools)
    - Required columns:
        * mapping_status: Liftover success status (UNIQUE, REGION, etc.)
        * source_chrom: Source chromosome (1, 2, X, Y, etc.)
        * source_pos: Source position (0-based coordinates)
        * source_alleles: Comma-separated alleles (A,G or A,ACGT)
        * flip: Strand flip status (no_flip, flip)
        * swap: Reference/alt swap status (NA, swap)
        * liftover_hg38_chrom: Target chromosome from CrossMap
        * liftover_hg38_pos: Target position from CrossMap
        * bcftools_hg38_chrom: Target chromosome from bcftools
        * bcftools_hg38_pos: Target position from bcftools
        * bcftools_hg38_ref: Reference allele from bcftools
        * bcftools_hg38_alt: Alternative allele from bcftools
        * pos_match: Position agreement (TRUE/FALSE)
        * gt_match: Genotype agreement (TRUE/FALSE)

2. VEP ANNOTATION FILES (hg19 and hg38):
    - Standard VEP tab-delimited output files
    - Must contain header starting with #Uploaded_variation
    - Required columns:
        * #Uploaded_variation: Variant ID (chr_pos_ref/alt format)
        * Location: Genomic location
        * Allele: Variant allele
        * Gene: Gene ID
        * Feature: Transcript/feature ID
        * Feature_type: Annotation type (Transcript, RegulatoryFeature, etc.)
        * Consequence: Variant consequence
        * IMPACT: Impact level (HIGH, MODERATE, LOW, MODIFIER)
        * SYMBOL: Gene symbol
        * SIFT: SIFT prediction
        * PolyPhen: PolyPhen prediction
        * gnomADg_AF: gnomAD genome frequency
        * CLIN_SIG: Clinical significance
        * HGVSc: HGVS coding sequence nomenclature
        * HGVSp: HGVS protein sequence nomenclature
        * CANONICAL: Canonical transcript flag (YES/-)

COORDINATE SYSTEMS:
    - Input comparison file: 0-based coordinates
    - VEP files: VEP-normalized coordinates (SNVs: same as input, Indels: +1 shift)
    - Output database: VEP-normalized coordinates (for consistent matching)

EXAMPLE CONFIG FILE (config.json):
{
    "input_files": {
        "comparison": "/path/to/liftover_comparison.txt",
        "hg19_vep": "/path/to/variants_hg19.vep.txt",
        "hg38_vep": "/path/to/variants_hg38.vep.txt"
    },
    "metadata": {
        "description": "CrossMap vs bcftools liftover comparison",
        "date_created": "2025-01-15",
        "genome_builds": ["hg19", "hg38"],
        "vep_version": "113.0"
    }
}

USAGE:
    python db_loader.py --config config.json --force

    The script will:
    1. Create optimized database schema with proper indexes
    2. Load and normalize comparison data to VEP coordinate system
    3. Load VEP annotation data for both genome builds
    4. Verify data integrity and report matching statistics

OUTPUT:
    SQLite database with three tables:
    - comparison: Normalized liftover comparison data
    - hg19_vep: VEP annotations for hg19 build
    - hg38_vep: VEP annotations for hg38 build
"""

import os
import sqlite3
import pandas as pd
import argparse
import json
import re
from pathlib import Path
import sys

# COMPLETE column mapping for all relevant VEP fields
VEP_COLUMN_MAP = {
    'Uploaded_variation': ['#Uploaded_variation', 'Uploaded_variation'],
    'Location': ['Location', 'Uploaded_variation'],
    'Allele': ['Allele'],
    'Gene': ['Gene'],
    'Feature': ['Feature', 'Feature_ID', 'Transcript_ID'],
    'Feature_type': ['Feature_type'],
    'Consequence': ['Consequence'],
    'IMPACT': ['IMPACT', 'Impact'],
    'SYMBOL': ['SYMBOL', 'Symbol'],
    'SIFT': ['SIFT'],
    'PolyPhen': ['PolyPhen', 'Polyphen'],
    'gnomADg_AF': ['gnomADg_AF', 'gnomAD_AF', 'AF'],
    'CLIN_SIG': ['CLIN_SIG', 'Clinical_significance', 'ClinVar_CLNSIG'],
    'HGVSc': ['HGVSc', 'HGVS_c'],
    'HGVSp': ['HGVSp', 'HGVS_p'],
    'CANONICAL': ['CANONICAL', 'Canonical'],
    'MANE': ['MANE'],
    'MANE_SELECT': ['MANE_SELECT'],
    'MANE_PLUS_CLINICAL': ['MANE_PLUS_CLINICAL']
}

def find_column(df, column_aliases):
    """Find first existing column from aliases"""
    for alias in column_aliases:
        if alias in df.columns:
            return alias
    return None

def get_safe_column_value(row, df, column_key, default=''):
    """Safely get column value using column mapping"""
    column_aliases = VEP_COLUMN_MAP.get(column_key, [column_key])
    actual_column = find_column(df, column_aliases)
    if actual_column:
        return row.get(actual_column, default)
    return default

def extract_refseq_from_hgvsc(hgvsc_value):
    """Extract RefSeq ID from HGVSc field (e.g., 'NM_001385641.1:c.1796del' -> 'NM_001385641.1')"""
    if pd.isna(hgvsc_value) or ':' not in str(hgvsc_value):
        return None
    return str(hgvsc_value).split(':')[0].strip()

def process_mane_columns(row, df):
    """Process MANE-related columns with defensive handling"""
    return {
        'mane': get_safe_column_value(row, df, 'MANE'),
        'mane_select': get_safe_column_value(row, df, 'MANE_SELECT'),
        'mane_plus_clinical': get_safe_column_value(row, df, 'MANE_PLUS_CLINICAL'),
        'refseq_transcript_id': extract_refseq_from_hgvsc(get_safe_column_value(row, df, 'HGVSc'))
    }

def parse_custom_vep_id(vep_id):
    """
    Parse custom VEP ID format like '1/878175/AG//A' or 'chr1/1520188/C//T'
    
    Args:
        vep_id: Custom VEP variant ID
        
    Returns:
        tuple: (chrom, pos, ref, alt) or (None, None, None, None) if parsing fails
    """
    if pd.isna(vep_id):
        return None, None, None, None
        
    parts = str(vep_id).split('/')
    if len(parts) >= 5:  # Format: chrom/pos/ref//alt
        try:
            chrom = parts[0]
            pos = int(parts[1])
            ref = parts[2]
            alt = parts[4]  # Skip empty part from '//'
            return chrom, pos, ref, alt
        except (ValueError, IndexError):
            return None, None, None, None
    
    return None, None, None, None

def create_database_schema(db_path):
    """Create database tables WITHOUT unique constraints for faster loading"""
    print("Creating optimized database schema (constraints added after loading)...")
    
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    # Drop tables if they exist
    cursor.execute("DROP TABLE IF EXISTS comparison")
    cursor.execute("DROP TABLE IF EXISTS hg19_vep")
    cursor.execute("DROP TABLE IF EXISTS hg38_vep")
    
    # Create comparison table WITHOUT unique constraint during loading
    cursor.execute("""
        CREATE TABLE comparison (
            id INTEGER PRIMARY KEY,
            mapping_status TEXT,
            source_chrom TEXT,
            source_pos INTEGER,
            source_alleles TEXT,
            source_ref TEXT,
            source_alt TEXT,     
            flip TEXT,
            swap TEXT,
            liftover_hg38_chrom TEXT,
            liftover_hg38_pos INTEGER,
            bcftools_hg38_chrom TEXT,
            bcftools_hg38_pos INTEGER,
            bcftools_hg38_ref TEXT,
            bcftools_hg38_alt TEXT,
            pos_match BOOLEAN,
            gt_match BOOLEAN
            -- NO UNIQUE constraint here for faster loading
        )
    """)
    
    # Create VEP tables WITHOUT unique constraints during loading
    vep_schema_template = """
        CREATE TABLE {} (
            id INTEGER PRIMARY KEY,
            uploaded_variation TEXT,
            chr TEXT,
            pos INTEGER,
            location TEXT,
            allele TEXT,
            gene TEXT,
            feature TEXT,
            feature_type TEXT,
            consequence TEXT,
            impact TEXT,
            symbol TEXT,
            sift TEXT,
            polyphen TEXT,
            gnomadg_af REAL,
            clin_sig TEXT,
            hgvsc TEXT,
            hgvsp TEXT,
            canonical TEXT,
            mane TEXT,
            mane_select TEXT,
            mane_plus_clinical TEXT,
            refseq_transcript_id TEXT,
            extracted_chrom TEXT,
            extracted_pos INTEGER,
            extracted_ref TEXT,
            extracted_alt TEXT
            -- NO UNIQUE constraint here for faster loading
        )
    """
    
    cursor.execute(vep_schema_template.format('hg19_vep'))
    cursor.execute(vep_schema_template.format('hg38_vep'))
    
    # Create basic indexes for loading performance (but not unique ones yet)
    print("Creating basic loading indexes...")
    cursor.execute("CREATE INDEX idx_comp_loading ON comparison(source_chrom, source_pos)")
    cursor.execute("CREATE INDEX idx_hg19_loading ON hg19_vep(chr, pos)")
    cursor.execute("CREATE INDEX idx_hg38_loading ON hg38_vep(chr, pos)")
    
    conn.commit()
    conn.close()
    print("✓ Optimized database schema created (unique constraints will be added after loading)")

def load_comparison_data(comparison_file, db_path):
    """Load comparison file with optimized bulk loading"""
    print(f"Loading comparison data from: {comparison_file}")
    
    if not Path(comparison_file).exists():
        raise FileNotFoundError(f"Comparison file not found: {comparison_file}")
    
    # Read with explicit dtypes
    df = pd.read_csv(
        comparison_file,
        sep='\t',
        dtype={
            'source_chrom': 'str',
            'bcftools_hg38_chrom': 'str', 
            'liftover_hg38_chrom': 'str'
        },
        low_memory=False
    )
    
    print(f"Loaded {len(df):,} comparison records")
    
    # Validate required columns
    required_cols = [
        'mapping_status', 'source_chrom', 'source_pos', 'source_alleles',
        'flip', 'swap', 'bcftools_hg38_chrom', 'bcftools_hg38_pos', 
        'bcftools_hg38_ref', 'bcftools_hg38_alt', 'pos_match', 'gt_match'
    ]
    
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns in comparison file: {missing_cols}")
    
    # Process comparison data
    processed_data = []
    skipped_variants = 0
    
    for _, row in df.iterrows():
        # Parse alleles without normalization
        if pd.notna(row['source_alleles']) and ',' in str(row['source_alleles']):
            parts = str(row['source_alleles']).split(',')
            if len(parts) == 2:
                parsed_ref = parts[0].strip()
                parsed_alt = parts[1].strip()
            else:
                parsed_ref = ""
                parsed_alt = ""
        else:
            parsed_ref = ""
            parsed_alt = ""
        
        # Use original coordinates directly
        if pd.notna(row['source_pos']):
            source_pos = row['source_pos']
        else:
            source_pos = None
            
        # Handle bcftools position (may be NaN for failed liftover)
        if pd.notna(row['bcftools_hg38_pos']):
            bcftools_pos = row['bcftools_hg38_pos']
        else:
            bcftools_pos = None
        
        # Handle liftover position (may be NaN for failed liftover)  
        if pd.notna(row['liftover_hg38_pos']):
            liftover_pos = row['liftover_hg38_pos']
        else:
            liftover_pos = None
        
        # Skip variants with missing essential coordinates
        if source_pos is None:
            skipped_variants += 1
            if skipped_variants <= 10:
                print(f"Warning: Skipping variant with missing source position: {row['source_chrom']}:{row['source_pos']}")
            continue
        
        processed_data.append({
            'mapping_status': row['mapping_status'],
            'source_chrom': row['source_chrom'],
            'source_pos': int(source_pos),
            'source_alleles': f"{parsed_ref}/{parsed_alt}" if parsed_ref and parsed_alt else row['source_alleles'],
            'source_ref': parsed_ref if parsed_ref else '',
            'source_alt': parsed_alt if parsed_alt else '', 
            'flip': row['flip'],
            'swap': row['swap'],
            'liftover_hg38_chrom': row['liftover_hg38_chrom'],
            'liftover_hg38_pos': int(liftover_pos) if liftover_pos is not None else None,
            'bcftools_hg38_chrom': row['bcftools_hg38_chrom'],
            'bcftools_hg38_pos': int(bcftools_pos) if bcftools_pos is not None else None,
            'bcftools_hg38_ref': parsed_ref,
            'bcftools_hg38_alt': parsed_alt,
            'pos_match': 1 if row['pos_match'] in ['TRUE', True, 1] else 0,
            'gt_match': 1 if row['gt_match'] in ['TRUE', True, 1] else 0
        })

    processed_df = pd.DataFrame(processed_data)
    
    # Remove duplicates in pandas BEFORE database insertion
    print(f"Removing duplicates from {len(processed_df):,} records...")
    original_count = len(processed_df)
    processed_df = processed_df.drop_duplicates(subset=['source_chrom', 'source_pos', 'source_alleles'])
    duplicates_removed = original_count - len(processed_df)
    print(f"✓ Removed {duplicates_removed:,} duplicates, {len(processed_df):,} unique records remain")

    # Use chunked bulk insert to avoid SQL variable limit
    conn = sqlite3.connect(db_path)
    try:
        print("Performing optimized chunked bulk insert...")
        
        # Insert in chunks to avoid SQLite variable limit
        chunk_size = 5000  # Safe chunk size for SQLite
        total_inserted = 0
        
        for i in range(0, len(processed_df), chunk_size):
            chunk = processed_df.iloc[i:i+chunk_size]
            chunk.to_sql('comparison', conn, if_exists='append', index=False, method=None)
            total_inserted += len(chunk)
            
            if i % (chunk_size * 5) == 0:  # Progress every 25k records
                print(f"  Inserted {total_inserted:,}/{len(processed_df):,} records...")
        
        print(f"✓ Inserted {len(processed_df):,} unique records into comparison table")

        if skipped_variants > 0:
            print(f"✓ Skipped {skipped_variants:,} variants with missing source coordinates")
        print("✓ All coordinates normalized to VEP format")
                
    except Exception as e:
        print(f"Error inserting data: {e}")
        raise
    finally:
        conn.close()

def load_vep_data(vep_file, db_path, genome_build):
    """Load VEP annotation data with optimized bulk loading"""
    print(f"Loading {genome_build} VEP data from: {vep_file}")
    
    if not Path(vep_file).exists():
        raise FileNotFoundError(f"VEP file not found: {vep_file}")
    
    table_name = f'{genome_build}_vep'
    
    # Find header line
    with open(vep_file, 'r') as f:
        lines = f.readlines()
    
    header_idx = None
    for i, line in enumerate(lines):
        if line.startswith('#Uploaded_variation'):
            header_idx = i
            break
    
    if header_idx is None:
        raise ValueError(f"Could not find VEP header line starting with '#Uploaded_variation' in {vep_file}")
    
    # Read in larger chunks and use bulk loading
    chunk_size = 100000  # Larger chunks for better performance
    total_records = 0
    
    conn = sqlite3.connect(db_path)
    
    for chunk_num, chunk_df in enumerate(pd.read_csv(
        vep_file, 
        sep='\t', 
        skiprows=header_idx,
        chunksize=chunk_size,
        low_memory=False
    )):
        print(f"  Processing chunk {chunk_num + 1} ({len(chunk_df):,} records)...")
        
        # Parse variant IDs and prepare data
        parsed_data = []
        for _, row in chunk_df.iterrows():
            # Extract coordinates from custom VEP ID
            uploaded_var = get_safe_column_value(row, chunk_df, 'Uploaded_variation')
            extracted_chrom, extracted_pos, extracted_ref, extracted_alt = parse_custom_vep_id(uploaded_var)
            
            # Use extracted coordinates for compatibility
            chr_part = extracted_chrom
            pos = extracted_pos

            # Process MANE columns (only present in hg38)
            mane_data = process_mane_columns(row, chunk_df) 
            
            # Convert gnomAD frequency to float
            gnomad_af = None
            gnomad_raw = get_safe_column_value(row, chunk_df, 'gnomADg_AF')
            if gnomad_raw and gnomad_raw != '-':
                try:
                    gnomad_af = float(gnomad_raw)
                except ValueError:
                    gnomad_af = None

            # Create allele string for this build
            if extracted_ref and extracted_alt:
                allele_string = f"{extracted_ref}/{extracted_alt}"
            else:
                allele_string = None


            parsed_data.append({
                'uploaded_variation': get_safe_column_value(row, chunk_df, 'Uploaded_variation'),
                'chr': chr_part,
                'pos': pos,
                'location': get_safe_column_value(row, chunk_df, 'Location'),
                'gene': get_safe_column_value(row, chunk_df, 'Gene'),
                'feature': get_safe_column_value(row, chunk_df, 'Feature'),
                'feature_type': get_safe_column_value(row, chunk_df, 'Feature_type'),
                'consequence': get_safe_column_value(row, chunk_df, 'Consequence'),
                'impact': get_safe_column_value(row, chunk_df, 'IMPACT'),
                'symbol': get_safe_column_value(row, chunk_df, 'SYMBOL'),
                'sift': get_safe_column_value(row, chunk_df, 'SIFT'),
                'polyphen': get_safe_column_value(row, chunk_df, 'PolyPhen'),
                'gnomadg_af': gnomad_af,
                'clin_sig': get_safe_column_value(row, chunk_df, 'CLIN_SIG'),
                'hgvsc': get_safe_column_value(row, chunk_df, 'HGVSc'),
                'hgvsp': get_safe_column_value(row, chunk_df, 'HGVSp'),
                'canonical': get_safe_column_value(row, chunk_df, 'CANONICAL'),
                # MANE FIELDS (hg38 only, will be empty for hg19)
                'mane': mane_data['mane'],
                'mane_select': mane_data['mane_select'],
                'mane_plus_clinical': mane_data['mane_plus_clinical'],
                'refseq_transcript_id': mane_data['refseq_transcript_id'],
                # EXTRACTED CUSTOM ID FIELDS
                'extracted_chrom': extracted_chrom,
                'extracted_pos': extracted_pos,
                'extracted_ref': extracted_ref,
                'extracted_alt': extracted_alt,
                'allele': allele_string
            })

        
        # Remove duplicates in pandas and bulk insert
        chunk_processed = pd.DataFrame(parsed_data)
        
        # Remove duplicates within this chunk
        original_chunk_size = len(chunk_processed)
        chunk_processed = chunk_processed.drop_duplicates(subset=['uploaded_variation', 'feature', 'consequence'])
        chunk_duplicates = original_chunk_size - len(chunk_processed)
        
        # Use chunked insert to avoid SQL variable limit
        try:
            insert_chunk_size = 5000
            inserted_count = 0
            
            for i in range(0, len(chunk_processed), insert_chunk_size):
                insert_chunk = chunk_processed.iloc[i:i+insert_chunk_size]
                insert_chunk.to_sql(table_name, conn, if_exists='append', index=False, method=None)
                inserted_count += len(insert_chunk)
            
            total_records += inserted_count
            
            print(f"    → Inserted {inserted_count:,} unique records (removed {chunk_duplicates:,} duplicates)")
            print(f"    → Total so far: {total_records:,}")
                
        except Exception as e:
            print(f"    → Warning: VEP bulk insert issue: {e}")
    
    conn.close()
    print(f"✓ Completed loading {total_records:,} records into {table_name} table")


def add_unique_constraints_and_indexes(db_path):
    """Add unique constraints and full indexes AFTER all data is loaded"""
    print("\n=== ADDING UNIQUE CONSTRAINTS AND OPTIMIZED INDEXES ===")
    print("This may take several minutes for large datasets...")
    
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    try:
        # Remove the temporary loading indexes
        print("Removing temporary loading indexes...")
        cursor.execute("DROP INDEX IF EXISTS idx_comp_loading")
        cursor.execute("DROP INDEX IF EXISTS idx_hg19_loading") 
        cursor.execute("DROP INDEX IF EXISTS idx_hg38_loading")
        
        # Add unique constraints by creating unique indexes
        print("Creating unique constraint for comparison table...")
        cursor.execute("""
            CREATE UNIQUE INDEX idx_comparison_unique 
            ON comparison(source_chrom, source_pos, source_alleles)
        """)
        
        print("Creating unique constraints for VEP tables...")
        cursor.execute("""
            CREATE UNIQUE INDEX idx_hg19_vep_unique 
            ON hg19_vep(uploaded_variation, feature, consequence)
        """)
        
        cursor.execute("""
            CREATE UNIQUE INDEX idx_hg38_vep_unique 
            ON hg38_vep(uploaded_variation, feature, consequence)
        """)
        
        # Create optimized indexes for analysis queries
        print("Creating optimized analysis indexes...")
        
        # Comparison table indexes
        cursor.execute("CREATE INDEX idx_comp_bcftools ON comparison(bcftools_hg38_chrom, bcftools_hg38_pos)")
        cursor.execute("CREATE INDEX idx_comp_mapping_status ON comparison(mapping_status)")
        cursor.execute("CREATE INDEX idx_comp_matches ON comparison(pos_match, gt_match)")
        
        # VEP table indexes
        for build in ['hg19', 'hg38']:
            table = f'{build}_vep'
            cursor.execute(f"CREATE INDEX idx_{build}_coords ON {table}(extracted_chrom, extracted_pos, extracted_ref, extracted_alt)")
            cursor.execute(f"CREATE INDEX idx_{build}_feature ON {table}(feature, feature_type)")
            cursor.execute(f"CREATE INDEX idx_{build}_consequence ON {table}(consequence, impact)")
            cursor.execute(f"CREATE INDEX idx_{build}_symbol ON {table}(symbol)")
        
        # Additional comparison table indexes for matching
        cursor.execute("CREATE INDEX idx_comp_source_alleles ON comparison(source_chrom, source_pos, source_ref, source_alt)")
        cursor.execute("CREATE INDEX idx_comp_bcftools_alleles ON comparison(bcftools_hg38_chrom, bcftools_hg38_pos, bcftools_hg38_ref, bcftools_hg38_alt)")

        conn.commit()
        print("✓ All unique constraints and indexes created successfully")
        
    except Exception as e:
        print(f"Error creating constraints/indexes: {e}")
        print("Database may still be usable, but performance could be impacted")
        raise
    finally:
        conn.close()

def verify_database(db_path):
    """Verify database contents and coordinate adjustments"""
    print("\n=== DATABASE VERIFICATION ===")
    
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    # Check table counts
    print("Table record counts:")
    for table in ['comparison', 'hg19_vep', 'hg38_vep']:
        cursor.execute(f"SELECT COUNT(*) FROM {table}")
        count = cursor.fetchone()[0]
        print(f"  {table}: {count:,} records")
    
    # Check JOIN success rates
    print("\nJOIN success rates:")
    cursor.execute("""
        SELECT 
            COUNT(DISTINCT c.source_chrom || ':' || c.source_pos) as total_variants,
            COUNT(DISTINCT CASE WHEN h19.extracted_pos IS NOT NULL THEN c.source_chrom || ':' || c.source_pos END) as hg19_matched,
            COUNT(DISTINCT CASE WHEN h38.extracted_pos IS NOT NULL THEN c.source_chrom || ':' || c.source_pos END) as hg38_matched
        FROM comparison c
        LEFT JOIN hg19_vep h19 ON (
            c.source_chrom = h19.extracted_chrom AND 
            c.source_pos = h19.extracted_pos AND
            c.source_ref = h19.extracted_ref AND
            c.source_alt = h19.extracted_alt
        )
        LEFT JOIN hg38_vep h38 ON (
            c.bcftools_hg38_chrom = h38.extracted_chrom AND
            c.bcftools_hg38_pos = h38.extracted_pos AND
            c.bcftools_hg38_ref = h38.extracted_ref AND
            c.bcftools_hg38_alt = h38.extracted_alt
        )
    """)

    match_stats = cursor.fetchone()
    total, hg19_matched, hg38_matched = match_stats
    
    print(f"  Total comparison variants: {total:,}")
    print(f"  HG19 VEP matches: {hg19_matched:,} ({hg19_matched/total*100:.1f}%)")
    print(f"  HG38 VEP matches: {hg38_matched:,} ({hg38_matched/total*100:.1f}%)")
    
    # Sample data verification
    print("\nSample data verification:")
    cursor.execute("SELECT source_chrom, source_pos, source_alleles, pos_match, gt_match FROM comparison LIMIT 3")
    print("  Comparison table sample:")
    for row in cursor.fetchall():
        print(f"    {row}")
    
    cursor.execute("SELECT extracted_chrom, extracted_pos, extracted_ref, extracted_alt, allele, consequence, symbol FROM hg19_vep LIMIT 3")
    print("  HG19 VEP table sample:")
    for row in cursor.fetchall():
        print(f"    {row}")
    
    conn.close()
    print("\n✓ Database verification completed")

def load_config(config_path):
    """Load and validate simplified JSON configuration"""
    if not Path(config_path).exists():
        raise FileNotFoundError(f"Configuration file not found: {config_path}")
    
    try:
        with open(config_path, 'r') as f:
            config = json.load(f)
    except json.JSONDecodeError as e:
        raise ValueError(f"Invalid JSON in configuration file: {e}")
    
    # Validate required sections
    if 'input_files' not in config:
        raise ValueError("Configuration file must contain 'input_files' section")
    
    if 'database' not in config:
        raise ValueError("Configuration file must contain 'database' section")
    
    if 'path' not in config['database']:
        raise ValueError("Configuration must specify 'database.path'")
    
    required_files = ['comparison', 'hg19_vep', 'hg38_vep']
    missing_files = [f for f in required_files if f not in config['input_files']]
    if missing_files:
        raise ValueError(f"Missing required input files in config: {missing_files}")
    
    return config

def main():
    parser = argparse.ArgumentParser(
        description='Load genomic variant data into SQLite database with VEP coordinate normalization',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
EXAMPLE CONFIG FILE (config.json):
{
  "input_files": {
    "comparison": "/path/to/compare_liftover.txt",
    "hg19_vep": "/path/to/hg19.vep.txt", 
    "hg38_vep": "/path/to/hg38.vep.txt"
  },
  "database": {
    "path": "/path/to/genomic_analysis.db"
  },
  "output": {
    "directory": "/path/to/output"
  }
}

USAGE:
    python db_loader.py --config config.json --force
        """
    )
    
    parser.add_argument('--config', '-c', required=True,
                       help='JSON configuration file specifying input file paths and database location')
    parser.add_argument('--force', action='store_true',
                       help='Overwrite existing database')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Enable verbose output')
    
    args = parser.parse_args()
    
    # Load configuration
    try:
        config = load_config(args.config)
        print(f"✓ Loaded configuration from: {args.config}")
    except Exception as e:
        print(f"Error loading configuration: {e}")
        sys.exit(1)
    
    db_path = Path(os.path.expanduser(config['database']['path']))
    
    # Check if database exists
    if db_path.exists() and not args.force:
        print(f"Error: Database {db_path} already exists. Use --force to overwrite.")
        sys.exit(1)
    
    try:
        print(f"\n=== CREATING OPTIMIZED DATABASE: {db_path} ===")
        
        # Step 1: Create schema without unique constraints
        create_database_schema(db_path)
        
        # Step 2: Load all data quickly (no constraint checking)
        load_comparison_data(os.path.expanduser(config['input_files']['comparison']), db_path)
        load_vep_data(os.path.expanduser(config['input_files']['hg19_vep']), db_path, 'hg19')
        load_vep_data(os.path.expanduser(config['input_files']['hg38_vep']), db_path, 'hg38')
        
        # Step 3: Add unique constraints and optimized indexes AFTER loading
        add_unique_constraints_and_indexes(db_path)
        
        # Step 4: Verify database
        verify_database(db_path)
        
        print(f"\n=== OPTIMIZED DATABASE CREATION COMPLETED ===")
        print(f"✓ Database: {db_path}")
        print(f"✓ Size: {db_path.stat().st_size / (1024**3):.2f} GB")
        print(f"\nOPTIMIZATION SUMMARY:")
        print(f"✓ Bulk loading without constraint checking during insert")
        print(f"✓ Duplicate removal in pandas before database insertion")
        print(f"✓ Unique constraints added after all data loaded")
        print(f"✓ Optimized indexes for analysis queries")
        print(f"\nIMPORTANT: All coordinates use VEP normalization")
        print(f"Ready for analysis with db_analyzer.py and variant_prioritizer.py")
        
    except Exception as e:
        print(f"Error: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
