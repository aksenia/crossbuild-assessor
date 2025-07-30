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

import sqlite3
import pandas as pd
import argparse
import json
import re
from pathlib import Path
import sys

def create_database_schema(db_path):
    """Create database tables with proper indexes for efficient analysis"""
    print("Creating database schema...")
    
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    # Drop tables if they exist
    cursor.execute("DROP TABLE IF EXISTS comparison")
    cursor.execute("DROP TABLE IF EXISTS hg19_vep")
    cursor.execute("DROP TABLE IF EXISTS hg38_vep")
    
    # Create comparison table
    cursor.execute("""
        CREATE TABLE comparison (
            id INTEGER PRIMARY KEY,
            mapping_status TEXT,
            source_chrom TEXT,
            source_pos INTEGER,
            source_alleles TEXT,
            flip TEXT,
            swap TEXT,
            liftover_hg38_chrom TEXT,
            liftover_hg38_pos INTEGER,
            bcftools_hg38_chrom TEXT,
            bcftools_hg38_pos INTEGER,
            bcftools_hg38_ref TEXT,
            bcftools_hg38_alt TEXT,
            pos_match BOOLEAN,
            gt_match BOOLEAN,
            UNIQUE(source_chrom, source_pos, source_alleles)
        )
    """)
    
    # Create VEP tables (identical structure for both builds)
    vep_schema = """
        CREATE TABLE {} (
            id INTEGER PRIMARY KEY,
            uploaded_variation TEXT,
            chr TEXT,
            pos INTEGER,
            pos_original INTEGER,
            ref_allele TEXT,
            alt_allele TEXT,
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
            UNIQUE(uploaded_variation, feature, consequence)
        )
    """
    
    cursor.execute(vep_schema.format('hg19_vep'))
    cursor.execute(vep_schema.format('hg38_vep'))
    
    # Create indexes for efficient joins and queries
    print("Creating database indexes...")
    
    # Comparison table indexes
    cursor.execute("CREATE INDEX idx_comp_source ON comparison(source_chrom, source_pos, source_alleles)")
    cursor.execute("CREATE INDEX idx_comp_bcftools ON comparison(bcftools_hg38_chrom, bcftools_hg38_pos)")
    cursor.execute("CREATE INDEX idx_comp_mapping_status ON comparison(mapping_status)")
    cursor.execute("CREATE INDEX idx_comp_matches ON comparison(pos_match, gt_match)")
    
    # VEP table indexes
    for build in ['hg19', 'hg38']:
        table = f'{build}_vep'
        cursor.execute(f"CREATE INDEX idx_{build}_coords ON {table}(chr, pos, ref_allele, alt_allele)")
        cursor.execute(f"CREATE INDEX idx_{build}_coords_orig ON {table}(chr, pos_original)")
        cursor.execute(f"CREATE INDEX idx_{build}_feature ON {table}(feature, feature_type)")
        cursor.execute(f"CREATE INDEX idx_{build}_consequence ON {table}(consequence, impact)")
        cursor.execute(f"CREATE INDEX idx_{build}_symbol ON {table}(symbol)")
    
    conn.commit()
    conn.close()
    print("Database schema created successfully")

def parse_comparison_alleles(alleles_str):
    """
    Parse comparison alleles and convert to VEP format
    
    Args:
        alleles_str: Comma-separated alleles (e.g., "A,AG" or "A,G")
    
    Returns:
        tuple: (vep_ref, vep_alt) in VEP notation
    
    Examples:
        "A,G" → ("A", "G")           # SNV
        "A,AG" → ("-", "G")          # Insertion
        "AG,A" → ("G", "-")          # Deletion
    """
    if pd.isna(alleles_str):
        return None, None
    
    parts = str(alleles_str).split(',')
    if len(parts) != 2:
        return None, None
    
    ref, alt = parts[0].strip(), parts[1].strip()
    
    # Convert to VEP notation for indels
    if len(ref) < len(alt) and alt.startswith(ref):
        # Insertion: "A,AG" → "-", "G" (VEP removes the common prefix)
        vep_ref = "-"
        vep_alt = alt[len(ref):]
    elif len(ref) > len(alt):
        # Deletion: handle various deletion formats
        if len(alt) == 0:
            # Pure deletion: "AGC," → "AGC", "-"
            vep_ref = ref
            vep_alt = "-"
        elif ref.startswith(alt):
            # Deletion with common prefix: "AGC,A" → "GC", "-"
            vep_ref = ref[len(alt):]
            vep_alt = "-"
        else:
            # Complex deletion
            vep_ref = ref
            vep_alt = alt
    else:
        # SNV: keep as-is
        vep_ref = ref
        vep_alt = alt
    
    return vep_ref, vep_alt

def adjust_comparison_coordinates(pos, ref, alt):
    """
    Adjust comparison coordinates to match VEP normalization
    
    VEP coordinate normalization:
    - SNVs: no change (VEP pos = original pos)
    - Indels: +1 (VEP pos = original pos + 1)
    
    Args:
        pos: Original position
        ref: Reference allele in VEP format
        alt: Alternative allele in VEP format
    
    Returns:
        int: VEP-normalized position
    """
    if ref == "-" or alt == "-":
        # Indel: VEP shifts position +1
        return pos + 1
    else:
        # SNV: VEP keeps same position
        return pos

def load_comparison_data(comparison_file, db_path):
    """Load comparison file into database with coordinate adjustment to VEP format"""
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
    
    # Convert coordinates and alleles to VEP format
    coordinate_adjustments = 0
    processed_data = []
    skipped_variants = 0
    
    for _, row in df.iterrows():
        # Parse alleles and convert to VEP format
        vep_ref, vep_alt = parse_comparison_alleles(row['source_alleles'])
        
        # Adjust coordinates to VEP normalization - handle NaN values
        if vep_ref is not None and pd.notna(row['source_pos']):
            vep_source_pos = adjust_comparison_coordinates(row['source_pos'], vep_ref, vep_alt)
            
            # Handle bcftools position (may be NaN for failed liftover)
            if pd.notna(row['bcftools_hg38_pos']):
                vep_bcftools_pos = adjust_comparison_coordinates(row['bcftools_hg38_pos'], vep_ref, vep_alt)
            else:
                vep_bcftools_pos = None
            
            # Handle liftover position (may be NaN for failed liftover)  
            if pd.notna(row['liftover_hg38_pos']):
                vep_liftover_pos = adjust_comparison_coordinates(row['liftover_hg38_pos'], vep_ref, vep_alt)
            else:
                vep_liftover_pos = None
            
            if vep_source_pos != row['source_pos']:
                coordinate_adjustments += 1
        else:
            # Keep original if parsing failed or source_pos is NaN
            vep_source_pos = row['source_pos'] if pd.notna(row['source_pos']) else None
            vep_bcftools_pos = row['bcftools_hg38_pos'] if pd.notna(row['bcftools_hg38_pos']) else None
            vep_liftover_pos = row['liftover_hg38_pos'] if pd.notna(row['liftover_hg38_pos']) else None
            vep_ref = ""
            vep_alt = ""
        
        # Skip variants with missing essential coordinates
        if vep_source_pos is None:
            skipped_variants += 1
            if skipped_variants <= 10:  # Only show first 10 warnings to avoid spam
                print(f"Warning: Skipping variant with missing source position: {row['source_chrom']}:{row['source_pos']}")
            continue
        
        processed_data.append({
            'mapping_status': row['mapping_status'],
            'source_chrom': row['source_chrom'],
            'source_pos': int(vep_source_pos),  # Safe now - we checked for None above
            'source_alleles': f"{vep_ref}/{vep_alt}" if vep_ref and vep_alt else row['source_alleles'],
            'flip': row['flip'],
            'swap': row['swap'],
            'liftover_hg38_chrom': row['liftover_hg38_chrom'],
            'liftover_hg38_pos': int(vep_liftover_pos) if vep_liftover_pos is not None else None,
            'bcftools_hg38_chrom': row['bcftools_hg38_chrom'],
            'bcftools_hg38_pos': int(vep_bcftools_pos) if vep_bcftools_pos is not None else None,
            'bcftools_hg38_ref': vep_ref,
            'bcftools_hg38_alt': vep_alt,
            'pos_match': 1 if row['pos_match'] in ['TRUE', True, 1] else 0,
            'gt_match': 1 if row['gt_match'] in ['TRUE', True, 1] else 0
        })
    
        processed_df = pd.DataFrame(processed_data)

        # Insert into database - HANDLE DUPLICATES
        conn = sqlite3.connect(db_path)
        try:
            cursor = conn.cursor()
            
            for _, row in processed_df.iterrows():
                cursor.execute("""
                    INSERT OR IGNORE INTO comparison (
                        mapping_status, source_chrom, source_pos, source_alleles,
                        flip, swap, liftover_hg38_chrom, liftover_hg38_pos,
                        bcftools_hg38_chrom, bcftools_hg38_pos, bcftools_hg38_ref,
                        bcftools_hg38_alt, pos_match, gt_match
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """, (
                    row['mapping_status'], row['source_chrom'], row['source_pos'],
                    row['source_alleles'], row['flip'], row['swap'],
                    row['liftover_hg38_chrom'], row['liftover_hg38_pos'],
                    row['bcftools_hg38_chrom'], row['bcftools_hg38_pos'],
                    row['bcftools_hg38_ref'], row['bcftools_hg38_alt'],
                    row['pos_match'], row['gt_match']
                ))
            
            conn.commit()
            unique_inserted = cursor.rowcount
            total_attempted = len(processed_df)
            duplicates_skipped = total_attempted - unique_inserted
            
            print(f"✓ Inserted {unique_inserted:,} unique records into comparison table")
            if duplicates_skipped > 0:
                print(f"✓ Skipped {duplicates_skipped:,} duplicate variants")
            if skipped_variants > 0:
                print(f"✓ Skipped {skipped_variants:,} variants with missing source coordinates")
            print(f"✓ Applied coordinate/allele adjustments to {coordinate_adjustments:,} variants")
            print("✓ All coordinates normalized to VEP format")
                
        except sqlite3.IntegrityError as e:
            print(f"Error inserting data: {e}")
        finally:
            conn.close()

def parse_vep_variant_id(variant_id):
    """
    Parse VEP variant ID into components
    
    Args:
        variant_id: VEP variant identifier (e.g., "1_69134_A/G" or "chr1_69134_A/G")
    
    Returns:
        tuple: (chromosome, position, ref_allele, alt_allele)
    """
    if pd.isna(variant_id):
        print(f"DEBUG: variant_id is NaN")
        return None, None, None, None
    
    print(f"DEBUG: Processing variant_id: '{variant_id}' (type: {type(variant_id)})")
    # Pattern: chr_pos_ref/alt
    pattern = r'^(.+?)_(\d+)_([^/]+)/(.+)$'
    match = re.match(pattern, str(variant_id))
    
    if match:
        chr_part = match.group(1).replace('chr', '')  # Remove chr prefix if present
        pos = int(match.group(2))
        ref = match.group(3)
        alt = match.group(4)
        print(f"DEBUG: Successfully parsed - chr:{chr_part}, pos:{pos}, ref:{ref}, alt:{alt}")
        return chr_part, pos, ref, alt
    else:
        print(f"DEBUG: REGEX FAILED for variant_id: '{variant_id}'")
        return None, None, None, None

def load_vep_data(vep_file, db_path, genome_build):
    """Load VEP annotation data into database"""
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
    
    # Read VEP file in chunks for memory efficiency
    chunk_size = 50000
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
        
        # Parse variant IDs (coordinates already in VEP format)
        parsed_data = []
        for _, row in chunk_df.iterrows():
            chr_part, pos, ref, alt = parse_vep_variant_id(row['#Uploaded_variation'])
            
            # Convert gnomAD frequency to float
            gnomad_af = None
            if 'gnomADg_AF' in row and pd.notna(row['gnomADg_AF']) and row['gnomADg_AF'] != '-':
                try:
                    gnomad_af = float(row['gnomADg_AF'])
                except ValueError:
                    gnomad_af = None
            
            parsed_data.append({
                'uploaded_variation': row['#Uploaded_variation'],
                'chr': chr_part,
                'pos': pos,
                'pos_original': pos,  # Same as pos since VEP coordinates are kept as-is
                'ref_allele': ref,
                'alt_allele': alt,
                'location': row.get('Location', ''),
                'allele': row.get('Allele', ''),
                'gene': row.get('Gene', ''),
                'feature': row.get('Feature', ''),
                'feature_type': row.get('Feature_type', ''),
                'consequence': row.get('Consequence', ''),
                'impact': row.get('IMPACT', ''),
                'symbol': row.get('SYMBOL', ''),
                'sift': row.get('SIFT', ''),
                'polyphen': row.get('PolyPhen', ''),
                'gnomadg_af': gnomad_af,
                'clin_sig': row.get('CLIN_SIG', '')
            })
        
        # Insert chunk into database - HANDLE DUPLICATES
        chunk_processed = pd.DataFrame(parsed_data)

        try:
            cursor = conn.cursor()
            inserted_count = 0
            
            for _, row in chunk_processed.iterrows():
                cursor.execute(f"""
                    INSERT OR IGNORE INTO {table_name} (
                        uploaded_variation, chr, pos, pos_original, ref_allele, alt_allele,
                        location, allele, gene, feature, feature_type, consequence,
                        impact, symbol, sift, polyphen, gnomadg_af, clin_sig
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """, (
                    row['uploaded_variation'], row['chr'], row['pos'], row['pos_original'],
                    row['ref_allele'], row['alt_allele'], row['location'], row['allele'],
                    row['gene'], row['feature'], row['feature_type'], row['consequence'],
                    row['impact'], row['symbol'], row['sift'], row['polyphen'],
                    row['gnomadg_af'], row['clin_sig']
                ))
                if cursor.rowcount > 0:
                    inserted_count += 1
            
            conn.commit()
            total_records += inserted_count
            duplicates_skipped = len(chunk_processed) - inserted_count
            
            print(f"    → Inserted {inserted_count:,} unique records (total: {total_records:,})")
            if duplicates_skipped > 0:
                print(f"    → Skipped {duplicates_skipped:,} duplicates in this chunk")
                
        except sqlite3.IntegrityError as e:
            print(f"    → Warning: VEP duplicate handling issue: {e}")
    
    conn.close()
    print(f"✓ Completed loading {total_records:,} records into {table_name} table")

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
            COUNT(DISTINCT CASE WHEN h19.pos IS NOT NULL THEN c.source_chrom || ':' || c.source_pos END) as hg19_matched,
            COUNT(DISTINCT CASE WHEN h38.pos IS NOT NULL THEN c.source_chrom || ':' || c.source_pos END) as hg38_matched
        FROM comparison c
        LEFT JOIN hg19_vep h19 ON (c.source_chrom = h19.chr AND c.source_pos = h19.pos)
        LEFT JOIN hg38_vep h38 ON (c.bcftools_hg38_chrom = h38.chr AND c.bcftools_hg38_pos = h38.pos)
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
    
    cursor.execute("SELECT chr, pos, ref_allele, alt_allele, consequence, symbol FROM hg19_vep LIMIT 3")
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
    
    db_path = Path(config['database']['path'])
    
    # Check if database exists
    if db_path.exists() and not args.force:
        print(f"Error: Database {db_path} already exists. Use --force to overwrite.")
        sys.exit(1)
    
    try:
        print(f"\n=== CREATING DATABASE: {db_path} ===")
        
        # Create database schema
        create_database_schema(db_path)
        
        # Load data from configured files
        load_comparison_data(config['input_files']['comparison'], db_path)
        load_vep_data(config['input_files']['hg19_vep'], db_path, 'hg19')
        load_vep_data(config['input_files']['hg38_vep'], db_path, 'hg38')
        
        # Verify database
        verify_database(db_path)
        
        print(f"\n=== DATABASE CREATION COMPLETED ===")
        print(f"✓ Database: {db_path}")
        print(f"✓ Size: {db_path.stat().st_size / (1024**3):.2f} GB")
        print(f"\nIMPORTANT: All coordinates use VEP normalization")
        print(f"  - SNVs: same position as input")
        print(f"  - Indels: input position + 1")
        print(f"\nReady for analysis with db_analyzer.py and variant_prioritizer.py")
        
    except Exception as e:
        print(f"Error: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
