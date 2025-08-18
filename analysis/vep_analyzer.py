"""
VEP Analyzer

Core VEP analysis engine for processing variant annotations and identifying
transcript-level discordances between genome builds.
"""

import pandas as pd
from utils.transcript_utils import extract_genotype_from_alleles
#from utils.data_utils import clean_string
from utils.clinical_utils import (
    is_pathogenic_clinical_significance,
    is_benign_clinical_significance, 
    parse_sift_prediction,
    parse_polyphen_prediction,
    normalize_clinical_significance
)
from config.constants import VEP_CONSEQUENCE_IMPACT
from utils.hgvs_utils import analyze_priority_transcript_hgvs

class VEPAnalyzer:
    """Analyzes VEP annotations to identify discordances between genome builds"""
    
    def __init__(self):
        """Initialize VEP analyzer"""
        pass
    
    def analyze_all_variants(self, conn):
        """
        Analyze all variants with comprehensive VEP processing
        
        Args:
            conn: Database connection
            
        Returns:
            pd.DataFrame: VEP analysis results for all variants
        """
        print("Calculating fresh VEP analysis...")
        
        # STEP 1: Get all unique variants (both concordant and discordant) - MEMORY EFFICIENT
        print("Step 1: Identifying all variants for VEP analysis...")
        variant_query = """
        SELECT DISTINCT 
            c.source_chrom,
            c.source_pos,
            c.source_ref,
            c.source_alt,
            c.bcftools_hg38_chrom,
            c.bcftools_hg38_pos,
            c.bcftools_hg38_ref, 
            c.bcftools_hg38_alt,
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
            chunk_analyses = self._process_variant_chunk(conn, chunk_variants)
            all_vep_analyses.extend(chunk_analyses)
            
            # Memory cleanup
            del chunk_variants
        
        print(f"Completed VEP analysis of {len(all_vep_analyses):,} variants")
        
        if len(all_vep_analyses) == 0:
            print("No variants found for analysis")
            return pd.DataFrame()
        
        return pd.DataFrame(all_vep_analyses)
    
    def _process_variant_chunk(self, conn, chunk_variants):
        """Process a chunk of variants with comprehensive VEP analysis"""
        vep_analyses = []
        
        for _, variant_row in chunk_variants.iterrows():

            chrom = variant_row['source_chrom']
            pos = variant_row['source_pos']
            hg38_chrom = variant_row['bcftools_hg38_chrom']
            hg38_pos = variant_row['bcftools_hg38_pos']
            
            # Get VEP annotations for this specific variant - TARGETED QUERIES
            hg19_query = """
            SELECT feature_type, consequence, impact, symbol, feature, sift, polyphen, gnomadg_af, clin_sig, hgvsc, hgvsp, 
                extracted_chrom, extracted_pos, extracted_ref, extracted_alt, allele,
                mane, mane_select, mane_plus_clinical, refseq_transcript_id,
            CASE WHEN CANONICAL = 'YES' THEN 1 ELSE 0 END as is_canonical
            FROM hg19_vep 
            WHERE extracted_chrom = ? AND extracted_pos = ? AND extracted_ref = ? AND extracted_alt = ?
            """
            
            hg38_query = """
            SELECT feature_type, consequence, impact, symbol, feature, sift, polyphen, gnomadg_af, clin_sig, hgvsc, hgvsp, 
                extracted_chrom, extracted_pos, extracted_ref, extracted_alt, allele,
                mane, mane_select, mane_plus_clinical, refseq_transcript_id,
            CASE WHEN CANONICAL = 'YES' THEN 1 ELSE 0 END as is_canonical
            FROM hg38_vep 
            WHERE extracted_chrom = ? AND extracted_pos = ? AND extracted_ref = ? AND extracted_alt = ?
            """
            
            hg19_annotations = pd.read_sql_query(hg19_query, conn, params=[chrom, pos, variant_row['source_ref'], variant_row['source_alt']])
            hg38_annotations = pd.read_sql_query(hg38_query, conn, params=[hg38_chrom, hg38_pos, variant_row['bcftools_hg38_ref'], variant_row['bcftools_hg38_alt']])
            
            # Apply comprehensive VEP analysis to this variant
            vep_analysis = self._analyze_single_variant(variant_row, hg19_annotations, hg38_annotations)
            
            if vep_analysis:
                vep_analyses.append(vep_analysis)
        
        return vep_analyses
    
    def _analyze_single_variant(self, variant_row, hg19_annotations, hg38_annotations): 
        """Apply comprehensive VEP analysis to a single variant (cached part) - CLEAN SET-BASED APPROACH"""
        from utils.data_utils import clean_string
        chrom = variant_row['source_chrom']
        pos = variant_row['source_pos']
        
        # Extract genotypes
        ref_allele, alt_allele = extract_genotype_from_alleles(variant_row['source_alleles'])
        
        # Focus on transcript annotations only for transcript analysis
        hg19_transcripts_df = hg19_annotations[hg19_annotations['feature_type'] == 'Transcript']
        hg38_transcripts_df = hg38_annotations[hg38_annotations['feature_type'] == 'Transcript']
        
        # Initialize analysis variables
        same_transcript_consequence_changes = 0
        problematic_transcripts_hg19 = []
        problematic_transcripts_hg38 = []
        
        if len(hg19_transcripts_df) > 0 or len(hg38_transcripts_df) > 0:
            # CLEAN SET-BASED ANALYSIS: Analyze all transcripts and consequences as sets
            from utils.data_utils import clean_string
            
            # Organize transcript data by build
            hg19_transcripts = {}
            hg38_transcripts = {}
            
            for _, row in hg19_transcripts_df.iterrows():
                hg19_transcripts[row['feature']] = {
                    'consequence': row['consequence'],
                    'impact': row['impact'],
                    'symbol': row['symbol'],
                    'feature_id': row['feature']
                }
            
            for _, row in hg38_transcripts_df.iterrows():
                hg38_transcripts[row['feature']] = {
                    'consequence': row['consequence'],
                    'impact': row['impact'],
                    'symbol': row['symbol'],
                    'feature_id': row['feature']
                }
            
            # Step 1: Same transcript ID analysis (for same_transcript_consequence_changes)
            hg19_tx_dict = {}  # full_transcript_id -> data
            hg38_tx_dict = {}  # full_transcript_id -> data
            
            # Use full transcript IDs directly
            for tx_id, data in hg19_transcripts.items():
                if tx_id:
                    hg19_tx_dict[tx_id] = data
            
            for tx_id, data in hg38_transcripts.items():
                if tx_id:
                    hg38_tx_dict[tx_id] = data
            
            # Find same transcript IDs (exact version match) with different consequences
            for tx_id in set(hg19_tx_dict.keys()) & set(hg38_tx_dict.keys()):
                hg19_data = hg19_tx_dict[tx_id]
                hg38_data = hg38_tx_dict[tx_id]
                
                # Use robust string comparison
                consequence_match = clean_string(hg19_data['consequence']) == clean_string(hg38_data['consequence'])
                
                if not consequence_match:
                    same_transcript_consequence_changes += 1
                    problematic_transcripts_hg19.append(f"{hg19_data['feature_id']}({hg19_data['consequence']})")
                    problematic_transcripts_hg38.append(f"{hg38_data['feature_id']}({hg38_data['consequence']})")
            
            # Step 2: Set-based transcript relationship analysis (exact version matching)
            all_hg19_transcript_ids = {tx_id for tx_id in hg19_transcripts.keys() if tx_id}
            all_hg38_transcript_ids = {tx_id for tx_id in hg38_transcripts.keys() if tx_id}
            
            # Transcript relationship (exact version matching)
            if len(all_hg19_transcript_ids) == 0 and len(all_hg38_transcript_ids) == 0:
                transcript_relationship = 'no_transcripts'
            elif all_hg19_transcript_ids == all_hg38_transcript_ids:
                transcript_relationship = 'matched'
            elif len(all_hg19_transcript_ids - all_hg38_transcript_ids) == 0:
                transcript_relationship = 'hg19_subset_of_hg38'
            elif len(all_hg38_transcript_ids - all_hg19_transcript_ids) == 0:
                transcript_relationship = 'hg38_subset_of_hg19'
            elif len(all_hg19_transcript_ids & all_hg38_transcript_ids) == 0:
                transcript_relationship = 'disjoint_transcripts'
            else:
                transcript_relationship = 'partial_overlap_transcripts'

            # Step 3: Set-based consequence relationship analysis
            all_hg19_consequences = set()
            all_hg38_consequences = set()

            # ROBUST consequence set building with comma-splitting
            for data in hg19_transcripts.values():
                consequence_field = data['consequence']
                if pd.notna(consequence_field) and str(consequence_field).strip():
                    # Split comma-separated consequences and clean each one
                    individual_consequences = [c.strip() for c in str(consequence_field).split(',')]
                    for cons in individual_consequences:
                        cleaned_cons = clean_string(cons)
                        if cleaned_cons:  # Only add non-empty cleaned consequences
                            all_hg19_consequences.add(cleaned_cons)

            for data in hg38_transcripts.values():
                consequence_field = data['consequence']
                if pd.notna(consequence_field) and str(consequence_field).strip():
                    # Split comma-separated consequences and clean each one
                    individual_consequences = [c.strip() for c in str(consequence_field).split(',')]
                    for cons in individual_consequences:
                        cleaned_cons = clean_string(cons)
                        if cleaned_cons:  # Only add non-empty cleaned consequences
                            all_hg38_consequences.add(cleaned_cons)
            
            # Calculate relationships step by step
            shared_consequences = all_hg19_consequences & all_hg38_consequences
            unique_hg19 = all_hg19_consequences - all_hg38_consequences  
            unique_hg38 = all_hg38_consequences - all_hg19_consequences

            # Determine relationship 
            if len(all_hg19_consequences) == 0 and len(all_hg38_consequences) == 0:
                consequence_relationship = 'no_consequences'
                unmatched_consequences = 0
            elif len(shared_consequences) > 0 and len(unique_hg19) == 0 and len(unique_hg38) == 0:
                consequence_relationship = 'matched'
                unmatched_consequences = 0
            elif len(shared_consequences) == 0:  # Truly disjoint
                consequence_relationship = 'disjoint_consequences'
                # Score based on significant consequences
                hg19_significant = any(VEP_CONSEQUENCE_IMPACT.get(cons, 'MODIFIER') in {'HIGH', 'MODERATE'} for cons in all_hg19_consequences)
                hg38_significant = any(VEP_CONSEQUENCE_IMPACT.get(cons, 'MODIFIER') in {'HIGH', 'MODERATE'} for cons in all_hg38_consequences)
                unmatched_consequences = 5 if (hg19_significant or hg38_significant) else 0
            elif len(unique_hg19) == 0:  # hg19 subset of hg38
                consequence_relationship = 'hg19_subset_of_hg38'
                unmatched_consequences = 0
            elif len(unique_hg38) == 0:  # hg38 subset of hg19 
                consequence_relationship = 'hg38_subset_of_hg19'
                unmatched_consequences = 0
            else:  # Partial overlap
                consequence_relationship = 'partial_overlap_consequences'
                # Score based on unique significant consequences
                hg19_significant = any(VEP_CONSEQUENCE_IMPACT.get(cons, 'MODIFIER') in {'HIGH', 'MODERATE'} for cons in unique_hg19)
                hg38_significant = any(VEP_CONSEQUENCE_IMPACT.get(cons, 'MODIFIER') in {'HIGH', 'MODERATE'} for cons in unique_hg38)
                unmatched_consequences = 1 if (hg19_significant or hg38_significant) else 0
            
            # Set other analysis variables for compatibility
            same_consequence_different_transcripts = 0  # Not used in new approach
            gene_changes = 0  # Will be calculated later
            impact_changes = 0  # Will be calculated later
            transcript_pairs_analyzed = len(hg19_transcripts) + len(hg38_transcripts)
        
        else:
            # No transcript data available
            same_transcript_consequence_changes = 0
            same_consequence_different_transcripts = 0
            unmatched_consequences = 0
            gene_changes = 0
            impact_changes = 0
            problematic_transcripts_hg19 = []
            problematic_transcripts_hg38 = []
            transcript_pairs_analyzed = 0
            transcript_relationship = 'no_transcripts'
            consequence_relationship = 'no_consequences'
        
        # Get representative VEP information with proper fallbacks
        # Priority: transcript data > any annotation data > empty string
        def get_best_gene(transcripts_df, annotations_df):
            """Get the best available gene symbol with proper fallback"""
            # First try transcript data
            if len(transcripts_df) > 0:
                gene_from_transcript = transcripts_df['symbol'].iloc[0]
                if pd.notna(gene_from_transcript) and gene_from_transcript != '' and gene_from_transcript != '-':
                    return gene_from_transcript
            
            # Fallback to any annotation data
            if len(annotations_df) > 0:
                gene_from_annotation = annotations_df['symbol'].iloc[0]
                if pd.notna(gene_from_annotation) and gene_from_annotation != '' and gene_from_annotation != '-':
                    return gene_from_annotation
            
            # Final fallback
            return ''
        
        def get_best_consequence(transcripts_df, annotations_df):
            """Get the best available consequence with proper fallback"""
            if len(transcripts_df) > 0:
                cons_from_transcript = transcripts_df['consequence'].iloc[0]
                if pd.notna(cons_from_transcript) and cons_from_transcript != '' and cons_from_transcript != '-':
                    return cons_from_transcript
            
            if len(annotations_df) > 0:
                cons_from_annotation = annotations_df['consequence'].iloc[0]
                if pd.notna(cons_from_annotation) and cons_from_annotation != '' and cons_from_annotation != '-':
                    return cons_from_annotation
            
            return ''
        
        def get_best_impact(annotations_df):
            """Get the best available impact"""
            if len(annotations_df) > 0:
                impact = annotations_df['impact'].iloc[0]
                if pd.notna(impact) and impact != '' and impact != '-':
                    return impact
            return ''
        
        # Apply the improved extraction
        hg19_gene = get_best_gene(hg19_transcripts_df, hg19_annotations)
        hg38_gene = get_best_gene(hg38_transcripts_df, hg38_annotations)
        hg19_consequence = get_best_consequence(hg19_transcripts_df, hg19_annotations)
        hg38_consequence = get_best_consequence(hg38_transcripts_df, hg38_annotations)
        hg19_impact = get_best_impact(hg19_annotations)
        hg38_impact = get_best_impact(hg38_annotations)
        
        # Calculate gene and impact changes using clean_string
        if clean_string(hg19_gene) != clean_string(hg38_gene):
            gene_changes = 1
        
        if clean_string(hg19_impact) != clean_string(hg38_impact):
            impact_changes = 1
        
        # Clinical significance and pathogenicity predictions
        def get_best_clinical_data(annotations_df, column):
            """Get clinical data with proper handling of missing values"""
            if len(annotations_df) > 0:
                value = annotations_df[column].iloc[0]
                if pd.notna(value) and value != '' and value != '-':
                    return value
            return ''
        
        hg19_clin_sig = get_best_clinical_data(hg19_annotations, 'clin_sig')
        hg38_clin_sig = get_best_clinical_data(hg38_annotations, 'clin_sig')
        
        # Normalize clinical significance
        hg19_clin_sig_normalized = normalize_clinical_significance(hg19_clin_sig)
        hg38_clin_sig_normalized = normalize_clinical_significance(hg38_clin_sig)
        hg19_sift = get_best_clinical_data(hg19_annotations, 'sift')
        hg38_sift = get_best_clinical_data(hg38_annotations, 'sift')
        hg19_polyphen = get_best_clinical_data(hg19_annotations, 'polyphen')
        hg38_polyphen = get_best_clinical_data(hg38_annotations, 'polyphen')
        
        # Parse pathogenicity predictions
        hg19_sift_pred, hg19_sift_score = parse_sift_prediction(hg19_sift)
        hg38_sift_pred, hg38_sift_score = parse_sift_prediction(hg38_sift)
        hg19_polyphen_pred, hg19_polyphen_score = parse_polyphen_prediction(hg19_polyphen)
        hg38_polyphen_pred, hg38_polyphen_score = parse_polyphen_prediction(hg38_polyphen)
        
        # Check for clinical significance changes (using normalized categories)
        hg19_is_pathogenic = (hg19_clin_sig_normalized == 'PATHOGENIC')
        hg38_is_pathogenic = (hg38_clin_sig_normalized == 'PATHOGENIC')
        hg19_is_benign = (hg19_clin_sig_normalized == 'BENIGN')
        hg38_is_benign = (hg38_clin_sig_normalized == 'BENIGN')
        
        # Create directional clinical significance change
        if hg19_clin_sig_normalized != hg38_clin_sig_normalized:
            clin_sig_change = f'{hg19_clin_sig_normalized}_TO_{hg38_clin_sig_normalized}'
        else:
            clin_sig_change = f'STABLE_{hg19_clin_sig_normalized}'
        
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
            'transcript_pairs_analyzed': transcript_pairs_analyzed,
            'same_transcript_consequence_changes': same_transcript_consequence_changes,
            'same_consequence_different_transcripts': same_consequence_different_transcripts,
            'unmatched_consequences': unmatched_consequences,
            'gene_changes': gene_changes,
            'impact_changes': impact_changes,
            'problematic_transcripts_hg19': '; '.join(problematic_transcripts_hg19) if problematic_transcripts_hg19 else '',
            'problematic_transcripts_hg38': '; '.join(problematic_transcripts_hg38) if problematic_transcripts_hg38 else '',
            'transcript_relationship': transcript_relationship,
            'consequence_relationship': consequence_relationship,
            
            # Representative VEP information (FIXED)
            'hg19_gene': hg19_gene,
            'hg38_gene': hg38_gene,
            'hg19_consequences_set': ', '.join(sorted(all_hg19_consequences)) if 'all_hg19_consequences' in locals() else '',
            'hg38_consequences_set': ', '.join(sorted(all_hg38_consequences)) if 'all_hg38_consequences' in locals() else '',
            'hg19_high_impact_consequences': ', '.join(sorted([cons for cons in (all_hg19_consequences if 'all_hg19_consequences' in locals() else set()) 
                                                  if VEP_CONSEQUENCE_IMPACT.get(cons, 'MODIFIER') in {'HIGH', 'MODERATE'}])),
            'hg38_high_impact_consequences': ', '.join(sorted([cons for cons in (all_hg38_consequences if 'all_hg38_consequences' in locals() else set()) 
                                                  if VEP_CONSEQUENCE_IMPACT.get(cons, 'MODIFIER') in {'HIGH', 'MODERATE'}])),

            'hg19_impact': hg19_impact,
            'hg38_impact': hg38_impact,
            'hg19_clin_sig': hg19_clin_sig,
            'hg38_clin_sig': hg38_clin_sig,
            'hg19_clin_sig_normalized': hg19_clin_sig_normalized,
            'hg38_clin_sig_normalized': hg38_clin_sig_normalized,
            'clin_sig_change': clin_sig_change,
            'hg19_gnomad_af': hg19_annotations['gnomadg_af'].iloc[0] if len(hg19_annotations) > 0 else None,
            'hg38_gnomad_af': hg38_annotations['gnomadg_af'].iloc[0] if len(hg38_annotations) > 0 else None,
            'hg19_sift': hg19_sift,
            'hg38_sift': hg38_sift,
            'sift_change': sift_change,
            'hg19_polyphen': hg19_polyphen,
            'hg38_polyphen': hg38_polyphen,
            'polyphen_change': polyphen_change,
            
            # Pathogenicity flags
            'hg19_is_pathogenic': hg19_is_pathogenic,
            'hg38_is_pathogenic': hg38_is_pathogenic,
            'hg19_is_benign': hg19_is_benign,
            'hg38_is_benign': hg38_is_benign
        }

        # MANE analysis 
        hg38_mane_flag, hg38_mane_transcript_id, hg38_mane_details = self._analyze_mane_annotations(hg38_annotations)
        
        # Check if MANE transcripts from hg38 are present in hg19 (exact version match)
        hg19_transcript_ids = set()
        for _, row in hg19_transcripts_df.iterrows():
            full_id = row['feature']
            if full_id:
                hg19_transcript_ids.add(full_id)
        
        # Find which MANE transcripts are present in hg19 (exact version match)
        hg19_mane_present = []
        if hg38_mane_transcript_id:
            if hg38_mane_transcript_id in hg19_transcript_ids:
                hg19_mane_present.append(hg38_mane_transcript_id)
        
        hg19_mane_transcript_id = hg19_mane_present[0] if hg19_mane_present else None

        # Format hg19 MANE details to match hg38 format
        if hg19_mane_present:
            # Use the same MANE type as hg38 for consistency
            if hg38_mane_flag == "MANE_Select":
                hg19_mane_details = f"MANE_Select:{hg19_mane_present[0]}"
            elif hg38_mane_flag == "MANE_Plus_Clinical":
                hg19_mane_details = f"MANE_Plus_Clinical:{hg19_mane_present[0]}"
            else:
                hg19_mane_details = hg19_mane_present[0]
        else:
            hg19_mane_details = "Not_Present"
        
        # Canonical transcript identification 
        hg19_canonical_transcript, hg38_canonical_transcript = self._identify_canonical_transcripts(hg19_transcripts_df, hg38_transcripts_df)

        # Priority transcript selection using MANE-first hierarchy
        transcript_crossbuild_status, priority_transcript_crossbuild = self._select_priority_transcripts(
            hg19_transcripts_df, hg38_transcripts_df, hg38_mane_flag, hg38_mane_transcript_id
        )

        # Analyze HGVSc concordance 
        priority_hgvs_analysis = analyze_priority_transcript_hgvs(
            hg19_transcripts_df, hg38_transcripts_df, transcript_crossbuild_status, priority_transcript_crossbuild
    )

        # Add HGVSc results to the return dictionary
        vep_analysis.update({
             # Basic canonical transcript identification
            'hg19_canonical_transcript': hg19_canonical_transcript,
            'hg38_canonical_transcript': hg38_canonical_transcript,
        
            # MANE information (hg38 source, hg19 presence check)
            'hg38_mane_flag': hg38_mane_flag,
            'hg38_mane_transcript_id': hg38_mane_transcript_id,
            'hg38_mane_details': hg38_mane_details,
            'hg19_mane_transcript_id': hg19_mane_transcript_id,
            'hg19_mane_details': hg19_mane_details,
            # Priority transcript selection using MANE-first hierarchy
            'transcript_crossbuild_status': transcript_crossbuild_status,
            'priority_transcript_crossbuild': priority_transcript_crossbuild,
            # Priority transcript HGVS analysis
            'priority_hgvsc_hg19': priority_hgvs_analysis['priority_hgvsc_hg19'],
            'priority_hgvsc_hg38': priority_hgvs_analysis['priority_hgvsc_hg38'],
            'priority_hgvsp_hg19': priority_hgvs_analysis['priority_hgvsp_hg19'],
            'priority_hgvsp_hg38': priority_hgvs_analysis['priority_hgvsp_hg38'],
            'priority_hgvsc_concordance': priority_hgvs_analysis['priority_hgvsc_concordance'],
            'priority_hgvsp_concordance': priority_hgvs_analysis['priority_hgvsp_concordance']
        })
  
        return vep_analysis
        
    def _analyze_mane_annotations(self, annotations_df):
        """
        Analyze all MANE annotations for a variant and determine priority
        
        Args:
            annotations_df: DataFrame with VEP annotations
            
        Returns:
            tuple: (mane_flag, mane_transcript_id, mane_details)
            - mane_flag: "MANE_Select" | "MANE_Plus_Clinical" | "Both" | "None"
            - mane_transcript_id: RefSeq ID of highest priority MANE transcript (Select > Clinical)
            - mane_details: String with all MANE transcripts found
        """
        if len(annotations_df) == 0:
            return "None", None, ""
        
        mane_select_transcripts = []
        mane_clinical_transcripts = []
        all_mane_details = []
        
        for _, row in annotations_df.iterrows():
            mane_field = row.get('mane', '')
            refseq_id = row.get('refseq_transcript_id', '')
            
            if pd.notna(mane_field) and str(mane_field).strip():
                mane_str = str(mane_field).strip()
                
                if "MANE_Select" in mane_str and refseq_id:
                        mane_select_transcripts.append(refseq_id)
                        all_mane_details.append(f"MANE_Select:{refseq_id}")
                elif "MANE_Plus_Clinical" in mane_str and refseq_id:
                    mane_clinical_transcripts.append(refseq_id)
                    all_mane_details.append(f"MANE_Plus_Clinical:{refseq_id}")
        
        # Determine flag and primary transcript ID
        has_select = bool(mane_select_transcripts)
        has_clinical = bool(mane_clinical_transcripts)
        
        if has_select and has_clinical:
            return "Both", mane_select_transcripts[0], "; ".join(all_mane_details)
        elif has_select:
            return "MANE_Select", mane_select_transcripts[0], "; ".join(all_mane_details)
        elif has_clinical:
            return "MANE_Plus_Clinical", mane_clinical_transcripts[0], "; ".join(all_mane_details)
        else:
            return "None", None, ""
        
    def _identify_canonical_transcripts(self, hg19_transcripts_df, hg38_transcripts_df):
        """
        Identify canonical transcripts for hg19 and hg38 builds
        
        Args:
            hg19_transcripts_df: hg19 transcript annotations
            hg38_transcripts_df: hg38 transcript annotations
            
        Returns:
            tuple: (hg19_canonical_transcript, hg38_canonical_transcript)
        """
        # Find canonical transcript for hg19
        hg19_canonical = ""
        canonical_hg19_rows = hg19_transcripts_df[hg19_transcripts_df['is_canonical'] == 1]
        if len(canonical_hg19_rows) > 0:
            hg19_canonical = canonical_hg19_rows.iloc[0]['feature']
        
        # Find canonical transcript for hg38
        hg38_canonical = ""
        canonical_hg38_rows = hg38_transcripts_df[hg38_transcripts_df['is_canonical'] == 1]
        if len(canonical_hg38_rows) > 0:
            hg38_canonical = canonical_hg38_rows.iloc[0]['feature']
        
        return hg19_canonical, hg38_canonical
    
    def _select_priority_transcripts(self, hg19_transcripts_df, hg38_transcripts_df, hg38_mane_flag, hg38_mane_transcript_id):
        """
        Select priority transcripts using MANE-first hierarchy with exact version matching
        
        Priority order: MANE Select > MANE Plus Clinical > Canonical fallback
        
        Args:
            hg19_transcripts_df: hg19 transcript annotations
            hg38_transcripts_df: hg38 transcript annotations
            hg38_mane_flag: MANE flag from hg38 ("MANE_Select", "MANE_Plus_Clinical", "Both", "None")
            hg38_mane_transcript_id: RefSeq ID of MANE transcript from hg38
            
        Returns:
            tuple: (transcript_crossbuild_status, priority_transcript_crossbuild)
        """
        # Build hg19 transcript set with exact versions
        hg19_transcript_ids = set()
        for _, row in hg19_transcripts_df.iterrows():
            full_id = row['feature']
            if full_id:
                hg19_transcript_ids.add(full_id)
        
        # Case 1: MANE transcript available in hg38
        if hg38_mane_flag in ["MANE_Select", "MANE_Plus_Clinical", "Both"] and hg38_mane_transcript_id:
            if hg38_mane_transcript_id in hg19_transcript_ids:
                # Exact match found - same transcript with same version
                if hg38_mane_flag == "MANE_Select" or hg38_mane_flag == "Both":
                    return "MANE_Select_Both_Builds", hg38_mane_transcript_id
                else:
                    return "MANE_Plus_Clinical_Both_Builds", hg38_mane_transcript_id
            else:
                return "MANE_hg38_Only", "NONE"
        
        # Case 2: No MANE - fallback to canonical with exact matching
        hg19_canonical = None
        hg38_canonical = None
        
        # Find canonical transcripts
        for _, row in hg19_transcripts_df.iterrows():
            if row.get('is_canonical', 0) == 1:
                hg19_canonical = row['feature']
                break
                
        for _, row in hg38_transcripts_df.iterrows():
            if row.get('is_canonical', 0) == 1:
                hg38_canonical = row['feature']
                break
        
        # Check if canonical transcripts match exactly
        if hg19_canonical and hg38_canonical and hg19_canonical == hg38_canonical:
            return "Canonical_Fallback_Both_Builds", hg19_canonical
        elif hg19_canonical or hg38_canonical:
            return "No_Matching_Transcripts", "NONE"
        else:
            return "No_Transcripts", "NONE"
