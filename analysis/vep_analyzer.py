"""
VEP Analyzer

Core VEP analysis engine for processing variant annotations and identifying
transcript-level discordances between genome builds.
"""

import pandas as pd
from utils.transcript_utils import normalize_transcript_id, extract_genotype_from_alleles
#from utils.data_utils import clean_string
from utils.clinical_utils import (
    is_pathogenic_clinical_significance,
    is_benign_clinical_significance, 
    parse_sift_prediction,
    parse_polyphen_prediction,
    normalize_clinical_significance
)
from config.constants import VEP_CONSEQUENCE_IMPACT

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
            hg19_tx_dict = {}  # transcript_base_id -> data
            hg38_tx_dict = {}  # transcript_base_id -> data
            
            # Normalize transcript IDs and group by base ID
            for tx_id, data in hg19_transcripts.items():
                base_id = normalize_transcript_id(tx_id)
                if base_id:
                    hg19_tx_dict[base_id] = data
            
            for tx_id, data in hg38_transcripts.items():
                base_id = normalize_transcript_id(tx_id)
                if base_id:
                    hg38_tx_dict[base_id] = data
            
            # Find same transcript IDs with different consequences
            for base_id in set(hg19_tx_dict.keys()) & set(hg38_tx_dict.keys()):
                hg19_data = hg19_tx_dict[base_id]
                hg38_data = hg38_tx_dict[base_id]
                
                # Use robust string comparison
                consequence_match = clean_string(hg19_data['consequence']) == clean_string(hg38_data['consequence'])
                
                if not consequence_match:
                    same_transcript_consequence_changes += 1
                    problematic_transcripts_hg19.append(f"{hg19_data['feature_id']}({hg19_data['consequence']})")
                    problematic_transcripts_hg38.append(f"{hg38_data['feature_id']}({hg38_data['consequence']})")
            
            # Step 2: Set-based transcript relationship analysis
            all_hg19_transcript_ids = {normalize_transcript_id(tx_id) for tx_id in hg19_transcripts.keys() if normalize_transcript_id(tx_id)}
            all_hg38_transcript_ids = {normalize_transcript_id(tx_id) for tx_id in hg38_transcripts.keys() if normalize_transcript_id(tx_id)}
            
            # Transcript relationship
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
            
            # Step 3: Set-based consequence relationship analysis - FIXED
            all_hg19_consequences = set()
            all_hg38_consequences = set()

            # Build sets more robustly
            for data in hg19_transcripts.values():
                cleaned_cons = clean_string(data['consequence'])
                if cleaned_cons:
                    all_hg19_consequences.add(cleaned_cons)

            for data in hg38_transcripts.values():
                cleaned_cons = clean_string(data['consequence'])
                if cleaned_cons:
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
        
        # FIXED: Get representative VEP information with proper fallbacks
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
  
        return vep_analysis