RESULTS_DIR = config['results_dir']
VEP_HG19_DIR = os.path.join(RESULTS_DIR, config['dirs']['vep_hg19'])
VEP_HG38_DIR = os.path.join(RESULTS_DIR, config['dirs']['vep_hg38'])

# VEP hg19

rule strip_id_hg19:
    input:
        vcf = rules.crossmap_liftover.output.vcf
    output:
        temp(f"{VEP_HG19_DIR}/{config['sample']}.noID.vcf")
    shell:
        """
        bcftools annotate -x ID {input.vcf} -Ov -o {output}
        """

rule vep_hg19:
    input:
        vcf = rules.strip_id_hg19.output
    output:
        txt = f"{VEP_HG19_DIR}/{config['sample']}.vep.hg19.txt"
    params:
        cache = config["vep_cache"]["hg19"],
        fasta = config["ref"]["hg19_fasta"]
    shell:
        """
        mkdir -p {VEP_HG19_DIR}
        vep \
            --dir_cache {params.cache} \
            --cache \
            --offline \
            --fasta {params.fasta} \
            --force_overwrite \
            --symbol \
            --refseq \
            --transcript_version \
            --hgvs \
            --hgvsg \
            --numbers \
            --domains \
            --regulatory \
            --canonical \
            --protein \
            --allele_number \
            --no_escape \
            --failed=1 \
            --exclude_predicted \
            --sift=b \
            --polyphen=b \
            --af_gnomadg \
	    --output_format tab \
	    --fields "Uploaded_variation,Location,Allele,Gene,Feature,Feature_type,Consequence,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,IMPACT,DISTANCE,STRAND,FLAGS,MINIMISED,SYMBOL,SYMBOL_SOURCE,HGNC_ID,CANONICAL,ENSP,REFSEQ_MATCH,REFSEQ_OFFSET,GIVEN_REF,USED_REF,BAM_EDIT,SIFT,PolyPhen,EXON,INTRON,DOMAINS,HGVSc,HGVSp,HGVS_OFFSET,HGVSg,gnomADg_AF,gnomADg_AFR_AF,gnomADg_AMI_AF,gnomADg_AMR_AF,gnomADg_ASJ_AF,gnomADg_EAS_AF,gnomADg_FIN_AF,gnomADg_MID_AF,gnomADg_NFE_AF,gnomADg_REMAINING_AF,gnomADg_SAS_AF,CLIN_SIG,SOMATIC,PHENO,BIOTYPE,MOTIF_NAME,MOTIF_POS,HIGH_INF_POS,MOTIF_SCORE_CHANGE,TRANSCRIPTION_FACTORS" \
            -i {input.vcf} \
            -o {output.txt}
        """

# VEP hg38
rule strip_id_hg38:
    input:
        vcf = rules.bcftools_liftover.output.vcf
    output:
        temp(f"{VEP_HG38_DIR}/{config['sample']}.noID.vcf")
    shell:
        """
        bcftools annotate -x ID {input.vcf} -Ov -o {output}
        """

rule vep_hg38:
    input:
        vcf = rules.strip_id_hg38.output
    output:
        txt = f"{VEP_HG38_DIR}/{config['sample']}.vep.hg38.txt"
    params:
        cache = config["vep_cache"]["hg38"],
        fasta = config["ref"]["hg38_fasta"]
    shell:
        """
        mkdir -p {VEP_HG38_DIR}
        vep \
            --dir_cache {params.cache} \
            --cache \
            --offline \
            --fasta {params.fasta} \
            --force_overwrite \
            --symbol \
            --mane \
            --refseq \
            --transcript_version \
            --hgvs \
            --hgvsg \
            --regulatory \
            --canonical \
            --protein \
            --no_escape \
            --failed=1 \
            --exclude_predicted \
            --sift=b \
            --polyphen=b \
            --af_gnomadg \
	    --output_format tab \
	    --fields "Uploaded_variation,Location,Allele,Gene,Feature,Feature_type,Consequence,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,IMPACT,DISTANCE,STRAND,FLAGS,MINIMISED,SYMBOL,SYMBOL_SOURCE,HGNC_ID,MANE,MANE_SELECT,MANE_PLUS_CLINICAL,CANONICAL,ENSP,REFSEQ_MATCH,REFSEQ_OFFSET,GIVEN_REF,USED_REF,BAM_EDIT,SIFT,PolyPhen,EXON,INTRON,DOMAINS,HGVSc,HGVSp,HGVS_OFFSET,HGVSg,gnomADg_AF,gnomADg_AFR_AF,gnomADg_AMI_AF,gnomADg_AMR_AF,gnomADg_ASJ_AF,gnomADg_EAS_AF,gnomADg_FIN_AF,gnomADg_MID_AF,gnomADg_NFE_AF,gnomADg_REMAINING_AF,gnomADg_SAS_AF,CLIN_SIG,SOMATIC,PHENO,BIOTYPE,MOTIF_NAME,MOTIF_POS,HIGH_INF_POS,MOTIF_SCORE_CHANGE,TRANSCRIPTION_FACTORS" \
            -i {input.vcf} \
            -o {output.txt}
        """
