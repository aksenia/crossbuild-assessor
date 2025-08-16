RESULTS_DIR = config['results_dir']
LIFTOVER_CRMP = os.path.join(RESULTS_DIR, config['dirs']['liftover_crossmap'])
LIFTOVER_BCFT = os.path.join(RESULTS_DIR, config['dirs']['liftover_bcftools'])

# CrossMap liftover
rule add_id_hg19:
    input:
        vcf = config["input_vcf"]
    output:
        temp(f"{LIFTOVER_CRMP}/{config['sample']}.ID.vcf")
    shell:
        """
        bcftools annotate --set-id '%CHROM/%POS/%REF//%ALT' {input.vcf} -Ov -o {output}
        """

rule crossmap_liftover:
    input:
        vcf = rules.add_id_hg19.output
    output:
        vcf = f"{LIFTOVER_CRMP}/{config['sample']}.vcf"
    params:
        script = config["tools"]["crossmap_script"],
        chain  = config["ref"]["liftover_chain"],
        fasta  = config["ref"]["hg38_fasta"]
    shell:
        """
        mkdir -p {LIFTOVER_CRMP}
        python {params.script} {input.vcf} {params.chain} {output.vcf} --output-bed
        """

# bcftools liftover
rule bcftools_liftover:
    input:
        vcf = rules.crossmap_liftover.output.vcf
    output:
        temp(f"{LIFTOVER_BCFT}/{config['sample']}.bt_noid.vcf")
    params:
        chain  = config["ref"]["liftover_chain"],
        srcfa  = config["ref"]["hg19_fasta"],
        dstfa  = config["ref"]["hg38_fasta"]
    shell:
        """
        mkdir -p {LIFTOVER_BCFT}
        bcftools +liftover -Ov {input.vcf} -o {output} \
            -- -c {params.chain} \
               --src-fasta-ref {params.srcfa} \
               --fasta-ref {params.dstfa} \
               --write-src \
               --reject {LIFTOVER_BCFT}/{config[sample]}_rejected.vcf \
               --write-reject
        """

rule add_id_hg38:
    input:
        vcf = rules.bcftools_liftover.output
    output:
        vcf = f"{LIFTOVER_BCFT}/{config['sample']}.bt.vcf"
    shell:
        """
        bcftools annotate --set-id '%CHROM/%POS/%REF//%ALT' {input.vcf} -Ov -o {output.vcf}
        """
