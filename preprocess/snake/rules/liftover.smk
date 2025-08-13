RESULTS_DIR = config['results_dir']
LIFTOVER_CRMP = os.path.join(RESULTS_DIR, config['dirs']['liftover_crossmap'])
LIFTOVER_BCFT = os.path.join(RESULTS_DIR, config['dirs']['liftover_bcftools'])

# CrossMap liftover
rule crossmap_liftover:
    input:
        vcf = config["input_vcf"]
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
        vcf = f"{LIFTOVER_BCFT}/{config['sample']}.bt.vcf"
    params:
        chain  = config["ref"]["liftover_chain"],
        srcfa  = config["ref"]["hg19_fasta"],
        dstfa  = config["ref"]["hg38_fasta"]
    shell:
        """
        mkdir -p {LIFTOVER_BCFT}
        bcftools +liftover -Ov {input.vcf} -o {output.vcf} \
            -- -c {params.chain} \
               --src-fasta-ref {params.srcfa} \
               --fasta-ref {params.dstfa} \
               --write-src \
               --reject {LIFTOVER_BCFT}/{config[sample]}_rejected.vcf \
               --write-reject
        """

