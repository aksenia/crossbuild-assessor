# CrossMap liftover
rule crossmap_liftover:
    input:
        vcf = config["input_vcf"]
    output:
        vcf = f"{config['dirs']['liftover_crossmap']}/{config['sample']}.vcf"
    params:
        script = config["tools"]["crossmap_script"],
        chain  = config["ref"]["liftover_chain"],
        fasta  = config["ref"]["hg38_fasta"]
    shell:
        """
        mkdir -p {config[dirs][liftover_crossmap]}
        python {params.script} {input.vcf} {output.vcf}
        """

# bcftools liftover
rule bcftools_liftover:
    input:
        vcf = rules.crossmap_liftover.output.vcf
    output:
        vcf = f"{config['dirs']['liftover_bcftools']}/{config['sample']}.bt.vcf"
    params:
        chain  = config["ref"]["liftover_chain"],
        srcfa  = config["ref"]["hg19_fasta"],
        dstfa  = config["ref"]["hg38_fasta"]
    shell:
        """
        mkdir -p {config[dirs][liftover_bcftools]}
        bcftools +liftover -Ov {input.vcf} -o {output.vcf} \
            -- -c {params.chain} \
               --src-fasta-ref {params.srcfa} \
               --fasta-ref {params.dstfa} \
               --write-src \
               --reject {config[dirs][liftover_bcftools]}/{config[sample]}_rejected.vcf \
               --write-reject
        """

