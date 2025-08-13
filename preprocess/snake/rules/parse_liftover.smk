rule parse_liftover:
    input:
        bcftools_vcf = rules.bcftools_liftover.output.vcf
    output:
        tsv = f"{config['dirs']['comparison']}/{config['sample']}_comparison.tsv"
    params:
        script = config["tools"]["parse_script"]
    shell:
        """
        mkdir -p {config[dirs][comparison]}
        python {params.script} {input.bcftools_vcf} {output.tsv}
        """

