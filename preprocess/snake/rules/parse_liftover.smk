RESULTS_DIR = config['results_dir']
COMPARISON_DIR = os.path.join(RESULTS_DIR, config['dirs']['comparison'])
rule parse_liftover:
    input:
        bcftools_vcf = rules.add_id_hg38.output.vcf
    output:
        tsv = f"{COMPARISON_DIR}/{config['sample']}_comparison.tsv"
    params:
        script = config["tools"]["parse_script"]
    shell:
        """
        mkdir -p {COMPARISON_DIR}
        python {params.script} {input.bcftools_vcf} {output.tsv}
        """

