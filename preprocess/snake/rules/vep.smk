# VEP hg19
rule vep_hg19:
    input:
        vcf = config["input_vcf"]
    output:
        txt = f"{config['dirs']['vep_hg19']}/{config['sample']}.vep.txt"
    params:
        cache = config["vep_cache"]["hg19"],
        fasta = config["ref"]["hg19_fasta"]
    container:
        "docker://your-container-image"
    shell:
        """
        mkdir -p {config[dirs][vep_hg19]}
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
            -i {input.vcf} \
            -o {output.txt}
        """

# VEP hg38
rule vep_hg38:
    input:
        vcf = rules.bcftools_liftover.output.vcf
    output:
        txt = f"{config['dirs']['vep_hg38']}/{config['sample']}.vep.txt"
    params:
        cache = config["vep_cache"]["hg38"],
        fasta = config["ref"]["hg38_fasta"]
    shell:
        """
        mkdir -p {config[dirs][vep_hg38]}
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
            -i {input.vcf} \
            -o {output.txt}
        """
