#!/usr/bin/env python3
import sys
import csv

if len(sys.argv) != 3:
    print(f"Usage: {sys.argv[0]} input.vcf output.tsv", file=sys.stderr)
    sys.exit(1)

vcf_file = sys.argv[1]
out_file = sys.argv[2]

with open(vcf_file) as vcf, open(out_file, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow([
        "mapping_status", "source_chrom", "source_pos", "source_alleles",
        "flip", "swap", "liftover_hg38_chrom", "liftover_hg38_pos",
        "bcftools_hg38_chrom", "bcftools_hg38_pos",
        "bcftools_hg38_ref", "bcftools_hg38_alt",
        "pos_match", "gt_match"
    ])

    for line in vcf:
        if line.startswith("#"):
            continue

        fields = line.strip().split("\t")
        bcftools_chr, bcftools_pos, _, bcftools_ref, bcftools_alt, _, _, info = fields[:8]

        # Initialize defaults
        map_status = "NA"
        src_chr = "NA"
        src_pos = "NA"
        src_alleles = "NA"
        flip = "no_flip"
        swap = "NA"
        liftover_chr = "NA"
        liftover_pos = "NA"

        # Parse INFO field like awk
        for item in info.split(";"):
            if item.startswith("hg38_map="):
                map_status = item.split("=", 1)[1]
            elif item.startswith("SRC_CHROM="):
                src_chr = item.split("=", 1)[1]
            elif item.startswith("SRC_POS="):
                src_pos = item.split("=", 1)[1]
            elif item.startswith("SRC_REF_ALT="):
                src_alleles = item.split("=", 1)[1]
            elif item == "FLIP":
                flip = "flip"
            elif item.startswith("SWAP="):
                swap = item.split("=", 1)[1]
            elif item.startswith("hg38_chr="):
                liftover_chr = item.split("=", 1)[1]
            elif item.startswith("hg38_start="):
                liftover_pos = item.split("=", 1)[1]

        # Remove "chr" prefix from bcftools_chr (like second awk)
        bcftools_chr_stripped = bcftools_chr.replace("chr", "")

        # pos_match logic
        pos_match = "TRUE" if (bcftools_chr_stripped == liftover_chr and bcftools_pos == liftover_pos) else "FALSE"

        # gt_match logic
        gt_field = f"{bcftools_ref},{bcftools_alt}"
        gt_match = "TRUE" if gt_field == src_alleles else "FALSE"

        writer.writerow([
            map_status, src_chr, src_pos, src_alleles,
            flip, swap, liftover_chr, liftover_pos,
            bcftools_chr, bcftools_pos, bcftools_ref, bcftools_alt,
            pos_match, gt_match
        ])

