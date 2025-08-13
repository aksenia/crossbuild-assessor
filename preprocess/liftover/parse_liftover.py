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
        "bcftools_hg38_chrom", "bcftools_hg38_pos", "bcftools_hg38_ref",
        "bcftools_hg38_alt", "pos_match", "gt_match"
    ])
    for line in vcf:
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        info_dict = {}
        for item in fields[7].split(";"):
            if "=" in item:
                k, v = item.split("=", 1)
                info_dict[k] = v
            else:
                info_dict[item] = True
        map_status = info_dict.get(

