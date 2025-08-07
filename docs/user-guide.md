# User guide

Complete usage guide for CrossBuild Assessor.

## Prerequisites

- Docker
- 8GB+ RAM for large datasets
- Required input files (see Data preparation)

## Installation

```bash
git clone https://github.com/your-org/crossbuild-assessor.git
cd crossbuild-assessor
docker build -t crossbuild-assessor .
```

## Data preparation

### 1. Liftover comparison file

Tab-separated file comparing CrossMap vs bcftools results.

**Required columns:**
```
mapping_status    source_chrom    source_pos    source_alleles
flip             swap            liftover_hg38_chrom    liftover_hg38_pos
bcftools_hg38_chrom    bcftools_hg38_pos    bcftools_hg38_ref
bcftools_hg38_alt    pos_match    gt_match
```

### 2. VEP annotation files

Standard VEP output for both genome builds:

- use "tab" option
- make sure "ID" field in your vcf is set to "." Then VEP will generate it by default as a first field in output and the assessor tool will parse it.

```bash
# For hg19
vep --input_file variants_hg19.vcf --output_file variants_hg19.vep.txt \
    --assembly GRCh37 --tab --cache --offline --sift b --polyphen b

# For hg38  
vep --input_file variants_hg38.vcf --output_file variants_hg38.vep.txt \
    --assembly GRCh38 --tab --cache --offline --sift b --polyphen b
```

### 3. Configuration file

```json
{
  "input_files": {
    "comparison": "/path/to/liftover_comparison.txt",
    "hg19_vep": "/path/to/variants_hg19.vep.txt",
    "hg38_vep": "/path/to/variants_hg38.vep.txt"
  },
  "database": {
    "path": "/path/to/genomic_analysis.db"
  }
}
```

## Pipeline workflow

### Step 1: Database creation

```bash
docker run -v /path/to/data:/data crossbuild-assessor \
  python db_loader.py --config /data/config.json
```

### Step 2: Quality control analysis

```bash
docker run -v /path/to/data:/data crossbuild-assessor \
  python db_analyzer.py --db-path /data/genomic_analysis.db --output-dir /data/qc/
```

### Step 3: Variant prioritization

```bash
docker run -v /path/to/data:/data crossbuild-assessor \
  python variant_prioritizer.py --db-path /data/genomic_analysis.db --output-dir /data/results/
```

### Step 4: Generate HTML report

```bash
docker run -v /path/to/data:/data crossbuild-assessor \
  python report_generator.py --input-dir /data/results/ --output /data/report.html
```

## Output files

**Database loader**: `genomic_analysis.db` - SQLite database  
**QC analyzer**: Concordance plots and statistics  
**Prioritizer**: `prioritized_variants.csv`, visualization plots, summary  
**Report generator**: `crossbuild_report.html` - Unified clinical dashboard

## Priority categories

- **CRITICAL**: Clinical interpretation changes, immediate review needed
- **HIGH**: Functionally significant changes, priority review  
- **MODERATE**: Prediction changes, standard review
- **LOW**: Technical issues, secondary review

## Usage tips

- Use `--force` to recalculate with new parameters
- Use `--no-plots` for faster processing
- Filter results with `--min-score` and `--max-variants`
- First run creates cache for faster subsequent analysis