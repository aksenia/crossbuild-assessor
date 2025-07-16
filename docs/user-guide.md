# User Guide

Complete guide to using CrossBuild Assessor for genomic variant analysis across genome builds.

## Prerequisites

### System requirements
- Docker
- 8GB+ RAM (for large datasets)
- 10GB+ storage space

### Required input files
1. **Liftover comparison file** (tab-separated)
2. **VEP annotation files** for both hg19 and hg38
3. **Configuration file** (JSON format)

## Installation

```bash
git clone https://github.com/your-org/crossbuild-assessor.git
cd crossbuild-assessor
docker build -t crossbuild-assessor .
```

## Data preparation

### 1. Liftover comparison file

Create a tab-separated file comparing liftover results between tools:

**Required columns:**
```
mapping_status    source_chrom    source_pos    source_alleles
flip             swap            liftover_hg38_chrom    liftover_hg38_pos
bcftools_hg38_chrom    bcftools_hg38_pos    bcftools_hg38_ref
bcftools_hg38_alt    pos_match    gt_match
```

### 2. VEP annotation files

Run VEP on your variants for both genome builds with standard parameters:

```bash
# For hg19
vep --input_file variants_hg19.vcf --output_file variants_hg19.vep.txt \
    --assembly GRCh37 --tab --cache --offline --sift b --polyphen b

# For hg38  
vep --input_file variants_hg38.vcf --output_file variants_hg38.vep.txt \
    --assembly GRCh38 --tab --cache --offline --sift b --polyphen b
```

### 3. Configuration file

Create a `config.json` file:

```json
{
  "input_files": {
    "comparison": "/path/to/liftover_comparison.txt",
    "hg19_vep": "/path/to/variants_hg19.vep.txt",
    "hg38_vep": "/path/to/variants_hg38.vep.txt"
  },
  "database": {
    "path": "/path/to/genomic_analysis.db"
  },
  "output": {
    "directory": "/path/to/output"
  }
}
```

## Pipeline workflow

### Step 1: Database creation

Load your data into an optimized SQLite database:

```bash
docker run -v /path/to/data:/data crossbuild-assessor \
  python db_loader.py --config /data/config.json
```

**Output:** Optimized SQLite database with indexed tables for fast analysis.

### Step 2: Liftover quality control

Analyze liftover performance between tools:

```bash
docker run -v /path/to/data:/data crossbuild-assessor \
  python db_analyzer.py --db-path /data/genomic_analysis.db --output-dir /data/qc_results/
```

**Output:** 
- Concordance plots and statistics
- Tool performance comparison
- Quality control reports

### Step 3: Variant prioritization

Generate ranked list of problematic variants:

```bash
docker run -v /path/to/data:/data crossbuild-assessor \
  python variant_prioritizer.py --db-path /data/genomic_analysis.db --output-dir /data/results/
```

**Output:**
- `prioritized_variants.csv` - Excel-compatible ranked variant list
- `variant_prioritization_plots.png` - 4-plot analysis dashboard
- `variant_prioritization_summary.txt` - Detailed statistics

## Output interpretation

### Priority categories

Variants are classified into five categories for clinical review:

- **CRITICAL** - Same transcript, different consequences (immediate review)
- **HIGH** - Clinical significance changes or major impact transitions
- **MODERATE** - Gene changes or pathogenicity prediction changes  
- **INVESTIGATE** - Unmatched consequences or minor transitions
- **LOW** - Benign variants with minimal impact

### Excel output columns

Key columns in `prioritized_variants.csv`:

- `Priority_Score` - Numerical score for ranking
- `Priority_Category` - Clinical review category
- `GT_hg19` / `GT_hg38` - Genotypes in both builds
- `Gene_hg19` / `Gene_hg38` - Gene symbols
- `Clinical_Significance_Change` - ClinVar changes between builds
- `Problematic_Transcripts_hg19/hg38` - Specific transcript issues

### Visualization dashboard

The 4-plot dashboard shows:

1. **Priority categories** - Distribution of variants by review level
2. **Clinical evidence** - Evidence types driving prioritization
3. **Discordance types** - Primary types of annotation differences
4. **Clinical transitions** - Dynamic vs static clinical evidence

## Advanced usage

### Custom scoring parameters

Modify scoring weights in `config/scoring_config.py`:

```python
BASE_SCORES = {
    'same_transcript_consequence_changes': 15,  # Increase weight
    'clinical_sig_benign_to_pathogenic': 12,   # Custom weight
    # ... other parameters
}
```

### Filtering options

Filter results by score or category:

```bash
# Only high-priority variants
python variant_prioritizer.py --min-score 10 --max-variants 1000

# Skip visualization for faster processing
python variant_prioritizer.py --no-plots
```

### Batch processing

Process multiple datasets:

```bash
for config in configs/*.json; do
  python db_loader.py --config "$config"
  python variant_prioritizer.py --db-path "$(jq -r .database.path "$config")"
done
```

## Troubleshooting

### Common issues

**Import errors**
- Ensure all config files are in the correct directory structure
- Check that file paths in config.json are absolute paths

**Memory issues**  
- Reduce chunk size in analysis settings
- Process smaller datasets or use a machine with more RAM

**Missing VEP data**
- Verify VEP files have the required columns
- Check that coordinate systems match between comparison and VEP files

**No results**
- Ensure input files have overlapping variants
- Check that bcftools coordinates are not NULL in comparison file

### Performance tips

- Use `--force` flag only when needed (ignores cache)
- Process QC analysis separately from prioritization
- Consider filtering input files to variants of interest
