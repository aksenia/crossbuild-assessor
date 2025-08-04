# crossbuild-assessor

**Tool for analyzing genomic variant annotation discordances between genome builds (hg19/hg38).**

## Quick overview

CrossBuild Assessor identifies variants that may affect clinical interpretation when lifted between reference genomes. the tool prioritizes clinical evidence changes over annotation differences.

**What it does:**
- Compares liftover tools (CrossMap vs bcftools) 
- Analyzes VEP annotations across genome builds
- Prioritizes variants by clinical impact
- Generates Excel reports and HTML summaries

---

## Quick start

```bash
git clone https://github.com/your-org/crossbuild-assessor.git
cd crossbuild-assessor
docker build -t crossbuild-assessor .

# 1. Load data
docker run -v /path/to/data:/data crossbuild-assessor \
  python db_loader.py --config /data/config.json

# 2. Quality control analysis  
docker run -v /path/to/data:/data crossbuild-assessor \
  python db_analyzer.py --db-path /data/genomic_analysis.db --output-dir /data/qc/

# 3. Variant prioritization
docker run -v /path/to/data:/data crossbuild-assessor \
  python variant_prioritizer.py --db-path /data/genomic_analysis.db --output-dir /data/results/

# 4. Generate HTML report
docker run -v /path/to/data:/data crossbuild-assessor \
  python report_generator.py --input-dir /data/results/ --output /data/crossbuild_report.html
```

## Pipeline components

| Script | Purpose | Key Output |
|--------|---------|------------|
| `db_loader.py` | Load liftover + VEP data | SQLite database |
| `db_analyzer.py` | Liftover quality control | Concordance analysis |
| `variant_prioritizer.py` | Clinical prioritization | Ranked CSV + plots |
| `report_generator.py` | HTML summary report | Clinical dashboard |

## Priority categories

- **CRITICAL**: Clinical interpretation changes (PATHOGENIC↔BENIGN, VUS→PATHOGENIC) 
- **HIGH**: Impact transitions, same transcript consequence changes
- **MODERATE**: Pathogenicity prediction changes (SIFT/PolyPhen)
- **LOW**: Technical issues, annotation differences

## Configuration

Create `config.json`:

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

## Input requirements

- **Liftover comparison**: Tab-separated CrossMap vs bcftools results
- **VEP annotations**: Standard VEP output for both builds
- **Memory**: 8GB+ RAM for large datasets

## Example output

Clinical prioritization CSV:

```text
Rank  Gene_hg19  Clinical_Significance_hg19  Clinical_Significance_hg38  Priority_Score  Priority_Category
1     BRCA1      PATHOGENIC                   BENIGN                      40000          CRITICAL
2     TP53       VUS                          PATHOGENIC                  30000          CRITICAL  
3     APOE       HIGH impact                  MODERATE impact             12000          HIGH
```

HTML report provides unified dashboard with embedded visualizations and clinical analysis.