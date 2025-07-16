# crossbuild-assessor

**Genomic variant analysis tool for evaluating liftover quality and identifying clinically relevant discordances between genome builds (hg19/hg38).**

## Quick overview

CrossBuild Assessor analyzes how genomic variants behave when lifted between reference genome builds. It identifies problematic variants that may affect clinical interpretation and provides quality control metrics for liftover operations.

**What it does:**
- Compares liftover tools (CrossMap vs bcftools) 
- Analyzes VEP annotations across genome builds
- Prioritizes variants needing clinical review
- Generates Excel-ready reports with actionable insights

**Who it's for:** Clinical genomics labs, researchers working with multi-build datasets, and teams transitioning between genome assemblies.

---

## Quick start

```bash
git clone https://github.com/your-org/crossbuild-assessor.git
cd crossbuild-assessor
docker build -t crossbuild-assessor .

# Run pipeline
docker run -v /path/to/data:/data crossbuild-assessor \
  python db_loader.py --config /data/config.json

docker run -v /path/to/data:/data crossbuild-assessor \
  python variant_prioritizer.py --db-path /data/genomic_analysis.db --output-dir /data/results/
```

## Pipeline components

| Script | Purpose | Output |
|--------|---------|---------|
| `db_loader.py` | Load liftover + VEP data into SQLite | Optimized database |
| `db_analyzer.py` | Liftover quality control analysis | QC plots & reports |
| `variant_prioritizer.py` | Clinical variant prioritization | Ranked Excel reports |

## Priority categories

- **CRITICAL**: Same transcript, different consequences
- **HIGH**: Clinically significant impact transitions OR clinical significance changes  
- **MODERATE**: Gene changes with clinical significance OR prediction changes
- **INVESTIGATE**: Unmatched consequences OR non-significant transitions
- **LOW**: Benign variants

## Documentation

- **[Installation](docs/installation.md)** - Setup and Docker usage
- **[User guide](docs/user-guide.md)** - Complete usage tutorial
- **[Configuration](docs/configuration.md)** - Scoring weights and settings
- **[Architecture](docs/architecture.md)** - Code structure overview
- **[Variant discrepancy scoring](docs/variant-discrepancy-scoring.md)** - Detailed scoring methodology

## Input requirements

- **Liftover comparison file**: Tab-separated CrossMap vs bcftools comparison
- **VEP annotation files**: Standard VEP output for both hg19 and hg38
- **Configuration file**: JSON format specifying file paths

See the [user guide](docs/user-guide.md) for detailed format requirements.

## Example output

Priority variants are ranked by clinical importance:

```
Rank  Chromosome  Position_hg19  Gene_hg19  Priority_Score  Priority_Category
1     1           69134          OR4F5      15.2           CRITICAL
2     1           865705         SAMD11     12.8           HIGH  
3     2           38513          FAM150A    8.4            MODERATE
```
