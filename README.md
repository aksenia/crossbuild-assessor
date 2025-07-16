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

## Build 

```
docker build -t crossbuild-assessor .

```


## Key Features

### 🔍 **Liftover Quality Control**
Compare CrossMap and bcftools performance with detailed concordance analysis and visualization plots.

### 🧬 **Discrepant variant prioritization** 
Evidence-driven scoring system that identifies variants with discordant annotations between genome builds, prioritizing those most likely to impact clinical interpretation.

### 📊 **VEP integration**
Comprehensive analysis of Variant Effect Predictor annotations, tracking consequence changes and transcript-level discordances.

### 💾 **Efficient processing**
Memory-optimized design with caching handles large datasets without performance issues.

### 📋 **Interpretation-ready output**
Excel-compatible reports with enhanced genotype tracking, problematic transcript identification, and priority scoring.

## Pipeline Components

| Script | Purpose | Output |
|--------|---------|---------|
| `db_loader.py` | Load liftover + VEP data into SQLite | Optimized database |
| `db_analyzer.py` | Liftover quality control analysis | QC plots & reports |
| `variant_prioritizer.py` | Clinical variant prioritization | Ranked Excel reports |

## Discrepancy scoring highlights

- **Smart filtering**: 90% score reduction for benign variants, 2x boost for pathogenic
- **Impact-weighted**: Follows VEP consequence hierarchy (HIGH > MODERATE > LOW > MODIFIER)  
- **Evidence-driven**: Prioritizes clinical significance changes and pathogenicity predictions
- **Transcript-focused**: Identifies same-transcript consequence changes (highest priority)

