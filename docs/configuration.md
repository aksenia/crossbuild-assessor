# Configuration

Configuration options for CrossBuild Assessor.

## Basic configuration

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

## Command-line options

### Database loader
```bash
python db_loader.py --config config.json [--force] [--verbose]
```

### Variant prioritizer
```bash
python variant_prioritizer.py \
  --db-path database.db \
  --output-dir results/ \
  [--max-variants 10000] \
  [--min-score 1] \
  [--no-plots] \
  [--force]
```

### Database analyzer
```bash
python db_analyzer.py --db-path database.db --output-dir qc_results/
```

### Report generator
```bash
python report_generator.py --input-dir results/ --output report.html
```

## Scoring customization

Modify `config/scoring_config.py`:

```python
BASE_SCORES = {
    'clinical_sig_benign_to_pathogenic': 20,   # Clinical changes
    'same_transcript_consequence_changes': 6,   # Transcript changes
    'sift_change': 5,                          # Prediction changes
    # ... other parameters
}
```

## Visualization customization

Modify `config/visualization_config.py`:

```python
PLOT_COLORS = {
    'priority': ['#d62728', '#ff7f0e', '#2ca02c', '#1f77b4'],
    'clinical': ['#8c564b', '#e377c2', '#7f7f7f', '#bcbd22']
}
```