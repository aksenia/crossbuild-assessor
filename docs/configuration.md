# Configuration

Configuration options for CrossBuild Assessor.

## Basic configuration

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

## Advanced configuration

### Scoring weights

Modify scoring in `config/scoring_config.py`:

```python
BASE_SCORES = {
    'same_transcript_consequence_changes': 10,  # CRITICAL variants
    'clinical_sig_benign_to_pathogenic': 10,   # Clinical changes
    'gene_changes_high_impact': 8,             # High impact gene changes
    # ... other parameters
}
```

### Clinical override

Adjust clinical evidence modifiers:

```python
CLINICAL_OVERRIDE = {
    'benign_reduction_factor': 0.1,   # 90% reduction for benign
    'pathogenic_boost_factor': 2.0    # 2x boost for pathogenic
}
```

### Visualization

Customize colors in `config/vizualisation_config.py`:

```python
PLOT_COLORS = {
    'priority': ['#d62728', '#ff7f0e', '#2ca02c', '#1f77b4', '#9467bd'],
    'clinical': ['#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
}
```

## Performance tuning

For large datasets, adjust memory settings in the code:

```python
chunk_size = 5000  # Reduce from default 10000 for limited memory
```
