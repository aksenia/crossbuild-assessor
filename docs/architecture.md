# Architecture

Modular architecture with clear separation of concerns.

## System overview

```text
CLI scripts (db_loader.py, variant_prioritizer.py, db_analyzer.py, report_generator.py)
                                      │
Core modules (analysis/, visualization/, utils/, config/)
                                      │
Data layer (SQLite database, cache files, output files)
```

## Module structure

### Configuration (`config/`)

```text
constants.py              # VEP consequence → impact mappings  
scoring_config.py         # Clinical scoring weights & rules
visualization_config.py   # Plot colors & layouts
```

### Analysis engine (`analysis/`)

```text
variant_processor.py     # Main orchestrator
vep_analyzer.py          # VEP annotation processing
scoring_engine.py        # Clinical evidence-driven scoring  
cache_manager.py         # Result caching system
```

### Utilities (`utils/`)

```text
clinical_utils.py       # Clinical significance parsing
transcript_utils.py     # Transcript ID normalization
data_utils.py           # Data transformation helpers
impact_utils.py         # Impact level calculations
```

### Visualization (`visualization/`)

```text
plot_generator.py        # 4-plot clinical dashboard
```

## Data flow

```text
Input files → db_loader.py → SQLite database
                                ↓
Database → VEPAnalyzer → Analysis results → CacheManager → Cache
                                ↓
Cache → ClinicalScorer → Scored variants → PrioritizationPlotter → Plots
```

## Design principles

**Separation of concerns**: Each module has single responsibility  
**Caching strategy**: VEP analysis cached, scores calculated fresh  
**Memory efficiency**: Chunked processing for large datasets  
**Clinical focus**: Evidence-first prioritization over annotation noise

## Key components

**VEPAnalyzer**: Processes VEP annotations, identifies transcript discordances  
**ClinicalScorer**: Evidence-driven scoring with clinical override  
**CacheManager**: Intelligent caching for performance  
**PrioritizationPlotter**: Clinical review visualizations