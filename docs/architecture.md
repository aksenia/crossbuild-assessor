# Architecture

CrossBuild Assessor uses a modular architecture with clear separation of concerns.

## System overview

```
CLI Scripts (db_loader.py, variant_prioritizer.py, db_analyzer.py)
                            │
Core Modules (analysis/, visualization/, utils/, config/)
                            │
Data Layer (SQLite Database, Cache Files, Output Files)
```

## Module structure

### Configuration (`config/`)
Centralized settings for scoring, visualization, and VEP mappings.

```
config/
├── constants.py              # VEP consequence → impact mappings
├── scoring_config.py         # Clinical scoring weights & rules
└── vizualisation_config.py   # Plot colors & layouts
```

### Utilities (`utils/`)
Pure functions with no dependencies or side effects.

```
utils/
├── clinical_utils.py     # Clinical significance & pathogenicity parsing
├── transcript_utils.py   # Transcript ID normalization & genotype extraction  
├── impact_utils.py       # Impact calculations & transitions
└── data_utils.py         # Data transformation helpers
```

### Analysis engine (`analysis/`)
Core variant analysis logic with VEP processing and clinical scoring.

```
analysis/
├── variant_processor.py     # Main orchestrator
├── vep_analyzer.py         # VEP annotation processing
├── scoring_engine.py       # Clinical evidence-driven scoring
└── cache_manager.py        # Intelligent caching system
```

### Visualization (`visualization/`)
Professional 4-plot dashboard for clinical review.

```
visualization/
└── plot_generator.py        # PrioritizationPlotter class
```

## Data flow

```
Input Files → db_loader.py → SQLite Database
                                ↓
Database → VEPAnalyzer → Raw Analysis → CacheManager → Cache
                                ↓
Cache → ClinicalScorer → Scored Variants → PrioritizationPlotter → Plots
```

## Design principles

### Separation of concerns
Each module has a single responsibility:
- **Analysis**: VEP processing and scoring
- **Visualization**: Plot generation  
- **Configuration**: Parameter management
- **Utilities**: Reusable functions

### Caching strategy
- **VEP analysis**: Cached (expensive to compute)
- **Scores**: Calculated fresh (easy to recalibrate)
- **Performance**: Dramatically faster subsequent runs

### Memory efficiency
- **Chunked processing**: Handle large datasets safely
- **Targeted queries**: Load only needed data
- **Clean-up**: Free memory after each chunk

## Key components

### VEPAnalyzer
Processes VEP annotations and identifies transcript-level discordances:
- Matches transcripts across builds using normalized IDs
- Identifies consequence changes within same transcripts
- Extracts clinical data (SIFT, PolyPhen, ClinVar)

### ClinicalScorer
Applies evidence-driven scoring with clinical focus:
- Magnitude-based impact scoring (HIGH→MODERATE more important than LOW→MODIFIER)
- Clinical evidence override (benign suppression, pathogenic boost)
- Conditional gene scoring (reduced weight for non-clinical changes)

### PrioritizationPlotter
Creates publication-ready 4-plot dashboard:
- Priority categories with concordant variants
- Clinical evidence with intelligent gap-filling
- Primary discordance types (aggregated)
- Clinical evidence transitions

This modular design makes the codebase maintainable, testable, and extensible while keeping the main scripts focused and readable., layouts
```

**Key Components**:
- **VEP Consequence Impact Mappings**: Official Ensembl hierarchy
- **Clinical Scoring Weights**: Evidence-driven prioritization rules
- **Visualization Settings**: Colorblind-friendly palettes, plot layouts

**Design Principles**:
- **Configurable**: Easy parameter tuning without code changes
- **Centralized**: Single source of truth for all constants
- **Hierarchical**: Logical grouping by functionality

### 2. Utilities Layer (`utils/`)

**Purpose**: Reusable pure functions with no side effects or dependencies.

```python
utils/
├── __init__.py           # Module exports
├── clinical_utils.py     # Clinical significance & pathogenicity parsing
├── transcript_utils.py   # Transcript ID normalization & genotype extraction
├── impact_utils.py       # Impact level calculations & transitions
└── data_utils.py         # Data transformation helpers
```

**Design Principles**:
- **Pure Functions**: No side effects, deterministic outputs
- **Single Responsibility**: Each function has one clear purpose
- **Testable**: Easy to unit test in isolation
- **Reusable**: Can be used across different modules

### 3. Analysis Engine (`analysis/`)

**Purpose**: Core variant analysis logic with sophisticated VEP processing and clinical scoring.

```python
analysis/
├── __init__.py              # Module exports
├── variant_processor.py     # Main orchestrator
├── vep_analyzer.py         # VEP annotation processing (~200 lines)
├── scoring_engine.py       # Clinical evidence-driven scoring (~150 lines)
└── cache_manager.py        # Intelligent caching system
```

**Component Details**:

#### VEPAnalyzer
- **Transcript Analysis**: Matches transcripts across builds using normalized IDs
- **Consequence Comparison**: Identifies same-transcript consequence changes
- **Clinical Data Extraction**: Parses SIFT, PolyPhen, ClinVar annotations
- **Memory Optimization**: Chunked processing for large datasets

#### ClinicalScorer  
- **Magnitude-Based Scoring**: Clinically significant impact transitions
- **Evidence-Driven Prioritization**: Clinical significance changes prioritized
- **Benign Suppression**: 90% score reduction for benign variants
- **Pathogenic Boost**: 2x score boost for pathogenic evidence

#### CacheManager
- **Intelligent Caching**: Stores VEP analysis, calculates scores fresh
- **Performance Optimization**: Dramatically faster subsequent runs
- **Flexible Recalibration**: Easy to tune scoring without re-analysis

### 4. Visualization System (`visualization/`)

**Purpose**: Professional plotting system for clinical review and quality assessment.

```python
visualization/
├── __init__.py              # Module exports
└── plot_generator.py        # PrioritizationPlotter class (~200 lines)
```

**Features**:
- **4-Plot Dashboard**: Priority categories, clinical evidence, discordance types, transitions
- **Clinical Focus**: Intelligent gap-filling for clinical evidence categorization
- **Professional Quality**: Publication-ready plots with proper statistics
- **Configurable**: Colors, styles, and layouts from config module

## Data Flow Architecture

### 1. Data Loading Pipeline

```
Liftover Files + VEP Files → db_loader.py → SQLite Database
                                ↓
                    Coordinate Normalization
                           ↓
                    Optimized Indexes
```

### 2. Analysis Pipeline

```
SQLite Database → VEPAnalyzer → Raw Analysis Results
                       ↓
                  CacheManager → Cached VEP Analysis
                       ↓
                ClinicalScorer → Scored Variants
                       ↓
              PrioritizationPlotter → Visualization
```

### 3. Caching Strategy

```
First Run:  Database → VEP Analysis → Cache → Scoring → Results
                                        ↓
Subsequent: Cache → Scoring → Results (Much Faster)
                    ↓
Parameter Tuning: Same Cache → New Scoring → Updated Results
```

## Design Patterns

### 1. Separation of Concerns
- **Analysis** vs **Visualization** vs **Configuration**
- Each module has a single, well-defined responsibility
- Clear interfaces between components

### 2. Dependency Injection
```python
# Components receive dependencies, don't create them
scorer = ClinicalScorer()  # Uses config internally
plotter = PrioritizationPlotter(PLOT_COLORS, FIGURE_CONFIG)
```

### 3. Factory Pattern
```python
# VariantProcessor orchestrates component creation
processor = VariantProcessor()
processor.process_all_variants(conn, cache_file, force_recalculate)
```

### 4. Caching Pattern
```python
# Intelligent cache management with fallback
if cache_manager.should_use_cache(force_recalculate):
    return cache_manager.load_cache()
else:
    results = analyzer.analyze_all_variants(conn)
    cache_manager.save_cache(results)
    return results
```

## Memory Management

### Chunked Processing
```python
# Process variants in memory-safe chunks
chunk_size = 10000
for chunk_idx in range(total_chunks):
    chunk_variants = variants_df.iloc[start_idx:end_idx]
    chunk_analyses = self._process_variant_chunk(conn, chunk_variants)
    # Memory cleanup after each chunk
```

### Efficient Database Queries
```python
# Targeted queries instead of loading all data
hg19_query = "SELECT ... FROM hg19_vep WHERE chr = ? AND pos = ?"
hg19_annotations = pd.read_sql_query(hg19_query, conn, params=[chrom, pos])
```

## Extensibility Points

### 1. New Scoring Algorithms
Add new scoring methods to `ClinicalScorer` without changing other components.

### 2. Additional Analysis Types
Create new analyzer classes implementing the same interface as `VEPAnalyzer`.

### 3. Custom Visualizations
Extend `PrioritizationPlotter` or create new plotting classes.

### 4. Different Data Sources
Modify `VEPAnalyzer` to read from different databases or file formats.

## Testing Strategy

### Unit Testing
```python
# Each module can be tested independently
def test_clinical_scorer():
    scorer = ClinicalScorer()
    mock_analysis = create_mock_vep_analysis()
    results = scorer.calculate_scores_from_analysis(mock_analysis)
    assert results['priority_score'].max() > 0
```

### Integration Testing
```python
# Test component interaction
def test_variant_processor():
    processor = VariantProcessor()
    results = processor.process_all_variants(test_conn, test_cache)
    assert len(results) > 0
```

### Performance Testing
```python
# Validate memory efficiency and speed
def test_large_dataset_processing():
    # Test with 100K+ variants
    processor = VariantProcessor()
    start_time = time.time()
    results = processor.process_all_variants(large_dataset_conn)
    processing_time = time.time() - start_time
    assert processing_time < MAX_ALLOWED_TIME
```

## Configuration Management

### Environment-Specific Settings
```python
# Development vs Production configurations
if ENVIRONMENT == "development":
    CHUNK_SIZE = 1000  # Smaller chunks for testing
else:
    CHUNK_SIZE = 10000  # Optimized for production
```

### Parameter Validation
```python
# Validate configuration on startup
def validate_scoring_config():
    assert all(score > 0 for score in BASE_SCORES.values())
    assert CLINICAL_OVERRIDE['benign_reduction_factor'] < 1.0
```

This architecture provides a robust foundation for genomic variant analysis while maintaining flexibility for future enhancements and clinical workflow integration.
