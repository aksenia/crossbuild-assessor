# Variant discrepancy scoring

CrossBuild Assessor uses a clinical evidence-first scoring system to prioritize annotation discrepancies between genome builds that are most likely to impact clinical interpretation.

## Scoring philosophy

### Clinical evidence first
The system prioritizes **clinical interpretation changes** over **VEP annotation differences**:

- **Highest priority**: Clinical significance changes (PATHOGENIC↔BENIGN, VUS→PATHOGENIC)
- **High priority**: Impact transitions and same transcript consequence changes
- **Lower priority**: Annotation differences and transcript mismatches

### Magnitude-based impact scoring
Transitions between impact levels are weighted by clinical significance:

- **HIGH ↔ MODERATE**: Clinically significant (+15 points)
- **HIGH ↔ LOW/MODIFIER**: Highly significant (+18-20 points)  
- **MODERATE ↔ LOW/MODIFIER**: Significant (+10-12 points)
- **LOW ↔ MODIFIER**: Annotation noise (+1 point)

## Priority categories

**CRITICAL** - Clinical interpretation changes OR high impact transitions (immediate review)  
**HIGH** - Functionally significant changes (priority review)  
**MODERATE** - Prediction changes and clinically relevant gene changes (standard review)  
**LOW** - Technical issues and annotation differences (secondary review)  
**INVESTIGATE** - Unclear cases requiring further investigation

## Scoring components

### 1. Clinical evidence scoring (Highest priority)

```
Clinical significance transitions (hg19→hg38):
  BENIGN → PATHOGENIC:     +20 points (CRITICAL)
  PATHOGENIC → BENIGN:     +18 points (CRITICAL)
  VUS → PATHOGENIC:        +15 points (CRITICAL)
  Other clinical changes:  +8 points (HIGH)
```

### 2. Impact transition scoring

```
High impact transitions (CRITICAL):
  HIGH ↔ MODERATE:    +15 points
  HIGH ↔ LOW:         +18 points  
  HIGH ↔ MODIFIER:    +20 points

Moderate impact transitions (HIGH):
  MODERATE ↔ LOW:     +10 points
  MODERATE ↔ MODIFIER: +12 points

Low priority transitions:
  LOW ↔ MODIFIER:     +1 point (annotation noise)
```

### 3. Functional impact scoring 

```
Same transcript consequence changes: +6 points 
```

### 4. Prediction and Gene scoring

```
Pathogenicity predictions (MODERATE):
  SIFT changes:       +5 points
  PolyPhen changes:   +5 points

Gene changes (conditional on clinical relevance):
  HIGH impact context:     +4 points per change
  MODERATE impact context: +3 points per change
  LOW impact context:      +2 points per change
  Non-clinical context:    +0.1 points (minimal weight)
```

### 5. Annotation differences

```
Unmatched consequences:                    +1 point 
Same consequence, different transcripts:   +0.5 points 
```

### 6. Technical issues

```
Position mismatch:        +3 points
Genotype mismatch:        +3 points  
Large position difference: +3 points (>100bp)
Medium position difference: +2 points (>10bp)
REF/ALT swap:             +2 points
```

## Evidence override

**Benign variant suppression**: 90% score reduction for LOW/MODIFIER impact + benign evidence  
**Pathogenic variant boost**: 2x score boost for HIGH impact or pathogenic evidence

## Priority category logic

### CRITICAL (Immediate review)
**Triggers**:
- Clinical interpretation changes: PATHOGENIC↔BENIGN, VUS→PATHOGENIC  
- High impact transitions involving HIGH impact level

**Examples**:
```
chr1:12345 C>T
ClinVar: benign (hg19) → pathogenic (hg38)
Score: +20 (clinical change) × 2.0 (pathogenic boost) = 40 points
```

### HIGH (Priority review)
**Triggers**:
- Moderate impact transitions (MODERATE↔LOW/MODIFIER)
- Other clinical significance changes
- Same transcript consequence changes

**Examples**:
```
chr2:67890 A>G
Impact: MODERATE (hg19) → LOW (hg38)
Score: +10 (impact transition) = 10 points
```

### MODERATE (Standard review)
**Triggers**:
- SIFT/PolyPhen prediction changes
- Gene changes with clinical relevance

**Examples**:
```
chr3:11111 G>A
SIFT: tolerated (hg19) → deleterious (hg38)
Score: +5 (SIFT change) = 5 points
```

### LOW (Secondary review)
**Triggers**:
- Technical liftover issues
- Annotation differences between builds
- Benign variants (with 90% score reduction)

**Examples**:
```
chr4:22222 T>C
Unmatched consequences + MODIFIER impact + benign evidence
Score: +1 (unmatched) × 0.1 (benign suppression) = 0.1 points
```

### INVESTIGATE (Further review)
**Triggers**:
- Unclear cases with insufficient clinical context
- Mixed evidence patterns

## Scoring Examples

### CRITICAL: Clinical interpretation change
```
Variant: chr1:12345 A>G
ClinVar: benign (hg19) → pathogenic (hg38)
Impact: MODERATE (both builds)

Scoring:
- Clinical significance change: +20 points
- Pathogenic boost: 2x multiplier
- Total: 40 points → CRITICAL
```

### CRITICAL: High impact transition  
```
Variant: chr2:67890 C>T
Impact: HIGH (hg19) → MODIFIER (hg38)
ENST00000123456: stop_gained → intron_variant

Scoring:
- Impact transition (HIGH→MODIFIER): +20 points
- Same transcript change: +6 points
- Total: 26 points → CRITICAL
```

### HIGH: Functional change
```
Variant: chr3:11111 G>A
Impact: MODERATE (hg19) → LOW (hg38)
ENST00000789012: missense_variant → synonymous_variant

Scoring:
- Impact transition (MODERATE→LOW): +10 points
- Same transcript change: +6 points
- Total: 16 points → HIGH
```

### LOW: Annotation noise (suppressed)
```
Variant: chr4:22222 T>C
Impact: MODIFIER (both builds)
ClinVar: benign (both builds)
Different transcript annotations only

Scoring:
- Unmatched consequences: +1 point
- Benign suppression: 0.1x multiplier
- Total: 0.1 points → LOW
```