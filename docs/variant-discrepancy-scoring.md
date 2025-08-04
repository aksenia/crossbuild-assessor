# Variant discrepancy scoring

Clinical evidence-first scoring system for prioritizing annotation discrepancies between genome builds.

## Scoring philosophy

**Clinical evidence first**: Clinical interpretation changes are prioritized over VEP annotation differences.

**Magnitude-based impact**: Transitions weighted by clinical significance:

- HIGH ↔ MODERATE: Clinically significant (+15 points)
- HIGH ↔ LOW/MODIFIER: Highly significant (+18-20 points)  
- MODERATE ↔ LOW/MODIFIER: Significant (+10-12 points)
- LOW ↔ MODIFIER: Annotation noise (+1 point)

## Priority categories

**CRITICAL**: Clinical interpretation changes OR high impact transitions  
**HIGH**: Functionally significant changes  
**MODERATE**: Prediction changes and clinically relevant gene changes  
**LOW**: Technical issues and annotation differences

## Scoring components

### Clinical evidence (Highest priority)

```text
Clinical significance transitions (hg19→hg38):
  BENIGN → PATHOGENIC:     +20 points (CRITICAL)
  PATHOGENIC → BENIGN:     +18 points (CRITICAL)
  VUS → PATHOGENIC:        +15 points (CRITICAL)
  Other clinical changes:  +8 points (HIGH)
```

### Impact transitions

```text
High impact transitions (CRITICAL):
  HIGH ↔ MODERATE:    +15 points
  HIGH ↔ LOW:         +18 points  
  HIGH ↔ MODIFIER:    +20 points

Moderate transitions (HIGH):
  MODERATE ↔ LOW:     +10 points
  MODERATE ↔ MODIFIER: +12 points

Low priority:
  LOW ↔ MODIFIER:     +1 point
```

### Functional impact

```text
Same transcript consequence changes: +6 points 
```

### Predictions and genes

```text
SIFT/PolyPhen changes: +5 points each
Gene changes: +0.1-4 points (conditional on clinical relevance)
```

### Technical issues

```text
Position mismatch:    +3 points
Genotype mismatch:    +3 points  
REF/ALT swap:         +2 points
```

## Evidence override

**Benign suppression**: 90% score reduction for LOW/MODIFIER impact + benign evidence  
**Pathogenic boost**: 2x score boost for HIGH impact or pathogenic evidence

## Example scoring

### CRITICAL: Clinical change

```text
Variant: chr1:12345 A>G
ClinVar: benign (hg19) → pathogenic (hg38)
Score: +20 × 2.0 (pathogenic boost) = 40 points → CRITICAL
```

### HIGH: Impact transition

```text
Variant: chr2:67890 C>T
Impact: HIGH (hg19) → MODIFIER (hg38)
Score: +20 (impact) + 6 (transcript) = 26 points → HIGH
```

### LOW: Annotation noise

```text
Variant: chr4:22222 T>C
Impact: MODIFIER (both), ClinVar: benign (both)
Score: +1 × 0.1 (benign suppression) = 0.1 points → LOW
```
