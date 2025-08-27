# Variant discrepancy scoring

Priority transcript-based scoring system for prioritizing annotation discrepancies between genome builds.

## Scoring philosophy

- HGVS concordance first: HGVS nomenclature mismatches on priority transcript drive prioritization.
- MANE-first transcript selection: Uses clinical standard transcript hierarchy:
  MANE Select → MANE Plus Clinical → Canonical → First available
- Same transcript compared between builds for consistency
- Clinical actionability: VUS and missing data de-prioritized (less clinically actionable).

## Priority categories

- CRITICAL: HGVS mismatches on priority transcript OR major clinical significance changes (Pathogenic↔Benign)
- MODERATE: Priority transcript unavailable OR moderate clinical changes OR serious functional differences
- LOW: Minor changes (VUS-Benign, prediction changes, gene symbol differences)
- CONCORDANT: Perfect priority transcript match with identical HGVS nomenclature
  
## Scoring components

HGVS Concordance (highest priority)

```text
Priority transcript HGVS mismatches:
  HGVSc mismatch:           +100 points (CRITICAL)
  HGVSp mismatch:           +50 points (CRITICAL)
  Missing HGVS data:        +5 points (minimal penalty)
Clinical Significance Changes
textMajor clinical changes (CRITICAL):
  Pathogenic ↔ Benign:      +90 points
  
Moderate changes (MODERATE):
  Pathogenic ↔ VUS:         +40 points
  
Minor changes (LOW):
  VUS ↔ Benign:            +25 points
  Likely ↔ Definite:       +20 points
  Missing clinical data:    +5 points
```

Transcript issues

```text
Priority transcript availability:
  MANE hg38 only:          +60 points (MODERATE)
  No matching transcripts:  +30 points (MODERATE)
  No transcript data:       +20 points (LOW)
```

Functional changes

```text
Consequence differences:
  Serious differences:      +35 points (MODERATE)
  Minor differences:        +20 points (LOW)
  
Annotation changes:
  Gene symbol changes:      +15 points (LOW)
  Impact level changes:     +15 points (LOW)
```

Technical issues

```text
Liftover problems:
  Position mismatch:        +20 points (LOW)
  Genotype mismatch:        +20 points (LOW)
  REF/ALT swap:            +10 points (LOW)
```

Prediction changes

```text
Pathogenicity predictions:
  SIFT changes:             +10 points (LOW)
  PolyPhen changes:         +10 points (LOW)
```

## Evidence override

- Benign suppression: 90% score reduction for LOW/MODIFIER impact + benign evidence
- Pathogenic boost: 2x score boost for HIGH impact or pathogenic evidence

## Bounded scoring

Multiple issues within one category cannot elevate to the next priority level, ensuring stable categorization.

## Example scoring

CRITICAL: HGVS Mismatch

```text
Variant: chr1:12345 A>G  
Priority transcript: NM_000059.4 (MANE Select)
HGVSc: c.1234G>A (hg19) vs c.1235G>A (hg38)
Score: +100 points → CRITICAL

```

MODERATE: Missing Priority Transcript

```text
Variant: chr2:67890 C>T
Transcript status: MANE_hg38_Only
Clinical: Pathogenic → VUS  
Score: +60 (transcript) + 40 (clinical) = 100 → CRITICAL

```

CONCORDANT: Perfect Match

```text
Variant: chr3:11111 T>C
Priority transcript: NM_012345.3 (MANE Select, both builds)
HGVSc: c.567T>C (identical)
Clinical: Pathogenic (stable)
Score: 0 points → CONCORDANT
```
