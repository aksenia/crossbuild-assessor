# Variant discrepancy scoring

CrossBuild Assessor uses an evidence-driven scoring system to prioritize annotation discrepancies between genome builds that are most likely to impact clinical interpretation.

## Scoring philosophy

### Clinical evidence over annotation differences
The system prioritizes **functional clinical impact** over **annotation model changes**:

- **High priority**: Same transcript with different consequences (functional impact)
- **Medium priority**: Clinical significance changes (interpretation impact)  
- **Low priority**: Gene symbol differences (often just synonyms)

### Magnitude-based impact scoring
Transitions between impact levels are weighted by clinical significance:

- **HIGH ↔ MODERATE**: Clinically significant (+10 points)
- **MODERATE ↔ LOW**: Moderately significant (+8 points)  
- **LOW ↔ MODIFIER**: Annotation noise (+1 point)

## Priority categories

**CRITICAL** - Same transcript, different consequences (immediate review)  
**HIGH** - Clinical significance changes or major impact transitions (priority review)  
**MODERATE** - Gene changes or pathogenicity prediction changes (standard review)  
**INVESTIGATE** - Unmatched consequences or minor transitions (secondary review)  
**LOW** - Benign variants with minimal impact (informational)

## Scoring components

### 1. Functional impact scoring

```
Same transcript consequence changes: +10 points (CRITICAL)

Impact transitions:
  HIGH ↔ MODERATE:    +10 points
  HIGH ↔ LOW:         +12 points  
  HIGH ↔ MODIFIER:    +15 points
  MODERATE ↔ LOW:     +8 points
  MODERATE ↔ MODIFIER: +10 points
  LOW ↔ MODIFIER:     +1 point (annotation noise)
```

### 2. Clinical evidence scoring

```
Clinical significance changes:
  Benign → Pathogenic:     +10 points
  Pathogenic → Benign:     +8 points
  Other clinical changes:  +5 points

Pathogenicity predictions:
  SIFT changes:       +5 points
  PolyPhen changes:   +5 points
```

### 3. Gene-level scoring

Gene symbol changes are scored conditionally based on clinical relevance:

- **Clinically relevant**: HIGH/MODERATE impact variants get full gene change scoring
- **Non-clinical**: LOW/MODIFIER variants get minimal weight (+0.1 points)

**Rationale**: Many gene symbol differences are synonyms and don't represent functional changes.

### 4. Evidence override

**Benign variant suppression**: 90% score reduction for LOW/MODIFIER impact + benign evidence  
**Pathogenic variant boost**: 2x score boost for HIGH impact or pathogenic evidence

## Scoring examples

### CRITICAL variant
```
chr1:12345 A>G
ENST00000123456: missense_variant (hg19) → stop_gained (hg38)
Score: +10 (consequence change) +10 (impact transition) = 20 points
```

### HIGH priority variant
```
chr2:67890 C>T  
ClinVar: benign (hg19) → pathogenic (hg38)
Score: +10 (clinical change) × 2.0 (pathogenic boost) = 20 points
```

### LOW priority variant (suppressed)
```
chr3:11111 G>A
MODIFIER impact + ClinVar benign
Score: +2 (minor changes) × 0.1 (benign suppression) = 0.2 points
```

This approach ensures focus on discrepancies most likely to impact patient care while filtering out annotation noise.# Clinical Scoring Methodology

CrossBuild Assessor uses a sophisticated, evidence-driven scoring system designed to prioritize variants most likely to impact clinical interpretation when transitioning between genome builds.

## Scoring Philosophy

### Clinical Evidence Over Annotation Differences
The scoring system prioritizes **functional clinical impact** over **annotation model changes**:

- **High Priority**: Same transcript with different consequences (functional impact)
- **Medium Priority**: Clinical significance changes (interpretation impact)  
- **Low Priority**: Gene symbol differences (often just synonyms)

### Magnitude-Based Impact Scoring
Transitions between impact levels are weighted by **clinical significance**:

- **HIGH ↔ MODERATE**: Clinically significant (score: +10)
- **MODERATE ↔ LOW**: Moderately significant (score: +8)  
- **LOW ↔ MODIFIER**: Annotation noise (score: +1)

## Priority Categories

### CRITICAL (Immediate Review Required)
**Definition**: Same transcript ID, different consequences between builds

**Clinical Significance**: These represent potential functional changes that could directly impact protein products.

**Examples**:
- ENST00000123456: `missense_variant` (hg19) → `stop_gained` (hg38)
- ENST00000789012: `synonymous_variant` (hg19) → `splice_region_variant` (hg38)

**Scoring**: Base score +10 points per occurrence

### HIGH (Priority Clinical Review)
**Definition**: Clinically significant impact transitions OR clinical significance changes

**Triggers**:
- Impact transitions involving HIGH or MODERATE levels
- Clinical significance changes: Benign ↔ Pathogenic

**Examples**:
- Impact change: HIGH → MODERATE (+10 points)
- ClinVar: `benign` → `pathogenic` (+10 points)
- Impact change: HIGH → LOW (+12 points)

### MODERATE (Standard Review)
**Definition**: Gene changes with clinical significance OR pathogenicity prediction changes

**Triggers**:
- Gene symbol changes with HIGH/MODERATE impact
- SIFT prediction changes: `tolerated` ↔ `deleterious`
- PolyPhen changes: `benign` ↔ `probably_damaging`

**Examples**:
- Gene change from BRCA1 to BRCA1P1 with MODERATE impact
- SIFT change: `tolerated` → `deleterious` (+5 points)

### INVESTIGATE (Secondary Review)
**Definition**: Unmatched consequences OR non-significant transitions

**Triggers**:
- Consequences present in one build but not the other
- LOW ↔ MODIFIER transitions (annotation updates)

**Examples**:
- Consequence in hg19 but missing in hg38
- Impact change: LOW → MODIFIER (+1 point)

### LOW (Informational)
**Definition**: Benign variants with minimal clinical impact

**Triggers**:
- LOW/MODIFIER impact + benign evidence
- 90% score reduction applied

**Examples**:
- MODIFIER impact + ClinVar `benign`
- LOW impact + SIFT `tolerated` + PolyPhen `benign`

## Scoring Components

### 1. Functional Impact Scoring

#### Same Transcript Consequence Changes (CRITICAL)
```
Score: +10 points per occurrence
Weight: Highest priority
```

**Rationale**: Same transcript with different consequences suggests functional changes that could directly impact protein products.

#### Impact Magnitude Transitions
```
HIGH ↔ MODERATE:    +10 points (clinically significant)
HIGH ↔ LOW:         +12 points (highly significant)  
HIGH ↔ MODIFIER:    +15 points (most significant)
MODERATE ↔ LOW:     +8 points (significant)
MODERATE ↔ MODIFIER: +10 points (significant)
LOW ↔ MODIFIER:     +1 point (annotation noise)
```

**Rationale**: Transitions involving HIGH or MODERATE impact are clinically relevant, while LOW ↔ MODIFIER transitions often represent annotation model updates rather than functional changes.

### 2. Clinical Evidence Scoring

#### Clinical Significance Changes
```
Benign → Pathogenic:     +10 points (critical for interpretation)
Pathogenic → Benign:     +8 points (important for interpretation)
Other clinical changes:  +5 points (noteworthy)
```

#### Pathogenicity Predictions
```
SIFT changes:       +5 points (deleterious ↔ tolerated)
PolyPhen changes:   +5 points (damaging ↔ benign)
```

### 3. Gene-Level Scoring (Conditional)

Gene symbol changes are scored **conditionally** based on clinical relevance:

#### Clinically Relevant Gene Changes
When variant has HIGH/MODERATE impact OR clinical significance changes:
```
HIGH impact genes:      +8 points per change
MODERATE impact genes:  +4 points per change  
LOW impact genes:       +2 points per change
Mixed impact:           +3 points per change
```

#### Non-Clinically Relevant Gene Changes
When variant has only LOW/MODIFIER impact with no clinical evidence:
```
Gene symbol differences: +0.1 points (minimal weight)
```

**Rationale**: Many gene symbol differences are synonyms (e.g., FAM123A vs C1orf456) and don't represent functional changes.

### 4. Technical Quality Scoring

#### Position and Genotype Issues
```
Position mismatch:        +3 points
Genotype mismatch:        +3 points  
Large position difference: +3 points (>100bp)
Medium position difference: +2 points (>10bp)
REF/ALT swap:             +2 points
```

#### Transcript Analysis Scoring
```
Unmatched consequences:                    +4 points
Same consequence, different transcripts:   +3 points
```

### 5. Clinical Evidence Override

#### Benign Variant Suppression
```
Conditions: LOW/MODIFIER impact + benign evidence
Effect: 90% score reduction (multiply by 0.1)
Benign evidence: ClinVar benign, SIFT tolerated, PolyPhen benign
```

**Rationale**: Reduces noise from annotation updates on variants with minimal clinical impact.

#### Pathogenic Variant Boost  
```
Conditions: HIGH impact OR pathogenic evidence
Effect: 2x score boost (multiply by 2.0)
Pathogenic evidence: ClinVar pathogenic, HIGH impact, SIFT deleterious, PolyPhen damaging
```

**Rationale**: Ensures clinically important variants receive attention even with minor annotation differences.

## Scoring Examples

### Example 1: CRITICAL Variant
```
Variant: chr1:12345 A>G
Transcript: ENST00000123456
hg19: missense_variant (MODERATE)
hg38: stop_gained (HIGH)

Scoring:
- Same transcript consequence change: +10 points
- Impact transition (MODERATE→HIGH): +10 points  
- HIGH impact bonus: +2 points
- Total: 22 points → CRITICAL
```

### Example 2: HIGH Priority Variant
```
Variant: chr2:67890 C>T  
ClinVar: benign (hg19) → pathogenic (hg38)
Impact: MODERATE (both builds)

Scoring:
- Clinical significance change: +10 points
- Has clinical data: +2 points
- Pathogenic boost: 2x multiplier
- Total: 24 points → HIGH
```

### Example 3: LOW Priority Variant (Suppressed)
```
Variant: chr3:11111 G>A
Impact: MODIFIER (both builds)  
ClinVar: benign (both builds)
Gene: GENE1 (hg19) → GENE1P (hg38)

Scoring:
- Gene change (non-clinical): +0.1 points
- Has clinical data: +2 points
- Benign suppression: 0.1x multiplier
- Total: 0.21 points → LOW
```

### Example 4: MODERATE Priority Variant
```
Variant: chr4:22222 T>C
Impact: HIGH (both builds)
SIFT: tolerated (hg19) → deleterious (hg38)

Scoring:  
- SIFT change: +5 points
- HIGH impact bonus: +2 points
- Has clinical data: +2 points
- Pathogenic boost: 2x multiplier
- Total: 18 points → MODERATE
```

## Clinical Workflow Integration

### Triage Recommendations

1. **CRITICAL variants**: Immediate manual review required
   - Same transcript with different consequences
   - Potential functional impact on protein products

2. **HIGH priority variants**: Priority review within clinical workflow
   - Clinical significance changes affecting interpretation
   - Significant impact level transitions

3. **MODERATE priority variants**: Standard review process
   - Pathogenicity prediction changes
   - Clinically relevant gene changes

4. **INVESTIGATE variants**: Secondary review or automated filtering
   - Unmatched consequences requiring investigation
   - Non-significant transitions

5. **LOW priority variants**: Informational only
   - Benign variants with minimal changes
   - Annotation updates with no clinical impact

### Quality Metrics

The scoring system provides quality metrics for liftover assessment:

- **Same transcript issues**: Indicate potential reference genome problems
- **Clinical significance changes**: Measure interpretation impact
- **Benign variant noise**: Quantifies annotation-only changes
- **Pathogenic variant coverage**: Ensures clinical variants are prioritized

This evidence-driven approach ensures clinical geneticists focus their time on variants most likely to impact patient care while filtering out noise from annotation model updates.
