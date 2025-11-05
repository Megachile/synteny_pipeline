# PHASE 1 REDESIGN PLAN: tblastn-First Locus Discovery

**Date:** 2025-01-04
**Reason:** Fix coordinate mismatch where BRAKER3 GFF annotations place genes at wrong locations (600kb error for DR_CM021340.1_a)

---

## Problem Identified

**Current workflow:**
```
Phase 1: Read GFF annotation → defines locus at annotated coordinates
Phase 2: Extract flanking from GFF coordinates
Phase 4: tblastn finds ACTUAL target location (600kb away!)
Phase 5: "Target too far from synteny" → REJECTED
```

**Example failure:**
- DR_CM021340.1_a GFF says ferritin at ~44M
- Actual genome sequence has ferritin at ~45M
- 600kb annotation error!
- Synteny block at 44M has NO ferritin
- tblastn finds ferritin at 45M
- Phase 5 rejects as "too far from synteny"

**Root cause:** Building locus around annotation coordinates instead of actual sequence coordinates.

---

## Solution: tblastn-First Approach

**User insight:** "it's more important that the target hit is precise than the flanking genes since their proximity is fuzzy across the 75 genomes anyway"

**New workflow:**
```
Phase 1 NEW: tblastn on genome → PRECISE target hit
             ↓
Phase 2: Extract flanking genes around THAT location
         (use BRAKER3/GFF if available)
         ↓
Phase 4: Finds same location as Phase 1 ✓
         ↓
Phase 5: Target IS in synteny block ✓
```

---

## Redesigned Phase 1 Workflow

### Step 1: Run tblastn on Reference Genome

**Input:**
- Query protein (BK ferritin, etc.)
- Reference genome FASTA (DR, TR, LB, etc.)

**Process:**
- Run tblastn against genome
- Filter hits by e-value, coverage, identity
- Select best hit as "seed" for locus

**Output:**
- Target coordinates: scaffold, start, end, strand
- Save as locus_seed.tsv

**Key decision:** What if multiple good hits?
- Create multiple loci (e.g., DR_CM021340.1_a, DR_CM021340.1_b)
- User can filter later based on synteny quality

### Step 2: Define Flanking Window

**Input:**
- Target coordinates from Step 1

**Process:**
- Define window: target ± 500kb (configurable)
- This captures typical synteny block sizes
- May be empty in gene-sparse regions

**Output:**
- Window coordinates: [target_start - 500000, target_end + 500000]

**Rationale:** 500kb is generous enough to capture synteny context but not so large we get unrelated genes

### Step 3: Query BRAKER3 GFF for Genes in Window

**Input:**
- Window coordinates
- BRAKER3 GFF3 file
- Genome proteome FASTA

**Process:**
- Parse GFF for gene features in window
- Filter to protein-coding genes only
- Extract gene IDs, coordinates, strand
- Match gene IDs to proteins in proteome

**Key handling:**
- If gene overlaps target hit → mark as TARGET
- Split into upstream/downstream based on strand
- Handle strand correctly (- strand reverses upstream/downstream)

**Output:**
- Flanking genes list with:
  - gene_id
  - protein_id
  - coordinates
  - position (upstream/downstream)
  - distance from target

**Edge cases to handle:**
- No GFF available → save target only, document "no flanking context"
- Sparse annotations → accept what's there, minimum 0 flanking genes
- Gene models overlap target → exclude or mark specially

### Step 4: Extract Flanking Protein Sequences

**Input:**
- Flanking gene list from Step 3
- Genome proteome FASTA

**Process:**
- Look up protein sequences by ID
- Label by position: U1, U2, ..., Un (upstream), D1, D2, ..., Dm (downstream)
- Order by distance from target

**Output:**
- flanking_proteins.faa with labeled sequences

**Labeling scheme:**
```
>g15427.t1|g15427 U5 [upstream protein 5]
MSLVRQNFHDEC...

>g15430.t1|g15430 U4 [upstream protein 4]
QGSVTLCPIES...

[TARGET would be here]

>g15435.t1|g15435 D1 [downstream protein 1]
PADQDWALQGL...
```

### Step 5: Validate and Save Locus

**Input:**
- Target coordinates
- Flanking proteins
- Window metadata

**Process:**
- Check minimum requirements:
  ✓ Target hit exists with good e-value
  ✓ At least 1 flanking gene found (warn if 0)
  ✓ Protein sequences extracted successfully

- Save locus definition:
  - locus_id
  - gene_family
  - reference_genome
  - target_coords (scaffold, start, end, strand)
  - num_flanking_upstream
  - num_flanking_downstream
  - window_size

**Output:**
- locus_definitions.tsv (updated)
- locus_coordinates.tsv (NEW - target coords for validation)
- flanking_proteins.faa
- locus_summary.txt (diagnostic info)

---

## File Structure

```
outputs/ferritin_phase1/
├── locus_definitions.tsv          # Existing format
├── locus_coordinates.tsv          # NEW: target locations
├── DR_CM021340.1_a/
│   ├── target_seed.tsv            # tblastn hit details
│   ├── window_definition.tsv      # Window coords used
│   ├── flanking_genes.tsv         # Genes found in window
│   ├── flanking_proteins.faa      # Extracted sequences
│   └── locus_summary.txt          # Human-readable report
```

---

## New File Format: locus_coordinates.tsv

```tsv
locus_id	genome	scaffold	start	end	strand	target_evalue	num_flanking_up	num_flanking_down
DR_CM021340.1_a	GCA_030998225.1	CM021340.1_RagTag	45041637	45045491	+	3.5e-38	8	12
```

**Purpose:**
- Enables Phase 2 validation (check synteny around THESE coords)
- Phase 6 can verify extraction matches Phase 1 location
- Debugging/diagnostics

---

## Edge Case Handling

### Case 1: No BRAKER3 GFF Available
```
- Run tblastn (still works!)
- Save target coordinates
- Set num_flanking_up = 0, num_flanking_down = 0
- Document: "No annotation available, synteny detection limited"
- Can still proceed to Phase 4 (will find same target)
```

### Case 2: Sparse Annotations in Window
```
- Accept whatever genes are found (even if 0)
- Document actual counts
- Synteny detection will adapt (need fewer matches)
- Better than wrong coordinates!
```

### Case 3: Target Overlaps Annotated Gene
```
Option A: Exclude that gene from flanking (RECOMMENDED)
Option B: Mark as "overlaps target" but include

Recommendation: Option A (cleaner)
```

### Case 4: Multiple Strong tblastn Hits
```
- Create separate loci for each hit above threshold
- Name: DR_scaffold_a, DR_scaffold_b, etc.
- Let synteny quality determine which are real paralogs
```

---

## Parameters

```python
# tblastn parameters
EVALUE_THRESHOLD = 1e-10     # Stricter than current?
MIN_COVERAGE = 0.5           # Hit must cover 50% of query
MIN_IDENTITY = 40            # Minimum % identity

# Window parameters
FLANKING_WINDOW = 500000     # 500kb each side
MIN_FLANKING_GENES = 0       # Accept even if sparse

# Gene filtering
EXCLUDE_TARGET_OVERLAP = True  # Don't use genes overlapping target
```

---

## Advantages

1. **Target precision:** tblastn on genome = ground truth
2. **Consistency:** Phase 1 and Phase 4 use same method
3. **No coordinate mismatches:** Everything built around actual target location
4. **Graceful degradation:** Works even if GFF sparse/missing
5. **Validation-ready:** Saves coordinates for Phase 2/6 cross-checking
6. **No augustus needed:** Uses existing BRAKER3 annotations

---

## Validation Strategy

After Phase 1 redesign, validate by:

1. Check DR_CM021340.1_a has target at ~45M (not 44M)
2. Synteny block in Phase 2 should be at ~45M
3. Phase 4 tblastn should find same location
4. Phase 6 extraction should succeed
5. No more `[empty]` for reference genomes at their own loci!

---

## Implementation Checklist

### Core Logic:
- [ ] Write tblastn wrapper function
- [ ] Implement hit filtering (e-value, coverage, identity)
- [ ] Window definition function (target ± 500kb)
- [ ] GFF parser for genes in window
- [ ] Upstream/downstream assignment logic (strand-aware)
- [ ] Protein sequence extraction
- [ ] Position labeling (U1, U2, D1, D2, etc.)

### File I/O:
- [ ] Save target_seed.tsv
- [ ] Save locus_coordinates.tsv (NEW)
- [ ] Save flanking_genes.tsv
- [ ] Save flanking_proteins.faa
- [ ] Update locus_definitions.tsv
- [ ] Generate locus_summary.txt

### Testing:
- [ ] Test on DR_CM021340.1_a (known problem case)
- [ ] Test on BK loci (NCBI annotations)
- [ ] Test on genome with no GFF
- [ ] Test on genome with sparse annotations
- [ ] Verify coordinates match Phase 4 results

### Integration:
- [ ] Update Phase 2 to use locus_coordinates.tsv
- [ ] Update Phase 6 validation to check coordinate matches
- [ ] Update WORKFLOW_BLUEPRINT.md
- [ ] Document in README

---

## Expected Impact

### Fixes:
- ✅ DR_CM021340.1_a will have target at correct location (45M)
- ✅ No more "reference genome missing target at own locus"
- ✅ Phase 4 and Phase 1 will be consistent
- ✅ Phase 5 won't reject targets as "too far"
- ✅ Phase 6 extraction will succeed

### May Reveal:
- Some current loci might disappear (no valid tblastn hit)
- Some genomes might gain new loci (multiple paralogs found)
- True extent of annotation quality issues

---

## Philosophy

**Build locus around where the target ACTUALLY IS, not where the annotation THINKS it is.**

Target precision is critical. Flanking gene context is useful but secondary.
