# BK Validation Analysis Summary

## Date: 2025-11-06

## Objective
Use BK (Belonocnema kinseyi, GCA_010883055.1) as ground truth to validate pipeline accuracy, since we have BK's own annotated genome and can check if the pipeline correctly recovers BK's target genes.

## Key Findings

### 1. Pipeline Completion Survey (32 gene families)
- **Phase 8 complete**: 17 families (53%)
- **Phase 7 complete**: 9 families (28%)
- **Phase 6 complete**: 1 family (3%)
- **Phase 3 complete**: 5 families (16%)

### 2. BK Target Recovery Rate
Across 17 families that completed Phase 8, analyzing 136 total BK loci:
- ✓ **Correct assignments**: 31 loci (22.8%)
- ↑ **More targets than expected**: 26 loci (19.1%)
- ↓ **Fewer targets than expected**: 3 loci (2.2%)
- ∅ **No targets assigned**: 76 loci (55.9%)

**Overall discrepancy rate: 77.2%**

### 3. Forensic Analysis Revealed Root Cause

**Phase 5 (synteny classification) is where BK targets are being lost.**

#### Evidence from MC117_BK_LB_batch (0/27 loci correct):

**Locus BK_chr5_c** (Expected: 1 target gene, Actual: 0):
- Phase 1: ✓ Locus defined (LOC117172861)
- Phase 4: ✓ **597 BK BLAST targets found** across all loci
- Phase 5: ✓ **0 BK targets classified as syntenic** for this locus ← THE PROBLEM
- Phase 8: Result = `[empty]`

**Working example** (natterin-4-like_MC59 - 5/5 correct):
- Phase 4: ✓ 34 BK targets found
- Phase 5: ✓ **1 BK target classified as syntenic** ← WORKS
- Phase 8: Result = 1 target ✓

**Pattern**: Of traced "no targets" loci:
- 0 lost targets in Phase 4 (BLAST working!)
- **6 lost targets in Phase 5** (classification failing!)
- 0 lost targets in Phase 8 (propagation correct)

### 4. Exonerate Refinement Attempt

Tested exonerate-based refinement to bypass Phase 5 classification:

**MC117 Results**:
- Original Phase 5: 43 targets, **3/27 loci** (11%)
- Exonerate refined: 312 targets, **27/27 loci** (100%) ✓

**Success**: Recovered all 27 loci
**Problem**: 312 targets vs 55 expected (567% overcalling)

### 5. Deep Dive: Why the Overcalling?

**Case study: BK_chr5_c**
- Phase 1 expected: 1 target gene
- Exonerate found: 12 targets
- Investigation revealed: Those 12 targets are **target genes from 12 different Phase 1 "loci"** (BK_chr5_a, c, d, e, f, i, m, p, s, t, w, {)

**All 12 are real MC117 family genes**, but they span **1.82 Mb** with:
- Average gap: **162.8 kb** between genes
- Largest gaps: **501 kb** and **800 kb**
- MC117 genes occupy only **1.4%** of the region

**This is a dispersed cluster, not a tight tandem array.**

### 6. Root Cause Identified

**Phase 1 correctly identified separate loci** (with 100-800 kb gaps between them).

**Phase 2/8 synteny blocks are too large** - creating mega-blocks spanning ~2 Mb that lump together what should be 12 separate loci.

**Phase 5 then assigns all targets** found in the mega-block to a single locus, rather than distributing them to the correct loci.

### 7. What We Tested

1. **Stricter Phase 2 clustering** (`max-gap-kb=100`): No change - still getting mega-blocks
2. **Stricter exonerate filtering** (score >75%, qcov >0.85): No significant improvement
3. **Various clustering parameters**: None prevented the mega-block problem

## Conclusion

**The discrepancy is in synteny block definition**, not in:
- Phase 1 locus calling (correctly identifies separate loci)
- Phase 4 BLAST (successfully finds targets)
- Phase 8 matrix generation (correctly propagates what Phase 5 provides)

**The issue**: Synteny blocks are spanning too large a region (>1 Mb), causing:
1. Multiple distinct loci to be grouped in one block
2. Phase 5 to incorrectly assign all targets to one locus
3. Other loci to be left empty

**Next steps should focus on**: Improving synteny block boundary detection to break at large gaps (>100kb) even when flanking protein conservation suggests continuity.

## Analysis Scripts Created
- `outputs/analyze_completion.py` - Survey pipeline completion
- `outputs/analyze_bk_assignments.py` - Check Phase 8 BK recovery
- `outputs/trace_bk_loci.py` - Forensic trace through all phases
- `outputs/BK_DISCREPANCY_FINDINGS.md` - Detailed findings document

## Key Metrics Files
- `outputs/pipeline_completion_analysis.json`
- `outputs/bk_assignment_analysis.json`
- `outputs/forensic_trace_output.txt`
