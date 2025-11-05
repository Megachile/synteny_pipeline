# Phase 1 tblastn-First Implementation Status

**Date:** 2025-11-05
**Status:** ✅ WORKING - tblastn-first successfully finds DR_CM021340.1 locus!

---

## Objective

Test **tblastn-first locus discovery for DR/TR (BRAKER3 genomes) ONLY** to fix spurious loci problem where DIAMOND found proteins but GFF coordinates were wrong (600kb error for DR_CM021340.1_a).

**Critical:** BK/LB (NCBI genomes) continue using old DIAMOND blastp approach (ground truth, don't touch).

---

## Problem Being Solved

### Old DIAMOND Approach:
- Used blastp against proteomes → found protein hits
- Built loci around GFF-annotated coordinates
- **Issue:** BRAKER3 GFF coordinates can be 600kb off from actual sequence location
- **Result:** Created spurious loci (e.g., DR_CM021340.1_a had synteny but no target found in Phase 4)

### New tblastn Approach (DR/TR only):
- Run tblastn on genome sequence FIRST → finds actual target location
- Build locus around THAT precise location (±500kb window)
- Parse BRAKER3 GFF for flanking genes within that window
- **Advantage:** Only creates loci where target actually exists on genome

---

## Implementation Status

### ✅ Phase 1 Script Working

**File:** `scripts/01_discover_paralogs_gene_level.py` (currently restored to OLD version)

**Test Results (manual run with LOC117167432):**

#### BK/LB Results (DIAMOND approach - unchanged):
```
BK: 1 locus with 14 flanking proteins ✓
LB: 1 locus with 79 flanking proteins ✓
```

#### TR Result (tblastn approach):
```
TR_CM036341.1_a at 9,533,067 with 141 flanking genes ✓
- Fixed: Scaffold name prefix matching (CM036341.1_Telenomus_remus... in GFF)
- Result: Successfully extracted 141 BRAKER3 flanking genes
```

#### DR Results (tblastn approach):
```
DR_CM021344.1_a at 34,254,843 with 79 flanking genes ✓
DR_CM021340.1_b at 45,041,643 with 50 flanking genes ✓ ← KEY FIX
```

**Critical Success:** DR_CM021340.1_b found at **45,041,643** (45M), not the incorrect ~44M location from GFF coordinates. This is the 600kb correction!

---

## Current Issue: Pipeline Integration

### Problem:
Pipeline jobs submitted but not completing Phase 1:

```bash
sbatch run_synteny_pipeline_master.slurm ferritin_MC102 ferritin_old_version_test
# Job 99321 - hung at Phase 1 input reading
```

### Evidence:
```
[2025-11-05 09:12:28] Input: inputs/ferritin_MC102*_test.tsv (auto-detected)
[2025-11-05 09:12:28] Output: outputs/ferritin_old_version_test/phase12_landmark/
# No further output
# No output files created in phase12_landmark/
```

### Expected vs Actual:
- **Manual script execution:** Works perfectly, creates all loci ✓
- **Pipeline execution:** Hangs at Phase 1, no output created ✗

### Why This Matters:
The Phase 1 changes shouldn't break downstream phases because:
1. Output format is identical (`locus_definitions.tsv` has same columns)
2. Flanking protein files have same U/D labeling format
3. Only difference is HOW targets are found (tblastn vs DIAMOND)

**Therefore:** Need to debug why pipeline integration is failing, not the Phase 1 logic itself.

---

## Validation Plan

Once pipeline runs complete, validate DR_CM021340.1_b through full pipeline:

### Phase 1: ✅ DONE
- tblastn finds target at 45,041,643
- Extracts 50 flanking genes (26U + 24D)
- Creates `DR_CM021340.1_b_flanking.faa`

### Phase 2: TODO
- Should detect synteny block around 45M (not 44M)
- Verify synteny coordinates match Phase 1 target location

### Phase 3: TODO
- Genome-wide tblastn should find same location (45M)

### Phase 4: TODO - CRITICAL TEST
- BLAST for target genes in synteny blocks
- **Key test:** Should find target at DR_CM021340.1_b (45M)
- **Old problem:** Found nothing at DR_CM021340.1_a (44M) because coordinates were wrong

### Phase 5: TODO
- Classify target as syntenic (should pass, target is AT the locus center)

### Phase 6: TODO
- Extract full gene structure with Exonerate
- Should succeed because target actually exists at these coordinates

### Phase 8: TODO
- DR should show target at DR_CM021340.1_b, not `[empty]`

---

## Next Steps

### Immediate (Debug Pipeline):
1. ✅ Confirmed Phase 1 script works in isolation
2. ❌ Pipeline integration failing - need to debug
3. Check:
   - Is input file being read correctly?
   - Are LOC IDs being parsed from TSV?
   - Is script being called with correct arguments?
   - Any errors in pipeline logs we're missing?

### Once Pipeline Runs:
1. Run full pipeline on ferritin with tblastn-first Phase 1
2. Validate DR_CM021340.1_b through all phases
3. Compare with old DR_CM021340.1_a results (should show target now, not `[empty]`)
4. If successful, document as fix for BRAKER3 coordinate mismatch issue

---

## Files Modified

### Active:
- `scripts/01_discover_paralogs_gene_level.py` - Currently OLD version (working)
- `scripts/01_discover_paralogs_gene_level_tblastn_first.py` - New tblastn version (isolated, not integrated)

### Backups:
- `scripts/01_discover_paralogs_gene_level.py.backup_20250104` - Original working version

### Test Outputs:
- `outputs/ferritin_tblastn_manual_test/` - Manual test with tblastn approach ✓
- `outputs/ferritin_test_manual_old/` - Manual test with old DIAMOND approach ✓

---

## Key Technical Details

### tblastn Parameters (for TR/DR):
```python
EVALUE_THRESHOLD = 1e-5      # Relaxed from 1e-10 (tblastn is more conservative)
MIN_COVERAGE = 0.25          # Relaxed from 0.5 (handles exon fragmentation)
MIN_IDENTITY = 35.0          # Relaxed from 40 (handles translation ambiguity)
WINDOW_SIZE_KB = 500         # ±500kb flanking window around tblastn hit
```

### GFF Scaffold Name Matching Fix:
```python
# Fixed for TR: GFF has long names like "CM036341.1_Telenomus_remus_isolate..."
# But tblastn finds hits on "CM036341.1"
# Solution: Prefix matching
if scaffold != window['scaffold'] and not scaffold.startswith(window['scaffold'] + '_'):
    continue
```

---

## Expected Outcome

If tblastn-first approach works correctly:

1. **DR_CM021340.1_b should have target through entire pipeline**
   - Phase 1: tblastn finds at 45M ✓
   - Phase 4: BLAST finds target at 45M (not empty)
   - Phase 6: Exonerate extracts successfully
   - Phase 8: Shows `[178I]` not `[empty]`

2. **No spurious loci created**
   - Only loci where tblastn actually found hits
   - No more "synteny exists but no target" situations

3. **BRAKER3 genomes usable for comparative genomics**
   - TR and DR loci are valid, not artifacts of bad coordinates
   - Can compare to BK/LB ground truth confidently

---

## Questions to Answer

1. **Why is pipeline hanging at Phase 1?**
   - Script works manually but not via SLURM
   - Need to check job logs, input file parsing, environment

2. **Does DR_CM021340.1_b pass validation through Phases 2-6?**
   - Will Phase 4 find target at 45M?
   - Will Phase 6 extraction succeed?

3. **Is tblastn creating the RIGHT loci and not spurious ones?**
   - Compare full pipeline results with old DIAMOND approach
   - Should see fewer empty loci, more successful extractions

---

## Status Summary

| Component | Status | Notes |
|-----------|--------|-------|
| Phase 1 script (manual) | ✅ Working | tblastn finds correct 45M location |
| TR flanking extraction | ✅ Fixed | Prefix matching for long scaffold names |
| DR flanking extraction | ✅ Working | 79 and 50 flanking genes extracted |
| Pipeline integration | ❌ Failing | Hangs at Phase 1, needs debugging |
| Full pipeline validation | ⏸️ Pending | Waiting on pipeline fix |

**Bottom line:** Phase 1 tblastn approach works correctly and finds the right loci. Need to fix pipeline integration issue so we can validate through Phases 2-6.
