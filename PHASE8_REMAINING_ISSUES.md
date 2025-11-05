# Phase 8 Remaining Issues

**Status:** 3 issues identified after initial fixes (2025-01-04)

---

## Issue #1: Downstream proteins missing from Phase 8a matrices

**Status:** Not started
**Priority:** HIGH
**Affects:** Phase 8a locus-specific matrices

**Problem:**
- Matrix only shows upstream columns (U14, U13, ..., U1)
- No downstream columns (D1, D2, D3, etc.) appear in output
- The script parses downstream proteins from the flanking .faa file but never adds them to the matrix

**Evidence:**
```bash
head -1 outputs/ferritin_phase3_real/08_matrices/BK_chr2_a_genome_swissprot_matrix.tsv
# Shows: ...U4, U3, U2, U1, TARGET
# Missing: D1, D2, D3, ...
```

**Diagnosis:**
- Script successfully parses downstream proteins (lines 181-186 in 08a_generate_locus_matrices.py)
- Script successfully identifies which downstream proteins are found in BLAST hits
- But downstream column creation code (lines 405-420) is never executed or has a bug

**Fix needed:**
- Verify downstream column creation loop actually runs
- Check if `downstream_proteins` list is populated correctly
- Debug column creation logic

---

## Issue #2: TR/DR loci have no SwissProt annotations

**Status:** Not started
**Priority:** MEDIUM
**Affects:** Phase 7 SwissProt annotation, Phase 8a matrices for TR/DR loci

**Problem:**
- TR (Torymus) and DR (Druon) locus flanking proteins show "no match" for ALL proteins
- BK and LB (Leptopilina) proteins have excellent SwissProt coverage
- This is because TR/DR proteins came from Braker3 gene predictions, not NCBI

**Evidence:**
```
# User observation: "not a single one of the TR or DR flanking genes has a swissprot match"
# BK proteins: XP_033209222.1 (NCBI RefSeq IDs)
# TR/DR proteins: Likely g12345.t1 format (Braker3 gene model IDs)
```

**Root Cause:**
1. Phase 1 extracted BK proteins from NCBI RefSeq (XP_ format IDs)
2. Phase 1 extracted TR/DR proteins from Braker3 predictions (different ID format)
3. Phase 7 SwissProt annotation blasts proteins against SwissProt database
4. SwissProt matching may only work for NCBI IDs or the Braker3 protein sequences don't match well

**Possible Solutions:**

### Option A: Re-run SwissProt annotation with TR/DR query proteins
- Modify Phase 7 to use TR/DR flanking proteins as queries
- Create TR/DR-specific SwissProt annotations
- Update Phase 8a to use TR/DR annotations as fallback when BK annotations missing

### Option B: Implement BK-based fallback (Issue #1 from original PHASE8_ISSUES.md)
- When target genome protein has no SwissProt match
- Check if the protein matched a BK query protein in Phase 2 BLAST
- Use the BK query protein's SwissProt annotation as fallback
- This would work for any protein that shares homology with BK proteins

### Option C: Extract TR/DR proteins from NCBI (if available)
- Check if TR/DR have NCBI RefSeq annotations available
- Re-extract proteins from NCBI annotations instead of Braker3
- Re-run Phase 2 synteny detection with NCBI-based proteins
- This would give XP_-format IDs that work with existing SwissProt pipeline

**Recommendation:** Option B (BK-based fallback) is most robust and was the original plan in PHASE8_ISSUES.md Issue #1. It will work for all loci regardless of protein ID format.

---

## Issue #3: Verify Phase 8b also needs downstream fixes

**Status:** To investigate
**Priority:** LOW
**Affects:** Phase 8b summary matrices (potentially)

**Question:**
- Does Phase 8b have the same downstream protein issue?
- Phase 8b creates summary matrices aggregating across loci
- Need to check if it has flanking protein columns at all

**Action:**
- Inspect Phase 8b output to see if it has flanking protein context
- If not, this may be by design (summary level doesn't need flanking detail)

---

## Completed Issues (from original PHASE8_ISSUES.md)

- ✅ **Issue #2:** Phase 8b showing `[0P]` entries - FIXED (filter length=0 targets)
- ✅ **Issue #3:** Phase 8a matrix cells empty - FIXED (protein ID matching)
- ✅ **Issue #4:** Phase 8a TARGET column showing `[0P]` - FIXED (extraction metadata + filter)
- ✅ **Issue #5:** Synteny block presence/absence ambiguous - FIXED (`[empty]` vs "No block")
