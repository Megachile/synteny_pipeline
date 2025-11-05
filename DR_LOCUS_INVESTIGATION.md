# Investigation: Why DR lacks targets at its own loci

**Date:** 2025-01-04
**Investigator:** Claude Code
**Issue:** DR_CM021340.1_a shows `[empty]` despite being DR's reference locus

---

## Summary

**Root Cause:** Locus was defined based on synteny block (flanking protein conservation), but the actual ferritin gene is **600kb away** from that synteny block. Phase 5 classification correctly rejected it as "not syntenic" because distance exceeds threshold.

**Key Finding:** This reveals a fundamental workflow design flaw - loci are defined by synteny WITHOUT validating that the target gene exists within the synteny block.

---

## Investigation Steps

### 1. Check Phase 4 BLAST Results

DR had 3 BLAST hits in Phase 4:

```bash
$ grep "GCA_030998225.1" outputs/ferritin_phase3_real/04_target_genes/all_target_loci.tsv

XP_033211301_targets_CM021340.1_RagTag_001  CM021340.1_RagTag  +  3  3  45041637  45045491  3.854  2  ferritin_MC102  3.50701e-38  GCA_030998225.1
XP_033208250_targets_CM021344.1_RagTag_001  CM021344.1_RagTag  +  -2  -2  34237318  34238103  0.785  2  ferritin_MC102  1.50254e-25  GCA_030998225.1
XP_033208250_targets_CM021344.1_RagTag_002  CM021344.1_RagTag  +  0  -3,-2  34254843  34255170  0.327  2  ferritin_MC102  5.05529e-33  GCA_030998225.1
```

**Key observation:** DR has ferritin hits, so Phase 4 target detection worked!

### 2. Check Phase 5 Classification

```bash
$ grep "GCA_030998225.1" outputs/ferritin_phase3_real/05_classified/all_targets_classified.tsv | grep "DR_"

# DR_CM021344.1_a: BOTH targets classified ✓
XP_033208250_targets_CM021344.1_RagTag_001  ...  GCA_030998225.1  ...  synteny  DR_CM021344.1_a
XP_033208250_targets_CM021344.1_RagTag_002  ...  GCA_030998225.1  ...  synteny  DR_CM021344.1_a

# DR_CM021340.1_a: NO targets for DR! ✗
# (search returns nothing)
```

**Finding:** DR_CM021344.1_a has targets, but DR_CM021340.1_a doesn't!

### 3. Compare Coordinates: DR_CM021340.1_a

**DR's synteny block:**
```
Scaffold: CM021340.1_RagTag
Start: 43,822,130
End: 44,409,502
Span: 587.4 kb
```

**DR's BLAST hit:**
```
Scaffold: CM021340.1_RagTag (SAME scaffold!)
Start: 45,041,637
End: 45,045,491
Span: 3.8 kb
```

**Distance between synteny block and target hit:** ~600kb

**Conclusion:** The target gene is NOT within the synteny block. Phase 5 classification correctly rejected it as "too far from synteny block."

### 4. Check Phase 8a Output

```bash
$ grep "GCA_030998225.1" DR_CM021340.1_a_genome_swissprot_matrix.tsv

synteny_pct: 100.0        # Perfect flanking conservation
num_proteins_found: 15    # All flanking proteins present
TARGET: ferritin_MC102 [empty]  # But no target gene!
```

### 5. Compare with DR_CM021344.1_a (working locus)

```bash
$ grep "GCA_030998225.1" DR_CM021344.1_a_genome_swissprot_matrix.tsv

synteny_pct: 100.0
num_proteins_found: 15
TARGET: ferritin_MC102 [178I]  # Has target! ✓
```

**DR_CM021344.1_a works correctly** - target is within synteny block.

---

## Root Cause Analysis

### How Loci Are Currently Defined (Phase 2)

1. User specifies reference genome coordinates (or Phase 2 auto-detects)
2. Extract flanking proteins around reference coordinates
3. BLAST flanking proteins against other genomes
4. Define synteny blocks based on flanking protein conservation
5. **Create locus WITHOUT validating target gene exists in synteny block**

### Why DR_CM021340.1_a Failed

```
Phase 2: Created locus based on synteny block (43.8-44.4M)
         Assumed ferritin exists within this block
         Never validated with exonerate extraction

Phase 4: Found ferritin 600kb away (45.0-45.0M)
         BLAST hit valid but distant from synteny

Phase 5: Classified hit as "too far from synteny block"
         Correctly rejected based on distance threshold
         DR ends up with no targets at this locus

Phase 8: Shows [empty] - synteny exists but no target
```

### The Core Problem

**Locus definition is based on synteny (flanking proteins), not target presence.**

This creates "phantom loci" where:
- ✓ Flanking proteins are conserved (synteny block exists)
- ✗ Target gene is missing or too distant
- Result: Reference genome has no target at "its own" locus

---

## Evidence This Is a Design Issue, Not a Bug

1. **Phase 5 classification is working correctly**
   - It properly filters targets >200kb from synteny blocks
   - The 600kb distance exceeds any reasonable threshold

2. **The locus definition process doesn't validate target presence**
   - No exonerate extraction in Phase 2
   - No check that reference genome has extractable gene
   - Synteny alone is insufficient to define a valid locus

3. **Other genomes successfully find targets at this locus**
   - Many genomes show targets at DR_CM021340.1_a
   - But DR itself (the reference!) doesn't

---

## Proposed Solution: Validate Reference Target in Phase 2

### Current Phase 2 Workflow:
```
1. Define coordinates (manual or automatic)
2. Extract flanking proteins
3. BLAST flanking proteins → synteny blocks
4. Create locus
```

### Proposed Phase 2 Workflow:
```
1. Define coordinates (manual or automatic)
2. Extract flanking proteins
3. BLAST flanking proteins → synteny blocks
4. **Run exonerate extraction on reference genome**
5. **ONLY create locus if:**
   - Synteny block exists AND
   - Target gene successfully extracted AND
   - Target is within acceptable distance of synteny block
6. Save extracted reference sequence for downstream
```

### Benefits:

1. **Fail fast** - Catch invalid loci in Phase 2, not Phase 8
2. **Quality control** - Only valid, extractable loci enter pipeline
3. **Honest reporting** - Every reference locus guaranteed to have reference sequence
4. **Earlier feedback** - User knows immediately if locus definition failed
5. **Reduced wasted computation** - Don't process invalid loci through entire pipeline

### Implementation Notes:

- Add exonerate extraction step after synteny block definition
- Set distance threshold (e.g., target must be within 200kb of synteny block)
- Log extraction failures with diagnostic information
- Create `locus_validation_report.tsv` showing which loci passed/failed
- Update locus_definitions.tsv to only include validated loci

---

## Additional Findings

### Protein Source: DR uses Braker3 predictions

```bash
$ head -5 DR_CM021340.1_a_flanking_filtered.faa

>g15427.t1|g15427 g15427.t1
MRVRLLWIEEQLAGIHTNFISSVSASLAIVHTTLAQEGYDYPRPSDPFISGIPSGTTSKP...
```

**Protein ID format:** `g15427.t1` (Braker3 gene models)

**Contrast with BK:**
```bash
$ head -5 BK_chr2_a_flanking_filtered.faa

>XP_033209112.1|LOC117167954 XP_033209112.1 probable elongator complex protein 2
```

**Protein ID format:** `XP_033209112.1` (NCBI RefSeq)

**Implication:** This explains why DR/TR proteins lack SwissProt annotations - they're Braker3 predictions, not curated NCBI proteins. Likely less complete/accurate annotations.

---

## Recommendations

### Immediate Actions:

1. **Document this as expected behavior**
   - DR_CM021340.1_a is not a valid locus
   - Synteny exists but reference lacks target gene
   - This is a locus definition issue, not an extraction failure

2. **Identify other phantom loci**
   - Check all reference loci for same pattern
   - Report which loci lack reference targets
   - May need to redefine or remove invalid loci

### Long-term Actions:

1. **Implement Phase 2 validation**
   - Add exonerate extraction before locus creation
   - Only create loci that pass validation
   - Save validated reference sequences

2. **Re-run Phase 2 with validation**
   - May reduce number of loci (some will fail validation)
   - But ensures all remaining loci are valid
   - Cleaner, more honest analysis

3. **Consider protein source implications**
   - Document which genomes use Braker3 vs NCBI proteins
   - May affect downstream annotation quality
   - Could seek NCBI RefSeq for DR/TR if available

---

## Conclusion

**DR_CM021340.1_a is not broken - it's an invalid locus that should never have been created.**

The synteny block exists (flanking proteins conserved) but the ferritin gene is 600kb away, far outside any reasonable synteny boundary. Phase 5 correctly rejected this as non-syntenic.

**The solution is workflow redesign, not bug fixing.** Phase 2 must validate that reference genomes have extractable target genes before creating loci.
