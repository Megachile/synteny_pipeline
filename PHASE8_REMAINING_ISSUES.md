# Phase 8 Remaining Issues

**Status:** 5 issues identified after initial fixes (2025-01-04)

---

## Issue #0: Empty brackets `[]` when target exists but extraction metadata missing

**Status:** Not started
**Priority:** CRITICAL
**Affects:** Phase 8a locus-specific matrices

**Problem:**
- Some genomes show `ferritin_MC102 []` instead of proper target or `[empty]`
- Empty brackets mean: targets exist in classification BUT extraction metadata lookup failed
- Example: Telenomus at TR_CM036338.1_a shows `[]` but has 100% synteny

**Evidence:**
```
# Telenomus at TR locus:
ferritin_MC102 []  # Empty brackets - bad!

# Targets file shows:
XP_033211301_targets_CM036338_001 ... GCA_020615435.1 ... synteny TR_CM036338.1_a

# But extraction directory is:
outputs/.../extracted_sequences/GCA_020615435.1/XP_033211301_targets_CM036338_002/
                                                                              ^^^
# Notice _002 vs _001 mismatch!
```

**Root Cause:**
- Phase 5 classification creates locus_name: `XP_033211301_targets_CM036338_001`
- Phase 6 extraction creates directory: `XP_033211301_targets_CM036338_002`
- Phase 8a tries to look up `_001` in metadata dict but only `_002` exists
- Lookup fails, target gets skipped, `parts` list is empty → `[]`

**Why this happens:**
- Phase 6 likely deduplicated tandems and renumbered remaining targets
- OR there's a bug in how locus names are generated in Phase 6

**Fix needed:**
1. Investigate Phase 6 extraction naming logic
2. Either: Fix Phase 6 to use consistent names with Phase 5
3. Or: Update Phase 8a to handle name mismatches (map by coordinates?)
4. Or: Re-run Phase 6 with fixed naming

**Impact:**
- Makes matrices misleading - looks like no target when one exists
- Different from `[empty]` which correctly means "synteny but no target found"

---

## Issue #0B: Reference loci exist but reference species has no target - workflow validation failure

**Status:** Not started
**Priority:** CRITICAL - WORKFLOW DESIGN ISSUE
**Affects:** Phase 2 locus definition, entire pipeline validity

**Problem:**
- DR loci exist (DR_CM021340.1_a, DR_CM021344.1_b) but DR has no targets at them
- Shows `[empty]` at its own reference locus
- **This shouldn't be possible** - locus definition implies target exists in reference genome
- How did these loci make it through Phase 2 if DR has no ferritin there?

**Evidence:**
```
# DR at DR_CM021340.1_a:
synteny_pct: 100.0     # Perfect synteny conservation
num_proteins_found: 15  # All flanking proteins found
TARGET: ferritin_MC102 [empty]  # But no target!!!

# DR at DR_CM021344.1_b:
synteny_pct: 100.0
TARGET: ferritin_MC102 [empty]

# No targets in classification:
grep GCA_030998225.1 all_targets_classified.tsv | grep "DR_CM021340.1_a\|DR_CM021344.1_b"
# Returns nothing!
```

**Why these loci shouldn't exist:**
1. Phase 4 target detection couldn't have been filtered:
   - Not tandem duplicates
   - No strength filter applied to syntenic hits
2. If DR has no ferritin at these locations, how did Phase 2 define them as "DR loci"?
3. Phase 2 should validate that reference genome actually HAS the gene before creating locus

**Root Cause Hypotheses:**

### Hypothesis A: Phase 6 extraction failed (user suggestion)
- Phase 2: synteny block defined (flanking proteins conserved)
- Phase 4: BLAST finds candidate target in DR genome
- Phase 5: classifies as syntenic
- **Phase 6: exonerate extraction FAILS** (couldn't extract actual sequence)
- Phase 8: shows [empty] because no extracted sequence exists

**Supporting evidence needed:**
- Check Phase 4 BLAST results - did DR have target hits?
- Check Phase 6 extraction logs - did exonerate fail for DR?
- Check if extracted_sequences/ has failed/empty DR entries

### Hypothesis B: Locus defined without validation
- Phase 2 creates locus based on:
  - User-specified reference coordinates?
  - Synteny block presence alone?
  - Assumed target exists without checking?
- Never validated that reference genome actually has extractable gene

**Proposed Solution: Immediate exonerate validation in Phase 2**

User suggestion: "maybe it would be worth trying that exonerate extraction on the candidate loci immediately in phase 2"

**Benefits:**
1. **Validate before defining** - don't create locus unless gene can be extracted
2. **Fail fast** - catch extraction issues early, not in Phase 8
3. **Quality control** - only valid, extractable loci enter the pipeline
4. **Honest loci** - every reference locus guaranteed to have reference sequence

**Implementation:**
```
Phase 2 (modified):
1. Define synteny block (flanking proteins)
2. Identify candidate target region (BLAST or coordinates)
3. **Run exonerate extraction on reference genome**
4. **ONLY create locus if extraction succeeds**
5. Save extracted reference sequence for downstream use
```

**Investigation needed:**
1. How were DR loci originally defined? (Manual coordinates? Automatic?)
2. Check Phase 4 results: did DR have BLAST hits?
3. Check Phase 6 logs: did exonerate fail for DR?
4. Check if this affects other loci (TR, LB, etc.)

**Priority justification:**
- If reference loci don't have reference sequences, the entire analysis is invalid
- Must understand and fix before trusting any locus definitions
- May require re-running Phase 2 with extraction validation

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

## Issue #2: TR/DR loci have no SwissProt annotations - Need to document root cause

**Status:** Investigation needed
**Priority:** MEDIUM
**Affects:** Phase 7 SwissProt annotation, Phase 8a matrices for TR/DR loci

**Problem:**
- TR and DR locus flanking proteins show "no match" for ALL proteins
- BK and LB proteins have excellent SwissProt coverage (55-90% identity)
- Need to understand WHY this happened

**Evidence:**
```
# User observation: "not a single one of the TR or DR flanking genes has a swissprot match"
# BK proteins show matches like: "Ferritin light chain (55%)", "Aurora kinase (79%)"
# TR/DR proteins: ALL show "no match"
```

**Hypothesis - Protein Source Difference:**
- BK/LB proteins likely from NCBI RefSeq (XP_033209222.1 format)
- TR/DR proteins possibly from Braker3 gene predictions (g12345.t1 format)
- Braker3 predictions may not match SwissProt as well as curated NCBI proteins

**Investigation needed:**
1. **Check protein IDs in flanking files:**
   - Look at BK_chr2_a_flanking_filtered.faa - what ID format?
   - Look at TR_CM036338.1_a_flanking_filtered.faa - what ID format?
   - Look at DR_CM021340.1_a_flanking_filtered.faa - what ID format?

2. **Check protein source:**
   - Where did TR/DR proteins come from? (Braker3? NCBI? Custom annotations?)
   - Are these genome assemblies less complete/annotated than BK?
   - Do TR/DR genomes have NCBI RefSeq annotations available?

3. **Check SwissProt BLAST results:**
   - Did Phase 7 actually BLAST TR/DR proteins?
   - Did they produce hits but with poor e-values/scores?
   - Or did BLAST fail entirely for these proteins?

**Documentation goal:**
- Understand and document the biological/technical reason for annotation differences
- NOT to copy BK annotations to TR/DR (would be dishonest)
- May reveal broader data quality issues or assembly completeness differences

**Potential outcomes:**
1. TR/DR proteins are real but divergent → document as limitation
2. Braker3 predictions are low quality → consider alternative protein sets
3. SwissProt BLAST parameters need adjustment → re-run Phase 7
4. TR/DR genomes lack quality annotations → document and move forward

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
