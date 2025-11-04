# Phase 8 Matrix Output Issues - Nov 4, 2025

## Summary

Phase 8a/8b successfully generate matrices, but have several output quality issues that need fixing.

## Issues Found

### 1. SwissProt Annotations Not Showing in Locus Matrices

**Problem**: Locus matrices show empty flanking protein columns despite 293,478 SwissProt annotations existing.

**Root Cause**: Phase 7 SwissProt annotation file has `bk_protein` column = "unknown" instead of actual BK protein IDs.

**Evidence**:
```bash
# SwissProt file has "unknown" in bk_protein column:
$ head -2 outputs/ferritin_phase3_real/07_swissprot/genome_specific_swissprot_annotations.tsv
locus   genome  bk_protein  query_id  ...
BK_chr2_a  GCA_037103525.1  unknown  GCA_037103525.1_XP_033209112.1|LOC117167954_...
```

**Expected**: `bk_protein` should be "XP_033209112.1" (extracted from query_id)

**Impact**: Phase 8a can't match SwissProt annotations to flanking proteins, leaving all U1-U14 columns empty.

**Fix Needed**: Phase 7 (scripts/07_swissprot_annotation.py) must extract BK protein ID from query_id and populate bk_protein column.

---

### 2. Redundant "All Genes Combined Summary" File

**Problem**: Phase 8b creates `all_genes_combined_summary.tsv` which is identical to single gene family summaries.

**File**: `outputs/ferritin_phase3_real/07_matrices/gene_type_summaries/all_genes_combined_summary.tsv`

**Impact**: Confusing for users, redundant for single gene families.

**Fix Needed**: Remove generation of "all genes combined" file from scripts/08b_generate_summary_matrices.py.

---

### 3. Unclear Flanking Protein Counts

**Question**: Are matrices showing ALL flanking proteins or the FILTERED set (15 proteins from MAX_FLANKING_FOR_SYNTENY)?

**Context**:
- Phase 1 extracts N flanking proteins per locus (variable, 14-85 in ferritin)
- Phase 3 uses config.MAX_FLANKING_FOR_SYNTENY = 15 for synteny detection
- Phase 8a shows U1-U14 columns in matrices

**Verification Needed**:
1. Check if Phase 8a uses all proteins from flanking.faa or filtered set
2. Check if column counts match MAX_FLANKING_FOR_SYNTENY setting
3. Clarify in WORKFLOW_BLUEPRINT what flanking proteins appear in matrices

---

### 4. Phase 6 Output Location Confusion

**Problem**: Phase 6 outputs ARE present but in unexpected location.

**Location**: `outputs/ferritin_phase3_real/04_target_genes/extracted_sequences/`
- 979 protein sequences successfully extracted
- Organized by genome/block_id/gene files

**Expected**: Based on config.py, should be in `STEP06_EXTRACTED` but actually using STEP04_TARGETS subdirectory.

**Impact**: Minor - outputs exist and work, but documentation/config doesn't match reality.

**Fix Needed**: Update WORKFLOW_BLUEPRINT.md and config.py to document actual output structure.

---

## Files Affected

**Scripts needing fixes**:
- `scripts/07_swissprot_annotation.py` - Fix bk_protein column population
- `scripts/08b_generate_summary_matrices.py` - Remove "all genes combined" file

**Documentation needing updates**:
- `WORKFLOW_BLUEPRINT.md` - Clarify flanking protein filtering, Phase 6 output locations
- `scripts/config.py` - Update output directory definitions to match reality

---

## Testing Plan

1. Fix Phase 7 bk_protein column
2. Rerun Phase 7 on ferritin data
3. Rerun Phase 8a to verify SwissProt annotations now appear
4. Fix Phase 8b to remove combined summary
5. Test on MC06 to ensure fixes work end-to-end

---

## Priority

- **P0**: Fix SwissProt annotation mapping (blocks matrix interpretation)
- **P1**: Remove redundant combined summary (user confusion)
- **P2**: Document flanking protein filtering (clarification)
- **P3**: Update Phase 6 output location docs (minor mismatch)
