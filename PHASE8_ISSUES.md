# Phase 8a Critical Issues - 2025-11-04

## Issue #1: SwissProt Fallback Not Working for TR/DR Loci

**Severity:** Critical
**Component:** `scripts/08a_generate_locus_matrices.py` lines 458-476

### Problem
Protein column headers for TR/DR loci still show raw protein IDs (`U13_protein LOC117183155`) instead of functional names from SwissProt annotations.

### Expected Behavior
Column headers should display functional names like `U9_Cilium assembly protein DZIP1L`

### Actual Behavior
```
U13_protein LOC117183155
U12_DNA replication complex GIN...
```

### Evidence
```bash
head -2 outputs/ferritin_phase3_real/08_matrices/BK_chr2_a_genome_swissprot_matrix.tsv | cut -f1,5,10-15
```

### Root Cause
SwissProt fallback logic added in commit 1cdb570 may not be matching `bk_protein_id` correctly for TR/DR locus proteins. The fallback only derives names from SwissProt annotations that exist, but TR/DR proteins may not have corresponding SwissProt matches.

### Files Affected
- `scripts/08a_generate_locus_matrices.py`
- `outputs/ferritin_phase3_real/07_swissprot_annotations/genome_specific_swissprot_annotations.tsv`

---

## Issue #2: All Flanking Proteins Labeled 'U' (Upstream) - Missing 'D' (Downstream)

**Severity:** Critical
**Component:** Phase 1/2 output format or Phase 8a parsing

### Problem
All flanking protein columns labeled with 'U' prefix (upstream). No 'D' prefix (downstream) columns appear.

### Expected Behavior
Phase 1 saves flanking proteins in genomic order with labels:
- `U1, U2, ..., Un` for upstream genes
- `D1, D2, ..., Dn` for downstream genes

### Actual Behavior
All columns show `U1, U2, U3, ...` with no downstream genes

### Root Cause Hypothesis
1. **Most Likely:** Ferritin was run with OLD pipeline (pre-Phase 1/2 redesign) which didn't save U/D labels
2. **Alternative:** Phase 1 flanking extraction bug not saving downstream genes
3. **Alternative:** Phase 8a parsing not recognizing D labels

### Evidence Needed
Check if ferritin has Phase 1/2 outputs:
```bash
ls outputs/ferritin_phase3_real/*/flanking*.faa
head outputs/ferritin_phase3_real/locus_definitions.tsv
```

Compare to MC6 v11:
```bash
head outputs/mc6_test_20251104_v11/phase12_landmark/locus_definitions.tsv
```

### Files Affected
- Phase 1: `scripts/01_discover_paralogs_gene_level.py`
- Phase 8a: `scripts/08a_generate_locus_matrices.py`

---

## Issue #3: Matrix Cells Completely Empty

**Severity:** Critical
**Component:** Phase 8a SwissProt annotation lookup

### Problem
Genome rows in the matrix have empty cells - no SwissProt annotations appearing.

### Expected Behavior
Each cell should contain SwissProt functional annotation for that protein in that genome (if synteny block exists).

### Actual Behavior
Matrix cells are blank/empty for most genomes.

### Root Cause Hypothesis
1. **Lookup mismatch:** `bk_protein_id` in SwissProt file doesn't match the protein IDs used for column headers
2. **Genome ID mismatch:** Genome identifiers in SwissProt annotations don't match matrix row genome_ids
3. **Missing data:** SwissProt annotations not loaded/indexed correctly

### Evidence Needed
```bash
# Check SwissProt file structure
head -5 outputs/ferritin_phase3_real/07_swissprot_annotations/genome_specific_swissprot_annotations.tsv

# Check what protein IDs are in SwissProt
cut -f3 outputs/ferritin_phase3_real/07_swissprot_annotations/genome_specific_swissprot_annotations.tsv | sort -u | head

# Check matrix for non-empty cells
grep -v "^genome_id" outputs/ferritin_phase3_real/08_matrices/BK_chr2_a_genome_swissprot_matrix.tsv | head -5
```

### Files Affected
- `scripts/08a_generate_locus_matrices.py` (SwissProt indexing/lookup logic)
- `outputs/ferritin_phase3_real/07_swissprot_annotations/genome_specific_swissprot_annotations.tsv`

---

## Issue #4: Target Locus Column Shows '0P' or Empty

**Severity:** High
**Component:** Phase 8a target gene parsing

### Problem
Target locus column (showing syntenic target gene status) displays `0P` or is empty instead of proper gene status/length.

### Expected Behavior
Should show target gene classification and protein count, e.g.:
- `ferritin_MC102 [3P]` - 3 proteins found
- `ferritin_MC102 [syntenic]` - synteny present

### Actual Behavior
```
ferritin_MC102 [0P]
ferritin_MC102 [not found]
```

Even when synteny blocks exist.

### Root Cause Hypothesis
1. **Targets file mismatch:** Target classification logic not matching genome/locus correctly
2. **Parsing bug:** Status/length extraction from targets file broken
3. **Missing data:** `all_targets_classified.tsv` doesn't contain expected data

### Evidence Needed
```bash
# Check targets file structure
head outputs/ferritin_phase3_real/05_classified/all_targets_classified.tsv

# Check for specific genome entries
grep "GCA_020615435.1" outputs/ferritin_phase3_real/05_classified/all_targets_classified.tsv | head
```

### Files Affected
- `scripts/08a_generate_locus_matrices.py` (target lookup logic)
- `outputs/ferritin_phase3_real/05_classified/all_targets_classified.tsv`

---

## Issue #5: Ferritin Run with OLD Pipeline

**Severity:** High
**Component:** Data compatibility

### Problem
Ferritin outputs were generated with the OLD pipeline (before Phase 1/2 redesign), missing critical metadata.

### Evidence
```bash
# Old format (ferritin)
head -3 outputs/ferritin_phase3_real/locus_definitions.tsv
# Shows: locus_id, gene_family (2 columns)

# New format (MC6)
head -3 outputs/mc6_test_20251104_v11/phase12_landmark/locus_definitions.tsv
# Shows: locus_id, gene_family, genome, chromosome, target_gene, upstream_count,
#        downstream_count, flanking_count, tandem_count, is_tandem, flanking_file (11 columns)
```

### Required Action
**Rerun ferritin through full master pipeline coordinator:**
```bash
sbatch run_synteny_pipeline_master.slurm ferritin ferritin_full_20251104
```

This will regenerate:
- Phase 1/2: Locus definitions with flanking protein files (U/D labeled)
- Phase 3a/3b: Synteny blocks with new data structure
- All downstream phases will work correctly

### Files Affected
- All ferritin outputs in `outputs/ferritin_phase3_real/`

---

## Recommended Fix Order

1. **FIRST:** Rerun ferritin with master coordinator (Issue #5)
   - This will provide properly formatted Phase 1/2 outputs
   - Will likely fix Issue #2 (U/D labeling)

2. **SECOND:** Debug SwissProt fallback (Issue #1)
   - Check protein ID matching logic
   - Verify TR/DR proteins exist in SwissProt annotations

3. **THIRD:** Debug matrix cell population (Issue #3)
   - Check genome ID matching
   - Verify SwissProt lookup logic

4. **FOURTH:** Fix target locus parsing (Issue #4)
   - Check targets file structure compatibility
   - Debug status/length extraction

---

## Test Data

**MC6 v11** (job 98501): Currently running with NEW pipeline - use this to verify fixes
**Ferritin**: Needs complete rerun with master coordinator

---

Generated: 2025-11-04 22:15 MST
