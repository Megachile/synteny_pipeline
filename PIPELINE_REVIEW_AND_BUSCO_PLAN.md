# Helixer Pipeline Review and BUSCO Integration Plan

**Date:** 2025-12-02
**Status:** Planning Document - No Script Changes

---

## 1. Pipeline Flow Overview

The Helixer pipeline performs synteny-based gene family analysis across cynipid wasp genomes using Helixer-predicted proteomes. The workflow proceeds through these phases:

```
Phase 1 → Phase 1b → Phase 4 → Phase 4b → Phase 2b → Phase 5 → Phase 5b → Phase 6 → Phase 7 → Phase 8a → Phase 8b
   ↓         ↓          ↓         ↓          ↓         ↓         ↓          ↓         ↓         ↓          ↓
 Locus    BK-LB      Target   Flanking   Synteny   Classify  Validate   Extract  SwissProt  Locus    Summary
Discovery Unify    Detection  Extract    Blocks    Targets   Novel     Sequences Annotate  Matrices  Matrix
```

### Phase Descriptions

| Phase | Script | Purpose | Key Inputs | Key Outputs |
|-------|--------|---------|------------|-------------|
| **1** | `01_phase1.py` | Extract target proteins & flanking genes from BK reference | LOC IDs, BK GFF | `locus_definitions.tsv`, `*_flanking.faa` |
| **1b** | `01b_merge_bk_lb_loci.py` | Find LB orthologs, unify BK+LB loci | `locus_definitions.tsv`, LB proteome | Updated `locus_definitions.tsv` with LB columns |
| **4** | `04_detect_targets_helixer.py` | DIAMOND search for targets in all genomes | `target_proteins.faa`, Helixer proteomes | `all_target_loci.tsv`, per-genome hits |
| **4b** | `04b_extract_flanking_helixer.py` | Extract flanking genes around detected targets | `all_target_loci.tsv`, Helixer GFFs | `flanking_genes.tsv`, `all_flanking_proteins.faa` |
| **2b** | `02b_detect_empty_blocks_helixer.py` | Search for synteny blocks via BK flanking genes | BK flanking proteins, Helixer proteomes | `phase2b_synteny_blocks.tsv`, `phase2b_flanking_details.tsv` |
| **5** | `05_classify_targets_helixer.py` | Classify targets as syntenic vs unplaceable | Phase 4/4b outputs, Phase 1 locus defs | `syntenic_targets.tsv`, `unplaceable_targets.tsv` |
| **5b** | `05b_validate_novel_loci.py` | Cluster unplaceables, validate as novel loci | `unplaceable_targets.tsv`, flanking data | `novel_loci_clustered.tsv` |
| **6** | `06_extract_sequences_helixer.py` | Extract protein sequences with length QC | Phase 5 targets, Helixer proteomes | `locus_sequences/`, length flags |
| **7** | `07_annotate_flanking_swissprot.py` | SwissProt annotation of flanking genes | Flanking FASTAs from 4b+2b | `flanking_swissprot.tsv` |
| **8a** | `08a_generate_locus_matrices_helixer.py` | Per-locus matrices with flanking presence | All prior phases | `*_locus_matrix.tsv` files |
| **8b** | `08b_generate_summary_matrices_helixer.py` | Gene family summary matrix | Phase 5/6 outputs | `summary_matrix.tsv` |

---

## 2. Input/Output Interface Issues

### 2.1 Column Name Inconsistencies

**Issue:** Several scripts normalize column names differently, creating fragile interfaces.

| Location | Pattern | Risk |
|----------|---------|------|
| `08a_generate_locus_matrices_helixer.py:197-200` | Renames `classification`→`placement`, `locus_id`→`assigned_to` | Fails if source column names change |
| `08b_generate_summary_matrices_helixer.py:140-143` | Same normalization | Duplicated logic |
| `05_classify_targets_helixer.py` | Outputs `classification` column | Other scripts expect `placement` |

**Recommendation:** Standardize column names across all scripts. Define a shared `column_names.py` constants file.

### 2.2 File Naming Evolution

**Issue:** Some scripts check for both old and new filename conventions.

| Script | Old Name | New Name |
|--------|----------|----------|
| Phase 7 | `empty_synteny_blocks.tsv` | `phase2b_synteny_blocks.tsv` |
| Phase 7 | `empty_block_details.tsv` | `phase2b_flanking_details.tsv` |
| Phase 8a | Same pattern | Same pattern |

**Recommendation:** Remove old filename fallbacks after verifying no outputs use old names. Current dual-check is appropriate during transition but adds maintenance burden.

### 2.3 Target ID Format Mismatch

**Issue:** Helixer proteomes use `.1` suffix on protein IDs, but some scripts don't account for this.

| Location | Behavior |
|----------|----------|
| `02b_detect_empty_blocks_helixer.py:200-205` | Strips `.1` suffix when loading proteomes |
| `04_detect_targets_helixer.py` | Preserves full IDs with suffix |
| `05_classify_targets_helixer.py` | Expects consistent format |

**Recommendation:** Standardize ID handling. Either always strip suffix at input or always preserve it. Current mixed approach risks join failures.

### 2.4 Optional Directory Handling

**Issue:** Some phase directories are optional, leading to silent failures or incomplete outputs.

| Script | Optional Input | Consequence if Missing |
|--------|----------------|------------------------|
| `06_extract_sequences_helixer.py` | `--phase5b-dir` | Novel loci sequences not extracted |
| `08a_generate_locus_matrices_helixer.py` | `--phase5b-dir` | Novel loci not in matrices |
| `08b_generate_summary_matrices_helixer.py` | `--phase5b-dir` | Novel loci not in summary |

**Recommendation:** Add explicit logging when optional inputs are missing. Consider `--require-novel-loci` flag for strict mode.

### 2.5 Metadata Loading Inconsistency

**Issue:** Phase 8a and 8b load target metadata differently.

| Script | Primary Source | Fallback |
|--------|----------------|----------|
| Phase 8a | `syntenic_targets.tsv` | None |
| Phase 8b | `target_proteins.faa` (parses headers) | `syntenic_targets.tsv` |

**Recommendation:** Use consistent metadata loading. The header-parsing approach is fragile.

---

## 3. Cleanup Recommendations

### 3.1 High Priority

1. **Standardize column names** - Create `helixer_pipeline/constants.py`:
   ```python
   # Standard column names
   COL_PLACEMENT = "placement"  # not "classification"
   COL_ASSIGNED_TO = "assigned_to"  # not "locus_id"
   COL_GENOME = "genome_id"
   COL_TARGET = "target_id"
   ```

2. **Remove deprecated filename checks** - After confirming all families use new names:
   - Remove `empty_synteny_blocks.tsv` fallback
   - Remove `empty_block_details.tsv` fallback

3. **Standardize target ID format** - Choose one approach:
   - Option A: Strip `.1` suffix immediately on load (simpler)
   - Option B: Preserve suffix everywhere (more accurate)

### 3.2 Medium Priority

4. **Add interface validation** - Each script should validate input columns exist:
   ```python
   REQUIRED_COLUMNS = ["genome_id", "target_id", "placement"]
   missing = set(REQUIRED_COLUMNS) - set(df.columns)
   if missing:
       raise ValueError(f"Missing columns: {missing}")
   ```

5. **Consolidate flanking file handling** - Phase 7 unions multiple sources; document expected inputs clearly.

6. **Add `--strict` mode** - Fail rather than warn when optional inputs missing.

### 3.3 Low Priority

7. **Refactor duplicate code** - Locus matrix creation logic appears in both 08a and 08b.

8. **Add type hints** - Improve maintainability of complex functions.

---

## 4. BUSCO Integration Plan

### 4.1 Current State

The pipeline uses `genome_quality_assessment.tsv` with a custom `mean_recovery` metric:

```
genome_id           species                 mean_recovery  n_observations  status
GCA_032370625.1     Bassettia_ligni        0.331058       92              best_kept
GCA_021234035.1     Andricus_foecundatrix  0.240000       5               best_kept
```

This metric appears to be gene family recovery percentage (0.24-0.35 range), **not BUSCO scores**.

### 4.2 BUSCO Score Integration

**Step 1: Generate BUSCO Scores**

Run BUSCO on all Helixer proteomes:
```bash
# Using insecta_odb10 database
busco -i proteome.faa -l insecta_odb10 -o busco_out -m proteins
```

**Step 2: Create Unified Quality Table**

New `genome_quality_comprehensive.tsv`:
```
genome_id           species                 busco_complete  busco_single  busco_dup  busco_frag  busco_missing  mean_recovery  status
GCA_032370625.1     Bassettia_ligni        0.876           0.845         0.031      0.067       0.057          0.331          best_kept
GCA_021234035.1     Andricus_foecundatrix  0.723           0.698         0.025      0.142       0.135          0.240          best_kept
```

**Step 3: Add Quality Loading to Scripts**

Create `helixer_pipeline/genome_quality.py`:
```python
def load_genome_quality(quality_file: str) -> Dict[str, GenomeQuality]:
    """Load genome quality metrics including BUSCO scores."""
    df = pd.read_csv(quality_file, sep='\t')
    return {
        row['genome_id']: GenomeQuality(
            busco_complete=row['busco_complete'],
            busco_fragmented=row['busco_frag'],
            busco_missing=row['busco_missing'],
            mean_recovery=row['mean_recovery'],
            confidence_weight=calculate_confidence(row)
        )
        for _, row in df.iterrows()
    }

def calculate_confidence(row) -> float:
    """
    Calculate confidence weight based on BUSCO completeness.
    Returns value between 0.0 (no confidence) and 1.0 (full confidence).
    """
    # Primary: BUSCO complete + single-copy
    busco_quality = row['busco_complete']

    # Penalize high duplication (assembly issues)
    if row['busco_dup'] > 0.1:
        busco_quality *= 0.9

    # Combine with recovery metric if available
    if 'mean_recovery' in row and pd.notna(row['mean_recovery']):
        combined = 0.7 * busco_quality + 0.3 * row['mean_recovery']
    else:
        combined = busco_quality

    return combined
```

### 4.3 BUSCO-Based Genome Exclusion

Before applying quality weighting to kept genomes, we first exclude genomes that are either too poor quality or redundant with better alternatives.

#### Exclusion Parameters

```python
HARD_FLOOR = 40.0         # Exclude genomes below this (always)
REDUNDANCY_GAP = 14.0     # Exclude if this many pts below group median
MIN_GOOD_IN_GROUP = 2     # Need this many good genomes to trigger redundancy exclusion
GOOD_THRESHOLD = 75.0     # What counts as "good" alternative
```

#### Decision Rules

```
For each genome with BUSCO score:

1. IF busco_complete < HARD_FLOOR (40%):
   → EXCLUDE: "below floor"

2. ELSE IF only genome in phylo_order:
   → KEEP: "sole representative" (even if marginal)

3. ELSE:
   - Compute group median for phylo_order
   - Count genomes with busco >= GOOD_THRESHOLD in group
   - IF n_good >= MIN_GOOD_IN_GROUP AND (median - busco) > REDUNDANCY_GAP:
     → EXCLUDE: "redundant"
   - ELSE:
     → KEEP: "passes filters"
```

#### Test Results (21 genomes, 2025-12-02)

| Phylo | Species | BUSCO | Decision | Reason |
|-------|---------|-------|----------|--------|
| 5 | Callaspidia notata | 14.9% | EXCLUDE | below 40% floor |
| 6 | Ganaspis sp | 50.7% | EXCLUDE | redundant (vs 87.8% median, 4 better) |
| 6 | Leptopilina clavipes | 73.0% | EXCLUDE | redundant (vs 87.8% median, 4 better) |
| 7 | Eschatocerus acacia | 40.4% | KEEP | sole representative |
| 16 | Andricus grossulariae | 32.2% | EXCLUDE | below 40% floor |
| 16 | Callirhytis spRG_2019_326 | 46.8% | KEEP | no good alternatives yet |

#### Implementation

**Script:** `helixer_pipeline/collect_busco_and_exclude.py`

**Inputs:**
- BUSCO outputs: `data/ragtag_output/*/busco_v6_hym/short_summary*.txt`
- Species mapping: `data/gca_to_species.tsv`

**Outputs:**
- Updated `data/gca_to_species.tsv` with columns:
  - `busco_complete`, `busco_single`, `busco_frag`, `busco_missing`
  - `busco_status` (keep/exclude/pending)
  - `busco_reason`

**Action:** Excluded genomes moved to `data/ragtag_output_excluded/`

```bash
# When ready (all Helixer + BUSCO complete):
python helixer_pipeline/collect_busco_and_exclude.py --dry-run   # Preview
python helixer_pipeline/collect_busco_and_exclude.py --execute   # Move excluded
```

### 4.4 BUSCO-Gated Novel Locus Proposal

**Problem:** Low-quality genomes have sparse annotations, meaning:
- Unplaceable hits may have few flanking genes
- Few flanking genes = weak anchoring signal
- Matches in other genomes could be coincidence rather than true synteny
- Finding "some" novel loci while unknowingly missing others is misleading

**Solution:** Only genomes above a BUSCO threshold can **propose** novel locus candidates.

```python
MIN_BUSCO_FOR_NOVEL_PROPOSAL = 60.0  # Genomes below this cannot introduce novel loci
```

**Rationale:**
- If a novel locus truly exists, a high-quality genome will propose it
- Low-quality genomes can still **receive** novel loci (be validated as targets)
- This prevents sparse annotations from introducing spurious candidates
- We don't lose real signal - just avoid false positives

**Effect (on current 21 genomes):**

| Can Propose | Cannot Propose |
|-------------|----------------|
| Neuroterus valhalla (62.6%) | Eschatocerus acacia (40.4%) |
| Melikaiella ostensackeni (63.3%) | Callirhytis spRG (46.8%) |
| All genomes ≥70% | Ganaspis sp (50.7%) |

**Implementation:** Add `--min-busco-for-proposal` flag to `05b_validate_novel_loci.py`

### 4.5 Phase 8 Display: Quality-Flagged Absences

**Problem:** Absences in low-quality genomes are likely false negatives. A "." in a 60% BUSCO genome doesn't mean "gene absent" - it means "gene not annotated."

**Solution:** Flag uncertain absences with `?` instead of `.`

```python
MIN_BUSCO_FOR_CONFIDENT_ABSENCE = 85.0  # Below this, absences are flagged as uncertain
```

**Display symbols:**

| Symbol | Meaning |
|--------|---------|
| `.` | Absent (genome ≥85% BUSCO, confident) |
| `?` | Absent (genome <85% BUSCO, may be false negative) |
| `*` | Target + synteny concordant (positive highlight) |
| `[X]` | Empty synteny block (flanking present, target absent, confident) |
| `[?]` | Empty synteny block (flanking present, target absent, low-quality genome) |

**Output changes:**
- **Phase 8a & 8b:** Add BUSCO% column on left side of matrix
- **Phase 8b:** Replace `.` with `?` for genomes below 85% BUSCO
- **Phase 8b:** Replace `[X]` with `[?]` for empty blocks in low-quality genomes

**Rationale:** Even 80% BUSCO means 20% of genes missing - absences can't be trusted. 85% is conservative threshold for "confident absence."

---

## 5. Summary: BUSCO Integration Decisions

After discussion, we decided on a **simple, targeted approach** rather than complex quality-weighting throughout the pipeline.

### 5.1 What We're NOT Doing

The following were considered but rejected as over-engineering:

- **Phase 5 threshold adjustment**: Current 15% flanking match threshold is already permissive; adjusting it based on genome quality adds complexity without clear benefit
- **Phase 6 length QC adjustment**: Low-quality genomes have missing genes, not bad genes; length ratios of found genes are still valid
- **Weighted genome counting**: Simple presence in N genomes is sufficient for validation

### 5.2 What We ARE Doing

| Section | Threshold | Action |
|---------|-----------|--------|
| 4.3 Genome Exclusion | 40% floor + redundancy | Physically exclude worst genomes |
| 4.4 Novel Loci Proposal | 60% minimum | Only good genomes can propose new loci |
| 4.5 Absence Flagging | 85% for confidence | Flag uncertain absences with `?` |

### 5.3 Rationale

**Presence we trust** - if we found a target or synteny block, it's real regardless of genome quality.

**Absence we doubt** - missing annotations in low-quality genomes could be false negatives, so flag them.

**Novel loci need good sources** - sparse annotations can't reliably anchor synteny patterns for new loci discovery.

---

## 6. Implementation Plan

### Already Complete
- `collect_busco_and_exclude.py` - Script written, ready to run when all BUSCO outputs available

### To Implement

1. **Add BUSCO to gca_to_species.tsv** (via `collect_busco_and_exclude.py --execute`)
2. **Update `05b_validate_novel_loci.py`** - Add `--min-busco-for-proposal` flag (default 60%)
3. **Update `08a_generate_locus_matrices_helixer.py`** - Add BUSCO% column
4. **Update `08b_generate_summary_matrices_helixer.py`** - Add BUSCO% column, implement `?` flagging for absences

---

## 7. Files to Modify

| File | Change |
|------|--------|
| `data/gca_to_species.tsv` | Add BUSCO columns via exclusion script |
| `05b_validate_novel_loci.py` | Add `--min-busco-for-proposal` flag |
| `08a_generate_locus_matrices_helixer.py` | Add BUSCO% column to output |
| `08b_generate_summary_matrices_helixer.py` | Add BUSCO% column, `?` flagging |

---

## 8. Next Steps

1. Wait for remaining Helixer + BUSCO jobs to complete
2. Run `collect_busco_and_exclude.py --execute` to apply exclusions
3. Implement Phase 5b and Phase 8 changes
4. Test on one family before full deployment

---

*Document updated 2025-12-02 after planning discussion.*
