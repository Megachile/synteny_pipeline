# Helixer Pipeline Phase 6-8 Fix Plan

## Problem Summary

The current Helixer Phase 6-7 and Phase 8 scripts produce outputs that don't match the original tblastn pipeline format:

| Aspect | Original Pipeline | My Helixer Version |
|--------|-------------------|-------------------|
| **Phase 6 structure** | `{genome}/{locus}/gene1_protein.fasta` | Single `target_proteins.faa` |
| **Phase 6 headers** | Metadata-rich: `length:169aa cov:95.1%` | Just gene IDs |
| **Phase 8a output** | 24+ per-locus matrices with U/D flanking columns | 2 simple TSV files |
| **Phase 8b output** | Family summary with synteny % per locus column | Non-existent |

## Root Cause

I rewrote the scripts from scratch instead of adapting the existing pipeline scripts to work with Helixer inputs. The original scripts expect specific input formats and produce specific output formats that downstream analyses depend on.

## Fix Strategy

**Adapt, don't rewrite.** Use the original Phase 6 and Phase 8 scripts as templates, modifying only the input parsing to handle Helixer's data format.

---

## Phase 6 Fix: `06_extract_sequences_helixer.py`

### Current Problem
- Outputs a single `target_proteins.faa` file
- No per-genome/per-locus organization
- No metadata in FASTA headers

### Required Output Structure
```
phase6_helixer_filtered/
├── Andricus_aestivalis_GCA/
│   ├── BK_chr1_a/
│   │   └── BK_chr1_a_gene1_protein.fasta
│   └── BK_chr2_c/
│       └── BK_chr2_c_gene1_protein.fasta
├── Andricus_quercusirregularis_GCA/
│   └── ...
└── all_extracted_genes.faa
```

### What to Adapt from Original

**From `06_extract_sequences.py`:**
1. Directory creation logic (lines 1431-1448):
   - Syntenic targets: `{output_dir}/{genome}/{assigned_to}/`
   - Unplaceable targets: `{output_dir}/{genome}/UNPLACED/{unique_tag}/`

2. FASTA header format (lines 967-1004):
   ```
   >gene1 {scaffold}:{start}-{end} strand:{strand} length:{length}aa query_len:{qlen}aa cov:{cov}%
   ```

3. Aggregation logic (lines 1526-1542):
   - Collect all `*_protein.fasta` files into `all_extracted_genes.faa`

### Key Differences for Helixer

1. **No Exonerate needed** - Helixer already has complete gene models
2. **Protein extraction** - Just fetch from `{genome}_helixer_proteins.faa`
3. **Coordinates** - Get from Phase 4 output (`all_target_loci.tsv`)
4. **Classification** - Get from Phase 5 (`syntenic_targets.tsv`, `unplaceable_targets.tsv`)

### Implementation Plan

```python
# Core logic for 06_extract_sequences_helixer.py

def extract_and_organize(phase4_dir, phase5_dir, helixer_dir, output_dir):
    # 1. Load Phase 5 classifications
    syntenic_df = load_syntenic_targets(phase5_dir)
    unplaceable_df = load_unplaceable_targets(phase5_dir)

    # 2. For each syntenic target:
    for target in syntenic_df:
        # Get assigned locus from Phase 5
        assigned_locus = target['assigned_locus']
        genome = target['genome']
        gene_id = target['target_gene_id']

        # Create directory: {output}/{genome}/{assigned_locus}/
        target_dir = output_dir / genome / assigned_locus
        target_dir.mkdir(parents=True, exist_ok=True)

        # Extract protein from Helixer proteome
        protein_seq = get_protein_from_helixer(helixer_dir, genome, gene_id)

        # Get coordinates from Phase 4
        coords = get_coordinates(phase4_dir, genome, gene_id)

        # Write with proper header
        fasta_path = target_dir / f"{assigned_locus}_gene1_protein.fasta"
        write_protein_fasta(fasta_path, gene_id, protein_seq, coords)

    # 3. For unplaceable targets:
    for target in unplaceable_df:
        # Similar but under UNPLACED/unique_tag
        ...

    # 4. Aggregate all proteins
    aggregate_all_proteins(output_dir)
```

---

## Phase 8a Fix: `08a_generate_locus_matrices_helixer.py`

### Current Problem
- Doesn't exist as proper script
- My Phase 8 script produces wrong format entirely

### Required Output Format
One TSV per locus: `{locus_id}_genome_swissprot_matrix.tsv`

Columns:
```
genome_id | species | phylo_order | synteny_pct | num_proteins_found | scaffold | strand | start | end | U14_... | U13_... | ... | U1_... | TARGET | D1_... | D2_... | ...
```

### What to Adapt from Original

**From `08a_generate_locus_matrices_core.py`:**

1. The `create_locus_matrix()` function (lines 227-553) is the core
2. Key inputs needed:
   - `locus_definitions.tsv` - defines BK/LB loci (Phase 1)
   - `synteny_blocks_filtered.tsv` - flanking gene matches (Phase 3)
   - `syntenic_targets.tsv` - target classifications (Phase 5)
   - `genome_specific_swissprot_annotations.tsv` - flanking gene annotations (Phase 7)
   - Extracted sequences directory - for length/status metadata (Phase 6)

3. The script reads:
   - Flanking proteins from Phase 1: `{locus}_flanking_dedup.faa`
   - SwissProt annotations from Phase 7
   - Target metadata from Phase 6 extraction

### Key Differences for Helixer

1. **Flanking genes** - Come from Phase 4b (`phase4b_helixer_filtered/`) instead of Phase 2
2. **Target info** - Come from Phase 5 Helixer (`phase5_helixer_filtered/`)
3. **SwissProt annotations** - Need to run Phase 7 for Helixer flanking genes

### Implementation Plan

The simplest fix is to **call the existing 08a script** with paths pointing to Helixer outputs:

```bash
python redesign_scripts/08a_generate_locus_matrices_core.py \
    --locus-defs outputs/$FAMILY/phase1_v2/locus_definitions.tsv \
    --synteny-dir outputs/$FAMILY/phase4b_helixer_filtered \
    --blocks outputs/$FAMILY/phase3_aggregated_v2/synteny_blocks_filtered.tsv \
    --targets outputs/$FAMILY/phase5_helixer_filtered/all_targets_classified.tsv \
    --swissprot outputs/$FAMILY/phase67_helixer_filtered/genome_specific_swissprot_annotations.tsv \
    --extracted-seqs outputs/$FAMILY/phase6_helixer_filtered \
    --species-map data/gca_to_species_order.tsv \
    --output-dir outputs/$FAMILY/phase8_helixer_filtered
```

BUT the original script expects specific column names and file structures. We need to:
1. Ensure Phase 5 Helixer outputs have same columns as original Phase 5
2. Ensure Phase 6 Helixer has same directory structure
3. Create proper `all_targets_classified.tsv` with both syntenic and unplaceable

---

## Phase 8b Fix: `08b_generate_summary_matrices_helixer.py`

### Required Output Format
One TSV per gene family: `{gene_family}_summary_matrix.tsv`

Columns:
```
genome_id | species | phylo_order | BK_chr1_a | BK_chr2_c | ... | {family}_unplaceable | total
```

Each cell shows: `{synteny_pct}% [{length1}; {length2}]` or `[empty]` or `[not found]`

### Implementation Plan

Same approach - call existing `08b_generate_summary_matrices_core.py` with Helixer paths:

```bash
python redesign_scripts/08b_generate_summary_matrices_core.py \
    --locus-defs outputs/$FAMILY/phase1_v2/locus_definitions.tsv \
    --blocks outputs/$FAMILY/phase3_aggregated_v2/synteny_blocks_filtered.tsv \
    --targets outputs/$FAMILY/phase5_helixer_filtered/all_targets_classified.tsv \
    --species-map data/gca_to_species_order.tsv \
    --extracted-seqs outputs/$FAMILY/phase6_helixer_filtered \
    --locus-matrices-dir outputs/$FAMILY/phase8_helixer_filtered \
    --output-dir outputs/$FAMILY/phase8_helixer_filtered
```

---

## Required Changes to Phase 5 Helixer

To make Phase 8 work, Phase 5 needs to output `all_targets_classified.tsv` with columns:
- `genome`
- `locus_name` (or `locus_id`)
- `parent_locus`
- `placement` ('synteny' or 'unplaceable')
- `assigned_to` (the BK/LB locus for syntenic targets)
- `scaffold`, `start`, `end`, `strand`
- `best_evalue`
- `gene_family`

---

## Implementation Order

1. **Fix Phase 5 output format** - Ensure `all_targets_classified.tsv` has required columns
2. **Rewrite Phase 6** - Create per-genome/per-locus structure with proper headers
3. **Test Phase 8a/8b** - Use existing scripts with Helixer paths
4. **Create wrapper SLURM script** - `run_helixer_phases_6_8_fixed.slurm`

---

## Files to Create/Modify

| File | Action | Purpose |
|------|--------|---------|
| `05_classify_targets_helixer.py` | Modify | Add `all_targets_classified.tsv` output |
| `06_extract_sequences_helixer.py` | Rewrite | Per-genome/per-locus organization |
| `08_helixer_summary_matrices.py` | Delete | Replace with calls to original scripts |
| `run_helixer_phases_6_8_array.slurm` | Modify | Call original 08a/08b with Helixer paths |

---

## Estimated Work

1. Phase 5 output fix: ~30 min
2. Phase 6 rewrite: ~2 hours (new script, adapting extraction logic)
3. Phase 8 integration: ~30 min (SLURM script modifications)
4. Testing: ~1 hour

Total: ~4 hours
