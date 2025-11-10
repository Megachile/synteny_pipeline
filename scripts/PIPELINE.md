# Synteny Scanning Pipeline (Hybrid Workflow)

This document describes the complete synteny scanning pipeline for gene family analysis across cynipid wasp genomes.

## Overview

The pipeline identifies syntenic gene family members across genomes by:
1. Defining loci in reference genomes (BK/LB) with flanking proteins
2. Detecting synteny blocks via flanking protein conservation
3. Finding target genes within syntenic regions
4. Classifying targets by confidence and extracting sequences
5. Annotating and generating comparative matrices

---

## Phase 1: Locus Discovery & Target Extraction

**Script:** `scripts/locus_discovery.py`

**Inputs:**
- Gene family name and LOC IDs
- Reference proteomes: `data/proteomes/`
- Reference genomes: `data/genomes/` (GFF files for coordinate extraction)

**Outputs:**
- `outputs/<family>/phase1/locus_definitions.tsv` - Locus metadata (scaffolds, coordinates, tandem clusters)
- `outputs/<family>/phase1/<locus>/flanking_proteins.faa` - Flanking protein sequences (U12-U1, D1-D12)
- `outputs/<family>/phase1/<locus>/target_proteins.faa` - Target gene sequences from BK/LB

**Key Features:**
- Deduplicates isoforms and tandem clusters
- Extracts 12 upstream + 12 downstream flanking proteins per locus
- Automatically extracts target protein sequences from reference genomes
- Generates flanking comparison between BK and LB loci

**Run:**
```bash
python scripts/locus_discovery.py --gene-family <name> --locs <LOC_IDs>
```

---

## Phase 2: Genome-wide Synteny Detection

**Script:** `scripts/synteny_detection.py`

**Inputs:**
- `outputs/<family>/phase1/locus_definitions.tsv`
- `outputs/<family>/phase1/<locus>/flanking_proteins.faa`
- `data/ragtag_dbs/*.nhr` - tBLASTn databases for all genomes

**Outputs:**
- `outputs/<family>/phase2_synteny/<locus>/blast_xml/*.xml` - Raw BLAST results
- `outputs/<family>/phase2_synteny/<locus>/<genome>_synteny_blocks.tsv` - Detected blocks
- `outputs/<family>/phase2_synteny/summary_synteny_by_locus.tsv` - Summary stats

**Parameters:**
- `--min-proteins 2` - Minimum unique flanking matches per block
- `--max-gap-kb 150` - Maximum gap between flanking hits
- `--max-block-span-kb 300` - Maximum total block span

**Critical:** Run at LOCUS level, not family level. Each locus processes 75 genomes independently.

**Run:**
```bash
# Example SLURM array
sbatch --array=0-3 run_phase2_locus_array.slurm <family>
```

---

## Phase 3: Aggregate and Filter Synteny Blocks

**Script:** `scripts/aggregate_and_filter_blocks.py`

**Inputs:**
- `outputs/<family>/phase2_synteny/<locus>/<genome>_synteny_blocks.tsv` (all loci)

**Outputs:**
- `outputs/<family>/phase3_filtered/all_synteny_blocks.tsv` - Combined blocks across loci
- `outputs/<family>/phase3_filtered/synteny_blocks_filtered.tsv` - Best blocks per genome/locus
- `outputs/<family>/phase3_filtered/summary_aggregate_filter.tsv` - Filter statistics

**Strategy:**
- `--keep-all-blocks` mode: Retains all blocks passing thresholds (not just best per genome)
- Enables recovery of multi-locus hits

**Run:**
```bash
python scripts/aggregate_and_filter_blocks.py \
    --synteny-dir outputs/<family>/phase2_synteny \
    --output-dir outputs/<family>/phase3_filtered \
    --keep-all-blocks
```

---

## Phase 4: Target Gene Detection

**Script:** `scripts/blast_targets.py`

**Inputs:**
- `outputs/<family>/phase1/<locus>/target_proteins.faa` (from all loci)
- `data/ragtag_dbs/*.nhr` - Genome tBLASTn databases

**Outputs:**
- `outputs/<family>/phase4_targets/all_target_loci.tsv` - All detected target hits
- `outputs/<family>/phase4_targets/combined_targets.faa` - Combined query sequences
- `outputs/<family>/phase4_targets/blast_xml/*.xml` - Raw BLAST results

**Run:**
```bash
python scripts/blast_targets.py \
    --locus-defs outputs/<family>/phase1/locus_definitions.tsv \
    --phase1-dir outputs/<family>/phase1 \
    --output-dir outputs/<family>/phase4_targets
```

---

## Phase 5: Classify Targets (3-Way Classification)

**Script:** `scripts/classify_targets.py`

**Inputs:**
- `outputs/<family>/phase4_targets/all_target_loci.tsv`
- `outputs/<family>/phase3_filtered/synteny_blocks_filtered.tsv`

**Outputs:**
- `outputs/<family>/phase5_classified/syntenic_targets.tsv` - Targets within synteny blocks
- `outputs/<family>/phase5_classified/unplaceable_targets.tsv` - Targets outside blocks
- `outputs/<family>/phase5_classified/all_targets_classified.tsv` - Combined with confidence
- `outputs/<family>/phase5_classified/dropped_targets.tsv` - Failed assignment
- `outputs/<family>/phase5_classified/debug_multilocus_candidates.tsv` - Ambiguous cases

**Classification System:**

1. **Syntenic:** Targets within synteny blocks
   - `single_unambiguous` - Clearly assigned to one locus
   - `multiple_ambiguous` - Could belong to multiple loci (NOT YET IMPLEMENTED)

2. **Unplaceable:** Targets outside synteny blocks
   - `near_block` - Within 200kb of a block (ambiguous/multilocus candidates)
   - `confidently_distinct` - >200kb from any block (truly distinct)

**Strategy:**
- Locus-first assignment prevents cross-locus contamination
- Distance-based confidence scoring for unplaceables

**Run:**
```bash
python scripts/classify_targets.py \
    --targets outputs/<family>/phase4_targets/all_target_loci.tsv \
    --blocks outputs/<family>/phase3_filtered/synteny_blocks_filtered.tsv \
    --output-dir outputs/<family>/phase5_classified \
    --keep-all-blocks
```

---

## Phase 6: Extract Gene Sequences with Exonerate

**Script:** `scripts/extract_sequences.py`

**Inputs:**
- `outputs/<family>/phase5_classified/syntenic_targets.tsv`
- `outputs/<family>/phase4_targets/combined_targets.faa` - Query proteins
- **Genome paths:**
  - Most genomes: `data/ragtag_output/<genome>/ragtag.scaffold.fasta`
  - **BK (Belonocnema_kinseyi_GCF):** Uses NCBI RefSeq with NC_ scaffolds
  - **LB genomes:** May use NCBI RefSeq with NW_ scaffolds

**Outputs:**
- `outputs/<family>/phase6_extracted/<genome>/<locus>/` - Per-locus extractions:
  - `*_gene1_cds.fasta` - CDS sequence
  - `*_gene1_protein.fasta` - Protein sequence
  - `*_gene1_genomic.fasta` - Genomic DNA
  - `*_exonerate_flank0.txt` - Exonerate alignment details
- `outputs/<family>/phase6_extracted/all_extracted_genes.faa` - Aggregated proteins

**Critical Requirements:**
- **Exonerate module:** Must load `module load exonerate/2.4.0-yl7q`
- **Genome scaffold matching:** Phase 5 targets must reference scaffolds that exist in genome files
- **Special handling for BK:** Use NCBI genome instead of RagTag assembly

**Quality Grades:**
- `intact` - Complete CDS, no stop codons
- `pseudogene` - Contains internal stop codons
- `fragment` - <90% query length

**Run:**
```bash
sbatch run_phase6_with_ncbi.slurm
```

---

## Phase 7: SwissProt Annotation

**Script:** `scripts/swissprot_annotation.py`

**Inputs:**
- `outputs/<family>/phase2_synteny/<locus>/<genome>_hits.tsv` - Flanking protein hits
- `outputs/<family>/phase3_filtered/synteny_blocks_filtered.tsv`
- SwissProt database

**Outputs:**
- `outputs/<family>/phase7_swissprot/<locus>_swissprot_annotations.tsv` - Per-locus annotations
- `outputs/<family>/phase7_swissprot/all_swissprot_combined.tsv` - Combined annotations

**Run:**
```bash
# Array job over loci
sbatch --array=0-N run_phase7_swissprot_array.slurm
```

---

## Phase 8a: Generate Locus Matrices

**Script:** `scripts/generate_locus_matrices.py`

**Inputs:**
- `outputs/<family>/phase1/locus_definitions.tsv`
- `outputs/<family>/phase2_synteny/` - Flanking protein hits
- `outputs/<family>/phase3_filtered/synteny_blocks_filtered.tsv`
- `outputs/<family>/phase5_classified/all_targets_classified.tsv`
- `outputs/<family>/phase6_extracted/graded_syntenic_targets.tsv` - Graded targets
- `outputs/<family>/phase7_swissprot/all_swissprot_combined.tsv`
- `outputs/<family>/phase6_extracted/` - Extracted sequences (for metadata)
- `data/gca_to_species.tsv` - Species mapping

**Outputs:**
- `outputs/<family>/phase8_matrices/<locus>_genome_swissprot_matrix.tsv` - Per-locus matrices

**Matrix Structure:**
- Rows: All genomes (sorted by phylogenetic order)
- Columns: Metadata (genome, species, phylo_order, synteny_pct, scaffold, coordinates) + flanking proteins (U12-U1, TARGET, D1-D12)
- Synteny %: Percentage of flanking proteins found
- Target notation: `[304I]` (length + status: I=intact, P=pseudogene, F=fragment)

**Run:**
```bash
python scripts/generate_locus_matrices.py \
    --locus-defs outputs/<family>/phase1/locus_definitions.tsv \
    --synteny-dir outputs/<family>/phase2_synteny \
    --blocks outputs/<family>/phase3_filtered/synteny_blocks_filtered.tsv \
    --graded-syntenic outputs/<family>/phase6_extracted/graded_syntenic_targets.tsv \
    --targets outputs/<family>/phase5_classified/all_targets_classified.tsv \
    --swissprot outputs/<family>/phase7_swissprot/all_swissprot_combined.tsv \
    --species-map data/gca_to_species.tsv \
    --extracted-seqs outputs/<family>/phase6_extracted \
    --output-dir outputs/<family>/phase8_matrices
```

---

## Phase 8b: Generate Summary Matrices

**Script:** `scripts/generate_summary_matrices.py`

**Inputs:**
- `outputs/<family>/phase1/locus_definitions.tsv`
- `outputs/<family>/phase3_filtered/synteny_blocks_filtered.tsv`
- `outputs/<family>/phase5_classified/all_targets_classified.tsv`
- `outputs/<family>/phase6_extracted/` - For metadata
- `data/gca_to_species.tsv`
- `outputs/<family>/phase8_matrices/<locus>_genome_swissprot_matrix.tsv` - For synteny percentages

**Outputs:**
- `outputs/<family>/phase8_matrices/<family>_summary_matrix.tsv` - Family-level summary

**Matrix Structure:**
- Rows: All genomes (sorted by phylogenetic order)
- Columns: Metadata + per-locus columns + `<family>_near_block` + `<family>_confidently_distinct`
- Locus columns show: `synteny% [grades]` (e.g., `30.3% [304I; 295P]`)
- Unplaceable columns show: `[hit; hit]` (count of unplaceable targets)

**Summary Statistics:**
- Syntenic: Targets within blocks (with grades)
- Ambiguous (near_block): Targets within 200kb of blocks
- Distinct (confidently_distinct): Targets >200kb from blocks

**Run:**
```bash
python scripts/generate_summary_matrices.py \
    --locus-defs outputs/<family>/phase1/locus_definitions.tsv \
    --blocks outputs/<family>/phase3_filtered/synteny_blocks_filtered.tsv \
    --targets outputs/<family>/phase5_classified/all_targets_classified.tsv \
    --species-map data/gca_to_species.tsv \
    --extracted-seqs outputs/<family>/phase6_extracted \
    --output-dir outputs/<family>/phase8_matrices
```

---

## Grading Step: Quality Assessment

**Script:** `scripts/grade_matches.py`

**Run between Phase 6 and Phase 8:**

```bash
python scripts/grade_matches.py \
    --syntenic outputs/<family>/phase5_classified/syntenic_targets.tsv \
    --extracted-dir outputs/<family>/phase6_extracted \
    --query-proteins outputs/<family>/phase4_targets/combined_targets.faa \
    --output outputs/<family>/phase6_extracted/graded_syntenic_targets.tsv
```

**Grades:**
- `intact` - Full-length, no stop codons
- `degraded_pseudogene` - Internal stop codons
- `degraded_fragment` - <90% query length
- `degraded_no_cds` - Exonerate extraction failed

---

## Critical Genome Handling Notes

### BK (Belonocnema kinseyi)
- **Genome code:** `Belonocnema_kinseyi_GCF` (used in targets/blocks)
- **Alternate code:** `GCA_010883055.1` (in species map)
- **Scaffold type:** NC_ (NCBI RefSeq chromosomes)
- **Path:** `/carc/scratch/projects/emartins/2016456/adam/genomes/belonocnema_kinseyi/GCF_010883055.1_B_treatae_v1_genomic.fna`
- **Why:** Phase 5 targets reference NC_ scaffolds, but RagTag assembly has CM_ scaffolds
- **Solution:** Use NCBI genome directly in Phase 6

### LB (Leptopilina boulardi)
- **Genome codes:** `GCA_003121605.1`, `GCA_011634795.1`, `GCA_015476485.1`, `GCA_019393585.1`, `GCA_032872485.1`
- **May require NCBI genomes** if RagTag assemblies have mismatched scaffolds

### Species Mapping
- **File:** `data/gca_to_species.tsv`
- **Format:** `accession<TAB>species<TAB>family<TAB>phylo_order`
- **Must include:** All genome codes used in pipeline, including alternate names

---

## Parallelization Strategy

| Phase | Parallelization | Method |
|-------|----------------|--------|
| 1 | Serial | Single family at a time |
| 2 | SLURM array | Per locus (4-6 loci per family) |
| 3 | Serial | Quick aggregation/filtering |
| 4 | Serial/Array | Can parallelize by genome if needed |
| 5 | Serial | Classification is fast |
| 6 | SLURM array | Per family (can do all families in one array) |
| 7 | SLURM array | Per locus |
| 8 | Serial | Matrix generation is fast |

---

## Key Thresholds and Parameters

| Parameter | Value | Purpose |
|-----------|-------|---------|
| `--min-proteins` | 2 | Minimum flanking matches per block |
| `--max-gap-kb` | 150 | Maximum gap between flanking hits |
| `--max-block-span-kb` | 300 | Maximum synteny block span |
| `--near-block-distance-kb` | 200 | Threshold for ambiguous unplaceables |
| `--intact-threshold` | 0.90 | Minimum query coverage for intact grade |
| `--unplaceable-evalue` | 1e-10 | E-value for unplaceable target search |

---

## Validation and Troubleshooting

### Common Issues

1. **Phase 6 "degraded_no_cds" for all targets:**
   - Check genome paths
   - Verify scaffold names match between Phase 5 targets and genome files
   - Ensure exonerate module is loaded

2. **Low BK/LB recovery in Phase 5:**
   - Use `--keep-all-blocks` in Phase 3
   - Verify locus-first classification in Phase 5

3. **Phase 8 showing 0% synteny:**
   - Regenerate locus matrices (Phase 8a) before summaries (Phase 8b)
   - Check that synteny blocks file has `strand` column

4. **Missing genomes in Phase 8:**
   - Add genome codes to `data/gca_to_species.tsv`
   - Regenerate both Phase 8a and 8b

### Output Validation

Check these files after each phase:
- Phase 1: `locus_definitions.tsv` should have all expected loci
- Phase 2: `summary_synteny_by_locus.tsv` should show blocks for most genomes
- Phase 3: `summary_aggregate_filter.tsv` should show retention rates
- Phase 5: `all_targets_classified.tsv` should have `confidence` column
- Phase 6: Count `*_gene1_cds.fasta` files for successful extractions
- Phase 8: Check synteny percentages and target grades in matrices

---

## Environment Requirements

**Conda Environment:** `busco_env`
- pandas, numpy, biopython
- matplotlib, seaborn (for plotting)

**Modules:**
- `exonerate/2.4.0-yl7q` (Phase 6 only)

**Databases:**
- tBLASTn DBs: `data/ragtag_dbs/`
- Reference genomes: `data/ragtag_output/` and special NCBI paths for BK/LB
- SwissProt: (path to be documented)

---

## Quick Start for New Gene Family

```bash
# 1. Locus discovery
python scripts/locus_discovery.py --gene-family my_family --locs LOC1,LOC2,LOC3

# 2. Synteny detection (array over loci)
sbatch --array=0-2 run_phase2_locus_array.slurm my_family

# 3. Aggregate and filter
python scripts/aggregate_and_filter_blocks.py \
    --synteny-dir outputs/my_family/phase2_synteny \
    --output-dir outputs/my_family/phase3_filtered \
    --keep-all-blocks

# 4. Target detection
python scripts/blast_targets.py \
    --locus-defs outputs/my_family/phase1/locus_definitions.tsv \
    --phase1-dir outputs/my_family/phase1 \
    --output-dir outputs/my_family/phase4_targets

# 5. Classify
python scripts/classify_targets.py \
    --targets outputs/my_family/phase4_targets/all_target_loci.tsv \
    --blocks outputs/my_family/phase3_filtered/synteny_blocks_filtered.tsv \
    --output-dir outputs/my_family/phase5_classified \
    --keep-all-blocks

# 6. Extract sequences
sbatch run_phase6_with_ncbi.slurm my_family

# 7. Grade
python scripts/grade_matches.py \
    --syntenic outputs/my_family/phase5_classified/syntenic_targets.tsv \
    --extracted-dir outputs/my_family/phase6_extracted \
    --query-proteins outputs/my_family/phase4_targets/combined_targets.faa \
    --output outputs/my_family/phase6_extracted/graded_syntenic_targets.tsv

# 8a. Locus matrices
python scripts/generate_locus_matrices.py \
    --locus-defs outputs/my_family/phase1/locus_definitions.tsv \
    --synteny-dir outputs/my_family/phase2_synteny \
    --blocks outputs/my_family/phase3_filtered/synteny_blocks_filtered.tsv \
    --graded-syntenic outputs/my_family/phase6_extracted/graded_syntenic_targets.tsv \
    --targets outputs/my_family/phase5_classified/all_targets_classified.tsv \
    --swissprot outputs/my_family/phase7_swissprot/all_swissprot_combined.tsv \
    --species-map data/gca_to_species.tsv \
    --extracted-seqs outputs/my_family/phase6_extracted \
    --output-dir outputs/my_family/phase8_matrices

# 8b. Summary matrix
python scripts/generate_summary_matrices.py \
    --locus-defs outputs/my_family/phase1/locus_definitions.tsv \
    --blocks outputs/my_family/phase3_filtered/synteny_blocks_filtered.tsv \
    --targets outputs/my_family/phase5_classified/all_targets_classified.tsv \
    --species-map data/gca_to_species.tsv \
    --extracted-seqs outputs/my_family/phase6_extracted \
    --output-dir outputs/my_family/phase8_matrices
```

---

Last updated: 2025-01-10
Pipeline version: Consolidated workflow with 3-way classification
