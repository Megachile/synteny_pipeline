# Synteny Scanning Hybrid Workflow - Blueprint

## Overview

This pipeline detects gene families across 75 wasp genomes using synteny-based comparative genomics. It identifies paralogs in landmark genomes (BK, LB, TR, DR), uses flanking proteins to detect syntenic regions in target genomes, then extracts and annotates target genes.

## Pipeline Phases

### Phase 1: Paralog Discovery (OPTIONAL - skip if using pre-built input)
**Script:** `scripts/01_discover_paralogs_gene_level.py`

**Purpose:** Find paralogs of a query LOC ID across 4 landmark genomes

**Inputs:**
- Query LOC ID or protein ID (e.g., `LOC117167432`)
- Landmark proteomes and GFFs in `data/proteomes/` and `data/genomes/`

**Outputs:**
- `<output_dir>/unique_loci.tsv` - Discovered loci across landmark genomes
- `<output_dir>/loc_<LOC>/flanking.faa` - Flanking proteins for each locus

**Parallelization:** Can run per-LOC in parallel

**Notes:** This phase generates the `locus_definitions.tsv` file that subsequent phases need. If you already have `locus_definitions.tsv`, skip this phase.

---

### Phase 2: Synteny Deduplication (OPTIONAL - deprecated)
**Status:** Currently skipped - deduplication happens in Phase 1

**Purpose:** Was intended to merge duplicate loci detected across landmarks

**Notes:** Phase 1 now handles this internally. This phase exists for backwards compatibility but is not used.

---

### Phase 3: Synteny Detection
**Script:** `scripts/02_synteny_detection.py`

**Purpose:** Use flanking proteins to find syntenic regions across 75 target genomes

**Inputs:**
- `SYNTENY_INPUT_FILE` env var → locus definitions TSV (columns: `locus_id`, `genome`, `chromosome`, `flanking_file`, `gene_family`)
- Flanking protein FASTA files (from Phase 1 or pre-built)
- RagTag genome databases in `data/ragtag_dbs/`

**Method:**
1. Run tBLASTn with flanking proteins against each target genome
2. Cluster BLAST hits by genomic proximity (500kb max gap)
3. Filter blocks requiring ≥5 query proteins matching

**Outputs:**
- `02_synteny_blocks/GCA_*/locus_<locus>_hits.tsv` - Raw BLAST hits per genome
- `02_synteny_blocks/GCA_*/locus_<locus>_blocks.tsv` - Synteny blocks per genome
- `02_synteny_blocks/all_synteny_blocks.tsv` - All blocks combined

**Parallelization:** Use `--locus <locus_id>` for SLURM arrays

**Key Parameters:**
- `FLANKING_BLAST_EVALUE = 1e-5`
- `SYNTENY_MAX_GAP_KB = 500` (max gap between proteins in same block)
- `MIN_PROTEINS_FOR_BLOCK = 5` (min query proteins for valid block)

---

### Phase 4: Target Gene Detection + Clustering
**Script:** `scripts/04_blast_targets.py`

**Purpose:** Find the actual target gene within each synteny block

**Inputs:**
- Synteny blocks from Phase 3
- BK query proteins from `data/reference/protein.faa`
- RagTag genome databases

**Method:**
1. Extract query protein from BK reference locus
2. Run tBLASTn for query protein against each genome
3. Cluster hits by reading frame and genomic proximity (10kb max gap, same reading frame)
4. Assign hit clusters to nearest synteny blocks

**Outputs:**
- `04_target_genes/combined_targets.faa` - All unique target queries
- `04_target_genes/GCA_*/target_hits.tsv` - Target BLAST hits
- `04_target_genes/GCA_*/target_clusters.tsv` - Clustered targets by frame
- `04_target_genes/all_targets.tsv` - All target clusters assigned to blocks

**Parallelization:** Use `--locus <locus_id>` for SLURM arrays

**Key Parameters:**
- `TARGET_BLAST_EVALUE = 1e-5`
- `TARGET_CLUSTER_GAP_KB = 10` (max gap for clustering hits)
- Reading frames: +1, +2, +3, -1, -2, -3 (clusters separately!)

**Critical Fix (Nov 3):** Added frame-aware clustering - targets in different reading frames are now treated as separate clusters

---

### Phase 5: Classify Targets
**Script:** `scripts/05_classify_targets.py`

**Purpose:** Classify each target as "syntenic" or "unplaceable"

**Inputs:**
- `04_target_genes/all_targets.tsv`
- `02_synteny_blocks/all_synteny_blocks.tsv`

**Method:**
- Check if target cluster overlaps or is near (≤500kb) a synteny block
- Label as "syntenic" if near block, "unplaceable" otherwise

**Outputs:**
- `05_classified/classified_targets.tsv` - All targets with classification
- `05_classified/syntenic_targets.tsv` - Filtered syntenic targets only

**Parallelization:** None (fast, runs on all targets at once)

---

### Phase 6: Extract Gene Structures
**Script:** `scripts/06_extract_sequences.py`
**Modules:** `scripts/extract_with_exonerate.py`, `scripts/exonerate_extract.py`

**Purpose:** Extract complete gene structures (CDS, genomic, protein) using Exonerate

**Inputs:**
- `05_classified/syntenic_targets.tsv`
- RagTag genome FASTA files in `data/ragtag_output/`
- BK query proteins from `data/reference/protein.faa`

**Method:**
1. For each target cluster, extract query protein
2. Run Exonerate protein2genome alignment with **adaptive windowing**:
   - Try window sizes: 0kb, 1kb, 5kb, 10kb, 15kb, 20kb, 50kb, 100kb, 200kb flanking
   - Stop when ≥90% of query protein is covered
   - Keep best result if none reach 90%
3. Parse Exonerate GFF output to extract gene features
4. Group exons into genes (fixed Nov 3 - was treating exons as separate genes!)
5. Extract CDS sequence, genomic sequence (with introns), and translate protein

**Outputs:**
- `04_target_genes/extracted_sequences/GCA_*/<block_id>/`
  - `<block>_gene<id>_cds.fasta` - CDS sequence
  - `<block>_gene<id>_genomic.fasta` - Full gene with introns
  - `<block>_gene<id>_protein.fasta` - Translated protein
  - `<block>_exonerate_flank<size>.txt` - Exonerate alignment details

**Parallelization:** None currently (processes all targets sequentially)

**Key Parameters:**
- Adaptive window sizes: 0, 1kb, 5kb, 10kb, 15kb, 20kb, 50kb, 100kb, 200kb flanking
- Completeness threshold: 90% query coverage
- Exonerate model: protein2genome

**Critical Fixes (Nov 3):**
- Fixed gene grouping bug (exons → single gene)
- Implemented adaptive windowing (handles genes with large introns)
- Added genomic + protein sequence outputs

---

### Phase 7: SwissProt Annotation
**Script:** `scripts/07_swissprot_annotation.py`

**Purpose:** Annotate all extracted hit sequences with SwissProt for functional context

**Inputs:**
- All `*_hitseq.faa` files from `02_synteny_blocks/GCA_*/`
- SwissProt DIAMOND database at `data/databases/swissprot.dmnd`

**Method:**
1. Collect all hit sequences from Phase 3 (flanking protein hits)
2. Run DIAMOND blastp against SwissProt
3. Keep best hit per protein

**Outputs:**
- `07_swissprot/all_hits.faa` - All protein sequences
- `07_swissprot/swissprot_results.tsv` - DIAMOND results
- `07_swissprot/annotations.tsv` - Best annotation per protein

**Parallelization:** None (single job processes all)

**Key Parameters:**
- `SWISSPROT_BLAST_EVALUE = 1e-5`
- `USE_DIAMOND = True` (DIAMOND is much faster than BLASTP)

**Status:** Completed for ferritin (293,478 proteins annotated, 100% match rate)

---

### Phase 8a: Generate Locus Matrices
**Script:** `scripts/08a_generate_locus_matrices.py`

**Purpose:** Create presence/absence matrices for each locus

**Inputs:**
- `02_synteny_blocks/all_synteny_blocks.tsv`
- `04_target_genes/all_targets.tsv`
- `05_classified/classified_targets.tsv`
- `07_swissprot/annotations.tsv`
- Extracted gene sequences from Phase 6

**Method:**
1. For each locus, create a matrix with:
   - Rows = genomes (sorted phylogenetically)
   - Columns = flanking proteins + target gene
2. Fill cells with:
   - Flanking: presence (1) or absence (0) + SwissProt annotation
   - Target: gene length + functional status (complete/partial/pseudogene)

**Outputs:**
- `07_matrices/locus_specific/<locus>_matrix.tsv` - One matrix per locus

**Parallelization:** None (fast, processes all loci)

---

### Phase 8b: Generate Summary Matrices
**Script:** `scripts/08b_generate_summary_matrices.py`

**Purpose:** Aggregate locus-level results into gene family summaries

**Inputs:**
- All locus matrices from Phase 8a
- Locus definitions

**Method:**
1. Aggregate all loci from same gene family
2. Count gene types per genome (functional, partial, pseudogene, absent)
3. Calculate summary statistics

**Outputs:**
- `07_matrices/gene_type_summaries/<gene_family>_summary.tsv`
- Overall statistics

**Parallelization:** None

---

## Input File Format

### For Phase 1 (paralog discovery):
Individual LOC IDs passed as command-line arguments

### For Phases 3+ (if skipping Phase 1):
**File:** `locus_definitions.tsv`

**Format:**
```tsv
locus_id	genome	chromosome	flanking_file	gene_family
BK_chr2_a	BK	NC_046658.1	outputs/phase12/loc_LOC117167432/BK_chr2_a_flanking.faa	ferritin_MC102
LB_scf7864_a	LB	NW_026137864.1	outputs/phase12/loc_LOC117167432/LB_scf7864_a_flanking.faa	ferritin_MC102
```

**Columns:**
- `locus_id`: Unique identifier for this locus (format: `<genome>_<chr>_<letter>`)
- `genome`: Landmark genome code (BK, LB, TR, DR)
- `chromosome`: Chromosome/scaffold ID
- `flanking_file`: Path to flanking proteins FASTA
- `gene_family`: Gene family name (e.g., `ferritin_MC102`)

### For batch processing (32 gene families):
**File:** `hpc_input_loci.tsv`

**Format:**
```tsv
gene_family	target_genes	num_locs
MC117	LOC117172449,LOC117172450,...	28
ferritin_MC102	LOC117167432	1
```

**Columns:**
- `gene_family`: Name of gene family
- `target_genes`: Comma-separated LOC IDs to discover
- `num_locs`: Expected number of loci (for validation)

---

## Output Organization

### Standard output structure (per gene family):
```
outputs/<gene_family>_<timestamp>/
├── 02_synteny_blocks/
│   ├── GCA_<accession>/
│   │   ├── locus_<id>_hits.tsv
│   │   ├── locus_<id>_blocks.tsv
│   │   └── <locus>_<qseqid>_hitseq.faa
│   └── all_synteny_blocks.tsv
├── 04_target_genes/
│   ├── combined_targets.faa
│   ├── GCA_<accession>/
│   │   ├── target_hits.tsv
│   │   └── target_clusters.tsv
│   ├── all_targets.tsv
│   └── extracted_sequences/
│       └── GCA_<accession>/
│           └── <block_id>/
│               ├── <block>_gene<id>_cds.fasta
│               ├── <block>_gene<id>_genomic.fasta
│               ├── <block>_gene<id>_protein.fasta
│               └── <block>_exonerate_flank<size>.txt
├── 05_classified/
│   ├── classified_targets.tsv
│   └── syntenic_targets.tsv
├── 07_swissprot/
│   ├── all_hits.faa
│   ├── swissprot_results.tsv
│   └── annotations.tsv
└── 07_matrices/
    ├── locus_specific/
    │   └── <locus>_matrix.tsv
    └── gene_type_summaries/
        └── <gene_family>_summary.tsv
```

---

## Key Configuration (config.py)

### Genome databases:
- `RAGTAG_DB_DIR` = `data/ragtag_dbs/` (tBLASTn databases for 75 genomes)
- `RAGTAG_FASTA_DIR` = `data/ragtag_output/` (Genome FASTA files)

### Reference data:
- `BK_GFF_FILE` = `data/reference/genomic.gff`
- `BK_PROTEINS_FILE` = `data/reference/protein.faa`

### BLAST parameters:
- `FLANKING_BLAST_EVALUE = 1e-5`
- `TARGET_BLAST_EVALUE = 1e-5`
- `SWISSPROT_BLAST_EVALUE = 1e-5`
- `BLAST_THREADS = 16`

### Synteny parameters:
- `SYNTENY_MAX_GAP_KB = 500` (max gap between proteins in block)
- `TARGET_CLUSTER_GAP_KB = 10` (max gap for target clustering)
- `MIN_PROTEINS_FOR_BLOCK = 5` (min query proteins for valid block)

### Environment variables:
- `SYNTENY_INPUT_FILE` - Path to locus_definitions.tsv
- `SYNTENY_OUTPUTS_DIR` - Base output directory

---

## Canonical Pipeline Execution

### Single gene family (e.g., ferritin):
```bash
# Set up environment
export SYNTENY_INPUT_FILE="outputs/ferritin_phase12/locus_definitions.tsv"
export SYNTENY_OUTPUTS_DIR="outputs/ferritin_phase3_real"

# Phase 3: Synteny detection
sbatch run_phase3.slurm

# Phase 4: Target gene detection
sbatch run_phase4.slurm

# Phase 6: Extract sequences (Phase 5 integrated into 4)
sbatch run_phase6_extract.slurm

# Phase 7: SwissProt annotation
sbatch run_phase7_swissprot.slurm

# Phase 8: Generate matrices
python scripts/08a_generate_locus_matrices.py
python scripts/08b_generate_summary_matrices.py
```

### Batch processing (32 gene families):
Use SLURM array jobs to parallelize across gene families (TO BE IMPLEMENTED)

---

## Script Naming Convention

**Current naming (historical):**
- `01_discover_paralogs_gene_level.py` - Phase 1
- `02_synteny_detection.py` - Phase 3 (!)
- `04_blast_targets.py` - Phase 4
- `05_classify_targets.py` - Phase 5 (integrated into 04)
- `06_extract_sequences.py` - Phase 6
- `07_swissprot_annotation.py` - Phase 7
- `08a_generate_locus_matrices.py` - Phase 8a
- `08b_generate_summary_matrices.py` - Phase 8b

**Note:** Numbers don't match phase numbers because:
- Phase 2 was deprecated
- Phase 3 kept old script name `02_`
- Phase 5 is integrated into Phase 4 (both use script `04_`)

**Supporting modules:**
- `config.py` - Pipeline configuration
- `exonerate_extract.py` - Exonerate utilities
- `extract_with_exonerate.py` - Exonerate wrapper with adaptive windowing

**Deprecated/unused:**
- `01_extract_proteins.py` - Old extraction method
- `03_filter_blocks.py` - Integrated into other phases
- `02b_aggregate_synteny_blocks.py` - Helper (may still be used?)
- `02_compare_synteny_gene_level.py` - Unknown purpose

---

## Recent Fixes (Nov 3, 2025)

### Phase 4 - Frame-aware clustering:
- **Bug:** Hits in different reading frames (+1, +2, +3, -1, -2, -3) were being merged
- **Fix:** Modified clustering to group by `(scaffold, strand, frame)` instead of `(scaffold, strand)`
- **Impact:** Prevents incorrect merging of overlapping genes in different frames

### Phase 6 - Gene grouping:
- **Bug:** Each exon was being treated as a separate gene
- **Fix:** Properly parse Exonerate GFF "gene" features and group CDS/exon features by gene_id
- **Impact:** Now correctly reports "1 gene with 3 exons" instead of "3 genes"

### Phase 6 - Adaptive windowing:
- **Bug:** Fixed 1kb flanking window was missing exons >1kb away from HSP
- **Fix:** Implemented progressive window expansion (0-200kb) with coverage-based stopping
- **Impact:** Can now extract genes with introns up to 200kb total span

### Phase 6 - Output buffering:
- **Bug:** Scripts appeared to hang with no output
- **Fix:** Added `flush=True` to all print statements
- **Impact:** Real-time progress visibility

### Phase 7 - Output buffering:
- **Bug:** Same as Phase 6
- **Fix:** Same as Phase 6
- **Impact:** Same as Phase 6

---

## Performance Notes

### Timing (ferritin test - 9 loci, 75 genomes):
- Phase 3: ~30 min
- Phase 4: ~1 hour
- Phase 6: ~42 hours (7 genomes in 4 hours = 6 hours/genome × 75 = 450 hours estimated, but adaptive windowing speeds up as it learns)
- Phase 7: <1 hour (293,478 proteins)

### Bottlenecks:
- Phase 6 is slowest (Exonerate with adaptive windowing)
- Phase 4 improved significantly with combined query optimization

### Parallelization opportunities:
- Phase 3: Per-locus SLURM arrays
- Phase 4: Per-locus SLURM arrays
- Phase 6: Could parallelize per-genome or per-target
- Batch mode: Per-gene-family SLURM arrays

---

## Known Issues

1. **Script numbering mismatch:** Scripts are numbered 01, 02, 04, 05, 06, 07, 08 but phases are 1, 3, 4, 5, 6, 7, 8
2. **Deprecated scripts:** Several old scripts in `scripts/` directory - unclear which are still used
3. **Phase 5 integration:** Phase 5 (classify) is called by Phase 4, not run standalone
4. **Window size limits:** Phase 6 maxes out at 200kb flanking - may miss genes with larger introns

---

## TODO

- [ ] Implement batch processing for 32 gene families
- [ ] Clarify/clean up deprecated scripts
- [ ] Consider renaming scripts to match phase numbers
- [ ] Add parallelization to Phase 6 (per-genome or per-target)
- [ ] Document Phase 1 input/output format more clearly
- [ ] Add validation/QC steps between phases
