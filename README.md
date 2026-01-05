# Helixer Proteome Pipeline

**The canonical synteny-scanning pipeline** for cynipid gene family analysis using Helixer ab initio gene predictions.

This folder is **self-contained** - all scripts needed to run the complete pipeline are here.

---

## Pipeline Flow

```
┌─────────────────────────────────────────────────────────────────────────┐
│                           PHASE 1: LOCUS DISCOVERY                       │
│  Input:  inputs/{family}_test.tsv (LOC IDs from NCBI)                   │
│  Output: outputs/{family}/phase1_v2/                                     │
│          ├── locus_definitions.tsv                                       │
│          ├── {locus_id}_flanking.faa                                     │
│          └── cluster_gene_coordinates.tsv                                │
└─────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────┐
│                       PHASE 1.5: BK-LB UNIFICATION                       │
│  Input:  phase1_v2/locus_definitions.tsv                                │
│  Output: Updated locus_definitions.tsv (merged BK-LB orthologs)         │
└─────────────────────────────────────────────────────────────────────────┘
                                    │
              ┌─────────────────────┴─────────────────────┐
              │                                           │
              ▼                                           ▼
┌──────────────────────────────┐       ┌──────────────────────────────────┐
│   PHASE 2b: FLANKING-FIRST   │       │    PHASE 4: TARGET-FIRST         │
│   SYNTENY BLOCK DETECTION    │       │    DIAMOND SEARCH                │
│                              │       │                                  │
│  Input:                      │       │  Input:                          │
│   - phase1_v2/ flanking FAA  │       │   - phase1_v2/ query proteins    │
│   - Helixer proteomes        │       │   - Helixer proteomes            │
│                              │       │                                  │
│  Output:                     │       │  Output:                         │
│   outputs/{family}/          │       │   outputs/{family}/              │
│   phase2b_helixer/           │       │   phase4_helixer_filtered/       │
│   ├── synteny_blocks.tsv     │       │   └── all_target_loci.tsv        │
│   └── flanking_details.tsv   │       │                                  │
└──────────────────────────────┘       └──────────────────────────────────┘
              │                                           │
              │                                           ▼
              │                        ┌──────────────────────────────────┐
              │                        │  PHASE 4b: EXTRACT FLANKING      │
              │                        │                                  │
              │                        │  Input:                          │
              │                        │   - phase4 targets               │
              │                        │   - Helixer GFF3 files           │
              │                        │                                  │
              │                        │  Output:                         │
              │                        │   outputs/{family}/              │
              │                        │   phase4b_helixer_filtered/      │
              │                        │   └── flanking_genes.tsv         │
              │                        └──────────────────────────────────┘
              │                                           │
              │                                           ▼
              │                        ┌──────────────────────────────────┐
              │                        │  PHASE 5: CLASSIFY TARGETS       │
              │                        │                                  │
              │                        │  Input:                          │
              │                        │   - phase4 targets               │
              │                        │   - phase4b flanking             │
              │                        │   - phase1 BK flanking           │
              │                        │                                  │
              │                        │  Output:                         │
              │                        │   outputs/{family}/              │
              │                        │   phase5_helixer_filtered/       │
              │                        │   ├── syntenic_targets.tsv       │
              │                        │   └── unplaceable_targets.tsv    │
              │                        └──────────────────────────────────┘
              │                                           │
              │                                           ▼
              │                        ┌──────────────────────────────────┐
              │                        │  PHASE 5b: VALIDATE NOVEL LOCI   │
              │                        │                                  │
              │                        │  Input:                          │
              │                        │   - phase5 unplaceables          │
              │                        │   - Helixer proteomes/GFF3       │
              │                        │   - BK/LB references             │
              │                        │                                  │
              │                        │  Output:                         │
              │                        │   outputs/{family}/              │
              │                        │   phase5b_helixer/               │
              │                        │   └── novel_loci_clustered.tsv   │
              │                        └──────────────────────────────────┘
              │                                           │
              └─────────────────────┬─────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────┐
│                      PHASE 6: EXTRACT SEQUENCES                          │
│  Input:  phase4, phase5, phase5b (syntenic + novel loci)                │
│  Output: outputs/{family}/phase6_helixer_filtered/                      │
│          ├── {genome}/{locus}/gene1_protein.fasta                       │
│          ├── NOVEL_{locus}/... (novel loci protein sequences)           │
│          └── all_extracted_genes.faa                                    │
└─────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────┐
│                   PHASE 7: SWISSPROT ANNOTATION                          │
│  Input:  phase4b flanking, phase5, phase5b (novel loci flanking)       │
│  Output: outputs/{family}/phase7_helixer_filtered/                      │
│          ├── flanking_matches_annotated.tsv                             │
│          └── novel_loci_flanking_annotated.tsv                          │
└─────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────┐
│                  PHASE 8a: LOCUS MATRICES                                │
│  Input:  phase1, phase5, phase5b, phase6, phase7                        │
│  Output: outputs/{family}/phase8_helixer_filtered/                      │
│          ├── {locus_id}_genome_swissprot_matrix.tsv (one per locus)     │
│          └── {family}_NOVEL_{locus}_matrix.tsv (novel loci matrices)    │
└─────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────────┐
│                  PHASE 8b: SUMMARY MATRIX                                │
│  Input:  phase1, phase2b, phase5, phase5b, phase6                       │
│  Output: outputs/{family}/phase8_helixer_filtered/                      │
│          └── {family}_summary_matrix.tsv                                │
└─────────────────────────────────────────────────────────────────────────┘
```

**Note**: Phase 5b novel loci are automatically integrated into Phases 6, 7, and 8a when `--phase5b-dir` is provided. Phase 9 is deprecated - novel loci detection now happens in Phase 5b.

---

## Quick Start

```bash
cd /carc/scratch/projects/emartins/2016456/adam/analysis/helixer_pipeline

# Run complete pipeline for all 33 families as array job:
sbatch run_full_pipeline.slurm

# Or run sequentially:
sbatch run_phase1.slurm                    # Phase 1 + 1.5
# Then:
sbatch run_helixer_phases_4_5_array.slurm  # Phases 4, 4b, 2b, 5, 5b
sbatch run_helixer_phases_6_8_array.slurm  # Phases 6, 7, 8a, 8b
```

---

## Script Reference

### Phase 1: Locus Discovery

**Two modes of operation:**

**Mode 1: Legacy (LOC ID lookup)** - Parses BK GFF to find gene coordinates from LOC IDs:
```bash
python 01_phase1.py \
    --loc-ids "LOC117167432,LOC117167433" \
    --gene-family "ferritin_MC102" \
    --output-dir "outputs/ferritin_MC102/phase1" \
    --genome-gff-dir "data/reference"
```

**Mode 2: Coordinates (pre-computed)** - Uses genomic coordinates directly, bypasses GFF lookup:
```bash
python 01_phase1.py \
    --coordinates-file "top50_effector_genomic_coordinates.tsv" \
    --gene-family "top50_effectors" \
    --output-dir "outputs/top50_effectors/phase1" \
    --genome-gff-dir "data/reference" \
    --n-flanking 20
```

The coordinates file should have columns: `subcluster_id`, `gene_id`, `chromosome`, `start`, `end`, `strand`

**Outputs:**
- `locus_definitions.tsv` - Locus metadata (chromosome, coordinates, flanking files)
- `{locus_id}_flanking.faa` - Flanking gene proteins for each locus
- `cluster_gene_coordinates.tsv` - Per-gene coordinates for tandem clusters

### Phase 1.5: BK-LB Unification

```bash
python 01b_merge_bk_lb_loci.py \
    --locus-defs "outputs/{family}/phase1/locus_definitions.tsv" \
    --lb-proteome-db "data/proteomes/GCF_019393585.1_Lboulardi_NCBI.dmnd" \
    --diamond-threads 4
```

**Outputs:**
- Updates `locus_definitions.tsv` in place (merges BK-LB orthologs)
- `bk_lb_locus_mapping.tsv` - Mapping of merged loci

### Phase 2b: Flanking-First Synteny Detection

```bash
python 02b_detect_empty_blocks_helixer.py \
    --family "ferritin_MC102" \
    --phase1-dir "outputs/ferritin_MC102/phase1_v2" \
    --phase4-targets "outputs/ferritin_MC102/phase4_helixer_filtered/all_target_loci.tsv" \
    --helixer-dir "data/ragtag_output" \
    --output-dir "outputs/ferritin_MC102/phase2b_helixer" \
    --threads 8
```

**Outputs:**
- `synteny_blocks.tsv` - Detected synteny blocks (concordant + empty)
- `flanking_details.tsv` - Individual flanking gene hits

### Phase 4: Target Detection (DIAMOND)

```bash
python 04_detect_targets_helixer.py \
    --family "ferritin_MC102" \
    --phase1-dir "outputs/ferritin_MC102/phase1" \
    --helixer-dir "data/ragtag_output" \
    --genome-list "helixer_genomes.tsv" \
    --output-dir "outputs/ferritin_MC102/phase4_helixer_filtered" \
    --threads 8 \
    --min-query-coverage 0.3
```

**Outputs:**
- `all_target_loci.tsv` - All DIAMOND hits passing 30% query coverage filter

### Phase 4b: Extract Flanking Genes

```bash
python 04b_extract_flanking_helixer.py \
    --family "ferritin_MC102" \
    --phase4-output "outputs/ferritin_MC102/phase4_helixer_filtered/all_target_loci.tsv" \
    --helixer-dir "data/ragtag_output" \
    --output-dir "outputs/ferritin_MC102/phase4b_helixer_filtered" \
    --n-flanking 10
```

**Outputs:**
- `flanking_genes.tsv` - Flanking gene IDs for each target
- `flanking_proteins.faa` - Combined flanking protein sequences

### Phase 5: Classify Targets

```bash
python 05_classify_targets_helixer.py \
    --family "ferritin_MC102" \
    --phase4-dir "outputs/ferritin_MC102/phase4_helixer_filtered" \
    --phase4b-dir "outputs/ferritin_MC102/phase4b_helixer_filtered" \
    --phase1-dir "outputs/ferritin_MC102/phase1" \
    --output-dir "outputs/ferritin_MC102/phase5_helixer_filtered" \
    --min-synteny 0.15 \
    --threads 8
```

**Outputs:**
- `syntenic_targets.tsv` - Targets assigned to BK/LB loci
- `unplaceable_targets.tsv` - Targets without synteny match

### Phase 5b: Validate Novel Loci

```bash
python 05b_validate_novel_loci.py \
    --family "ferritin_MC102" \
    --phase5-dir "outputs/ferritin_MC102/phase5_helixer_filtered" \
    --phase1-dir "outputs/ferritin_MC102/phase1" \
    --helixer-dir "data/ragtag_output" \
    --output-dir "outputs/ferritin_MC102/phase5b_helixer" \
    --threads 8
```

**Outputs:**
- `novel_loci_clustered.tsv` - Validated novel loci with cross-genome support

### Phase 6: Extract Sequences

```bash
python 06_extract_sequences_helixer.py \
    --family "ferritin_MC102" \
    --phase1-dir "outputs/ferritin_MC102/phase1" \
    --phase4-dir "outputs/ferritin_MC102/phase4_helixer_filtered" \
    --phase5-dir "outputs/ferritin_MC102/phase5_helixer_filtered" \
    --helixer-dir "data/ragtag_output" \
    --output-dir "outputs/ferritin_MC102/phase6_helixer_filtered" \
    --phase5b-dir "outputs/ferritin_MC102/phase5b_helixer"  # Optional: novel loci
```

**Outputs:**
- `{genome}/{locus}/gene1_protein.fasta` - Per-locus protein sequences
- `all_extracted_genes.faa` - Combined FASTA
- `NOVEL_{locus_name}/` - Novel loci protein sequences (if `--phase5b-dir` provided)

### Phase 7: Annotate Flanking (NR annotations from headers)

Uses NR annotations already present in Helixer proteome FASTA headers instead of running a separate DIAMOND search. This simplified approach parses annotations like "Cytochrome P450 3A31" directly from the FASTA description lines.

```bash
python 07_annotate_flanking.py \
    --family "ferritin_MC102" \
    --phase4b-dir "outputs/ferritin_MC102/phase4b_helixer_filtered" \
    --phase5-dir "outputs/ferritin_MC102/phase5_helixer_filtered" \
    --helixer-dir "data/ragtag_output" \
    --output-dir "outputs/ferritin_MC102/phase7_helixer_filtered" \
    --phase2b-dir "outputs/ferritin_MC102/phase2b_helixer"  # Optional
    --phase5b-dir "outputs/ferritin_MC102/phase5b_helixer"  # Optional: novel loci
```

**Outputs:**
- `flanking_annotations.tsv` - NR annotations for flanking genes + chromosome info
- `flanking_matches_annotated.tsv` - Phase 5 matches with annotations
- `phase2b_flanking_annotated.tsv` - Phase 2b flanking with annotations (if `--phase2b-dir` provided)

### Phase 8a: Locus Matrices

```bash
python 08a_generate_locus_matrices_helixer.py \
    --family "ferritin_MC102" \
    --phase1-dir "outputs/ferritin_MC102/phase1" \
    --phase5-dir "outputs/ferritin_MC102/phase5_helixer_filtered" \
    --phase6-dir "outputs/ferritin_MC102/phase6_helixer_filtered" \
    --phase7-dir "outputs/ferritin_MC102/phase7_helixer_filtered" \
    --species-map "data/gca_to_species.tsv" \
    --output-dir "outputs/ferritin_MC102/phase8_helixer_filtered" \
    --phase5b-dir "outputs/ferritin_MC102/phase5b_helixer"  # Optional: novel loci
```

**Outputs:**
- `{locus_id}_genome_swissprot_matrix.tsv` - One matrix per locus
- `{family}_NOVEL_{locus_name}_matrix.tsv` - Novel loci matrices (if `--phase5b-dir` provided)

### Phase 8b: Summary Matrix

```bash
python 08b_generate_summary_matrices_helixer.py \
    --family "ferritin_MC102" \
    --phase1-dir "outputs/ferritin_MC102/phase1" \
    --phase2b-dir "outputs/ferritin_MC102/phase2b_helixer" \
    --phase5-dir "outputs/ferritin_MC102/phase5_helixer_filtered" \
    --phase5b-dir "outputs/ferritin_MC102/phase5b_helixer" \
    --phase6-dir "outputs/ferritin_MC102/phase6_helixer_filtered" \
    --species-map "data/gca_to_species.tsv" \
    --output-dir "outputs/ferritin_MC102/phase8_helixer_filtered"
```

**Outputs:**
- `{family}_summary_matrix.tsv` - Family-wide summary across all genomes

---

## Required External Inputs

| Path | Description |
|------|-------------|
| `inputs/{family}_test.tsv` | LOC IDs for each family (col 2) |
| `data/ragtag_output/{genome}/` | Helixer proteomes (`*_helixer_proteins.faa`) and GFF3 |
| `data/gca_to_species.tsv` | GCA ID → species name mapping |
| `data/reference/` | BK/LB genome GFF files |
| `data/proteomes/` | BK/LB DIAMOND databases |
| `helixer_genomes.tsv` | List of Helixer-annotated genomes |

**Optional** (used for BK flanking SwissProt annotation):
| `/carc/scratch/projects/emartins/2016456/adam/databases/uniprot_sprot.dmnd` | SwissProt DB |

---

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--min-query-coverage` | 0.30 | Minimum query coverage for DIAMOND hits |
| `--min-synteny` | 0.15 | Minimum synteny score (3/20 flanking match) |
| `--n-flanking` | 10 | Flanking genes per side |
| `--evalue` | 1e-5 | DIAMOND e-value threshold |

---

## SLURM Scripts

| Script | Phases | Time | Memory |
|--------|--------|------|--------|
| `run_phase1.slurm` | 1 + 1.5 | 1 hr | 8G |
| `run_helixer_phases_4_5_array.slurm` | 4 + 4b + 5 | 1 hr | 16G |
| `run_helixer_phases_6_8_array.slurm` | 6 + 7 + 8a + 8b | 30 min | 8G |
| `run_full_pipeline.slurm` | All phases | 3 hr | 16G |

---

## 33 Gene Families

Array indices 0-32:
```
alaserpin_MC35, ApoD_MC13, beta_lectin_MC147, CAP_MC28, carboxypeptidase_MC26,
der_f_21_MC2, endoglucanase_Z_MC6, ester_hydrolase_MC19, extensin-like_MC192,
ferritin_MC102, Glutenin_MC21, LRR_15_MC91, LRR_4_MC10, MC117, MC174, MC18,
MC211, MC258, MC3, natterin-4-like_MC59, neurexin_MC14, OS-D_MC11,
pectin_lyase_PL, peptidoglycan-recognition_protein_SA-like___peptidoglycan-recognition_protein_LC-like_MC125,
rhamnogalacturonate_lyase-like_MC47, ribonuclease_MC33, sprouty_MC112,
thioredoxin-2-like_MC108, TIL_MC9, venom_acid_phosphatase_Acph-1-like_transcript_variant_X1_MC63,
venom_R-like_MC1, venom_serine_carboxypeptidase_MC203, venom_serine_protease_34-like_MC175
```
