# Helixer Proteome Pipeline

Synteny-based target detection using Helixer ab initio gene predictions instead of tblastn/exonerate.

## Pipeline Overview

```
Phase 4:   DIAMOND blastp against Helixer proteomes (30% query coverage filter)
Phase 4b:  Extract flanking genes from Helixer GFF3
Phase 5:   Classify targets by flanking gene synteny to BK/LB loci
Phase 6:   Extract and organize target proteins (per-genome/per-locus)
Phase 8a:  Generate per-locus matrices (synteny + targets)
Phase 8b:  Generate family summary matrix
```

## Output Format

The pipeline produces outputs **matching the original tblastn pipeline format**:

### Phase 6 (Sequence Extraction)
```
outputs/{family}/phase6_helixer_filtered/
├── {genome}/
│   ├── {assigned_locus}/
│   │   └── {locus}_gene1_protein.fasta
│   └── UNPLACED/
│       └── {scaffold}_{start}_{end}/
│           └── {tag}_gene1_protein.fasta
└── all_extracted_genes.faa
```

### Phase 8 (Summary Matrices)
```
outputs/{family}/phase8_helixer_filtered/
├── {locus_id}_genome_swissprot_matrix.tsv  (24 files, one per locus)
└── {family}_summary_matrix.tsv             (family-wide summary)
```

## Usage

From `hybrid_workflow/` directory:

```bash
# Run Phases 4-5 (target detection + classification)
sbatch helixer_pipeline/run_helixer_phases_4_5_array.slurm

# Run Phases 6-8 (extraction + summaries)
sbatch helixer_pipeline/run_helixer_phases_6_8_array.slurm
```

## Required Inputs

- `outputs/{family}/phase1_v2/` - BK/LB locus definitions and query proteins
- `data/ragtag_output/{genome}/` - Helixer proteomes and GFF3 files
- `data/gca_to_species.tsv` - Species mapping with phylogenetic order
- `helixer_genomes.tsv` - List of genomes with Helixer annotations

## Key Parameters

- **Query coverage filter**: 30% minimum (removes partial domain matches)
- **Flanking genes**: 10 on each side
- **Synteny threshold**: 15% flanking gene matches

## Scripts

| Script | Phase | Description |
|--------|-------|-------------|
| `04_detect_targets_helixer.py` | 4 | DIAMOND blastp + deduplication |
| `04b_extract_flanking_helixer.py` | 4b | Parse GFF3 for flanking genes |
| `05_classify_targets_helixer.py` | 5 | Synteny scoring via BK proteome |
| `06_extract_sequences_helixer.py` | 6 | Per-genome/per-locus protein extraction |
| `08a_generate_locus_matrices_helixer.py` | 8a | Per-locus matrices with U/D flanking |
| `08b_generate_summary_matrices_helixer.py` | 8b | Family summary with synteny % |

## Comparison with tblastn Pipeline

| Aspect | tblastn Pipeline | Helixer Pipeline |
|--------|------------------|------------------|
| Gene finding | tblastn + Exonerate | Helixer ab initio |
| Target extraction | Exonerate models | Helixer proteins |
| Output format | **Identical** | **Identical** |
| Flanking genes | Phase 2 tblastn | Phase 4b GFF3 parse |
