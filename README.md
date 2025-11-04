# Hybrid Synteny Scanning Workflow

**Status**: In development
**Date**: 2025-11-03

## Overview

This hybrid workflow combines the best approaches from multiple methods:

1. **Landmark discovery** (Phases 1-2) from gene_level_workflow - uses NCBI proteomes for BK/LB
2. **Genome-level scanning** (Phase 3+) from original working pipeline - uses 75 RagTag genomes
3. **Exonerate extraction** - replaces BLAST HSPs with proper gene structures

## Architecture

### Phase 1-2: Landmark Discovery
- `01_discover_paralogs_gene_level.py` - Find paralogs in BK, LB (NCBI proteomes)
- `02_compare_synteny_gene_level.py` - Deduplicate loci by synteny comparison

### Phase 3+: Genome-Level Synteny Scanning
- `01_extract_proteins.py` - Extract flanking proteins for each locus
- `02_synteny_detection.py` - **tBLASTn + Exonerate extraction** (KEY CHANGE)
- `03_filter_blocks.py` - Quality filter synteny blocks
- `04_blast_targets.py` - Search for target genes in blocks
- `05_classify_targets.py` - Classify gene functional status
- `06_swissprot_annotation.py` - Functional annotation
- `07a_generate_locus_matrices.py` - Per-locus presence/absence matrices
- `07b_generate_summary_matrices.py` - Gene family summary matrices

### Support Modules
- `exonerate_extract.py` - Exonerate extraction functions
- `config.py` - Configuration (USE_PROTEOME_MODE = False)

## Key Innovation: Exonerate Instead of BLAST HSPs

**Old method** (from failed proteome approach):
```
tBLASTn → extract HSP sequences → use directly
```

**New method** (this hybrid approach):
```
tBLASTn → get coordinates → Exonerate extraction → full gene structures with introns/exons
```

**Benefits**:
- Full CDS sequences (not fragments)
- Proper intron/exon structure
- Detects pseudogenes and frameshifts
- Better quality for phylogenetics
- Expected 70-90% detection rate (vs 30% with BRAKER3 proteomes)

## Development Tracking

Use Beads to track progress:
```bash
cd hybrid_workflow
bd list                    # Show all issues
bd show hybrid_workflow-N  # Show issue details
bd start hybrid_workflow-N # Mark issue as in progress
bd close hybrid_workflow-N # Mark issue as completed
```

## Dependencies

**Environment**:
```bash
micromamba activate trinity_new_env  # Python, BLAST
module load exonerate/2.4.0-yl7q     # Exonerate
```

**Data Requirements**:
- `data/reference/genomic.gff` - BK annotations
- `data/reference/protein.faa` - BK proteins
- `data/ragtag_dbs/*.n{hr,in,sq}` - 75 genome BLAST databases
- `data/gca_to_species.tsv` - Species mapping
- `data/databases/swissprot.dmnd` - SwissProt for annotation

## Testing Plan

1. Test on ferritin_MC102 (simple 1-locus case)
2. Verify Phases 1-2 find paralogs correctly
3. Verify Phase 3+ with Exonerate extracts gene structures
4. Compare results with failed proteome approach
5. Scale to all 32 gene families via SLURM array

## Current Status

**Completed**:
- ✓ Directory structure created
- ✓ All scripts copied into place
- ✓ Beads tracking initialized
- ✓ Development roadmap created

**In Progress**:
- See `bd list` for current task status

---

**Pipeline Version**: 1.0 (Hybrid: Landmark + Genome + Exonerate)
**Last Updated**: 2025-11-03
