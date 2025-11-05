# Pipeline End-to-End Automation - Obstacles & Solutions

**Goal**: Run entire pipeline automatically with just a locus input
**Date**: 2025-11-04
**Status**: Identified 6 major obstacles

---

## OBSTACLE 1: Outdated Master Coordinator Script

**File**: `run_ferritin_full_pipeline.slurm`

**Issues**:
1. References old script names that don't exist:
   - Line 170: `scripts/06_swissprot_annotation.py` (doesn't exist)
   - Line 182: `scripts/07a_generate_locus_matrices.py` (should be `08a`)
   - Line 192: `scripts/07b_generate_summary_matrices.py` (should be `08b`)

2. Hardcoded for ferritin only - not parameterized by locus

3. Missing Phase 3 filter step (only calls 02_synteny_detection.py, not 03_filter_blocks.py)

**Solution**: Rewrite coordinator with correct phase numbering and parameterization

---

## OBSTACLE 2: Phase 7 Signature Changed

**Old signature** (in `run_phase7_8_ferritin.slurm`):
```bash
python3 scripts/07_swissprot_annotation.py \
  --extracted-dir outputs/.../06_extracted_sequences \
  --swissprot-db data/databases/uniprot_sprot.fasta \
  --output-dir ... \
  --threads 8
```

**NEW signature** (current working version):
```bash
python scripts/07_swissprot_annotation.py \
  --synteny-dir outputs/.../02_synteny_blocks \
  --filtered-blocks outputs/.../03_filtered_blocks/synteny_blocks_filtered.tsv \
  --swissprot-db data/reference/swissprot.dmnd \
  --output-dir ... \
  --evalue 1e-5 \
  --threads 16
```

**Why it changed**:
- Phase 7 now annotates ONLY proteins in synteny blocks (not all extracted sequences)
- Uses DIAMOND database (.dmnd) for 100-1000x speedup
- Requires filtered blocks file to know which hits to annotate

**Solution**: Update all coordinator scripts to use new signature

---

## OBSTACLE 3: Missing Phase 3 Aggregation

**Current Phase 3 consists of TWO scripts**:
1. `02_synteny_detection.py` - Runs tBLASTn and clusters hits into blocks
2. `03_filter_blocks.py` - Selects best block per genome per locus

**Additional script**: `02b_aggregate_synteny_blocks.py` (purpose unclear)

**What's missing**: The master coordinator only calls `02_synteny_detection.py` but not `03_filter_blocks.py`

**Impact**: Phase 3 produces raw blocks but not the filtered set needed for downstream phases

**Solution**:
```bash
# Phase 3a: Detect synteny blocks
python scripts/02_synteny_detection.py \
  --locus-defs outputs/.../locus_definitions.tsv \
  --genome-db-dir data/blast_dbs \
  --output-dir outputs/.../02_synteny_blocks \
  ...

# Phase 3b: Filter to best block per genome
python scripts/03_filter_blocks.py \
  --input outputs/.../02_synteny_blocks/all_synteny_blocks.tsv \
  --output-dir outputs/.../03_filtered_blocks
```

---

## OBSTACLE 4: Phase 8a Requires Additional Inputs

**Phase 8a (matrix generation) requires**:
```bash
--locus-defs          # From Phase 2: locus_definitions.tsv
--synteny-dir         # From Phase 3: 02_synteny_blocks/
--blocks              # From Phase 3: 03_filtered_blocks/synteny_blocks_filtered.tsv
--targets             # From Phase 5: 05_classified/all_targets_classified.tsv (?)
--swissprot           # From Phase 7: 07_swissprot_annotations/genome_specific_swissprot_annotations.tsv
--reference-proteins  # FIXED: data/reference/protein.faa (BK proteins)
--species-map         # FIXED: data/gca_to_species.tsv
--output-dir          # New: 08_matrices/
```

**Status**:
- ✓ `reference-proteins` exists: `data/reference/protein.faa`
- ✓ `species-map` exists: `data/gca_to_species.tsv`
- ❓ `--targets` argument name might be wrong (check Phase 5 output filename)

**Solution**: Verify Phase 5 outputs `all_targets_classified.tsv` and use correct path

---

## OBSTACLE 5: Phase Numbering Confusion

**Historical evolution**:
```
OLD numbering:                NEW numbering (current):
Phase 1: Paralog discovery    Phase 1: Paralog discovery
Phase 2: Synteny dedup         Phase 2: Synteny dedup
Phase 3: Synteny detection     Phase 3: Synteny detection + Filter blocks
Phase 4: Target BLAST          Phase 4: Target BLAST
Phase 5: Classification        Phase 5: Classification
Phase 6: SwissProt annotation  Phase 6: Extract sequences (Exonerate)
Phase 7a/7b: Matrices          Phase 7: SwissProt annotation
                               Phase 8a/8b: Matrices
```

**Impact**: Multiple scripts (coordinators, documentation) reference old phase numbers

**Solution**: Standardize on new numbering throughout

---

## OBSTACLE 6: Missing Dependency Chain Documentation

**Current pipeline phases**:
```
Phase 1: 01_discover_paralogs_gene_level.py
  └─> outputs/{output_dir}/unique_loci.tsv

Phase 2: 02_compare_synteny_gene_level.py
  └─> outputs/{output_dir}/locus_definitions.tsv
  └─> outputs/{output_dir}/synteny_groups.json

Phase 3a: 02_synteny_detection.py
  ├─ Inputs: locus_definitions.tsv, data/blast_dbs/*.nhr
  └─> outputs/{output_dir}/02_synteny_blocks/{locus_id}/flanking_blast_all.tsv
  └─> outputs/{output_dir}/02_synteny_blocks/{locus_id}/hit_sequences/*.fasta

Phase 3b: 03_filter_blocks.py
  ├─ Inputs: 02_synteny_blocks/all_synteny_blocks.tsv
  └─> outputs/{output_dir}/03_filtered_blocks/synteny_blocks_filtered.tsv

Phase 4: 04_blast_targets.py
  ├─ Inputs: 03_filtered_blocks/synteny_blocks_filtered.tsv, data/blast_dbs/*
  └─> outputs/{output_dir}/04_target_genes/all_target_loci.tsv

Phase 5: 05_classify_targets.py
  ├─ Inputs: 04_target_genes/all_target_loci.tsv, 03_filtered_blocks/synteny_blocks_filtered.tsv
  └─> outputs/{output_dir}/05_classified/syntenic_targets.tsv
  └─> outputs/{output_dir}/05_classified/unplaceable_targets.tsv

Phase 6: 06_extract_sequences.py
  ├─ Inputs: 05_classified/syntenic_targets.tsv, 05_classified/unplaceable_targets.tsv
  ├─ Requires: Exonerate module, genome FASTA files
  └─> outputs/{output_dir}/06_extracted_sequences/{genome}/{locus}/*.faa

Phase 7: 07_swissprot_annotation.py
  ├─ Inputs: 02_synteny_blocks/, 03_filtered_blocks/synteny_blocks_filtered.tsv
  ├─ Requires: data/reference/swissprot.dmnd
  └─> outputs/{output_dir}/07_swissprot_annotations/genome_specific_swissprot_annotations.tsv

Phase 8a: 08a_generate_locus_matrices.py
  ├─ Inputs: All above + data/reference/protein.faa + data/gca_to_species.tsv
  └─> outputs/{output_dir}/08_matrices/{locus}_genome_swissprot_matrix.tsv

Phase 8b: 08b_generate_summary_matrices.py
  ├─ Inputs: 08_matrices/{locus}_*.tsv
  └─> outputs/{output_dir}/08_matrices/summary_matrix.tsv
```

---

## SOLUTION: Create New Master Coordinator

**Requirements**:
1. Parameterized by locus/gene family name
2. Correct script names and arguments
3. Proper error handling and logging
4. Dependency chain validation
5. Output verification between phases

**Template structure**:
```bash
#!/bin/bash
# Master coordinator: run_synteny_pipeline.slurm
# Usage: sbatch run_synteny_pipeline.slurm <gene_family> <output_base>

GENE_FAMILY=$1
OUTPUT_BASE=$2

# Phase 1-2: Discover and deduplicate
# Phase 3a: Synteny detection
# Phase 3b: Filter blocks
# Phase 4: Target BLAST
# Phase 5: Classification
# Phase 6: Extract sequences
# Phase 7: SwissProt annotation (NEW signature)
# Phase 8a: Generate matrices
# Phase 8b: Generate summaries
```

---

## IMMEDIATE ACTION ITEMS

1. **Update Phase 7 in all coordinator scripts** to use new signature
2. **Add Phase 3b (filter_blocks.py)** to pipeline
3. **Create new master coordinator** with correct parameterization
4. **Test end-to-end** on ferritin data
5. **Document** all required input files and their locations

---

## FILES TO UPDATE

- `run_ferritin_full_pipeline.slurm` → Rewrite completely
- `run_phase7_8_ferritin.slurm` → Update Phase 7 args
- Create: `run_synteny_pipeline_master.slurm` (new parameterized version)
- Update: Pipeline documentation with correct phase numbers
