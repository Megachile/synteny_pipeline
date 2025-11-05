# Synteny Pipeline Master Coordinator - Usage Guide

## Quick Start

### Run entire pipeline with one command:
```bash
sbatch run_synteny_pipeline_master.slurm ferritin ferritin_test_run1
```

Arguments:
- `ferritin` - Gene family name (looks for `inputs/ferritin*.tsv`)
- `ferritin_test_run1` - Output directory name (creates `outputs/ferritin_test_run1/`)

---

## What It Does

The master coordinator runs all 8 phases automatically:

**Phase 1**: Discover paralogs in landmark genomes (BK, LB, TR, DR)
- Input: `inputs/ferritin*.tsv`
- Output: `outputs/{run}/phase12_landmark/unique_loci.tsv`

**Phase 2**: Deduplicate by synteny (group loci with same flanking genes)
- Output: `outputs/{run}/phase12_landmark/locus_definitions.tsv`

**Phase 3a**: Genome-wide synteny detection (tBLASTn + clustering)
- Output: `outputs/{run}/02_synteny_blocks/`

**Phase 3b**: Filter to best block per genome
- Output: `outputs/{run}/03_filtered_blocks/synteny_blocks_filtered.tsv`

**Phase 4**: BLAST for target genes
- Output: `outputs/{run}/04_target_genes/all_target_loci.tsv`

**Phase 5**: Classify as syntenic/tandem/unplaceable
- Output: `outputs/{run}/05_classified/syntenic_targets.tsv`

**Phase 6**: Extract full gene structures with Exonerate
- Output: `outputs/{run}/06_extracted_sequences/`

**Phase 7**: SwissProt annotation of flanking proteins (ONLY in synteny blocks)
- Output: `outputs/{run}/07_swissprot_annotations/genome_specific_swissprot_annotations.tsv`

**Phase 8a**: Generate locus-specific matrices
- Output: `outputs/{run}/08_matrices/{locus}_genome_swissprot_matrix.tsv`

**Phase 8b**: Generate summary matrices
- Output: `outputs/{run}/08_matrices/{gene}_summary.tsv`

---

## Required Files (Already Present)

✓ Inputs:
- `inputs/{gene_family}*.tsv` - Landmark genome loci

✓ Databases:
- `data/blast_dbs/*.nhr` - Genome BLAST databases
- `data/reference/swissprot.dmnd` - SwissProt DIAMOND database
- `data/reference/protein.faa` - BK reference proteins
- `data/gca_to_species.tsv` - Species mapping

✓ Genome sequences:
- `data/ragtag_output/` - Assembled genome FASTA files

---

## Monitoring Progress

**Check job status:**
```bash
squeue -u $USER
```

**Watch progress log (live):**
```bash
tail -f logs/pipeline_ferritin_progress.log
```

**Check phase-specific logs:**
```bash
ls logs/phase*_ferritin.log
```

**Check overall output:**
```bash
tail -50 logs/pipeline_synteny_pipeline_<jobid>.out
```

---

## Expected Runtime

- **Phase 1-2**: 5-10 minutes (landmark discovery)
- **Phase 3**: 1-2 hours (tBLASTn genome-wide)
- **Phase 4**: 20-30 minutes (target BLAST)
- **Phase 5**: 5 minutes (classification)
- **Phase 6**: 2-4 hours (Exonerate extraction)
- **Phase 7**: 10-15 minutes (DIAMOND SwissProt)
- **Phase 8**: 5 minutes (matrix generation)

**Total**: ~4-8 hours for ~75 genomes

---

## Output Structure

```
outputs/ferritin_test_run1/
├── phase12_landmark/
│   ├── unique_loci.tsv
│   ├── locus_definitions.tsv
│   └── synteny_groups.json
├── 02_synteny_blocks/
│   ├── {locus_id}/
│   │   ├── flanking_blast_all.tsv
│   │   └── hit_sequences/
│   └── all_synteny_blocks.tsv
├── 03_filtered_blocks/
│   └── synteny_blocks_filtered.tsv
├── 04_target_genes/
│   └── all_target_loci.tsv
├── 05_classified/
│   ├── syntenic_targets.tsv
│   └── unplaceable_targets.tsv
├── 06_extracted_sequences/
│   └── {genome}/
│       └── {locus}/*.faa
├── 07_swissprot_annotations/
│   └── genome_specific_swissprot_annotations.tsv
├── 08_matrices/
│   ├── {locus}_genome_swissprot_matrix.tsv
│   └── {gene}_summary.tsv
└── locus_definitions.tsv (copy)
```

---

## Troubleshooting

**If a phase fails:**
1. Check the error log: `logs/pipeline_synteny_pipeline_<jobid>.err`
2. Check phase-specific log: `logs/phase{N}_{gene}.log`
3. The pipeline will stop at the failed phase
4. Fix the issue and either:
   - Restart from beginning: `sbatch run_synteny_pipeline_master.slurm ...`
   - OR manually run remaining phases

**Common issues:**
- **"Input file not found"**: Check `inputs/{gene_family}*.tsv` exists
- **"Module not found"**: Check `module load exonerate/2.4.0-yl7q` works
- **"BLAST database not found"**: Check `data/blast_dbs/` has `.nhr` files
- **Phase 6 slow**: Exonerate is CPU-intensive, use more cores if available

---

## Testing Individual Phases

You can also run phases individually for testing:

```bash
# Phase 3a only
sbatch run_phase3_direct.slurm

# Phase 7 only
sbatch outputs/ferritin_phase3_real/run_phase7_swissprot.slurm

# Phase 8 only (requires Phase 7 complete)
python scripts/08a_generate_locus_matrices.py --locus-defs ... --output-dir ...
```

---

## Next Steps After Completion

1. **Check detection rate**:
   ```bash
   wc -l outputs/{run}/03_filtered_blocks/synteny_blocks_filtered.tsv
   ```
   Expected: ~50-70 blocks for ferritin (70-90% of 75 genomes)

2. **View locus matrices**:
   ```bash
   head outputs/{run}/08_matrices/*_genome_swissprot_matrix.tsv
   ```

3. **Analyze summary matrix**:
   ```bash
   column -t outputs/{run}/08_matrices/*_summary.tsv | less -S
   ```

4. **Check for specific genome**:
   ```bash
   grep "GCA_010883055.1" outputs/{run}/08_matrices/*_summary.tsv
   ```

---

## Advanced: Running Multiple Gene Families

```bash
# Ferritin
sbatch run_synteny_pipeline_master.slurm ferritin ferritin_run1

# Another gene family (if prepared)
sbatch run_synteny_pipeline_master.slurm geneX geneX_run1
```

Each runs independently in its own output directory.
