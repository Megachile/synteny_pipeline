# Overnight Pipeline Test - 2025-11-04

## Goal
Test the complete end-to-end pipeline on ferritin using the new master coordinator.

## Command to Run

```bash
# Wait for current Phase 7 (job 98091) to complete first
squeue -u $USER

# Then launch full pipeline test
sbatch run_synteny_pipeline_master.slurm ferritin ferritin_full_test_20251104
```

## What Will Happen

The pipeline will run **all 8 phases automatically** in sequence:
1. Phase 1-2: Discover paralogs (5-10 min)
2. Phase 3: Genome-wide synteny detection (1-2 hrs)
3. Phase 4: Target BLAST (20-30 min)
4. Phase 5: Classification (5 min)
5. Phase 6: Exonerate extraction (2-4 hrs)
6. Phase 7: SwissProt annotation (10-15 min)
7. Phase 8: Matrix generation (5 min)

**Expected total runtime**: 4-8 hours

## Output Location

```
outputs/ferritin_full_test_20251104/
├── phase12_landmark/         # Landmark discovery
├── 02_synteny_blocks/        # tBLASTn results
├── 03_filtered_blocks/       # Best blocks
├── 04_target_genes/          # Target genes found
├── 05_classified/            # Syntenic/tandem classification
├── 06_extracted_sequences/   # Full gene structures
├── 07_swissprot_annotations/ # Protein annotations
└── 08_matrices/              # Final matrices
```

## Monitoring

**Check job status:**
```bash
squeue -u $USER
```

**Watch progress live:**
```bash
tail -f logs/pipeline_ferritin_progress.log
```

**Check logs:**
```bash
ls -lh logs/phase*_ferritin.log
ls -lh logs/pipeline_synteny_pipeline_*.out
```

## What to Check in the Morning

1. **Job completed successfully?**
   ```bash
   tail -50 logs/pipeline_synteny_pipeline_*.out
   ```
   Should see: "PIPELINE COMPLETE - SUCCESS!"

2. **How many genomes detected?**
   ```bash
   wc -l outputs/ferritin_full_test_20251104/03_filtered_blocks/synteny_blocks_filtered.tsv
   ```
   Expected: ~50-70 blocks (70-90% of 75 genomes)

3. **Final matrices created?**
   ```bash
   ls -lh outputs/ferritin_full_test_20251104/08_matrices/
   ```
   Should see: `*_genome_swissprot_matrix.tsv` and summary files

4. **Check a specific locus matrix:**
   ```bash
   head -20 outputs/ferritin_full_test_20251104/08_matrices/BK_chr2_a_genome_swissprot_matrix.tsv | column -t
   ```

5. **Any errors?**
   ```bash
   cat logs/pipeline_synteny_pipeline_*.err
   ```

## Success Criteria

- ✓ All 8 phases complete
- ✓ ≥50 synteny blocks detected (≥70% genome coverage)
- ✓ Locus matrices created with SwissProt annotations
- ✓ Gene order matches landmark genome (BK)

## If Something Fails

1. Check which phase failed: `tail logs/pipeline_ferritin_progress.log`
2. Check phase-specific log: `cat logs/phase{N}_ferritin.log`
3. Check error log: `cat logs/pipeline_synteny_pipeline_*.err`

The pipeline stops at the first failure, so you can fix and continue from there.

## Alternative: Test Individual Phases First

If you want to be cautious, test phases individually:

```bash
# Just Phase 1-3 (synteny detection)
sbatch run_phase3_direct.slurm

# Then Phase 4-5 (targets + classification)
sbatch run_phase4.slurm

# Then Phase 6-8 (extraction + annotation + matrices)
# ... etc
```

But the master coordinator should handle all this automatically!

---

## Notes

- Current Phase 7 job (98091) is running - wait for it to finish first
- The new master coordinator fixes all the issues we found:
  - ✓ Correct Phase 7 signature (uses DIAMOND, filters to synteny blocks)
  - ✓ Includes Phase 3b (filter_blocks.py)
  - ✓ Parameterized by gene family
  - ✓ Proper error handling
  - ✓ Saves bk_protein_id for Phase 8 matching
