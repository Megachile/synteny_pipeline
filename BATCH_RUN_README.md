# Batch Processing: 32 Gene Families

## Overview

Process all 32 gene families through Phases 1-5 in parallel using SLURM array jobs.

**Phases included:**
- Phase 1: Paralog discovery with tandem LOC deduplication
- Phase 3: Synteny detection via tBLASTn
- Phase 4: Target gene detection with frame-aware clustering
- Phase 5: Target classification (syntenic vs unplaceable)

**Phases excluded:**
- Phase 6: Exonerate extraction (just fixed, will run separately)
- Phase 7: SwissProt annotation (depends on Phase 6)
- Phase 8: Matrix generation (depends on Phase 6)

## Input

`hpc_input_loci.tsv` - 32 gene families with BK LOC IDs
- Families range from 1 LOC (ferritin_MC102) to 28 LOCs (MC117)
- LOCs are comma-separated per family
- Phase 1 detects and groups tandem duplications automatically

## Output Structure

```
outputs/batch_run/
├── MC117/
│   ├── 01_paralogs/
│   │   ├── unique_loci.tsv
│   │   ├── input_tandem_clusters.json
│   │   └── *_flanking.faa
│   ├── 02_synteny_blocks/
│   ├── 03_filtered_blocks/
│   ├── 04_target_genes/
│   └── 05_classified/
├── der_f_21_MC2/
│   └── ...
├── ApoD_MC13/
│   └── ...
... (32 families total)
```

## Running

### Submit array job:
```bash
sbatch run_batch_array.slurm
```

This launches 32 parallel tasks (one per gene family).

### Monitor progress:
```bash
# Check running jobs
squeue -u $USER

# Check specific task output
tail -f logs/batch_JOBID_TASKID.out

# Count completed tasks
ls outputs/batch_run/*/05_classified/ 2>/dev/null | wc -l
```

## Resource Allocation

- **Time**: 8 hours per family
- **Memory**: 16GB per task
- **CPUs**: 4 per task
- **Total resources**: 32 tasks × 16GB = 512GB concurrent (if all run at once)

## Expected Completion

- **Best case**: All 32 families complete in ~8 hours (if sufficient nodes available)
- **Typical**: Staggered completion as resources become available
- **Check completion**: All families should have `05_classified/` directories

## Troubleshooting

**Check failed tasks:**
```bash
grep -l "ERROR\|failed" logs/batch_*.err
```

**Rerun individual family:**
```bash
bash process_gene_family.sh "FAMILY_NAME" "LOC1,LOC2,LOC3"
```

**Common issues:**
- Missing BLAST databases (check data/proteomes/*.dmnd)
- Environment not activated (trinity_new_env required)
- Insufficient memory (increase --mem if tasks fail)
