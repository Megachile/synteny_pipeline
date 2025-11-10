# Phase 2 Batch Organization - CRITICAL README

## ⚠️ CRITICAL: Phase 2 MUST be run at LOCUS level, NOT family level

### Why Locus-Level Processing is Essential

**Family-level processing WILL FAIL because:**
- Large families have 38+ loci × 75 genomes = 2,850+ BLAST searches
- Jobs timeout even with 24-hour limits
- Cannot get compute allocation (stuck in Priority queue)
- Corrupt XML files when jobs timeout mid-write

**Locus-level processing succeeds because:**
- Each job handles 1 locus × 75 genomes = ~75 BLAST searches
- Completes within 3-hour window
- Gets immediate compute allocation
- Clean failure boundaries (one locus at a time)

## 📋 Batch Organization Strategy

### Step 1: Create Locus List
```bash
# Generate list of all loci from Phase 1 outputs
python create_locus_list.py
# Creates: all_loci.txt with 293 loci
```

Format of all_loci.txt:
```
family_name<TAB>locus_id<TAB>path_to_locus_definitions.tsv
```

### Step 2: Submit SLURM Array Job
```bash
# Submit array job with one task per locus
sbatch run_phase2_locus_array.slurm
```

Key SLURM parameters:
- **Array size**: 1-293 (one per locus)
- **CPUs**: 16 (for tBLASTn threading)
- **Memory**: 20GB (sufficient for one locus)
- **Time**: 3 hours (most complete in 1-2 hours)

### Step 3: Monitor Progress
```bash
# Check running jobs
squeue -u $USER | grep phase2_locus

# Count completed loci
ls outputs/*/phase2_synteny/*/blast_xml/*.xml | wc -l

# Check for incomplete XML files
find outputs/*/phase2_synteny -name "*.xml" -exec sh -c \
  'tail -1 "$1" | grep -q "</BlastOutput>" || echo "$1: Incomplete"' sh {} \;
```

## 🔄 Resume Strategy

The pipeline automatically handles resumption:

1. **XML Validation**: Checks for `</BlastOutput>` closing tag
2. **Auto-cleanup**: Deletes incomplete XML files
3. **Skip completed**: Only runs missing BLAST searches

```python
# In synteny_detection.py:
if output_xml.exists() and is_valid_blast_xml(output_xml):
    print(f"Using existing BLAST results")
else:
    if output_xml.exists():
        print(f"Removing incomplete BLAST results...")
        output_xml.unlink()
    # Run tBLASTn...
```

## 📊 Resource Calculation

### Per Locus:
- **Flanking proteins**: 50-80 (typically ~60)
- **Target genomes**: 75
- **BLAST searches**: ~60 × 75 = 4,500 per locus
- **Time per BLAST**: 1-3 seconds
- **Total time**: 1-2 hours per locus

### Per Family (DO NOT USE):
- **Loci**: 1-38 (average ~9)
- **BLAST searches**: up to 38 × 75 = 2,850 per family
- **Time**: 24+ hours for large families
- **Result**: Timeouts and corrupt files

## 🚨 Common Pitfalls to Avoid

1. **DO NOT** run Phase 2 at family level
2. **DO NOT** use less than 3-hour time limit
3. **DO NOT** trust XML files without validation
4. **DO NOT** skip the cleanup step after timeouts

## 📝 Example SLURM Script (Correct Way)

```bash
#!/bin/bash
#SBATCH --job-name=phase2_locus
#SBATCH --array=1-293
#SBATCH --cpus-per-task=16
#SBATCH --mem=20G
#SBATCH --time=03:00:00

# Get locus info for this array task
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" all_loci.txt)
family=$(echo "$line" | cut -f1)
locus_id=$(echo "$line" | cut -f2)
locus_defs=$(echo "$line" | cut -f3)

# Run for SPECIFIC LOCUS (not whole family!)
python scripts/synteny_detection.py \
    --locus-defs "$locus_defs" \
    --genome-db-dir data/ragtag_dbs \
    --output-dir "outputs/${family}/phase2_synteny" \
    --locus "$locus_id" \  # CRITICAL: Specify individual locus!
    --evalue 1e-5 \
    --max-targets 50 \
    --threads 16
```

## 📈 Expected Output Structure

```
outputs/
├── family_name/
│   ├── phase1/
│   │   └── locus_definitions.tsv
│   └── phase2_synteny/
│       ├── BK_chr1_a/
│       │   ├── blast_xml/
│       │   │   ├── genome1.xml
│       │   │   ├── genome2.xml
│       │   │   └── ... (75 XML files)
│       │   ├── flanking_blast_all.tsv
│       │   └── hit_sequences/
│       ├── BK_chr2_b/
│       │   └── ... (same structure)
│       └── combined_synteny_blocks.tsv
```

## ✅ Verification Checklist

- [ ] Created all_loci.txt with 293 entries
- [ ] SLURM script uses --locus parameter
- [ ] Array size matches locus count (1-293)
- [ ] Time limit ≥ 3 hours
- [ ] XML validation enabled in synteny_detection.py
- [ ] Monitoring for incomplete XML files
- [ ] Resume capability tested

## 🔥 Emergency Recovery

If jobs fail with corrupt XML files:
```bash
# 1. Find and remove incomplete XML files
./cleanup_incomplete_xml.sh

# 2. Resubmit failed tasks only
failed_tasks=$(sacct -j $JOBID --format=JobID,State | \
  grep TIMEOUT | cut -d_ -f2)
sbatch --array=$failed_tasks run_phase2_locus_array.slurm
```

---
**Last Updated**: November 2025
**Critical Change**: Switched from family-level to locus-level processing after 9,494 XML files corrupted in family-level timeout