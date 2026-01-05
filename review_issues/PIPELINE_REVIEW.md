# Helixer Pipeline Review

**Date:** 2026-01-04
**Reviewer:** Claude Code
**Scope:** Code quality, consistency, organization, SLURM optimization

---

## 1. TEMP CODE / SHORTCUTS

### 1.1 SLURM Script Path Mismatch
**File:** `run_full_pipeline.slurm:30`
```bash
cd /carc/scratch/projects/emartins/2016456/adam/synteny_scanning/hpc_deployment_package/hybrid_workflow
```
**Issue:** The SLURM script changes to a completely different directory (`hpc_deployment_package/hybrid_workflow`) instead of the pipeline's actual location (`analysis/helixer_pipeline`). This suggests the pipeline was copied/moved and paths weren't updated.

**Impact:** Script won't find the Python scripts if run from current location.

---

### 1.2 Phase 7 Script Name Mismatch
**File:** `run_full_pipeline.slurm:251`
```bash
python helixer_pipeline/07_annotate_flanking_swissprot.py \
```
**Issue:** SLURM calls `07_annotate_flanking_swissprot.py` but the active script is `07_annotate_flanking.py`. There's a deprecated version in `deprecated/07_annotate_flanking_swissprot.py`.

**Impact:** SLURM job will fail with "file not found".

---

### 1.3 Phase 9 Input File Mismatch
**File:** `09_detect_novel_loci_helixer.py:47-55`
```python
def load_orphan_targets(phase5b_dir: Path) -> pd.DataFrame:
    refined_file = phase5b_dir / "unplaceable_refined.tsv"
    ...
    orphans = df[df['refined_class'] == 'orphan'].copy()
```
**Issue:** Phase 9 expects `unplaceable_refined.tsv` with a `refined_class` column containing "orphan" values. But Phase 5b (`05b_validate_novel_loci.py`) produces:
- `novel_candidates.tsv`
- `novel_locus_validation.tsv`
- `novel_loci_clustered.tsv`

There's no `unplaceable_refined.tsv` output and no `refined_class` column.

**Impact:** Phase 9 will always return "No orphans to analyze" and do nothing.

---

### 1.4 Symlink Dependency
**File:** `data_paths.py:30-31`
```python
DATA_DIR = PIPELINE_DIR / "data"
```
**Issue:** Pipeline expects a `data` symlink in the pipeline directory pointing to `pancynipoid_helixer`. This symlink may not exist.

**Impact:** Scripts will fall back to hardcoded path, but this is fragile.

---

### 1.5 Test Run Debris
**Location:** `test_run/endoglucanase_Z_MC6/` and `test_run/endo_test/`

**Issue:** Test run folders exist with partial outputs. These should either be cleaned up or documented as examples.

---

## 2. APPARENT MISTAKES

### 2.1 SLURM Argument Mismatch for Phase 2b
**File:** `run_full_pipeline.slurm:179-185`
```bash
python helixer_pipeline/02b_detect_empty_blocks_helixer.py \
    --family "$FAMILY" \
    --phase1-dir "$PHASE1_DIR" \
    --phase4-dir "$PHASE4_DIR" \   # WRONG!
```
**Script expects:** `--phase4-targets` (Path to TSV file, not directory)

**Actual script signature:** `02b_detect_empty_blocks_helixer.py:487-489`
```python
parser.add_argument(
    "--phase4-targets", type=Path, default=None,
    help="Phase 4 targets file (to exclude regions with hits)"
)
```

**Impact:** Phase 2b won't receive target information for concordance checking.

---

### 2.2 Missing Phase 5b Arguments in SLURM
**File:** `run_full_pipeline.slurm:214-221`
```bash
python helixer_pipeline/05b_validate_novel_loci.py \
    --family "$FAMILY" \
    --phase5-dir "$PHASE5_DIR" \
    --phase4b-dir "$PHASE4B_DIR" \    # Script expects --phase1-dir, not --phase4b-dir
    --helixer-dir "$HELIXER_DIR" \
    --bk-proteome "$BK_PROTEOME_DB" \  # Script expects --phase1-dir for BK reference
```

**Actual script signature:** `05b_validate_novel_loci.py:386-392`
```python
parser.add_argument("--phase1-dir", type=Path, required=True,
    help="Phase 1 directory (for BK/LB reference proteomes)")
```

**Impact:** Phase 5b won't know where reference proteomes are.

---

### 2.3 Phase 9 Wrong Argument in SLURM
**File:** `run_full_pipeline.slurm:329-334`
```bash
python helixer_pipeline/09_detect_novel_loci_helixer.py \
    ...
    --helixer-dir "$HELIXER_DIR" \  # Script doesn't take this argument!
```

**Actual script signature:** `09_detect_novel_loci_helixer.py:183-192`
```python
parser.add_argument("--phase5b-dir", type=Path, required=True, ...)
parser.add_argument("--phase4b-dir", type=Path, required=True, ...)
parser.add_argument("--phase5-dir", type=Path, required=True, ...)  # For diamond_vs_bk.tsv
parser.add_argument("--output-dir", type=Path, required=True, ...)
# NO --helixer-dir argument exists!
```

**Impact:** Script will fail with "unrecognized arguments".

---

### 2.4 Documentation vs Reality: Phase 7
**File:** `README.md` and `IMPLEMENTATION_PLAN.md`

Documentation claims Phase 7 runs DIAMOND against SwissProt. Actual `07_annotate_flanking.py` just parses NR annotations from FASTA headers (no DIAMOND search).

```python
# From 07_annotate_flanking.py docstring:
"""
SIMPLIFIED VERSION: Uses annotations already present in Helixer proteome FASTA headers
(from DIAMOND NR search) instead of running a separate SwissProt search.
"""
```

**Impact:** Documentation is misleading about what Phase 7 does.

---

### 2.5 Broken Phase 5b → Phase 9 Pipeline
**Chain:** Phase 5b outputs → Phase 9 inputs

Phase 5b outputs:
- `novel_candidates.tsv` (columns: target_gene_id, genome, scaffold, start, end, best_bitscore, ...)
- `novel_locus_validation.tsv`
- `novel_loci_clustered.tsv`

Phase 9 expects:
- `unplaceable_refined.tsv` with column `refined_class == 'orphan'`

**Impact:** Phase 9 is completely non-functional as currently wired.

---

## 3. DOCUMENTATION UPDATES NEEDED

### 3.1 README.md
- Phase numbering is confusing (jumps 01 → 01b → 02b → 04)
- Missing description of what Phase 7 actually does now
- SLURM example paths are wrong

### 3.2 IMPLEMENTATION_PLAN.md
- Still describes old SwissProt DIAMOND workflow for Phase 7
- Doesn't reflect Phase 9 being broken
- Phase dependencies diagram is incomplete

---

## 4. ORGANIZATION PROPOSAL FOR FAST SLURM

### Current State
Scripts use inconsistent numbering:
```
01_phase1.py
01b_merge_bk_lb_loci.py
02b_detect_empty_blocks_helixer.py
04_detect_targets_helixer.py
04b_extract_flanking_helixer.py
05_classify_targets_helixer.py
05b_validate_novel_loci.py
06_extract_sequences_helixer.py
07_annotate_flanking.py
08a_generate_locus_matrices_helixer.py
08b_generate_summary_matrices_helixer.py
09_detect_novel_loci_helixer.py
```

**Problems:**
- Phases 2, 3 are skipped (confusing)
- "b" suffix inconsistent (01b vs 02b vs 04b vs 05b)
- Phase 9 overlaps conceptually with Phase 5b

### Proposed Renumbering
```
Phase 1: Reference locus discovery
  01_discover_loci.py              (was 01_phase1.py)
  02_merge_reference_loci.py       (was 01b_merge_bk_lb_loci.py)

Phase 2: Target detection
  03_detect_targets.py             (was 04_detect_targets_helixer.py)
  04_extract_flanking.py           (was 04b_extract_flanking_helixer.py)

Phase 3: Synteny classification
  05_detect_synteny_blocks.py      (was 02b_detect_empty_blocks_helixer.py)
  06_classify_targets.py           (was 05_classify_targets_helixer.py)
  07_validate_novel_loci.py        (was 05b_validate_novel_loci.py)

Phase 4: Output generation
  08_extract_sequences.py          (was 06_extract_sequences_helixer.py)
  09_annotate_flanking.py          (was 07_annotate_flanking.py)
  10_generate_matrices.py          (merge 08a + 08b)
```

### For Fast SLURM Array Execution

**Single entry point script:**
```bash
#!/bin/bash
# run_family.sh FAMILY_NAME
FAMILY=$1
python run_pipeline.py --family "$FAMILY" --config pipeline_config.yaml
```

**Python wrapper with timing:**
```python
# run_pipeline.py
import time
from pathlib import Path

def run_phase(name, func, *args, **kwargs):
    start = time.time()
    result = func(*args, **kwargs)
    elapsed = time.time() - start
    print(f"[{name}] completed in {elapsed:.1f}s")
    return result

# Each phase as a function call with timing
```

**SLURM array:**
```bash
#SBATCH --array=0-32
#SBATCH --time=00:30:00   # 30 min should be plenty per family
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4

FAMILIES=(...)
./run_family.sh ${FAMILIES[$SLURM_ARRAY_TASK_ID]}
```

---

## 5. TIME BENCHMARKING

### Expected Timing (per family)
Based on operations:
- Phase 1 (LOC parsing, FASTA extraction): **<10s**
- Phase 2 (BK-LB merge, DIAMOND): **30-60s**
- Phase 3 (Target DIAMOND vs all proteomes): **1-3 min** (parallelizable)
- Phase 4 (Flanking extraction): **<30s**
- Phase 5 (Synteny DIAMOND + clustering): **1-2 min**
- Phase 6 (Classification): **<30s**
- Phase 7 (Sequence extraction): **<30s**
- Phase 8 (Header parsing): **<10s**
- Phase 9 (Matrix generation): **<30s**

**Total: ~5-8 minutes per family** if everything works correctly.

### Recommendation
Add timing wrapper to each phase. Example:
```python
# In each script main():
import time
start = time.time()
# ... do work ...
print(f"Phase completed in {time.time() - start:.1f}s")
```

---

## 6. PRIORITY FIXES

### Critical (Pipeline Won't Run)
1. Fix SLURM `cd` path (line 30)
2. Fix Phase 7 script name in SLURM
3. Fix Phase 2b argument (`--phase4-targets` not `--phase4-dir`)
4. Fix Phase 5b missing `--phase1-dir` argument
5. ~~Fix or remove Phase 9 (broken input expectations)~~ **DONE: Deprecated 2026-01-04**

### High (Incorrect Results)
1. Phase 2b concordance checking won't work without correct target file
2. Phase 5b reference proteome lookup will fail

### Medium (Documentation)
1. Update README with actual Phase 7 behavior
2. Fix IMPLEMENTATION_PLAN phase descriptions
3. Document expected runtime per phase

### High (Architectural Improvement)

1. **Phase 1: Use genomic coordinates directly instead of LOC lookup** ✅ **DONE 2026-01-04**

   Implemented `--coordinates-file` parameter for Phase 1. Now supports two modes:
   - Legacy mode: `--loc-ids` (unchanged)
   - Coordinates mode: `--coordinates-file` with pre-computed TSV

   Example:
   ```bash
   python 01_phase1.py \
       --coordinates-file top50_effector_genomic_coordinates.tsv \
       --gene-family top50_effectors \
       --output-dir outputs/phase1
   ```

### High (Missing Features)

1. **Extract target NR annotations** ✅ **DONE 2026-01-04**

   Modified `06_extract_sequences_helixer.py` to:
   - Parse NR annotations from Helixer FASTA headers
   - Add `helixer_annotation` column to `length_qc.tsv`
   - Display annotations in verbose output

   Now the `length_qc.tsv` includes a `helixer_annotation` column for each target.

2. **Extract CDS nucleotide sequences** - Required for positive selection (dN/dS) analysis. Currently only protein sequences are extracted.

   **Implementation:** Parse Helixer GFF3 for exon/CDS coordinates, extract from genome assembly, splice together handling strand.

3. **Extract gene structure (exons + introns)** - Exon/intron architecture is phylogenetically informative and useful for gene structure evolution analysis.

   **Implementation:** Record exon and intron coordinates, calculate lengths and phases. Output as `*_structure.tsv` per target with columns:
   ```
   target_id  genome  scaffold  strand  gene_start  gene_end  n_exons  exon_lengths  n_introns  intron_lengths  intron_phases  total_cds_length
   ```

### Low (Cleanup)
1. Remove test_run/ debris or document as examples
2. Consider renumbering for clarity
3. Add timing instrumentation
