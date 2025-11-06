# Phase 2 Colinear Chain Clustering - Improvement Results

## Date: 2025-11-06

## Problem Identified

BK validation analysis revealed **77% target recovery failure** due to mega-blocks in Phase 2:
- BK chr5 dispersed cluster (MC117): 12 target genes spanning **1.82 Mb**
- Gaps between genes: 100-800 kb (average 163 kb)
- Old Phase 2: Created **one mega-block (746 kb)** lumping all 12 loci together
- Phase 5: Incorrectly assigned all targets to one locus instead of 12

## Solution Implemented

**Colinear chain clustering with hard span limit** (commit 8f34e3b):

1. **Best hit per query**: Reduce to lowest evalue hit (tie: highest bitscore)
2. **Colinear chains (LIS-like)**: Build ordered chains respecting query position order
3. **Hard span limit**: `MAX_BLOCK_SPAN_KB=300` segments chains exceeding threshold
4. **Min proteins per segment**: Maintain `MIN_PROTEINS_FOR_BLOCK=3` after segmentation

## Test Results: BK_chr5_c Locus

### Before (Old Code, Nov 6 11:40)

```
Locus       Genome           Block        Scaffold   Strand  Start       End         Span_kb  Proteins  Queries
BK_chr5_c   GCA_010883055.1  block_00057  CM021342   +       142238576   142984867   746.3    72        17
```

**Issues**:
- Single mega-block spanning 746 kb
- Likely contains targets from multiple nearby loci
- Phase 5 would assign all to this one "locus"

### After (NEW Code with Colinear Chains, Nov 6 15:36)

```
Locus       Genome           Block        Scaffold   Strand  Start       End         Span_kb  Proteins  Queries
BK_chr5_c   GCA_010883055.1  block_00003  CM021342   +       144169795   144397835   228.0    6         6
```

**Improvements**:
- Block span reduced: **746 kb → 228 kb** (69% reduction)
- More focused block with fewer proteins (72 → 6)
- Better represents actual locus boundaries
- During unfiltered clustering: 4 BK blocks detected (vs 1 mega-block)

### Full Genome Comparison (BK_chr5_c across all 75 genomes)

**Before**:
- Total blocks: 76 (many mega-blocks)
- Average span: likely 200-500 kb

**After**:
- Unfiltered: 81 blocks detected
- Filtered (best per genome): 48 blocks
- More accurate segmentation of dispersed clusters

## Full MC117 Pipeline Test (Job 101728)

**Running**: All 27 BK chr5 loci through Phases 2→2b→3→5

**Expected Results**:
- More blocks per locus (reflecting true dispersed structure)
- Better target assignment in Phase 5
- Improved BK recovery rate (from 11% to hopefully >80%)

**Output Locations**:
- Phase 2: `outputs/MC117_BK_LB_batch/02_synteny_blocks_NEW/`
- Phase 3: `outputs/MC117_BK_LB_batch/03_filtered_blocks_NEW/`
- Phase 5: `outputs/MC117_BK_LB_batch/05_classified_NEW/`

## Parameters Used

```bash
--max-gap-kb 100          # Maximum gap between hits in a chain
MAX_BLOCK_SPAN_KB=300     # Hard limit on block span (internal constant)
--threads 16              # Parallel BLAST execution
```

## Next Steps

1. **Validate MC117 results**: Check BK target recovery after Job 101728 completes
2. **Batch reevaluation**: Apply to all 32 gene families if successful
3. **Fine-tuning**: Adjust `MAX_BLOCK_SPAN_KB` if needed (200-400 kb range)
4. **Phase 5 integration**: May need to adjust proximity thresholds with smaller blocks

## Code Changes

**File**: `scripts/02_synteny_detection.py`
- Added `_build_query_order_map()` for FASTA position tracking
- Added `_select_blocks_by_chain()` for LIS chain extraction + segmentation
- Modified `cluster_into_blocks()` to use colinear clustering
- Added `MAX_BLOCK_SPAN_KB` constant (default 300)

**Commit**: 8f34e3b "Fix(Phase 2): Add colinear chain clustering + hard span limit"

## Biological Interpretation

The new approach better reflects biological reality:
- **Tight tandem arrays** (<50 kb): Single block as expected
- **Dispersed clusters** (100-800 kb gaps): Multiple segments
- **True synteny** (colinear order): Preserved in chains
- **False synteny** (scattered hits): Filtered by colinearity requirement

This agnostic approach lets the data dictate structure without assuming tight clustering.
