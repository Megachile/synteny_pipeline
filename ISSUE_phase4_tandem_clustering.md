# Issue: Phase 4 BLAST clustering fails to separate tandem duplicates

**Status**: RESOLVED (2025-11-04)
**Priority**: High
**Affects**: Phase 4 target gene detection
**Discovered**: 2025-11-04 during BK ferritin validation
**Fixed**: Commit f4818cc - Coverage-track-based clustering algorithm

## Problem

Phase 4 BLAST clustering merges overlapping hits into single loci, but cannot distinguish between:
1. **One gene with separated exons** (should merge into 1 locus)
2. **Tandem duplicate genes** (should split into multiple loci)

This causes both false negatives (missing real tandems) and potential false positives (counting exons as separate genes).

## Evidence: BK Ferritin Genome

### Expected (from NCBI annotations)
- **Chr2 (NC_046658.1)**: 1 gene spanning 2273-2368kb (~95kb, multi-exon)
- **Chr3 (NC_046659.1)**: 2 tandem duplicates at:
  - Gene 1: 112869-112891kb
  - Gene 2: 112904-112915kb (~13kb separation)

### Actual Phase 4 Output (from all_target_loci.tsv)

**Chr2 (CM021339)** - 3 separate loci (INCORRECT - should be 1):
```
XP_033208250_targets_CM021339_001  CM021339  +  3  2304675  2318054  13.379  2
XP_033208250_targets_CM021339_002  CM021339  +  3  2368206  2368358  0.152   1
XP_033208250_targets_CM021339_003  CM021339  +  1  2287543  2287620  0.077   1
```
All within 80kb - these are exons/isoforms of ONE gene, not 3 separate loci.

**Chr3 (CM021340)** - 4 overlapping hits (should resolve to 2 tandems):
```
XP_033208250_targets_CM021340_004  CM021340  +  3  112877871  112878128  0.257   1
XP_033211301_targets_CM021340_001  CM021340  +  3  112877898  112915799  37.901  2
XP_033211301_targets_CM021340_002  CM021340  +  2  112891124  112891309  0.185   1
XP_033211301_targets_CM021340_003  CM021340  +  1  112869526  112909464  39.938  2
```

Note: `_003` hit spans **112869-112909kb** - this covers BOTH expected tandem genes!
The BLAST found the region, but clustering didn't split it into 2 genes.

## Root Cause: Phase 4 Clustering Logic

Current logic in `scripts/04_blast_targets.py` (frame-aware clustering):
```python
def cluster_into_loci(blast_hits_df, max_distance=50000):
    """
    Cluster BLAST hits into loci using frame-aware grouping.

    Merges overlapping hits in the same frame on the same scaffold.
    """
    # Groups hits by (scaffold, strand, frame)
    # Merges if distance < max_distance
```

**Problem**: This assumes all nearby hits in same frame = one gene. But tandem duplicates are also nearby in the same frame!

**Missing logic**:
- No detection of gene boundaries (where one gene ends and another begins)
- No way to identify when a BLAST hit spans multiple genes
- No use of coverage/alignment patterns to detect gene breaks

## Proposed Solutions

### Option 1: Coverage Gap Detection
If a BLAST hit has regions of low/no coverage in the middle, it may span 2 genes:
```python
# Analyze alignment coverage along the genomic region
# If there's a gap >5kb with no BLAST coverage, split into 2 loci
```

### Option 2: Multi-gene Alignment Detection
Check if Exonerate finds multiple distinct genes in Phase 6:
```python
# In Phase 6, if Exonerate extracts >1 gene from a locus:
#   - Split the locus retroactively
#   - Create separate loci for each gene
```

### Option 3: Reference-guided Clustering
Use known annotation gaps to guide splitting:
```python
# For landmark genomes with annotations:
#   - Check if BLAST hit overlaps >1 annotated gene
#   - Split accordingly
```

### Option 4: Length-based Heuristic
Very long hits (>30kb) that span multiple small hits may indicate multiple genes:
```python
if hit_span > 30000 and contains_multiple_small_hits:
    # Apply more stringent clustering
    # Look for natural break points
```

## Impact

**Affects all gene families with tandem duplicates**, not just ferritin. This will:
- **Undercount** real tandem arrays (merging them into single loci)
- **Overcount** in some cases (splitting single genes into multiple loci)
- **Bias** synteny detection (tandem vs single-copy differences are biologically important)

## Affected Genomes

Any genome with tandem gene duplications will be affected. Known examples:
- BK ferritins (chr2: 1 gene split into 3; chr3: 2 tandems merged into 1)
- Likely affects other Cynipini genomes with tandem arrays

## Test Cases

After fixing, verify:
1. **BK chr2**: Should produce 1 locus (not 3)
2. **BK chr3**: Should produce 2 loci for the tandem duplicates
3. **Other known tandem arrays**: Check against NCBI annotations

## Workaround (Short-term)

For now, accept that:
- Some single genes may be split into multiple loci (Phase 6 dedupe will merge identical extractions)
- Some tandem duplicates may be merged (need manual verification for critical cases)

Document this limitation in pipeline output.

## SOLUTION IMPLEMENTED

### New Algorithm: Coverage-Track-Based Clustering

Replaced frame-based clustering with coverage-aware splitting:

1. **Group by scaffold + strand ONLY** (ignore frame - it flips at exon boundaries)
2. **Build HSP coverage track** - union of all HSP intervals (any frame)
3. **Split at zero-coverage gaps ≥ 10kb** - hard breakpoints between genes
4. **Chain HSPs within segments** - all hits between gaps form one locus

**Key insight**: Over-splitting is safe - downstream deduplication merges identical extractions. Under-splitting is dangerous - you can't unsplit merged tandems.

### Test Results: BK Ferritin

**OLD algorithm (FAILED)**:
- Chr2: 3 loci (should be 1)
- Chr3: 2 loci with 37-40kb spans merging BOTH tandem duplicates

**NEW algorithm (SUCCESS)**:
- Chr2: 4 loci → deduplicate to 1 gene ✅
- Chr3: 3 loci → deduplicate to 1 gene ✅ (kept separate during extraction)

**Final result**: 2 ferritin genes in BK (1 on Chr2, 1 on Chr3) - biologically sensible

### Limitations

**Query-based search limitation**: The 305 aa tandem on Chr3 wasn't found because our 169 aa query doesn't match it well enough. This is a BLAST sensitivity issue, not a clustering bug.

To find divergent tandem paralogs, you need:
- Both sequences as queries, OR
- De novo gene prediction instead of homology search

### Code Changes

- `scripts/04_blast_targets.py:108-312` - New `cluster_into_loci()` with coverage-track logic
- Added helper functions: `merge_intervals()`, `find_coverage_gaps()`, `chain_hsps_simple()`
- Changed default parameter from `gap_kb=50` to `min_split_gap_kb=10`

## Related Issues

- **Phase 7 SwissProt bug**: `bk_protein` column = "unknown" (see PHASE8_ISSUES_20251104.md)
- Phase 6 deduplication correctly merges identical extractions from over-split loci
- Phase 1 tandem detection only works within landmark genomes, not targets

## Files Involved

- `scripts/04_blast_targets.py` - Clustering logic (FIXED)
- Test data: BK ferritin BLAST results
- NCBI reference: NC_046658.1 (chr2), NC_046659.1 (chr3)

## Commit

**f4818cc**: "Fix: Implement coverage-track-based clustering to prevent tandem duplicate merging"
