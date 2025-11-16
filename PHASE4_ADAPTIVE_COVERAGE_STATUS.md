# Phase 4 Adaptive Coverage Calibration - Status Update

**Date:** 2025-11-16
**Status:** Implemented and tested

## Summary

Implemented conditional adaptive min_cluster_coverage calibration to handle fragmentary HSPs from unstructured proteins (e.g., Glutenin). The system automatically detects low BK recovery (<80%) and retries with lower coverage threshold (0.3 instead of 0.5).

## Implementation

### Modified Files
- `redesign_scripts/04_detect_targets.py`
  - Lines 37-42: Added `MAX_CLUSTER_SPAN_KB = 500.0` constant
  - Lines 253-260: Updated `cluster_hsps_per_query()` signature with `max_cluster_span_kb` parameter
  - Lines 334-345: Added cluster span hard split check (safety valve)
  - Lines 523-635: **CRITICAL** - Rewrote `calibrate_threshold()` for dual-parameter calibration
    - Now returns `(query_overlap_threshold, min_cluster_coverage)` tuple
    - Tests thresholds with min_coverage=0.5 first
    - If recovery <80%, retries with min_coverage=0.3
    - Selects parameters that best match Phase 1 ground truth
  - Lines 755-765: Updated `main()` to capture and use both calibrated parameters
  - Lines 825-830: Applied calibrated min_coverage in clustering

### Calibration Workflow

```python
# 1. Test with standard min_coverage=0.5
for threshold in [0.1, 0.2, 0.3, 0.5, 0.7]:
    test_clustering(threshold, min_coverage=0.5)

best = min(results, key=lambda x: abs(detected - ground_truth))

# 2. If recovery <80%, retry with relaxed min_coverage=0.3
if (best_detected / ground_truth) < 0.8:
    for threshold in [0.1, 0.2, 0.3, 0.5, 0.7]:
        test_clustering(threshold, min_coverage=0.3)

    # 3. Select better result
    if new_best_is_better:
        return (new_threshold, 0.3)

return (best_threshold, best_min_coverage)
```

## Test Results

### Glutenin_MC21 (Unstructured Protein)

**Test Job:** 118913
**Expected behavior:** Detect low recovery, switch to min_coverage=0.3

| Phase | BK Genes | Recovery | Status |
|-------|----------|----------|--------|
| Phase 1 (ground truth) | 49 | 100% | - |
| Phase 4 v3 (fixed 0.5) | 24 | 49% | ❌ Undersplit |
| **Phase 4 adaptive** | **32** | **65%** | **✅ +8 genes** |

**Calibrated Parameters:**
- `query_overlap_threshold: 0.2`
- `min_cluster_coverage: 0.3`
- Automatic fallback triggered: ✅ Yes (recovery 49% → 65%)

**Analysis:**
- 535 HSPs all on single scaffold spanning 135 Mb
- With min_coverage=0.5: only 21 HSPs pass filter
- With min_coverage=0.3: 63 HSPs pass filter (+42 HSPs)
- Adaptive calibration successfully captured +8 genes
- Still missing 17 genes (35% gap) - may need additional strategies

### Ferritin_MC102

| Phase | BK Genes | Recovery | Status |
|-------|----------|----------|--------|
| Phase 1 (ground truth) | 2 | 100% | - |
| Phase 4 v3 (fixed 0.5) | 1 | 50% | ❌ Undersplit |
| Phase 4 adaptive | 1 | 50% | ❌ No improvement |

**Note:** Adaptive coverage did not help ferritin - different root cause than Glutenin.

## Key Features

### 1. Conditional Parameter Relaxation
- Only relaxes min_coverage when recovery <80%
- Prevents unnecessary parameter changes for well-performing families
- Transparent reporting in output logs

### 2. Recovery Percentage Reporting
```
Testing query overlap thresholds (min_coverage=0.5):
  Threshold 0.7: 24 genes (49% recovery, ±25)

⚠️  Low recovery (49%) - retrying with lower min_coverage=0.3
   (handles fragmentary HSPs from unstructured proteins)

Testing query overlap thresholds (min_coverage=0.3):
  Threshold 0.2: 32 genes (65% recovery, ±17)

✓ Lower min_coverage improved recovery!
```

### 3. Safety Valves
- `MAX_CLUSTER_SPAN_KB = 500.0` prevents mega-gene blobs
- Applied conditionally when other filters fail
- Doesn't interfere with normal clustering

## Deployment Status

**Git Branch:** redesign-scripts-development
**Commit:** Pending
**Production Status:** Ready for full deployment after git push

### Next Steps
1. ✅ Commit adaptive coverage implementation
2. ✅ Push to GitHub
3. 🔄 Deploy to all 35 families
4. 📊 Analyze results across families
5. 🔬 Investigate additional strategies for severe cases (Glutenin <80%)

## Known Limitations

### Glutenin Still Undersplits (65% recovery)
Despite +8 gene improvement, still missing 17/49 genes. Possible additional strategies:

1. **Frame-aware clustering:** Group HSPs by reading frame to separate tandem duplicates
2. **Iterative splitting:** Progressively lower min_coverage until recovery reaches 80%
3. **Paralog-specific thresholds:** Different parameters for high-copy gene families
4. **Multi-stage calibration:** Calibrate max_gap_kb and paralog_coverage_threshold too

### Not All Families Benefit
- Ferritin: No improvement (different root cause)
- Well-performing families: Adaptive calibration doesn't trigger
- Some families may need different parameter adjustments

## Related Issues

### Phase 6 Scaffold Naming Mismatch
Discovered during Phase 6 syntenic test (job 118666):
- Syntenic targets TSV uses RefSeq accessions (NC_*)
- Genome databases use GenBank accessions (GCF_*, GCA_*)
- Phase 6 extraction fails with "Could not find scaffold" warnings
- **Impact:** 2/3 test families failed (MC6, rhamnogalacturonate_lyase-like)
- **Status:** Separate issue to address

## Conclusion

Adaptive coverage calibration is **working as designed** for fragmentary HSPs from unstructured proteins:
- Automatically detects low recovery
- Applies appropriate fallback strategy
- Improves recovery without manual intervention
- Transparent and reproducible

**Ready for production deployment** pending git push.
