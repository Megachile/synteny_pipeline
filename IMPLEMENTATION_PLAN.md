# Helixer Pipeline Enhancement - Implementation Plan

## Executive Summary

The current Helixer pipeline finds target proteins via DIAMOND, then checks flanking genes only around those hits. This misses three important cases:

1. **Empty synteny blocks** - Conserved flanking genes exist but target is absent (deletion/pseudogenization)
2. **Unplaceable tandems** - Some "unplaceables" are actually secondary copies near a syntenic hit
3. **Novel loci** - Genuinely new loci with their own conserved flanking pattern

This plan addresses all three while preventing an explosion of false positives from fragmented genomes.

---

## Problem 1: Empty Synteny Block Detection

### What We Lost

Original tblastn Phase 2/3 ran tblastn on FLANKING GENES against ALL genomes, finding conserved synteny blocks even without target hits. This detects:
- Gene losses (target deleted but flanking conserved)
- Pseudogenes (target degenerated)
- Annotation gaps (target exists but wasn't predicted)

### Proposed Solution

**New script: `02b_detect_empty_blocks_helixer.py`**

Search for BK flanking gene orthologs in ALL Helixer proteomes, then look for conserved clusters that DON'T have a target hit nearby.

### Implementation Details

```
Input:
- Phase 1 flanking proteins (BK_*_flanking.faa)
- All Helixer proteomes
- Phase 4 target hits (to exclude regions already covered)

Process:
1. DIAMOND blastp: BK flanking proteins → ALL Helixer proteomes
2. For each genome, cluster flanking gene hits by scaffold/position
3. Score clusters: (n_BK_flanking_genes_found / n_total_BK_flanking_genes)
4. Filter: require ≥3 BK flanking genes within max_span_kb
5. Check if cluster overlaps with existing target hit → if yes, skip (already covered)
6. If no target overlap → "empty synteny block"

Output:
- empty_synteny_blocks.tsv: locus_id, genome, scaffold, start, end, n_flanking_matched, score
- empty_block_details.tsv: individual flanking gene hits
```

### Key Parameters

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| min_flanking_genes | 3 | Minimum for synteny signal |
| max_span_kb | 500 | Derived from Phase 1 locus spans |
| evalue | 1e-5 | Same as Phase 4 |
| min_score | 0.15 | Same as Phase 5 threshold (3/20) |

### Pitfalls and Mitigations

1. **False positives from scattered hits**
   - Mitigation: Require clustering (not just count) - flanking genes must be physically clustered
   - Test: Check if detected empty blocks are on same scaffold as known targets

2. **Tandem flanking genes creating multiple hits**
   - Mitigation: Deduplicate by BK gene ID before scoring
   - Test: Count unique BK flanking genes, not total hits

3. **Partial matches on short scaffolds**
   - Mitigation: Track scaffold length; require hits on both sides of expected target position
   - Test: Flag scaffolds < 50kb as "low_confidence_scaffold"

### Tests

```
TEST 1: Known empty locus
- Use endoglucanase_Z_MC6 (known to have many empty blocks)
- Compare to tblastn Phase 2/3 combined_synteny_blocks.tsv
- Verify script detects similar block locations

TEST 2: Coverage check
- Empty blocks + target blocks should approximate tblastn synteny blocks
- Count should be similar across genomes

TEST 3: BK/LB empty blocks for novel loci
- If we discover novel loci (Problem 3), BK/LB WILL have empty blocks for those
- Initially expect 0 for known loci; may increase after novel loci detection
```

---

## Problem 2: Unplaceable Duplicate Classification

### Current Behavior

Phase 5 classifies targets as:
- **syntenic**: flanking genes match a reference locus
- **unplaceable**: no flanking gene match above threshold

But "unplaceable" is a catch-all that includes:
- **Tandems**: Secondary copies physically near a syntenic hit
- **Partial duplicates**: Match a known locus but not as well as "best" hit
- **Orphans**: No synteny signal at all
- **Novel**: Genuinely new loci (addressed in Problem 3)

### Proposed Solution

**Modify: `05_classify_targets_helixer.py`** (add subclassification)
**New script: `05b_refine_unplaceables_helixer.py`**

### Implementation Details

**Phase 5 modifications** - Add columns to unplaceable_targets.tsv:
- `best_locus_match`: Locus with highest (non-qualifying) score
- `best_locus_score`: Score for that locus
- `n_flanking_hits`: How many flanking genes had ANY BK hit
- `reason`: Why unplaceable (no_hits, below_threshold, multiple_weak)

**Phase 5b: Refine unplaceables**

```
Input:
- Phase 5 syntenic_targets.tsv
- Phase 5 unplaceable_targets.tsv
- Phase 4 target coordinates

Process:
1. For each unplaceable target:
   a. Check distance to nearest SYNTENIC target on same scaffold
   b. If within tandem_distance_kb → classify as "tandem_of_{locus}"
   c. Else if best_locus_score > 0 → classify as "weak_match_{locus}"
   d. Else → classify as "orphan"

Output:
- unplaceable_refined.tsv: adds tandem_of, weak_match, orphan classifications
- tandem_candidates.tsv: targets near syntenic hits (potential tandems)
```

### Key Parameters

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| tandem_distance_kb | 50 | Within ~10-15 genes of syntenic target |
| weak_match_threshold | 0.05 | At least 1/20 flanking genes match |

### Pitfalls and Mitigations

1. **Tandem clusters with inverted order**
   - Mitigation: Check for nearby syntenic target regardless of strand
   - Test: Manually review tandem candidates to check strand patterns

2. **Distant tandems on same scaffold**
   - Mitigation: Use generous tandem_distance_kb for initial flag, refine later
   - Alternative: Also check if unplaceable is between two syntenic targets for same locus

3. **Weak matches that are false positives**
   - Mitigation: Track which BK flanking genes matched - if only 1-2, might be false
   - Test: Look for functional annotation of matched flanking genes

### Tests

```
TEST 1: Known tandem family
- Use a family known to have tandems (e.g., venom serine protease)
- Verify tandems are correctly identified

TEST 2: Counts should be consistent
- len(syntenic) + len(tandem) + len(weak_match) + len(orphan) == len(all_targets)

TEST 3: Tandem distance distribution
- Plot histogram of tandem distances
- Should show clear bimodal (close tandems vs distant)
```

---

## Problem 3: Novel Locus Detection

### The Challenge

Some unplaceables might represent genuinely new loci - gene duplications to new genomic locations that are CONSERVED across species but weren't in BK/LB Phase 1.

### Why This Could Explode

If we naively call every unplaceable with conserved flanking a "novel locus":
- 61 genomes × potential false positives = explosion
- Fragmented genomes create spurious flanking patterns
- Need reciprocal evidence across multiple genomes

### Proposed Solution

**New script: `09_detect_novel_loci_helixer.py`**

Key insight: A real novel locus should have CONSERVED flanking genes across multiple genomes. We can detect this by doing the Phase 1/1.5 analysis in reverse.

### Implementation Details

```
Input:
- Phase 5 orphan targets (from 05b)
- All Helixer proteomes + GFF3s
- Phase 4b flanking gene data

Process:
1. For each orphan target, get its Helixer flanking genes
2. DIAMOND each orphan's flanking genes against ALL other Helixer proteomes
3. Find "flanking gene clusters" - sets of flanking genes that co-occur across genomes
4. Score each orphan: n_genomes_with_similar_flanking / n_genomes_with_target
5. Filter: require similar flanking in ≥3 genomes

Clustering:
- Group orphans by their flanking gene "signature"
- If multiple orphans share >50% of flanking gene orthologs → same novel locus

Output:
- novel_locus_candidates.tsv: novel_locus_id, n_genomes, signature_genes
- novel_locus_members.tsv: target_id, genome, novel_locus_id, flanking_score
```

### Critical Filters to Prevent Explosion

| Filter | Threshold | Rationale |
|--------|-----------|-----------|
| min_genomes | 3 | Minimum for conservation signal |
| min_flanking_score | 0.20 | Higher than Phase 5 threshold |
| max_novel_loci | 50 | Hard cap per family |
| scaffold_length_min_kb | 100 | Exclude very short scaffolds |
| flanking_both_sides | True | Must have flanking on BOTH sides |

### Additional Quality Filters

```python
def is_high_quality_novel_candidate(target, flanking_data, genome_quality):
    """
    Require multiple quality signals to call a novel locus.
    """
    # Must have flanking genes on both sides
    n_upstream = len([f for f in flanking_data if f['position'].startswith('upstream')])
    n_downstream = len([f for f in flanking_data if f['position'].startswith('downstream')])
    if n_upstream < 2 or n_downstream < 2:
        return False

    # Scaffold must be reasonably long (not short fragment)
    if target['scaffold_length_kb'] < 100:
        return False

    # Genome must have reasonable BUSCO completeness
    if genome_quality.get(target['genome'], {}).get('busco_complete', 0) < 70:
        return False

    return True
```

### Workflow

```
Phase 5 output
    ↓
┌───────────────────┐
│ Phase 5b: Refine  │
│ - tandem          │
│ - weak_match      │
│ - orphan          │
└───────────────────┘
    ↓ (orphans only)
┌───────────────────┐
│ Phase 9: Novel    │
│ Locus Detection   │
│ - reciprocal      │
│ - multi-genome    │
│ - quality filters │
└───────────────────┘
    ↓
novel_loci.tsv
```

### Pitfalls and Mitigations

1. **Explosion from fragmented genomes**
   - Mitigation: Require ≥3 genomes with similar flanking
   - Mitigation: Require flanking on BOTH sides
   - Mitigation: Use BUSCO scores to weight genome quality

2. **False clusters from common genes**
   - Mitigation: Exclude "promiscuous" flanking genes that appear in >50% of orphans
   - Mitigation: Weight flanking genes by specificity

3. **Computational cost**
   - All-vs-all flanking comparison could be expensive
   - Mitigation: Pre-cluster orphans by scaffold before full analysis
   - Mitigation: Use sequence clustering (cd-hit) to reduce redundancy

4. **Overfitting to sequencing artifacts**
   - Mitigation: Require novel locus to be in ≥2 different clades
   - Mitigation: Down-weight genomes from same species

### Tests

```
TEST 1: No novel loci in reference genomes
- BK and LB targets should NOT be called as novel
- They should all be syntenic to known loci

TEST 2: Known novel locus recovery (if any)
- If we have prior knowledge of a novel locus, verify detection

TEST 3: Quality filter effectiveness
- Compare novel locus calls WITH and WITHOUT quality filters
- Quality filters should reduce count by >50%

TEST 4: Flanking gene specificity
- For each novel locus, check if flanking genes are specific or ubiquitous
- Flag if >50% of flanking genes appear in >20% of all targets
```

---

## Implementation Order

### Phase A: Foundation (modify existing)

1. **Modify Phase 5** - Add detailed output for unplaceables
   - Add `best_locus_match`, `best_locus_score`, `n_flanking_hits`, `reason` columns
   - Estimated: 1-2 hours

### Phase B: Unplaceable Refinement

2. **New Phase 5b** - Refine unplaceables
   - Tandem detection based on proximity
   - Weak match classification
   - Estimated: 2-3 hours

3. **Test Phase 5b** on ferritin
   - Verify tandem detection
   - Check classification counts

### Phase C: Empty Blocks

4. **New Phase 2b** - Detect empty synteny blocks
   - DIAMOND flanking → all proteomes
   - Cluster and score
   - Estimated: 3-4 hours

5. **Test Phase 2b** on ferritin
   - Compare to tblastn Phase 2/3 results
   - Verify no false positives in BK/LB

### Phase D: Novel Loci (most complex)

6. **New Phase 9** - Detect novel loci
   - Reciprocal flanking analysis
   - Multi-genome requirement
   - Quality filters
   - Estimated: 4-6 hours

7. **Test Phase 9** extensively
   - Run on multiple families
   - Tune filters to prevent explosion
   - Manual review of top candidates

### Phase E: Integration

8. **Modify Phase 8** - Include new classifications in matrices
   - Add empty block column
   - Add refined unplaceable classification
   - Add novel locus membership
   - Estimated: 2-3 hours

---

## File Structure After Implementation

```
helixer_pipeline/
├── 02b_detect_empty_blocks_helixer.py    # NEW
├── 04_detect_targets_helixer.py          # existing
├── 04b_extract_flanking_helixer.py       # existing
├── 05_classify_targets_helixer.py        # MODIFY
├── 05b_refine_unplaceables_helixer.py    # NEW
├── 06_extract_sequences_helixer.py       # existing
├── 07_annotate_swissprot_helixer.py      # existing
├── 08a_generate_locus_matrices.py        # MODIFY
├── 08b_generate_summary_matrices.py      # MODIFY
├── 09_detect_novel_loci_helixer.py       # NEW
├── run_helixer_full_pipeline.slurm       # MODIFY
└── IMPLEMENTATION_PLAN.md                # THIS FILE
```

---

## Expected Outputs After Full Implementation

For a family like ferritin_MC102:

```
outputs/ferritin_MC102/
├── phase2b_empty_blocks/
│   ├── empty_synteny_blocks.tsv       # Genomes with flanking but no target
│   └── empty_block_details.tsv        # Flanking gene hits for empty blocks
├── phase5_helixer/
│   ├── syntenic_targets.tsv           # Targets with clear locus assignment
│   ├── unplaceable_targets.tsv        # Targets without synteny (ENHANCED)
│   └── all_targets_classified.tsv
├── phase5b_refined/
│   ├── unplaceable_refined.tsv        # With tandem/weak_match/orphan
│   ├── tandem_candidates.tsv          # Potential tandem duplicates
│   └── classification_summary.tsv
├── phase9_novel/
│   ├── novel_locus_candidates.tsv     # Putative new loci
│   ├── novel_locus_members.tsv        # Targets assigned to novel loci
│   └── flanking_signatures.tsv        # Flanking gene patterns
└── phase8_matrices/
    ├── {locus}_matrix.tsv             # Per-locus with empty block column
    └── {family}_summary.tsv           # Family summary with all classifications
```

---

## Risk Assessment

| Risk | Likelihood | Impact | Mitigation |
|------|------------|--------|------------|
| Novel loci explosion | High | High | Multiple quality filters, hard cap |
| False empty blocks | Medium | Medium | Clustering requirement, scaffold length filter |
| Tandem false positives | Low | Low | Distance threshold conservative |
| Computational cost | Medium | Medium | Incremental processing, caching |
| Integration complexity | Medium | Medium | Careful column naming, version tracking |

---

## Success Criteria

1. **Empty blocks**: Detect ≥1 empty synteny block in genomes known to lack ferritin
2. **Tandem detection**: Correctly classify known tandem families
3. **Novel loci**: Generate <20 novel locus candidates per family (not explosion)
4. **Quality**: All classifications traceable back to flanking gene evidence
5. **Performance**: Full pipeline runs in <2 hours per family on SLURM

---

## Questions to Resolve Before Implementation

1. **Genome quality data**: Do we have BUSCO scores for all Helixer proteomes? (Beads task helixer_pipeline-2)
2. **Known tandems**: Which families have known tandems we can test against?
3. **Empty block validation**: Do we have ground truth for empty blocks from tblastn pipeline?
4. **Novel loci priors**: Any prior knowledge of cynipid-specific novel gene loci?

---

## Next Steps

1. Review this plan with user
2. Start with Phase A (modify Phase 5) - lowest risk, immediate value
3. Test on ferritin_MC102 after each phase
4. Iterate on filters based on results
