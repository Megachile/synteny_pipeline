# Session Notes – 2025‑11‑06

## What We Investigated

- Drove the hybrid synteny workflow end‑to‑end with a focus on BK (positive control) and MC117.
- Added two post‑Phase‑6 steps to improve specificity without capping:
  - 06c: Exonerate‑guided block re‑selection (`scripts/06c_reselect_blocks_with_exonerate.py`)
  - 6b: Strict Exonerate refinement (`scripts/06b_refine_targets_with_exonerate.py`)
- Tightened Phase 2 clustering to prevent mega‑blocks.

## Key Results (MC117, BK)

- Original Phase 5: 43 targets, 3/27 loci (11%) — under‑recovery.
- Exonerate loose (pre‑strict): 2074 targets, 27/27 loci — too permissive.
- Exonerate strict (score ≥ 40%): 390 targets, 27/27 loci.
- Exonerate stricter defaults (score ≥ 60%, qcov ≥ 60%, cluster gap 3000bp): 312 targets, 27/27 loci.

Observation: 6b recovers all 27 BK loci but still ~12 targets per “locus” in the BK chr5 region — the targets correspond to many nearby BK loci packed within ~1.8 Mb.

## Biological Context (BK chr5, MC117)

- Measured gaps between consecutive MC117 target genes on BK chr5:
  - Average gap ~163 kb; two large gaps ~501 kb and ~800 kb; total span ~1.82 Mb.
  - Target genes occupy ~1.4% of the region.
- Conclusion: This is a dispersed cluster, not a tight tandem array.
  - Phase 1 is reasonable to keep these as separate loci.
  - Phase 8 “mega‑block” behavior lumped many dispersed loci when blocks spanned the entire region.

## Why Phase 2 Tighten Didn’t Move BK Counts

- Even with `--max-gap-kb=100` and chain logic, if a block’s bounds still cover a broad “dispersed cluster,” 6b will find many high‑confidence targets inside it (because they are real and inside the block).
- The root issue was mega‑blocks spanning 0.5–2 Mb; Phase 2 now includes:
  - Colinearity by query order
  - Hard span limit per block (`MAX_BLOCK_SPAN_KB`, default 300)
  - One best hit per query prior to chaining
  - (Planned) split at zero‑coverage troughs

## Current Controls (No Capping, Agnostic)

- 06c: Reselects one best block per (locus, genome) by Exonerate overlap.
- 6b strict:
  - Filters: min score %, min query coverage (from Exonerate “Query range” and query length), overlap with selected block
  - Dedup: one best per query per locus
  - Clustering: group placements within cluster‑gap and keep one best per cluster
  - No tandem capping anywhere (including BK) to remain fully agnostic.

## Why We Still See ~12 Targets Per Locus

- Those are real BK target genes dispersed across ~1.8 Mb in the same broad block region.
- Strict 6b prunes duplicates but retains distinct genes (clusters) per locus because the block is still large enough to encompass many genes.

## Next Steps (To Do Next Session)

1) Phase 2 segmentation at gaps and density troughs
   - Build coverage track over locus region from tBLASTn hits and split blocks where gap between adjacent hits exceeds `--gap-split-kb` (e.g., 50–100 kb).
   - Maintain colinearity and span limit; keep segments with ≥ `MIN_PROTEINS_FOR_BLOCK`.

2) 6b thresholds tuning (BK stress test)
   - Try `--min-score 75`, `--min-query-cov 0.85`, `--cluster-gap-bp 10000` on MC117.
   - Optional: require ≥50–70% of aligned query span to lie inside the block (filters edge alignments).

3) Optional Phase 4 query set refinement (agnostic)
   - Build query FASTA per locus from Phase 1 cluster members only (for BK), rather than family‑wide. Keeps discovery agnostic in non‑BK but reduces BK pre‑inflation of queries.

4) Batch reevaluation after Phase 2 fix
   - Re-run 02→02b→06c→6b strict across finished families and regenerate `outputs/refined_bk_recovery_summary_exo.tsv`.

## Files & Scripts We Touched

- Overview: `PIPELINE_OVERVIEW.md`
- Phase 2: `scripts/02_synteny_detection.py` (added colinearity chains, span limit)
- Phase 3: `scripts/03_filter_blocks.py` (scoring hooks)
- Exo reselect: `scripts/06c_reselect_blocks_with_exonerate.py`
- Exo refine (strict): `scripts/06b_refine_targets_with_exonerate.py`
- Batch: `scripts/run_exo_refinement_all.py`, `scripts/run_exo_refinement_all.slurm`

## Open Questions

- How aggressive should `MAX_BLOCK_SPAN_KB` and gap‑split be to respect the dispersed nature of MC117?
- Should we add a “min‑inside‑block fraction” in 6b to filter edge placements?
- Do we prefer segmenting the chr5 region into ~12 smaller blocks (Phase 2) or adjusting Phase 4 queries to avoid family‑wide inflation?

