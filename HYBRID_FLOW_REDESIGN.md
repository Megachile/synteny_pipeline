# Hybrid Synteny Workflow — Redesign Plan (MC6, MC117 pilots)

## Motivation

- Tandem-rich regions frequently fragment into multiple “loci” or multiple short synteny blocks, which propagates ambiguity downstream.
- We want locus-first reasoning: merge long tandem strings into a single locus envelope early, use robust flanking anchors to detect synteny in other genomes, then use target hits to choose the winning block within that envelope.
- Evidence-of-absence should only be called when there is a strong synteny signal (good envelope) and no target hits in a reasonable nearby window.

## High-Level Principles

- Treat tandem strings as one locus if they are contiguous at the gene-level, even with interspersed non-family genes.
- Carry a locus-scale parameter through all phases to adapt thresholds to local gene density (“tight” vs “wide” regions).
- Prefer blocks on the expected chromosome (RagTag-aligned reference chromosome) when available.
- Keep all blocks through detection; reduce only at the point we have target evidence to decide.

## New Locus-Scale Parameter

- Name: `locus_scale` (float), propagated across phases.
- Computation (Phase 1):
  - Gather per-gene genomic positions for all flanking anchors (U/I/D) and compute:
    - Median inter-anchor spacing (bp) and its MAD.
    - Locus span (envelope) = [min(anchor), max(anchor)].
    - Local gene density = (number of anchors + interspersed non-family genes) / span.
  - Normalize against cohort medians to produce `locus_scale`:
    - Example: `locus_scale = clamp(0.5, 2.0, median_spacing / cohort_median_spacing)`
    - Tight region → `locus_scale < 1.0`; wide region → `> 1.0`.
- Use sites:
  - Phase 1: dynamic locus separation threshold.
  - Phase 2: set `max-gap-kb`, `max-span-kb`, and optionally `max-targets` as functions of `locus_scale`.
  - Phase 4: set target split gap (`min_split_gap_kb`) as a function of `locus_scale`.
  - Phase 5: dynamic nearby window adds a `locus_scale` term and caps.

## RagTag Chromosome Leverage

- Inputs: `data/ragtag_output/<genome>/*` and `data/ragtag_dbs/*.nhr`.
- Expected chromosome:
  - Parse locus_id for BK/LB (e.g., `BK_chr2_a` → chr2). Store `expected_chromosome="chr2"` in `locus_definitions.tsv`.
  - Maintain a per-genome alias map from RagTag headers to reference chromosomes (use `normalize_scaffold` + regex on header tokens like `CM/NC/NW` and `_RagTag` suffixes). Persist as `ragtag_chrom_alias.tsv` per run.
- Scoring bonuses:
  - Phase 2 (optional): annotate blocks with `chrom_match` if scaffold aliases to expected chromosome.
  - Phase 3: include a chromosome-match bonus in block scoring (soft bias, not a hard filter).
  - Phase 5: prefer blocks within the chosen locus whose scaffold aliases to the expected chromosome (tie-breaker below overlap/proximity and target-support).

## Phase-by-Phase Design

### Phase 1 — Locus Discovery (BK/LB)

- Widen separation adaptively:
  - Replace fixed `MAX_LOCUS_GAP_KB` with `MAX_LOCUS_GAP_KB * locus_scale` (capped)
  - Merge long tandem strings into a single locus envelope unless an inter-anchor distance exceeds `k * locus_scale`.
- Outputs:
  - `locus_definitions.tsv` adds columns:
    - `flanking_file` → set to deduplicated file (`*_flanking_dedup.faa`).
    - `locus_span_start`, `locus_span_end`, `locus_span_kb`.
    - `locus_scale` (0.5–2.0 typical).
    - `expected_chromosome` (from `locus_id`), plus `bk_chr`, `lb_chr` if both sources exist.
- Keep targets separate: targets are not part of anchors; we use U/I/D anchors as flankings.

### Phase 2 — Synteny Detection (tBLASTn)

- Inputs: deduplicated flanking file preferred (already patched) and `locus_scale`.
- Parameters:
  - `max-gap-kb = base_gap_kb * locus_scale` (e.g., base 300–500).
  - `max-span-kb = base_span_kb * locus_scale` (e.g., base 600–1200 with cap).
  - `max-targets`: 1000 default; optionally `ceil(1000 * locus_scale)`.
  - `merge-adjacent = true` with `merge-gap-kb` and `merge-span-kb` scaled by `locus_scale`.
- Behavior:
  - Keep all blocks (≥2–3 anchors) and include `chrom_match` annotation using RagTag alias map.
  - Emit per-locus TSVs; aggregated summaries include span and matches distributions.

### Phase 3 — Aggregate + Threshold (No Best-Per-Genome)

- Add a `--threshold-only` mode:
  - Write `all_synteny_blocks.tsv` (unchanged) and `all_synteny_blocks_ge{min}.tsv` (blocks with `num_query_matches >= min`).
  - Keep current `synteny_blocks_filtered.tsv` for strict mode.
- Scoring (when used): add soft bonus for `chrom_match` and slight density weighting scaled by `locus_scale`.

### Phase 4 — Target Detection

- Reduce over-splitting:
  - `min_split_gap_kb = base_split_kb * locus_scale` (e.g., base 30–50 kb; cap).
  - Group by scaffold+strand only (already done) and allow mixed frames.
- Optionally bias search to envelope first: when scanning a genome, prioritize scaffolds that host Phase 2 blocks for that locus.

### Phase 5 — Locus-First Classification

- Input blocks: use `all_synteny_blocks_ge{min}.tsv` (not pre-reduced best).
- Locus envelopes:
  - Build merged segments per locus (`segment_gap_kb` scaled by `locus_scale`).
  - Assign targets by (overlap → proximity) to the locus envelope first, then to the best block within the chosen locus.
- Block choice within locus:
  - Primary: overlap in bp with target chain.
  - Secondary: number of targets supported by the block.
  - Tertiary: `num_query_matches` and `chrom_match` bonus.
- Evidence-of-absence:
  - Conditions: good locus envelope in genome, zero targets within dynamic window `max(min_kb, frac * span, locus_scale_term)` capped by upper bound.
  - Optional relaxed tBLASTn retry on the envelope scaffold before final call.
- Diagnostics: per-locus placement table and histogram of min distances from targets to nearest block (kept).

## Implementation Phases (MC6, MC117)

1. Phase 1 updates
   - Compute and emit `locus_scale`, spans, expected chromosomes.
   - Point `flanking_file` to `*_flanking_dedup.faa`.
   - Parameterize separation with `locus_scale`.

2. Phase 2 parameterization
   - Read `locus_scale` and adjust gap/span/merge windows.
   - Keep-all blocks; annotate `chrom_match` using alias map.

3. Phase 3 threshold-only output
   - Add `--threshold-only` and emit `all_synteny_blocks_ge{min}.tsv`.

4. Phase 4 scaling
   - Use `locus_scale` to set `min_split_gap_kb`; add optional envelope-first priority.

5. Phase 5 locus-first + chromosome bias
   - In grouped mode, prefer blocks with `chrom_match` when ties.
   - Accept `all_synteny_blocks_ge{min}.tsv` as the blocks input.

6. A/B tests (MC6, MC117)
   - Strict mode (current): best-per-genome Phase 3 + overlap-only.
   - Locus-first mode: threshold-only Phase 3 + grouped + dynamic nearby.
   - Metrics: synteny spans, placed targets, ambiguous placements, absence calls.

## Data Artifacts & Compatibility

- New columns in `locus_definitions.tsv`: `locus_scale`, spans, `expected_chromosome`.
- New `ragtag_chrom_alias.tsv`: per-genome mapping of scaffolds to reference-style chromosomes.
- New `all_synteny_blocks_ge{min}.tsv`: threshold-only blocks for classification.

## Open Questions

- How to best calibrate `locus_scale` cohort normalization (per family vs global)?
- Should chromosome bonus propagate into Phase 2 chain segmentation (likely “no”, keep as annotation)?
- When RagTag names don’t carry reference chromosome hints, do we infer via synteny to BK/LB only?

## Notes on Current State

- Phase 2 already prefers deduplicated flanking when present (patched).
- Group-by-locus and dynamic nearby are implemented in Phase 5 and can be leveraged immediately once we pass non-reduced blocks.

