# Hybrid Workflow – Pipeline Overview (Current Version)

This document summarizes the end‑to‑end flow of the hybrid synteny pipeline as currently implemented, the role of each phase/script, key inputs/outputs, and major parameters. It reflects the latest additions: Exonerate‑guided block re‑selection (06c) and strict Exonerate refinement (6b).

## High‑Level Phases

- Phase 1 (Landmarks & Flanking)
  - Script: `scripts/01_discover_paralogs_gene_level.py`
  - Purpose: Discover landmarks (BK/LB), map proteins→genes (LOC), deduplicate isoforms/tandems, extract ordered flanking proteins (U/D) per locus.
  - Inputs: Reference GFF/FAA (BK/LB), BK/LB DIAMOND DBs
  - Outputs: `outputs/<FAMILY>/phase12_landmark/locus_definitions.tsv`, per-locus flanking FASTA files
  - Key Params: genome‑specific flanking distance caps; max flanking gene count

- Phase 2 (Synteny Detection)
  - Script: `scripts/02_synteny_detection.py`
  - Purpose: tBLASTn of flanking proteins against all target genomes; cluster hits into synteny blocks.
  - Inputs: locus_definitions.tsv (flanking FASTA), genome BLAST DBs
  - Outputs (per locus): `02_synteny_blocks/<LOC>_synteny_blocks.tsv`, `all_synteny_blocks.tsv`
  - Current Logic (updated):
    - Reduce to one best hit per query (lowest evalue, tie: highest bitscore).
    - Build colinear chains by query order vs genome coordinates; segment chains by hard span limit (`MAX_BLOCK_SPAN_KB`) and min hits (`MIN_PROTEINS_FOR_BLOCK`).
  - Important CLI params: `--max-gap-kb` (legacy), `MAX_BLOCK_SPAN_KB` (source constant), min hits.

- Phase 2b (Aggregate)
  - Script: `scripts/02b_aggregate_synteny_blocks.py`
  - Purpose: Combine per‑locus block TSVs to `02_synteny_blocks/all_synteny_blocks.tsv`

- Phase 3 (Best Block Filter)
  - Script: `scripts/03_filter_blocks.py`
  - Purpose: Pick one best block per (locus, genome). Scoring can penalize span and fragmented scaffolds, and reward density.
  - Output: `03_filtered_blocks/synteny_blocks_filtered.tsv` (+ per‑locus files)

- Phase 4 (Target BLAST)
  - Script: `scripts/04_blast_targets.py`
  - Purpose: Combined multi‑query tBLASTn of target proteins; cluster target hits into loci via coverage tracks.
  - Output: `04_target_genes/all_target_loci.tsv` and per‑genome target XMLs
  - Note: The query set controls how many distinct targets can be discovered; future refinement may use per‑locus cluster members only.

- Phase 5 (Classification)
  - Script: `scripts/05_classify_targets.py`
  - Purpose: Classify targets as syntenic/unplaceable based on proximity/overlap with filtered blocks.
  - Outputs: `05_classified/all_targets_classified.tsv`, `syntenic_targets.tsv`, `unplaceable_targets.tsv`, `dropped_targets.tsv`
  - Important: No capping; proximity defaults tightened (overlap-first); dropped targets reported (reason).

- Phase 6 (Exonerate Extraction)
  - Scripts: `scripts/exonerate_extract.py`, `scripts/extract_with_exonerate.py`
  - Purpose: Extract gene structures with Exonerate (protein2genome) across hits; generate gene/cds sequences and GFF.
  - Output tree: `06_extracted_sequences/<GENOME>/<QUERY>_*`
  - Note: We do not re-run Exonerate in 6b/6c; we reuse these outputs.

- Phase 06c (Block Re‑selection with Exonerate) – NEW
  - Script: `scripts/06c_reselect_blocks_with_exonerate.py`
  - Purpose: Re‑select best block per (locus, genome) using overlap with Exonerate placements (fallback to num_query_matches).
  - Inputs: `02_synteny_blocks/all_synteny_blocks.tsv`, `05_classified_exonerate_refined/refined_targets.tsv`
  - Output: `03_filtered_blocks_exo/synteny_blocks_filtered.tsv`

- Phase 6b (Exonerate Refinement, Strict) – UPDATED
  - Script: `scripts/06b_refine_targets_with_exonerate.py`
  - Purpose: Refine target placements using Exonerate outputs (no re-run). Strict filters:
    - min score % (default 60; tune 70–75)
    - min query coverage (default 0.6; tune 0.8–0.9)
    - require block overlap (nearby optional)
    - dedup: one best per query_id per (genome,locus);
    - cluster by genomic proximity (`--cluster-gap-bp`, default 3000bp) and keep one best per cluster
  - Inputs: `03_filtered_blocks(_exo)/synteny_blocks_filtered.tsv` and `04_target_genes/all_target_loci.tsv` (for query gating)
  - Outputs: `05_classified_exonerate_refined(_exo)/refined_targets.tsv`, `refined_summary.tsv`
  - Important: No capping. Uses Phase 4 query gating per genome.

- Phase 7 (SwissProt Annotation)
  - Script: `scripts/07_swissprot_annotation.py`
  - Purpose: Annotate flanking proteins that are actual block members.
  - Inputs: `03_filtered_blocks/synteny_blocks_filtered.tsv`, `02_synteny_blocks/*/flanking_blast_all.tsv`

- Phase 8 (Matrices)
  - Scripts: `scripts/08a_generate_locus_matrices.py`, `scripts/08b_generate_summary_matrices.py`
  - Purpose: Build locus‑specific and gene‑type summary matrices.

## Orchestration & Batch

- Full array: `scripts/run_batch_all_families.slurm` (env activation inside SLURM)
- Exo reselect + refine (across finished families):
  - `scripts/run_exo_refinement_all.py` + `scripts/run_exo_refinement_all.slurm`

## Current Investigation Notes (Design Knobs)

- BK positive‑control (MC117):
  - 27/27 loci recovered (good) but inflated targets per locus due to large, dispersed region on chr5.
  - 6b strict reduces false positives; combining with 06c helps anchor blocks; Phase 2 chain/segment constraints reduce mega‑blocks.

- Key parameters to tune:
  - 02: `MAX_BLOCK_SPAN_KB`, gap split logic, colinearity, min hits
  - 6b: `--min-score`, `--min-query-cov`, `--cluster-gap-bp`, `--proximity-kb`

## Inputs/Outputs Quick Table

- Phase 1 → `phase12_landmark/locus_definitions.tsv`
- Phase 2 → `02_synteny_blocks/*_synteny_blocks.tsv`, `all_synteny_blocks.tsv`
- Phase 3 → `03_filtered_blocks/synteny_blocks_filtered.tsv`
- Phase 4 → `04_target_genes/all_target_loci.tsv`
- Phase 5 → `05_classified/*_targets.tsv`, `dropped_targets.tsv`
- Phase 6 → `06_extracted_sequences/<GENOME>/<QUERY>_*`
- Phase 06c → `03_filtered_blocks_exo/synteny_blocks_filtered.tsv`
- Phase 6b → `05_classified_exonerate_refined(_exo)/refined_targets.tsv`
- Phase 7 → `07_swissprot_annotations/`
- Phase 8 → `08_matrices/`

## Environment

- Activate inside SLURM scripts only. Baseline: `busco_env` (BLAST, pandas, biopython) and Exonerate module.

