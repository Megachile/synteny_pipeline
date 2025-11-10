Canonical Pipeline (Hybrid Workflow)

This document lists each phase (clean names), inputs/outputs, and parallelization. Legacy numbered scripts/dirs are deprecated and moved to scripts/deprecated.

01 — Locus Discovery & Flanking (BK/LB) + Landmark Comparison
- Script: scripts/01_phase1_combined.py
- Input: gene family + LOC IDs
- Output: outputs/<family>/phase1/locus_definitions.tsv and per‑locus flanking FASTA; synteny_comparisons.tsv
- Notes: deduplicates isoforms and tandem clusters; only BK/LB landmarks; includes BK↔LB flanking comparison

02 — (Included in phase 1)
- Landmark flanking comparison runs inside phase 1 combined. No separate step.

03 — Genome‑wide Synteny Search (per‑locus tBLASTn)
- Script: scripts/synteny_detection.py
- Input: phase1/locus_definitions.tsv, data/ragtag_dbs/*.nhr
- Output: outputs/<family>/phase2_synteny/<locus>/blast_xml/*.xml
- **CRITICAL**: Run at LOCUS level, NOT family level
  - Each locus = 1 SLURM array task (75 genomes × ~50 flanking proteins)
  - Use --locus parameter to specify individual locus
  - See run_phase2_locus_array.slurm for proper setup
- Parallel: SLURM array over individual loci (not families!)
- Threshold: keep best block per genome; exclude only blocks with <3 unique flanking matches

04 — Aggregate Synteny Blocks
- Script: scripts/aggregate_and_filter_blocks.py (combined aggregation + filtering)
- Input: synteny_blocks/*_synteny_blocks.tsv
- Output: filtered_blocks/all_synteny_blocks.tsv (written alongside filtered TSV)

05 — Filter to Best Block per Genome (combined in step 04)
- Script: scripts/aggregate_and_filter_blocks.py
- Input: aggregated all_synteny_blocks
- Output: filtered_blocks/synteny_blocks_filtered.tsv (+ summary)

06 — Extract Target Proteins (BK/LB)
- Script: scripts/03c_extract_target_proteins.py
- Input: phase1/locus_definitions.tsv, data/proteomes, data/genomes (GFF)
- Output: reference_targets/<locus>/<locus>_targets.faa

07 — Target Gene Detection (tBLASTn)
- Script: scripts/04_blast_targets.py
- Input: reference_targets/, data/ragtag_dbs/*.nhr
- Output: target_genes/all_target_loci.tsv, target_genes/combined_targets.faa

08 — Classify Targets (Syntenic vs Unplaceable)
- Script: scripts/05_classify_targets.py
- Input: target_genes/all_target_loci.tsv, filtered_blocks/synteny_blocks_filtered.tsv
- Output: classified/syntenic_targets.tsv, unplaceable_targets.tsv, all_targets_classified.tsv

09 — Exonerate Extraction (Gene Structures)
- Script: scripts/06_extract_sequences.py
- Input: classified/*, target_genes/combined_targets.faa, data/ragtag_output
- Output: extracted_sequences/<genome>/<locus>/*, extracted_sequences/all_extracted_genes.faa

10 — Grade Match Quality
- Script: scripts/07_grade_matches.py
- Input: classified/syntenic_targets.tsv, extracted_sequences/, target_genes/combined_targets.faa
- Output: classified/syntenic_targets_graded.tsv

11 — SwissProt Annotation (per‑locus)
- Script: scripts/07_swissprot_annotation.py
- Input: synteny_blocks/, filtered_blocks/synteny_blocks_filtered.tsv, SwissProt DB
- Output: swissprot_annotations/<locus>_*.tsv, swissprot_annotations/genome_specific_swissprot_annotations.tsv
- Parallel: SLURM array over loci

12 — Matrices (Locus + Summary)
- Scripts: scripts/08a_generate_locus_matrices.py, scripts/08b_generate_summary_matrices.py
- Input: locus_definitions.tsv, filtered_blocks, classified, extracted_sequences, species map, SwissProt
- Output: matrices/*_locus_matrix.tsv, matrices/summary_matrix.tsv

Validation Helpers (TBD additions)
- Phase 1: scripts/summarize_phase1.py → summary_phase1.tsv (flags: missing/zero/low flanking, count mismatches)
- Phase 03: summarize synteny detection (per-locus blocks, zeros, span outliers) — TBD
- Phase 04/05: summarize aggregate/filter (per-locus kept/dropped breakdown) — TBD
- Phase 07: summarize target detection distribution (per-genome counts, zeros) — TBD
- Phase 08: summarize placement (syntenic vs unplaceable per genome/locus) — TBD
- Phase 09/10: summarize Exonerate + grades (status distributions) — TBD
- Phase 11: summarize SwissProt hit rates per locus — TBD

Conventions & Thresholds
- Scaffold normalization: scripts/synteny_utils.py::normalize_scaffold
- Synteny block inclusion: drops only blocks with <3 unique flanking matches (applied at aggregation/filter); otherwise picks best per genome
- Landmarks limited to BK/LB; BRAKER3 TR/DR excluded
