# Hybrid Workflow Redesign Scripts (Pilot)

This folder hosts new versions of the pipeline scripts implementing the redesigned flow described in `../HYBRID_FLOW_REDESIGN.md`.

We will implement the phases incrementally. For now, Phase 1 is implemented here as `01_phase1.py`.

## Phases (planned)

- 01_phase1.py — Locus discovery with deduplicated flanking preference, locus envelope, `locus_scale`, and expected chromosome metadata.
- 02_synteny_detection.py — Synteny detection with scaled gap/span, keep-all blocks, merge-adjacent, chromosome annotation.
- 03_aggregate_threshold_only.py — Aggregate per-locus blocks and emit threshold-only file (no best-per-genome reduction).
- 04_detect_targets.py — Target detection with scaled split-gaps and envelope-aware prioritization.
- 05_classify_targets.py — Locus-first classification; choose block within locus using target evidence; evidence-of-absence rules.

We’ll add the subsequent scripts after validating Phase 1 on MC6 and MC117.

## Running Phase 1 (pilot)

Example:

```
python redesign_scripts/01_phase1.py \
  --loc-ids LOC12345,LOC67890 \
  --gene-family MC6 \
  --output-dir outputs/MC6/phase1_v2
```

This internally runs the existing discovery implementation, then post-processes outputs to:
- Prefer deduplicated flanking in `locus_definitions.tsv`.
- Add `expected_chromosome`, `locus_span_*`, and `locus_scale` columns.

Outputs are written to the directory you pass via `--output-dir`.
