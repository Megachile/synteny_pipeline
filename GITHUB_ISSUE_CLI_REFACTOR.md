# GitHub Issue: Refactor pipeline scripts to use explicit CLI arguments

**Title**: Refactor pipeline scripts to use explicit CLI arguments

**Labels**: enhancement, refactor

**Body**:

## Problem

Currently, all scripts use `config.py` with environment variables to determine input/output paths. This makes it hard to:
- Test scripts on specific output directories (e.g., ferritin vs mc6)
- Run multiple analyses in parallel
- Understand data flow between phases
- Debug failures (paths are implicit)

## Solution

Phase 8a has been refactored (commit a335926) to use explicit command-line arguments:

```bash
python 08a_generate_locus_matrices.py \
    --locus-defs outputs/ferritin_phase3_real/locus_definitions.tsv \
    --synteny-dir outputs/ferritin_phase3_real/02_synteny_blocks \
    --blocks outputs/ferritin_phase3_real/03_filtered_blocks/synteny_blocks_filtered.tsv \
    --targets outputs/ferritin_phase3_real/05_classified/all_targets_classified.tsv \
    --swissprot outputs/ferritin_phase3_real/07_swissprot/genome_specific_swissprot_annotations.tsv \
    --reference-proteins data/reference/protein.faa \
    --species-map data/gca_to_species.tsv \
    --output-dir outputs/ferritin_phase3_real/07_matrices/locus_specific
```

## Benefits

- **Testable**: Can run on any output directory without changing config
- **Composable**: Clear data flow between phases
- **Debuggable**: Explicit paths visible in command
- **Standalone**: No environment variables needed
- **Parallel**: Can run multiple analyses simultaneously

## Scripts to Refactor

- [ ] 01_discover_paralogs_gene_level.py
- [ ] 02_synteny_detection.py
- [ ] 03_filter_blocks.py
- [ ] 04_blast_targets.py
- [ ] 05_classify_targets.py
- [ ] 06_extract_sequences.py
- [ ] 07_swissprot_annotation.py
- [x] 08a_generate_locus_matrices.py ✓
- [ ] 08b_generate_summary_matrices.py

## Template

See `scripts/08a_generate_locus_matrices.py` for the pattern:
- Use `argparse` with `Path` types
- Make all inputs/outputs explicit arguments
- Keep `config.py` only for database paths (SwissProt, RagTag, etc)
- Print input summary at start for transparency

## Related

- Beads: hybrid_workflow-25
- Fixes issue from hybrid_workflow-24 (column headers)

## Instructions

To create this GitHub issue, run:
```bash
cd /carc/scratch/projects/emartins/2016456/adam/synteny_scanning/hpc_deployment_package/hybrid_workflow
# Install gh CLI if needed, then:
gh issue create --title "Refactor pipeline scripts to use explicit CLI arguments" --body-file GITHUB_ISSUE_CLI_REFACTOR.md --label enhancement,refactor
```

Or create manually at: https://github.com/Megachile/synteny_scanning/issues/new
