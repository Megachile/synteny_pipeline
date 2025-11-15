#!/usr/bin/env python3
"""
Phase 8b (redesign wrapper): Generate gene-type summary matrices using existing generator.

Resolves redesigned outputs and calls scripts/generate_summary_matrices.py.
"""

from __future__ import annotations

from pathlib import Path
import argparse
import subprocess
import sys


def parse_args():
    p = argparse.ArgumentParser(description='Phase 8b (redesign): summary matrices wrapper')
    p.add_argument('--family', required=True, help='Gene family name (folder under outputs)')
    p.add_argument('--species-map', type=Path, default=Path('data/gca_to_species.tsv'))
    p.add_argument('--out-dir', type=Path, help='Override output directory (default outputs/<family>/phase8b_v2)')
    return p.parse_args()


def main():
    args = parse_args()

    base = Path('outputs') / args.family
    phase1 = base / 'phase1_v2'
    phase3 = base / 'phase3_v2'
    phase5 = base / 'phase5_v2'
    phase6 = base / 'phase6_extracted_v2'
    # Write Phase 8b outputs into the same unified Phase 8 folder by default
    out_dir = args.out_dir if args.out_dir else (base / 'phase8_v2')
    out_dir.mkdir(parents=True, exist_ok=True)

    locus_defs = phase1 / 'locus_definitions.tsv'
    blocks = phase5 / 'selected_synteny_blocks.tsv'
    if not blocks.exists():
        cand = phase3 / 'synteny_blocks_filtered.tsv'
        blocks = cand if cand.exists() else blocks
    targets = phase5 / 'all_targets_classified.tsv'

    if not locus_defs.exists():
        print(f'ERROR: Missing locus_definitions.tsv: {locus_defs}', file=sys.stderr)
        sys.exit(1)
    if not targets.exists():
        print(f'ERROR: Missing all_targets_classified.tsv: {targets}', file=sys.stderr)
        sys.exit(1)

    cmd = [
        sys.executable,
        str(Path(__file__).resolve().parent / '08b_generate_summary_matrices_core.py'),
        '--locus-defs', str(locus_defs),
        '--blocks', str(blocks),
        '--targets', str(targets),
        '--species-map', str(args.species_map),
        '--extracted-seqs', str(phase6),
        '--output-dir', str(out_dir),
    ]
    # Tell generator where to find locus matrices
    # Prefer the same output directory (unified Phase 8 folder),
    # but fall back to legacy phase8a_v2 if needed.
    lm_dir = out_dir
    if not lm_dir.exists():
        legacy = base / 'phase8a_v2'
        lm_dir = legacy if legacy.exists() else lm_dir
    if lm_dir.exists():
        cmd.extend(['--locus-matrices-dir', str(lm_dir)])

    print('Running:', ' '.join(cmd))
    subprocess.run(cmd, check=True)
    print('Phase 8b complete:', out_dir)


if __name__ == '__main__':
    main()
