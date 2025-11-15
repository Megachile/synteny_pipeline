#!/usr/bin/env python3
"""
Phase 3 (redesign): Aggregate per-locus synteny blocks and emit a threshold-only file.

Outputs:
  - all_synteny_blocks.tsv                 (concatenation of per-locus *_synteny_blocks.tsv)
  - all_synteny_blocks_ge{min}.tsv         (blocks with num_query_matches >= min, or fallback count)

Best-per-genome reduction is intentionally skipped here — classification (Phase 5)
will choose the winning block within a locus using target evidence.
"""

from __future__ import annotations

from pathlib import Path
import argparse
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description='Phase 3 (redesign): aggregate synteny blocks, threshold only')
    p.add_argument('--synteny-dir', type=Path, required=True, help='Directory with per-locus *_synteny_blocks.tsv')
    p.add_argument('--output-dir', type=Path, required=True, help='Directory to write aggregated outputs')
    p.add_argument('--min-proteins', type=int, default=3, help='Minimum unique matches per block (default 3)')
    return p.parse_args()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Collect per-locus block files (exclude combined files)
    per_locus = [p for p in args.synteny_dir.glob('*_synteny_blocks.tsv')
                 if p.name not in ('all_synteny_blocks.tsv', 'combined_synteny_blocks.tsv')]
    if not per_locus:
        print(f'ERROR: no per-locus synteny files in {args.synteny_dir}')
        return

    frames = []
    for f in sorted(per_locus):
        try:
            df = pd.read_csv(f, sep='\t')
            if not df.empty:
                frames.append(df)
        except Exception:
            continue

    if not frames:
        print('ERROR: could not read any per-locus synteny TSVs')
        return

    all_df = pd.concat(frames, ignore_index=True)
    all_out = args.output_dir / 'all_synteny_blocks.tsv'
    all_df.to_csv(all_out, sep='\t', index=False)
    print(f'Aggregated: {len(all_df)} blocks -> {all_out}')

    # Choose count column
    count_col = None
    for c in ('num_query_matches', 'num_target_proteins', 'num_proteins'):
        if c in all_df.columns:
            count_col = c
            break
    if count_col is None:
        # Fallback — keep everything since we can't compute counts
        thr_df = all_df.copy()
        print('WARNING: no count column found; skipping thresholding')
    else:
        thr_df = all_df[all_df[count_col] >= args.min_proteins].copy()
        print(f'Threshold ≥{args.min_proteins}: {len(thr_df)} blocks remain (from {len(all_df)})')

    thr_out = args.output_dir / f'all_synteny_blocks_ge{args.min_proteins}.tsv'
    thr_df.to_csv(thr_out, sep='\t', index=False)
    print(f'Wrote threshold-only blocks: {thr_out}')


if __name__ == '__main__':
    main()

