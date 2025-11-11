#!/usr/bin/env python3
"""
Phase 1 (redesign): Locus discovery + locus envelope metadata + dedup flanking preference.

This wraps the existing discovery implementation to produce flanking/targets, then
post-processes locus_definitions.tsv to:
  - Prefer deduplicated flanking files ("*_flanking_dedup.faa")
  - Add expected chromosome from locus_id (e.g., BK_chr2_x → chr2)
  - Add locus envelope coordinates from debug flanking files when available
  - Compute a per-locus scale factor (locus_scale) based on span vs cohort median

Usage:
    python 01_phase1.py \
        --loc-ids LOC117167432,LOC117167433 \
        --gene-family ferritin_MC102 \
        --output-dir outputs/ferritin_MC102/phase1_v2
"""

from __future__ import annotations

from pathlib import Path
import sys
import argparse
import pandas as pd


def _import_legacy_scripts():
    here = Path(__file__).resolve().parent
    scripts_dir = here.parent / 'scripts'
    sys.path.insert(0, str(scripts_dir))
    return scripts_dir


def run_legacy_discovery(loc_ids: str, gene_family: str, out_dir: Path) -> None:
    """Invoke existing phase1_discovery_impl in-process, writing to out_dir.

    This mirrors scripts/01_phase1_combined.py's approach.
    """
    _import_legacy_scripts()
    import phase1_discovery_impl as _disc  # type: ignore
    old_argv = sys.argv[:]
    try:
        sys.argv = [
            _disc.__file__,
            '--loc-ids', loc_ids,
            '--output-dir', str(out_dir),
            '--gene-family', gene_family,
        ]
        if hasattr(_disc, 'main'):
            _disc.main()
        else:
            raise RuntimeError('phase1_discovery_impl has no main()')
    finally:
        sys.argv = old_argv


def _expected_chr_from_locus_id(locus_id: str) -> str:
    # Examples: BK_chr2_a → chr2; LB_scf7864_a → scf7864
    try:
        parts = locus_id.split('_')
        if len(parts) >= 3 and parts[1].startswith('chr'):
            return parts[1]
        if len(parts) >= 3 and parts[1].startswith('scf'):
            return parts[1]
        return ''
    except Exception:
        return ''


def _read_debug_boundaries(phase1_dir: Path) -> pd.DataFrame:
    rows = []
    for dbg in phase1_dir.glob('debug_flanking_*.tsv'):
        try:
            df = pd.read_csv(dbg, sep='\t')
            # Expect columns: boundary_start, boundary_end, locus_id
            if {'boundary_start', 'boundary_end', 'locus_id'} <= set(df.columns):
                rows.append(df[['locus_id', 'boundary_start', 'boundary_end']])
        except Exception:
            continue
    if not rows:
        return pd.DataFrame(columns=['locus_id', 'boundary_start', 'boundary_end'])
    merged = pd.concat(rows, ignore_index=True)
    # Some files may duplicate; keep first per locus_id
    merged = merged.drop_duplicates(subset=['locus_id'], keep='first')
    return merged


def _compute_locus_scale(df: pd.DataFrame) -> pd.Series:
    # Compute per-genome median span (kb), then ratio per locus, clamped to [0.5, 2.0]
    if 'genome' not in df.columns or 'locus_span_kb' not in df.columns:
        return pd.Series([1.0] * len(df))
    med_by_genome = df.groupby('genome')['locus_span_kb'].median().to_dict()
    def scale(row):
        base = med_by_genome.get(row['genome'], max(1.0, df['locus_span_kb'].median()))
        if base <= 0:
            base = 1.0
        r = float(row['locus_span_kb']) / float(base)
        return max(0.5, min(2.0, r))
    return df.apply(scale, axis=1)


def postprocess_phase1(phase1_dir: Path) -> Path:
    locus_defs = phase1_dir / 'locus_definitions.tsv'
    if not locus_defs.exists():
        raise FileNotFoundError(f'Expected locus_definitions.tsv at {locus_defs}')

    df = pd.read_csv(locus_defs, sep='\t')

    # Prefer deduplicated flanking file
    def pick_dedup(path_str: str) -> str:
        p = Path(path_str)
        dedup = p.with_name(f"{p.stem}_dedup{p.suffix}")
        return str(dedup if dedup.exists() else p)

    if 'flanking_file' in df.columns:
        df['flanking_file'] = df['flanking_file'].astype(str).apply(pick_dedup)

    # Add expected chromosome from locus_id
    df['expected_chromosome'] = df['locus_id'].astype(str).apply(_expected_chr_from_locus_id)

    # Join locus envelope boundaries from debug files
    dbg = _read_debug_boundaries(phase1_dir)
    if not dbg.empty:
        df = df.merge(dbg, on='locus_id', how='left')
        df.rename(columns={'boundary_start': 'locus_span_start', 'boundary_end': 'locus_span_end'}, inplace=True)
    else:
        # Fallback to target extents
        if {'target_start', 'target_end'} <= set(df.columns):
            df['locus_span_start'] = df['target_start']
            df['locus_span_end'] = df['target_end']

    # Compute span kb if possible
    if {'locus_span_start', 'locus_span_end'} <= set(df.columns):
        try:
            df['locus_span_kb'] = ((df['locus_span_end'].astype(float) - df['locus_span_start'].astype(float)).clip(lower=0)) / 1000.0
        except Exception:
            df['locus_span_kb'] = 0.0
    else:
        df['locus_span_kb'] = 0.0

    # Compute locus_scale
    df['locus_scale'] = _compute_locus_scale(df)

    # Write updated locus_definitions.tsv (in-place)
    df.to_csv(locus_defs, sep='\t', index=False)
    return locus_defs


def parse_args():
    p = argparse.ArgumentParser(description='Phase 1 (redesign): locus discovery + envelope + scaling')
    p.add_argument('--loc-ids', required=True, help='Comma-separated LOC IDs (e.g., LOC...,LOC...)')
    p.add_argument('--gene-family', required=True, help='Gene family label')
    p.add_argument('--output-dir', required=True, type=Path, help='Output directory for phase1_v2 outputs')
    return p.parse_args()


def main():
    args = parse_args()
    out_dir = args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    print('=' * 80)
    print('PHASE 1 (REDESIGN): DISCOVERY + ENVELOPE + DEDUP PREFERENCE')
    print('=' * 80)
    print(f'Output dir: {out_dir}')

    # 1) Run legacy discovery into the requested directory
    print('\n[1] Running legacy discovery...')
    run_legacy_discovery(args.loc_ids, args.gene_family, out_dir)

    # 2) Post-process locus definitions
    print('\n[2] Post-processing locus definitions...')
    locus_defs = postprocess_phase1(out_dir)
    print(f'  Updated: {locus_defs}')

    # 3) Summary
    print('\n[3] Summary (preview):')
    try:
        df = pd.read_csv(locus_defs, sep='\t')
        cols = [c for c in ['locus_id','genome','expected_chromosome','locus_span_kb','locus_scale'] if c in df.columns]
        if cols:
            print(df[cols].to_string(index=False))
    except Exception:
        pass

    print('\nDone.')


if __name__ == '__main__':
    main()

