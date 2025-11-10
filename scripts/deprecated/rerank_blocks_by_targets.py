#!/usr/bin/env python3
"""
Rerank synteny blocks per (locus, genome) using overlap with detected targets.

Inputs (within a family root dir):
  - 02_synteny_blocks/<locus>_synteny_blocks.tsv (candidates)
  - 04_target_genes/all_target_loci.tsv

Outputs:
  - 03_filtered_blocks_reranked/synteny_blocks_filtered.tsv
  - 03_filtered_blocks_reranked/<locus>_synteny_blocks_filtered.tsv

Notes:
  - No BLAST is run; this only uses existing outputs.
  - Score for each candidate block is primary: total bp overlap with any targets
    from the same genome on the same scaffold; tie-breakers: num_query_matches desc,
    then smaller span_kb (ascending).
"""

from pathlib import Path
import argparse
import pandas as pd


def _load_candidates(fam_dir: Path) -> pd.DataFrame:
    cand_rows = []
    blocks_dir = fam_dir / '02_synteny_blocks'
    if not blocks_dir.exists():
        return pd.DataFrame()

    for tsv in blocks_dir.glob('*_synteny_blocks.tsv'):
        if tsv.name == 'all_synteny_blocks.tsv':
            # Prefer per-locus split files to preserve locus_id
            continue
        try:
            df = pd.read_csv(tsv, sep='\t')
            if len(df) == 0:
                continue
            # Ensure required columns
            exp = {'locus_id', 'genome', 'scaffold', 'start', 'end'}
            if not exp.issubset(set(df.columns)):
                continue
            cand_rows.append(df)
        except Exception:
            continue

    if cand_rows:
        return pd.concat(cand_rows, ignore_index=True)
    return pd.DataFrame()


def _load_targets(fam_dir: Path) -> pd.DataFrame:
    tsv = fam_dir / '04_target_genes' / 'all_target_loci.tsv'
    if not tsv.exists():
        return pd.DataFrame()
    try:
        df = pd.read_csv(tsv, sep='\t')
        exp = {'genome', 'scaffold', 'start', 'end', 'gene_family', 'locus_name'}
        if not exp.issubset(set(df.columns)):
            return pd.DataFrame()
        # Standardize numeric types
        df['start'] = df['start'].astype(int)
        df['end'] = df['end'].astype(int)
        return df
    except Exception:
        return pd.DataFrame()


def _interval_overlap(a_start, a_end, b_start, b_end) -> int:
    s = max(a_start, b_start)
    e = min(a_end, b_end)
    return max(0, e - s + 1)


def rerank_family(fam_dir: Path) -> Path:
    fam_dir = fam_dir.resolve()
    candidates = _load_candidates(fam_dir)
    targets = _load_targets(fam_dir)

    out_dir = fam_dir / '03_filtered_blocks_reranked'
    out_dir.mkdir(exist_ok=True, parents=True)

    if candidates.empty or targets.empty:
        # Create an empty file for consistent downstream
        empty = out_dir / 'synteny_blocks_filtered.tsv'
        if not empty.exists():
            pd.DataFrame(columns=['locus_id','genome','block_id','scaffold','strand','start','end','span_kb','num_target_proteins','num_query_matches']).to_csv(empty, sep='\t', index=False)
        return empty

    # Normalize types
    for col in ('start', 'end'):
        candidates[col] = candidates[col].astype(int)
    if 'span_kb' not in candidates.columns:
        candidates['span_kb'] = (candidates['end'] - candidates['start']) / 1000.0

    # Build quick index for targets by (genome, scaffold)
    t_index = {}
    for (g, s), df in targets.groupby(['genome', 'scaffold']):
        # Keep as list of intervals to loop (targets per scaffold usually modest)
        t_index[(g, s)] = df[['start', 'end']].values.tolist()

    rows = []
    # Score per (locus, genome, block_id)
    for (locus_id, genome), group in candidates.groupby(['locus_id', 'genome']):
        best = None
        best_key = None
        for _, row in group.iterrows():
            intervals = t_index.get((genome, row['scaffold']), [])
            # Total overlap with any target interval on same scaffold
            ov = 0
            for ts, te in intervals:
                ov += _interval_overlap(int(row['start']), int(row['end']), int(ts), int(te))

            # Tie-breakers: overlap desc, num_query_matches desc, span_kb asc
            # Missing num_query_matches -> treat as 0
            nq = row.get('num_query_matches', 0)
            key = (ov, int(nq), -float(row['span_kb']))
            if (best_key is None) or (key > best_key):
                best_key = key
                best = row

        if best is not None:
            rows.append(best)

    if not rows:
        empty = out_dir / 'synteny_blocks_filtered.tsv'
        if not empty.exists():
            pd.DataFrame(columns=['locus_id','genome','block_id','scaffold','strand','start','end','span_kb','num_target_proteins','num_query_matches']).to_csv(empty, sep='\t', index=False)
        return empty

    filtered = pd.DataFrame(rows)
    # Save combined
    combined = out_dir / 'synteny_blocks_filtered.tsv'
    filtered.to_csv(combined, sep='\t', index=False)

    # Save per-locus
    for locus_id, df in filtered.groupby('locus_id'):
        (out_dir / f'{locus_id}_synteny_blocks_filtered.tsv').write_text(df.to_csv(sep='\t', index=False))

    return combined


def find_families(outputs_root: Path):
    families = []
    for p in outputs_root.iterdir():
        if not p.is_dir():
            continue
        if (p / '04_target_genes' / 'all_target_loci.tsv').exists() and (p / '02_synteny_blocks').exists():
            families.append(p)
    return sorted(families)


def main():
    ap = argparse.ArgumentParser(description='Rerank synteny blocks by overlap with targets (no BLAST).')
    ap.add_argument('--family-dir', type=Path, help='Path to a single family root (outputs/<FAM>)')
    ap.add_argument('--outputs-root', type=Path, default=Path('outputs'), help='Root outputs dir containing families')
    args = ap.parse_args()

    if args.family_dir:
        fams = [args.family_dir]
    else:
        fams = find_families(args.outputs_root)

    print(f'Found {len(fams)} families to rerank')
    results = []
    for fam in fams:
        print(f'Reranking: {fam.name}')
        out = rerank_family(fam)
        results.append({'family': fam.name, 'output': str(out)})

    if results:
        if len(fams) == 1:
            out = fams[0] / 'rerank_summary.tsv'
            pd.DataFrame(results).to_csv(out, sep='\t', index=False)
            print(f'Summary: {out}')
        else:
            pd.DataFrame(results).to_csv('rerank_summary.tsv', sep='\t', index=False)
            print('Summary: rerank_summary.tsv')


if __name__ == '__main__':
    main()
