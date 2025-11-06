#!/usr/bin/env python3
"""
Phase 06c: Reselect best synteny blocks using Exonerate placements.

Uses per-(locus, genome) candidate blocks from Phase 2 and refines the
selection using Exonerate gene placements (from Phase 6b). Falls back to
num_query_matches when no Exonerate placement exists.

Inputs (per family):
  - 02_synteny_blocks/all_synteny_blocks.tsv
  - 05_classified_exonerate_refined/refined_targets.tsv

Output:
  - 03_filtered_blocks_exo/synteny_blocks_filtered.tsv
"""

from pathlib import Path
import argparse
import pandas as pd


def norm_scaffold(s: str) -> str:
    if pd.isna(s) or s is None:
        return ''
    s = str(s).strip()
    s = s.replace('RagTag', '').replace('_RagTag', '')
    s = s.split()[0]
    base = s.split(':', 1)[0]
    return base.split('.')[0]


def load_candidates(cand_path: Path) -> pd.DataFrame:
    df = pd.read_csv(cand_path, sep='\t')
    for col in ['start', 'end']:
        if col in df.columns:
            df[col] = df[col].astype(int)
    if 'num_query_matches' not in df.columns:
        if 'num_target_proteins' in df.columns:
            df['num_query_matches'] = df['num_target_proteins']
        else:
            df['num_query_matches'] = 0
    df['scf_norm'] = df['scaffold'].map(norm_scaffold)
    return df


def load_exo(exo_path: Path) -> pd.DataFrame:
    if not exo_path.exists():
        return pd.DataFrame(columns=['genome','locus_id','scaffold','start','end','placement','scf_norm'])
    df = pd.read_csv(exo_path, sep='\t')
    for col in ['start','end']:
        if col in df.columns:
            df[col] = df[col].astype(int)
    df['scf_norm'] = df['scaffold'].map(norm_scaffold)
    keep = ['genome','locus_id','scf_norm','start','end','placement']
    for k in list(df.columns):
        if k not in keep:
            df.drop(columns=k, inplace=True)
    return df


def overlap(a_start, a_end, b_start, b_end) -> int:
    if b_end < a_start or b_start > a_end:
        return 0
    return max(0, min(a_end, b_end) - max(a_start, b_start) + 1)


def select_blocks(cands: pd.DataFrame, exo_df: pd.DataFrame) -> pd.DataFrame:
    out_rows = []
    for (locus, genome), gdf in cands.groupby(['locus_id','genome']):
        ehits = exo_df[(exo_df['locus_id']==locus) & (exo_df['genome']==genome)]
        inblock = ehits[ehits['placement']=='synteny']
        use = inblock if len(inblock)>0 else ehits

        best_idx = None
        best_key = None
        for idx, row in gdf.iterrows():
            score = 0
            if len(use)>0:
                u = use[use['scf_norm']==row['scf_norm']]
                for _, r in u.iterrows():
                    score += overlap(int(row['start']), int(row['end']), int(r['start']), int(r['end']))
            q = int(row.get('num_query_matches', 0))
            span = int(row['end']) - int(row['start'])
            key = (score, q, -span)
            if best_idx is None or key > best_key:
                best_idx = idx
                best_key = key
        out_rows.append(gdf.loc[best_idx].to_dict())
    return pd.DataFrame(out_rows)


def main():
    ap = argparse.ArgumentParser(description='Phase 06c: Reselect blocks using Exonerate placements')
    ap.add_argument('--family-dir', type=Path, required=True)
    ap.add_argument('--candidates', type=Path)
    ap.add_argument('--exonerate-refined', type=Path)
    ap.add_argument('--output-dir', type=Path)
    args = ap.parse_args()

    fam_dir = args.family_dir
    cand_path = args.candidates or (fam_dir/'02_synteny_blocks'/'all_synteny_blocks.tsv')
    exo_path = args.exonerate_refined or (fam_dir/'05_classified_exonerate_refined'/'refined_targets.tsv')
    out_dir = args.output_dir or (fam_dir/'03_filtered_blocks_exo')
    out_dir.mkdir(parents=True, exist_ok=True)

    print('='*80)
    print('PHASE 06c: BLOCK RE-SELECTION USING EXONERATE')
    print('='*80)
    print('Family:', fam_dir.name)
    print('Candidates:', cand_path)
    print('Exonerate refined:', exo_path)
    print('Output dir:', out_dir)

    if not cand_path.exists():
        raise SystemExit(f'Candidates not found: {cand_path}')

    cand_df = load_candidates(cand_path)
    exo_df = load_exo(exo_path)
    sel = select_blocks(cand_df, exo_df)
    out = out_dir/'synteny_blocks_filtered.tsv'
    sel.to_csv(out, sep='\t', index=False)
    print('Wrote', out, 'rows:', len(sel))

if __name__ == '__main__':
    main()

