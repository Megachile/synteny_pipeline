#!/usr/bin/env python3
"""
Step 06b: Refine syntenic targets using existing Exonerate outputs (no re-run).

Specificity safeguards:
- Filter exonerate placements by minimum score percent (default 60%).
- Optionally restrict to query_ids present in 04_target_genes/all_target_loci.tsv per genome.
- Require overlap with the selected block (overlap-only by default; --proximity-kb enables nearby).
- Deduplicate: keep one best placement per query_id per (genome, locus),
  then cluster placements by genomic proximity per (genome,locus) and keep one best per cluster.
- No tandem capping: keep all distinct cluster representatives.
"""

from pathlib import Path
import argparse
import re
import pandas as pd
from collections import defaultdict


def norm_scaffold(s: str) -> str:
    if s is None:
        return ''
    s = str(s).strip()
    s = s.replace('RagTag', '').replace('_RagTag', '')
    s = s.split()[0]
    base = s.split(':', 1)[0]
    return base.split('.')[0]


def parse_exonerate_alignment(exo_path: Path):
    data = {
        'query_id': None,
        'target_seqid': None,
        'target_region_start': None,
        'target_region_end': None,
        'target_hit_start': None,
        'target_hit_end': None,
        'strand': None,
        'score_percent': None,
        'query_len': None,
        'q_ab': None,
        'q_ae': None,
    }
    with open(exo_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('Query:') and not data['query_id']:
                m = re.search(r'Query:\s+([^\s]+)', line)
                if m:
                    data['query_id'] = m.group(1)
            if line.startswith('# Query:'):
                # e.g. "# Query: XP_... (221 bp)"
                mlen = re.search(r'\((\d+)\s*bp\)', line)
                if mlen:
                    data['query_len'] = int(mlen.group(1))
            if 'Query range:' in line:
                mm = re.findall(r'(\d+)', line)
                if len(mm) >= 2:
                    data['q_ab'] = int(mm[0]); data['q_ae'] = int(mm[1])
            if line.startswith('Target:') and not data['target_seqid']:
                m = re.search(r'Target:\s+([^\s\[]+)', line)
                if m:
                    data['target_seqid'] = m.group(1)
                sm = re.search(r'strand:([+-])', line)
                if sm:
                    data['strand'] = sm.group(1)
                if data['target_seqid'] and ':' in data['target_seqid']:
                    base, rng = data['target_seqid'].split(':', 1)
                    rm = re.match(r'(\d+)-(\d+)', rng)
                    if rm:
                        data['target_region_start'] = int(rm.group(1))
                        data['target_region_end'] = int(rm.group(2))
            if 'Percent:' in line and data['score_percent'] is None:
                pm = re.search(r'Percent:\s*([0-9.]+)', line)
                if pm:
                    data['score_percent'] = float(pm.group(1))
            if line.startswith('# Target range:') or line.startswith('Target range:'):
                mm = re.findall(r'(\d+)', line)
                if len(mm) >= 2:
                    data['target_hit_start'] = int(mm[0])
                    data['target_hit_end'] = int(mm[1])
    return data


def to_absolute_coords(seqid: str, region_start: int, hit_start: int, hit_end: int):
    start = min(hit_start, hit_end)
    end = max(hit_start, hit_end)
    if seqid and ':' in seqid and region_start is not None:
        abs_start = region_start + start - 1
        abs_end = region_start + end - 1
        base = seqid.split(':', 1)[0]
        return base, abs_start, abs_end
    else:
        return seqid, start, end


def load_blocks(blocks_path: Path):
    df = pd.read_csv(blocks_path, sep='\t')
    df['start'] = df['start'].astype(int)
    df['end'] = df['end'].astype(int)
    df['scf_norm'] = df['scaffold'].astype(str).map(norm_scaffold)
    return df


def overlap(a_start, a_end, b_start, b_end) -> int:
    if b_end < a_start or b_start > a_end:
        return 0
    return max(0, min(a_end, b_end) - max(a_start, b_start) + 1)


def main():
    ap = argparse.ArgumentParser(description='Refine targets using Exonerate placements (no capping, stricter)')
    ap.add_argument('--family-dir', type=Path, required=True)
    ap.add_argument('--blocks', type=Path, required=True)
    ap.add_argument('--output-dir', type=Path, required=True)
    ap.add_argument('--proximity-kb', type=int, default=0, help='Allow nearby (kb) if no overlap (default 0)')
    ap.add_argument('--min-score', type=float, default=60.0, help='Minimum Exonerate score percent (default 60)')
    ap.add_argument('--min-query-cov', type=float, default=0.6, help='Minimum query coverage fraction (default 0.6)')
    ap.add_argument('--cluster-gap-bp', type=int, default=3000, help='Cluster placements within this bp and keep best per cluster (default 3000)')
    ap.add_argument('--targets-file', type=Path, help='Optional: 04_target_genes/all_target_loci.tsv to restrict queries per genome')
    args = ap.parse_args()

    fam_dir = args.family_dir
    blocks_path = args.blocks
    out_dir = args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    print('=' * 80)
    print('PHASE 06b: EXONERATE-REFINED TARGETS (STRICTER)')
    print('=' * 80)
    print('Family:', fam_dir.name)
    print('Blocks:', blocks_path)
    print('Output:', out_dir)
    print('Min score %:', args.min_score)
    print('Min query cov:', args.min_query_cov)
    print('Proximity kb:', args.proximity_kb)
    print('Cluster gap (bp):', args.cluster_gap_bp)

    blocks_df = load_blocks(blocks_path)
    proximity_bp = max(0, args.proximity_kb) * 1000
    min_score = float(args.min_score)
    min_qcov = float(args.min_query_cov)
    cluster_gap = int(args.cluster_gap_bp)

    # Allowed query ids from Phase 4 per genome (optional)
    allowed_queries_by_genome = defaultdict(set)
    tf_path = args.targets_file or (fam_dir / '04_target_genes' / 'all_target_loci.tsv')
    if tf_path.exists():
        try:
            tdf = pd.read_csv(tf_path, sep='\t')
            if 'genome' in tdf.columns and 'query_id' in tdf.columns:
                for g, gdf in tdf.groupby('genome'):
                    allowed_queries_by_genome[str(g)].update(set(gdf['query_id'].astype(str)))
        except Exception:
            pass

    refined = []
    exo_dir = fam_dir / '06_extracted_sequences'
    if not exo_dir.exists():
        print('No exonerate directory:', exo_dir)
    else:
        for genome_dir in sorted(exo_dir.iterdir()):
            if not genome_dir.is_dir():
                continue
            genome = genome_dir.name
            candidates = []
            for sub in genome_dir.iterdir():
                if not sub.is_dir():
                    continue
                for exo_file in sorted(sub.glob('*_exonerate_flank*.txt')):
                    info = parse_exonerate_alignment(exo_file)
                    if not info.get('target_seqid') or info.get('target_hit_start') is None:
                        continue
                    sc = info.get('score_percent') or 0.0
                    if sc < min_score:
                        continue
                    qid = str(info.get('query_id',''))
                    if allowed_queries_by_genome[str(genome)] and qid and (qid not in allowed_queries_by_genome[str(genome)]):
                        continue
                    # Query coverage filter
                    qlen = info.get('query_len')
                    qab = info.get('q_ab'); qae = info.get('q_ae')
                    if qlen and qab and qae:
                        qcov = (abs(qae - qab) + 1) / float(qlen)
                        if qcov < min_qcov:
                            continue
                    base, abs_start, abs_end = to_absolute_coords(info['target_seqid'], info.get('target_region_start'), info.get('target_hit_start'), info.get('target_hit_end'))
                    candidates.append({
                        'genome': genome,
                        'query_id': qid,
                        'scaffold': base,
                        'scf_norm': norm_scaffold(base),
                        'start': int(abs_start),
                        'end': int(abs_end),
                        'strand': info.get('strand','.'),
                        'score_percent': float(sc),
                        'source_file': str(exo_file)
                    })

            if not candidates:
                continue
            cand_df = pd.DataFrame(candidates)
            # Deduplicate multiple flank runs: keep best per (query_id, scaffold, start, end)
            cand_df.sort_values(['query_id','score_percent'], ascending=[True, False], inplace=True)
            cand_df = cand_df.drop_duplicates(subset=['query_id','scaffold','start','end'], keep='first')

            # For each locus in this genome, match candidates to the best block and cluster
            loci = blocks_df[blocks_df['genome']==genome]['locus_id'].unique().tolist()
            for locus_id in loci:
                blk = blocks_df[(blocks_df['genome']==genome) & (blocks_df['locus_id']==locus_id)]
                if blk.empty:
                    continue
                blk = blk.iloc[0]
                sub = cand_df[cand_df['scf_norm']==blk['scf_norm']].copy()
                if sub.empty:
                    continue
                def ovp(row):
                    s = int(row['start']); e = int(row['end'])
                    if not (e < int(blk['start']) or s > int(blk['end'])):
                        return max(0, min(e, int(blk['end'])) - max(s, int(blk['start'])) + 1)
                    if proximity_bp > 0:
                        dist = int(blk['start']) - e if e < int(blk['start']) else s - int(blk['end'])
                        return -dist if dist <= proximity_bp else 0
                    return 0
                sub['overlap_bp'] = sub.apply(ovp, axis=1)
                sub = sub[sub['overlap_bp'] != 0]
                if sub.empty:
                    continue
                # Per query_id: keep best
                best_per_q = []
                for qid, qdf in sub.groupby('query_id'):
                    qdf = qdf.sort_values(by=['overlap_bp','score_percent'], ascending=[False, False])
                    best_per_q.append(qdf.iloc[0])
                if not best_per_q:
                    continue
                sub = pd.DataFrame(best_per_q).sort_values(by=['overlap_bp','score_percent'], ascending=[False, False])
                # Cluster by genomic proximity and keep best per cluster
                clustered = []
                used = []
                for _, r in sub.iterrows():
                    s = int(r['start']); e = int(r['end'])
                    placed = False
                    for i, (cs, ce, idx_best) in enumerate(used):
                        # if overlaps or within gap
                        if not (e < cs - cluster_gap or s > ce + cluster_gap):
                            # already kept best per cluster; since sub sorted, first is best; skip rest
                            placed = True
                            break
                    if not placed:
                        used.append((s, e, len(clustered)))
                        clustered.append(r)
                for _, rr in pd.DataFrame(clustered).iterrows():
                    refined.append({
                        'genome': genome,
                        'locus_id': locus_id,
                        'scaffold': rr['scaffold'],
                        'start': int(rr['start']),
                        'end': int(rr['end']),
                        'strand': rr['strand'],
                        'score_percent': rr['score_percent'],
                        'placement': 'synteny' if rr['overlap_bp']>0 else 'nearby',
                        'source': 'exonerate',
                        'query_id': rr['query_id'],
                        'block_scaffold': blk['scaffold'],
                        'block_start': int(blk['start']),
                        'block_end': int(blk['end']),
                        'source_file': rr['source_file']
                    })

    if not refined:
        print('No refined targets found')
        return

    ref_df = pd.DataFrame(refined)
    out_file = out_dir / 'refined_targets.tsv'
    ref_df.to_csv(out_file, sep='\t', index=False)
    print('Refined targets written:', out_file, f'({len(ref_df)})')

    # Summary
    summary = []
    for (loc, g), grp in ref_df.groupby(['locus_id','genome']):
        summary.append({'locus_id': loc, 'genome': g, 'in_block': (grp['placement']=='synteny').sum(), 'nearby': (grp['placement']=='nearby').sum()})
    sum_df = pd.DataFrame(summary)
    sum_file = out_dir / 'refined_summary.tsv'
    sum_df.to_csv(sum_file, sep='\t', index=False)
    print('Summary written:', sum_file)


if __name__ == '__main__':
    main()
