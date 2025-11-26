#!/usr/bin/env python3
"""
Phase 5 (redesign): Locus-first classification using threshold-only blocks.

Inputs:
  --targets  (from Phase 4_v2/all_target_loci.tsv)
  --blocks   (from Phase 3_v2/all_synteny_blocks_ge{min}.tsv)
  --locus-defs (optional) to read locus_scale (scales segment merge gap and dynamic nearby)

Outputs:
  - all_targets_classified.tsv
  - syntenic_targets.tsv
  - unplaceable_targets.tsv
"""

from __future__ import annotations

from pathlib import Path
import argparse
import pandas as pd
from collections import defaultdict

# -----------------------------------------------------------------------------
# Helper: compute per-target query-locus exclusivity
# -----------------------------------------------------------------------------
def compute_query_exclusivity(targets_df: pd.DataFrame, overlap_only: bool = True, gap_kb: int = 0) -> pd.DataFrame:
    """
    For each target row, determine the set of parent_locus values from ALL Phase 4
    loci that overlap (or are within gap_kb) on the same (genome, scaffold).

    Adds columns:
      - query_hit_loci: comma-joined sorted set of loci with hits at that spot
      - query_hit_count: size of that set
      - query_exclusive: True if query_hit_count == 1, else False
    """
    by_scaf = defaultdict(list)
    for i, t in targets_df.iterrows():
        try:
            g = str(t['genome']); s = str(t['scaffold'])
            a = int(t['start']); b = int(t['end'])
            loc = str(t.get('parent_locus', t.get('locus_id', '')))
        except Exception:
            continue
        by_scaf[(g, s)].append((a, b, loc))

    gap_bp = int(gap_kb) * 1000
    hit_loci = []
    for _, t in targets_df.iterrows():
        try:
            g = str(t['genome']); s = str(t['scaffold'])
            a = int(t['start']); b = int(t['end'])
        except Exception:
            hit_loci.append(([], 0, False))
            continue
        loci = set()
        for (x1, x2, loc) in by_scaf.get((g, s), []):
            # overlap or within gap_bp
            if not overlap_only:
                # distance between [a,b] and [x1,x2]
                if b < x1:
                    d = x1 - b
                elif a > x2:
                    d = a - x2
                else:
                    d = 0
                if d <= gap_bp:
                    loci.add(loc)
            else:
                if not (b < x1 or a > x2):
                    loci.add(loc)
        loci_list = sorted(loci)
        hit_loci.append((loci_list, len(loci_list), len(loci_list) == 1))

    targets_df = targets_df.copy()
    targets_df['query_hit_loci'] = [''] * len(targets_df)
    targets_df['query_hit_count'] = 0
    targets_df['query_exclusive'] = False
    for i, (lst, cnt, ex) in enumerate(hit_loci):
        targets_df.at[i, 'query_hit_loci'] = ','.join(lst)
        targets_df.at[i, 'query_hit_count'] = cnt
        targets_df.at[i, 'query_exclusive'] = ex
    return targets_df


def build_locus_segment_index(blocks_df: pd.DataFrame, merge_gap_kb: int) -> dict:
    idx = {}
    pad_bp = merge_gap_kb * 1000
    for _, b in blocks_df.iterrows():
        key = (b['genome'], b['scaffold'])
        li = idx.setdefault(key, {})
        loc = b['locus_id']
        bucket = li.setdefault(loc, {'segments': [], 'blocks': []})
        bucket['blocks'].append({'start': int(b['start']), 'end': int(b['end'])})
    # merge contiguous blocks per locus
    for key, loci in idx.items():
        for loc, data in loci.items():
            arr = sorted(data['blocks'], key=lambda x: x['start'])
            segs = []
            cur = None
            for bl in arr:
                if cur is None:
                    cur = {'start': bl['start'], 'end': bl['end']}
                    continue
                gap = bl['start'] - cur['end'] if bl['start'] > cur['end'] else 0
                if gap <= pad_bp:
                    cur['end'] = max(cur['end'], bl['end'])
                else:
                    segs.append((cur['start'], cur['end']))
                    cur = {'start': bl['start'], 'end': bl['end']}
            if cur is not None:
                segs.append((cur['start'], cur['end']))
            data['segments'] = segs
    return idx


def classify_targets(targets_df: pd.DataFrame, blocks_df: pd.DataFrame,
                     locus_defs: pd.DataFrame | None,
                     pad_kb: int, dynamic_nearby: bool,
                     nearby_frac_span: float, nearby_upper_kb: int,
                     base_segment_gap_kb: int,
                     rescue_same_scaffold: bool = False,
                     rescue_upper_kb: int = 350,
                     rescue_min_evalue: float = 1e-50) -> pd.DataFrame:
    # Optional per-locus scaling
    scale_map = {}
    if locus_defs is not None and 'locus_id' in locus_defs.columns and 'locus_scale' in locus_defs.columns:
        scale_map = {str(r['locus_id']): float(r['locus_scale']) for _, r in locus_defs.iterrows()}

    # Build per-locus segment envelopes by (genome, scaffold)
    # Merge gap scaled by locus_scale
    def seg_gap_kb_for(locus_id: str) -> int:
        s = float(scale_map.get(locus_id, 1.0))
        return max(20, min(300, int(round(base_segment_gap_kb * s))))

    # Index blocks by scaffold for direct selection later
    by_scaf = defaultdict(list)
    for _, b in blocks_df.iterrows():
        by_scaf[(b['genome'], b['scaffold'])].append(b.to_dict())

    # Build locus segment index (per-locus scaling)
    # We need one index per merge gap; to keep it simple, build a base index and use it directly.
    # For finer scaling, we could rebuild per-locus, but this is acceptable for initial redesign.
    locus_index = build_locus_segment_index(blocks_df, merge_gap_kb=base_segment_gap_kb)

    rows = []
    pad_bp = pad_kb * 1000
    for _, t in targets_df.iterrows():
        genome = t['genome']; scf = t['scaffold']; s = int(t['start']); e = int(t['end'])
        locus_id = str(t.get('parent_locus', t.get('locus_id', '')))
        seg_gap_kb = seg_gap_kb_for(locus_id)

        # Candidate blocks selection with query exclusivity logic:
        # - If this target's BLAST support is exclusive to a single parent_locus
        #   (only one locus produced hits at this spot), restrict to that locus.
        # - If non-exclusive (conflicting across multiple loci), consider all blocks
        #   on this scaffold to avoid biasing by query locus.
        exclusive = bool(t.get('query_exclusive', False))
        if exclusive:
            candidate_blocks = [b for b in by_scaf.get((genome, scf), []) if str(b.get('locus_id')) == locus_id]
            if not candidate_blocks:
                candidate_blocks = by_scaf.get((genome, scf), [])
        else:
            candidate_blocks = by_scaf.get((genome, scf), [])

        chosen = None
        best_overlap = -1
        best_dist = None

        # Precompute envelope segments for this locus on this scaffold
        env = None
        if (genome, scf) in locus_index and locus_id in locus_index[(genome, scf)]:
            env = [(a, b) for (a, b) in locus_index[(genome, scf)][locus_id]['segments']]

        # Choose block by overlap, then by nearest distance within dynamic window
        for b in candidate_blocks:
            bs, be = int(b['start']), int(b['end'])
            ov = max(0, min(e, be) - max(s, bs) + 1)
            if ov > 0:
                if ov > best_overlap:
                    chosen = b; best_overlap = ov; best_dist = 0
                else:
                    # tie-breaker by center distance
                    t_c = (s + e) // 2
                    b_c = (bs + be) // 2
                    ch_c = ((int(chosen['start']) + int(chosen['end'])) // 2) if chosen else None
                    if best_overlap == ov and ch_c is not None and abs(t_c - b_c) < abs(t_c - ch_c):
                        chosen = b; best_dist = 0
                continue
            # proximity only if enabled
            if dynamic_nearby or pad_kb > 0:
                if e < bs:
                    d = bs - e
                else:
                    d = s - be
                # dynamic threshold
                span_kb = float(b.get('span_kb', (be - bs) / 1000.0))
                thr_kb = nearby_frac_span * span_kb if dynamic_nearby else 0
                thr_kb = max(thr_kb, 0)
                thr_kb = min(max(thr_kb, 30), nearby_upper_kb) if dynamic_nearby else nearby_upper_kb
                thr_bp = thr_kb * 1000 + pad_bp
                if d <= thr_bp:
                    if best_dist is None or d < best_dist:
                        chosen = b; best_dist = d

        if chosen is not None:
            placement = 'synteny'
            assigned_to = chosen['locus_id']
            chosen_block_id = chosen.get('block_id', '')
            chosen_block_scaf = chosen.get('scaffold', '')
            chosen_block_start = chosen.get('start', '')
            chosen_block_end = chosen.get('end', '')
        else:
            placement = 'unplaceable'
            assigned_to = f"{t.get('locus_id', t.get('parent_locus','unknown'))}_unplaceable"
            chosen_block_id = ''
            chosen_block_scaf = ''
            chosen_block_start = ''
            chosen_block_end = ''

        # Rescue pass: if unplaceable but on same scaffold as parent locus blocks and within rescue distance,
        # assign to parent locus and mark rescue flag.
        rescue_flag = False
        rescue_dist = ''
        if placement == 'unplaceable' and rescue_same_scaffold:
            parent_loc = locus_id
            parent_blocks = [b for b in by_scaf.get((genome, scf), []) if str(b.get('locus_id')) == parent_loc]
            if parent_blocks:
                # distance to locus envelope segments if available; otherwise to closest block
                def interval_distance(a1, a2, b1, b2):
                    if a2 < b1:
                        return b1 - a2
                    if b2 < a1:
                        return a1 - b2
                    return 0
                # compute distance to nearest block
                best_b = None; best_d = None
                for b in parent_blocks:
                    bs, be = int(b['start']), int(b['end'])
                    d = interval_distance(s, e, bs, be)
                    if best_d is None or d < best_d:
                        best_d = d; best_b = b
                try:
                    ev = float(t.get('best_evalue', 'inf'))
                except Exception:
                    ev = float('inf')
                if best_d is not None and best_d <= rescue_upper_kb * 1000 and ev <= rescue_min_evalue:
                    placement = 'synteny'
                    assigned_to = parent_loc
                    chosen_block_id = best_b.get('block_id', '')
                    chosen_block_scaf = best_b.get('scaffold', '')
                    chosen_block_start = best_b.get('start', '')
                    chosen_block_end = best_b.get('end', '')
                    rescue_flag = True
                    rescue_dist = best_d

        out = t.to_dict()
        out['placement'] = placement
        out['assigned_to'] = assigned_to
        out['assigned_block_id'] = chosen_block_id
        out['assigned_block_scaffold'] = chosen_block_scaf
        out['assigned_block_start'] = chosen_block_start
        out['assigned_block_end'] = chosen_block_end
        out['rescue'] = rescue_flag
        out['rescue_distance_bp'] = rescue_dist
        rows.append(out)

    return pd.DataFrame(rows)


def parse_args():
    p = argparse.ArgumentParser(description='Phase 5 (redesign): locus-first classification')
    p.add_argument('--targets', type=Path, required=True, help='Phase 4_v2 all_target_loci.tsv')
    p.add_argument('--blocks', type=Path, required=True, help='Phase 3_v2 all_synteny_blocks_ge{min}.tsv')
    p.add_argument('--output-dir', type=Path, required=True)
    p.add_argument('--locus-defs', type=Path, help='Optional Phase1_v2 locus_definitions.tsv (for locus_scale)')
    p.add_argument('--pad-kb', type=int, default=20)
    p.add_argument('--dynamic-nearby', action='store_true')
    p.add_argument('--nearby-frac-span', type=float, default=0.3)
    p.add_argument('--nearby-upper-kb', type=int, default=200)
    p.add_argument('--segment-gap-kb', type=int, default=50)
    p.add_argument('--rescue-same-scaffold', action='store_true', help='Rescue unplaceables near parent locus on same scaffold')
    p.add_argument('--rescue-upper-kb', type=int, default=350, help='Max rescue distance (kb)')
    p.add_argument('--rescue-min-evalue', type=float, default=1e-50, help='Min strength (evalue <=) to allow rescue')
    return p.parse_args()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    targets_df = pd.read_csv(args.targets, sep='\t')
    blocks_df = pd.read_csv(args.blocks, sep='\t')
    locus_defs = pd.read_csv(args.locus_defs, sep='\t') if args.locus_defs and Path(args.locus_defs).exists() else None

    # Compute per-target query exclusivity across ALL loci on a scaffold.
    # If exclusive (only one parent_locus hits a spot), we can safely restrict
    # candidate blocks to that locus. If not exclusive (conflict), we avoid
    # using the query locus to bias assignment, and consider all locus blocks
    # on the scaffold.
    targets_df = compute_query_exclusivity(targets_df, overlap_only=True, gap_kb=0)

    classified = classify_targets(targets_df, blocks_df, locus_defs,
                                  pad_kb=args.pad_kb,
                                  dynamic_nearby=args.dynamic_nearby,
                                  nearby_frac_span=args.nearby_frac_span,
                                  nearby_upper_kb=args.nearby_upper_kb,
                                  base_segment_gap_kb=args.segment_gap_kb,
                                  rescue_same_scaffold=args.rescue_same_scaffold,
                                  rescue_upper_kb=args.rescue_upper_kb,
                                  rescue_min_evalue=args.rescue_min_evalue)

    # Add locus_name column for compatibility with legacy matrix scripts
    classified['locus_name'] = classified['locus_id']

    # Add gene_family column for compatibility with legacy matrix scripts
    if locus_defs is not None and 'gene_family' in locus_defs.columns:
        family_map = {str(r['locus_id']): str(r['gene_family']) for _, r in locus_defs.iterrows()}
        classified['gene_family'] = classified['locus_id'].astype(str).map(family_map)
    else:
        # Fallback: infer from output directory path or use 'unknown'
        classified['gene_family'] = 'unknown'

    all_out = args.output_dir / 'all_targets_classified.tsv'
    classified.to_csv(all_out, sep='\t', index=False)

    syn = classified[classified['placement'] == 'synteny']
    unp = classified[classified['placement'] == 'unplaceable']
    syn.to_csv(args.output_dir / 'syntenic_targets.tsv', sep='\t', index=False)
    unp.to_csv(args.output_dir / 'unplaceable_targets.tsv', sep='\t', index=False)

    print(f'Classified: {len(classified)} | synteny={len(syn)} | unplaceable={len(unp)}')
    print(f'Wrote: {all_out}')

    # Derive selected synteny blocks per (locus, genome) from assigned targets
    # Pick the block_id with the most target assignments; ties broken by larger span_kb in blocks_df
    # IMPORTANT: Use 'assigned_to' (where target was placed) not 'locus_id' (parent_locus) for block lookup,
    # because assigned_block_id belongs to the assigned_to locus, not necessarily the parent_locus.
    sel_rows = []
    if not syn.empty and {'assigned_block_id','locus_id','genome','assigned_to'} <= set(syn.columns):
        # Count assignments per block - group by assigned_to (where target was placed)
        counts = syn.groupby(['assigned_to','genome','assigned_block_id']).size().reset_index(name='n_targets')
        counts = counts.rename(columns={'assigned_to': 'locus_id'})  # Rename for consistency
        # For each (locus, genome), choose the best
        for (locus_id, genome), sub in counts.groupby(['locus_id','genome']):
            # Ignore rows with empty block id
            sub = sub[sub['assigned_block_id'].astype(str) != '']
            if sub.empty:
                continue
            # Join span_kb from blocks_df for tie-breaker
            bcol = ['locus_id','genome','block_id','scaffold','strand','start','end','span_kb','num_query_matches','num_target_proteins']
            bdf = blocks_df[bcol].copy() if set(bcol) <= set(blocks_df.columns) else blocks_df.copy()
            merged = sub.merge(bdf, left_on=['locus_id','genome','assigned_block_id'], right_on=['locus_id','genome','block_id'], how='left')
            merged['span_kb'] = pd.to_numeric(merged.get('span_kb', 0), errors='coerce').fillna(0)
            merged = merged.sort_values(by=['n_targets','span_kb'], ascending=[False, False])
            best = merged.iloc[0]
            sel_rows.append({
                'locus_id': best['locus_id'],
                'genome': best['genome'],
                'block_id': best['assigned_block_id'],
                'scaffold': best.get('scaffold',''),
                'strand': best.get('strand',''),
                'start': int(best.get('start', 0)) if pd.notnull(best.get('start', 0)) else 0,
                'end': int(best.get('end', 0)) if pd.notnull(best.get('end', 0)) else 0,
                'span_kb': float(best.get('span_kb', 0)),
                'num_query_matches': int(best.get('num_query_matches', 0)) if pd.notnull(best.get('num_query_matches', 0)) else 0,
                'num_target_proteins': int(best.get('num_target_proteins', 0)) if pd.notnull(best.get('num_target_proteins', 0)) else 0,
                'n_assigned_targets': int(best['n_targets']),
            })

    if sel_rows:
        sel_df = pd.DataFrame(sel_rows)
    else:
        sel_df = pd.DataFrame(columns=['locus_id','genome','block_id','scaffold','strand','start','end','span_kb','num_query_matches','num_target_proteins','n_assigned_targets'])

    # Fallback selection for (locus, genome) pairs with no target-supported block:
    # choose the strongest block by base_count, then density, then smaller span.
    # This ensures we still have a block to annotate when targets are absent.
    if not blocks_df.empty:
        have = set((r['locus_id'], r['genome']) for _, r in sel_df.iterrows()) if not sel_df.empty else set()
        all_pairs = set((str(r['locus_id']), str(r['genome'])) for _, r in blocks_df.iterrows())
        missing_pairs = all_pairs - have

        def base_count_row(r):
            if 'num_query_matches' in r:
                return r['num_query_matches']
            if 'num_target_proteins' in r:
                return r['num_target_proteins']
            return r.get('num_proteins', 0)

        add_rows = []
        for (locus_id, genome) in missing_pairs:
            cand = blocks_df[(blocks_df['locus_id'] == locus_id) & (blocks_df['genome'] == genome)].copy()
            if cand.empty:
                continue
            # Ensure span_kb is present
            if 'span_kb' not in cand.columns and {'start','end'} <= set(cand.columns):
                cand['span_kb'] = (cand['end'].astype(float) - cand['start'].astype(float)) / 1000.0
            cand['__base__'] = cand.apply(lambda r: base_count_row(r), axis=1)
            # density = base_count / max(1kb, span_kb)
            cand['__density__'] = cand.apply(lambda r: (r['__base__'] / max(1.0, float(r.get('span_kb', 0) or 0))), axis=1)
            # Sort by base desc, density desc, span asc
            cand = cand.sort_values(by=['__base__','__density__','span_kb'], ascending=[False, False, True])
            b = cand.iloc[0]
            add_rows.append({
                'locus_id': locus_id,
                'genome': genome,
                'block_id': b.get('block_id',''),
                'scaffold': b.get('scaffold',''),
                'strand': b.get('strand',''),
                'start': int(b.get('start', 0)) if pd.notnull(b.get('start', 0)) else 0,
                'end': int(b.get('end', 0)) if pd.notnull(b.get('end', 0)) else 0,
                'span_kb': float(b.get('span_kb', 0) or 0),
                'num_query_matches': int(b.get('num_query_matches', 0) or 0),
                'num_target_proteins': int(b.get('num_target_proteins', 0) or 0),
                'n_assigned_targets': 0,
            })

        if add_rows:
            add_df = pd.DataFrame(add_rows)
            sel_df = pd.concat([sel_df, add_df], ignore_index=True) if not sel_df.empty else add_df

    sel_out = args.output_dir / 'selected_synteny_blocks.tsv'
    if not sel_df.empty:
        sel_df.to_csv(sel_out, sep='\t', index=False)
        print(f'Selected synteny blocks for SwissProt: {sel_out} ({len(sel_df)})')
    else:
        print('No selected synteny blocks could be derived (no blocks found).')


if __name__ == '__main__':
    main()
