#!/usr/bin/env python3
"""
Phase 2b: Re-cluster existing BLAST hits with updated merge parameters.

Reads flanking_blast_all.tsv files and re-clusters using flanking_span_kb
from locus_definitions.tsv to properly merge adjacent blocks.

Usage:
    python 02b_recluster_blocks.py \
        --locus-defs outputs/FAMILY/phase1_v2/locus_definitions.tsv \
        --synteny-dir outputs/FAMILY/phase2_synteny_v2
"""

from pathlib import Path
import argparse
import pandas as pd
from collections import defaultdict


def cluster_hits_into_blocks(hits_df: pd.DataFrame, max_gap_kb: int, max_span_kb: int, min_proteins: int) -> list[dict]:
    """
    Cluster BLAST hits into synteny blocks based on coordinate proximity.
    """
    if hits_df.empty:
        return []

    blocks = []

    # Group by genome, scaffold, strand
    for (genome, scaffold, strand), group in hits_df.groupby(['genome', 'sseqid', 'strand']):
        group = group.sort_values('coord_start')

        current_block = None
        block_id = 1

        for _, hit in group.iterrows():
            if current_block is None:
                current_block = {
                    'genome': genome,
                    'scaffold': scaffold,
                    'strand': strand,
                    'start': hit['coord_start'],
                    'end': hit['coord_end'],
                    'hits': [hit.to_dict()],
                    'block_id': f'block_{block_id:05d}'
                }
            else:
                gap = hit['coord_start'] - current_block['end']
                new_span = hit['coord_end'] - current_block['start']

                if gap <= max_gap_kb * 1000 and new_span <= max_span_kb * 1000:
                    # Extend current block
                    current_block['end'] = max(current_block['end'], hit['coord_end'])
                    current_block['hits'].append(hit.to_dict())
                else:
                    # Save current block and start new one
                    if len(set(h['qseqid'] for h in current_block['hits'])) >= min_proteins:
                        blocks.append(current_block)
                    block_id += 1
                    current_block = {
                        'genome': genome,
                        'scaffold': scaffold,
                        'strand': strand,
                        'start': hit['coord_start'],
                        'end': hit['coord_end'],
                        'hits': [hit.to_dict()],
                        'block_id': f'block_{block_id:05d}'
                    }

        # Don't forget last block
        if current_block and len(set(h['qseqid'] for h in current_block['hits'])) >= min_proteins:
            blocks.append(current_block)

    return blocks


def merge_adjacent_blocks(blocks: list[dict], merge_gap_kb: int, merge_span_kb: int) -> list[dict]:
    """
    Merge adjacent blocks on the same scaffold/strand if within merge thresholds.
    """
    if not blocks:
        return blocks

    # Group by scaffold and strand
    grouped = defaultdict(list)
    for b in blocks:
        key = (b['genome'], b['scaffold'], b['strand'])
        grouped[key].append(b)

    merged_all = []

    for key, arr in grouped.items():
        arr.sort(key=lambda x: x['start'])

        current = None
        for b in arr:
            if current is None:
                current = {**b, 'hits': list(b['hits'])}
            else:
                gap = b['start'] - current['end']
                span = max(current['end'], b['end']) - min(current['start'], b['start'])

                if gap <= merge_gap_kb * 1000 and span <= merge_span_kb * 1000:
                    # Merge
                    current['start'] = min(current['start'], b['start'])
                    current['end'] = max(current['end'], b['end'])
                    current['hits'].extend(b['hits'])
                else:
                    merged_all.append(current)
                    current = {**b, 'hits': list(b['hits'])}

        if current:
            merged_all.append(current)

    return merged_all


def recluster_locus(locus_id: str, synteny_dir: Path, flanking_span_kb: float,
                    base_max_gap_kb: int = 300, base_max_span_kb: int = 800,
                    min_proteins: int = 3) -> pd.DataFrame:
    """
    Re-cluster BLAST hits for a single locus using updated parameters.
    """
    locus_dir = synteny_dir / locus_id
    hits_file = locus_dir / 'flanking_blast_all.tsv'

    if not hits_file.exists():
        print(f"  Warning: {hits_file} not found, skipping")
        return pd.DataFrame()

    # Read existing hits
    hits_df = pd.read_csv(hits_file, sep='\t')
    if hits_df.empty:
        return pd.DataFrame()

    # Calculate parameters from flanking_span_kb
    if flanking_span_kb and flanking_span_kb > 0:
        max_gap_kb = int(flanking_span_kb / 2.5)
        max_span_kb = int(flanking_span_kb * 1.1)
        merge_gap_kb = int(flanking_span_kb / 5)
        merge_span_kb = int(flanking_span_kb * 1.3)
        param_source = "data-driven"
    else:
        max_gap_kb = base_max_gap_kb
        max_span_kb = base_max_span_kb
        merge_gap_kb = 50
        merge_span_kb = 800
        param_source = "default"

    print(f"  {locus_id}: {param_source} merge_gap={merge_gap_kb}kb merge_span={merge_span_kb}kb")

    # Cluster hits
    blocks = cluster_hits_into_blocks(hits_df, max_gap_kb, max_span_kb, min_proteins)

    # Merge adjacent blocks
    blocks = merge_adjacent_blocks(blocks, merge_gap_kb, merge_span_kb)

    # Convert to DataFrame
    rows = []
    for b in blocks:
        unique_queries = set(h['qseqid'] for h in b['hits'])
        rows.append({
            'locus_id': locus_id,
            'genome': b['genome'],
            'block_id': b['block_id'],
            'scaffold': b['scaffold'],
            'strand': b['strand'],
            'start': int(b['start']),
            'end': int(b['end']),
            'width_kb': round((b['end'] - b['start']) / 1000, 1),
            'num_target_proteins': len(unique_queries),
            'num_query_matches': len(unique_queries),
        })

    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description='Re-cluster BLAST hits with updated merge parameters')
    parser.add_argument('--locus-defs', type=Path, required=True)
    parser.add_argument('--synteny-dir', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, help='Output dir (default: same as synteny-dir)')
    args = parser.parse_args()

    output_dir = args.output_dir or args.synteny_dir

    # Read locus definitions
    locus_df = pd.read_csv(args.locus_defs, sep='\t')
    print(f"Loaded {len(locus_df)} locus definitions")

    all_blocks = []

    for _, row in locus_df.iterrows():
        locus_id = row['locus_id']
        flanking_span_kb = row.get('flanking_span_kb', None)

        blocks_df = recluster_locus(locus_id, args.synteny_dir, flanking_span_kb)

        if not blocks_df.empty:
            # Save per-locus file
            out_file = output_dir / f'{locus_id}_synteny_blocks.tsv'
            blocks_df.to_csv(out_file, sep='\t', index=False)
            print(f"    Wrote {len(blocks_df)} blocks to {out_file.name}")
            all_blocks.append(blocks_df)

    # Save combined file
    if all_blocks:
        combined = pd.concat(all_blocks, ignore_index=True)
        combined_file = output_dir / 'combined_synteny_blocks.tsv'
        combined.to_csv(combined_file, sep='\t', index=False)
        print(f"\nCombined: {len(combined)} blocks -> {combined_file}")

    print("\nDone! Now run Phase 3 to re-aggregate.")


if __name__ == '__main__':
    main()
