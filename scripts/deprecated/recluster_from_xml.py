#!/usr/bin/env python3
"""
Recluster synteny blocks from existing tBLASTn XML (no new BLAST), then optionally rerank by target overlap.

Inputs (within a family root dir):
  - 02_synteny_blocks/<locus>/blast_xml/*.xml (per-genome BLAST XML saved by Phase 2)
  - 04_target_genes/all_target_loci.tsv

Parameters:
  --max-gap-kb {100,200,300}
  --max-span-kb {300,500,800}
  --min-proteins {2,3}

Outputs:
  - 02_synteny_blocks_tuned/<grid_key>/<locus>_synteny_blocks.tsv
  - 02_synteny_blocks_tuned/<grid_key>/all_synteny_blocks.tsv
  - 03_filtered_blocks_tuned_<grid_key>/synteny_blocks_filtered.tsv (reranked by target overlap)
"""

from pathlib import Path
import argparse
import pandas as pd
import sys

# Import helpers from existing Phase 2 script
sys.path.insert(0, str(Path(__file__).parent))
from synteny_detection import parse_blast_xml, _build_query_order_map, cluster_into_blocks
from synteny_utils import normalize_scaffold
from rerank_blocks_by_targets import rerank_family


def recluster_locus(locus_dir: Path, locus_id: str, grid_key: str, max_gap_kb: int, max_span_kb: int, min_proteins: int) -> pd.DataFrame:
    xml_dir = locus_dir / 'blast_xml'
    if not xml_dir.exists():
        return pd.DataFrame()

    # Build dummy query order from original flanking file to preserve intended order (if present)
    flanking_fa = locus_dir.parent.parent / 'phase12_landmark' / f'{locus_id}_flanking.faa'
    qmap = _build_query_order_map(flanking_fa) if flanking_fa.exists() else {}

    # Collect hits per genome XML and cluster them
    blocks = []
    for xml in sorted(xml_dir.glob('*.xml')):
        genome_name = xml.stem
        hits = parse_blast_xml(xml)
        if not hits:
            continue
        # Cluster into blocks with provided parameters
        # monkey-patch min proteins by filtering afterward
        b = cluster_into_blocks(hits, qmap, max_gap_kb=max_gap_kb, max_span_kb=max_span_kb)
        for blk in b:
            if blk.get('num_query_matches', 0) < min_proteins:
                continue
            blk['genome'] = genome_name
            blk['locus_id'] = locus_id
            blk['span_kb'] = round((blk['end'] - blk['start']) / 1000, 1)
            blk['scaffold'] = normalize_scaffold(blk['scaffold'])
            blocks.append(blk)

    if not blocks:
        return pd.DataFrame()
    return pd.DataFrame(blocks)


def recluster_family(fam_dir: Path, max_gap_kb: int, max_span_kb: int, min_proteins: int) -> Path:
    fam_dir = fam_dir.resolve()
    grid_key = f"gap{max_gap_kb}_span{max_span_kb}_min{min_proteins}"
    out02 = fam_dir / '02_synteny_blocks_tuned' / grid_key
    out02.mkdir(parents=True, exist_ok=True)

    loci_files = list((fam_dir / '02_synteny_blocks').glob('*_synteny_blocks.tsv'))
    locus_ids = sorted({p.stem.replace('_synteny_blocks','') for p in loci_files if p.name != 'all_synteny_blocks.tsv'})
    all_blocks = []
    for locus_id in locus_ids:
        locus_dir = fam_dir / '02_synteny_blocks' / locus_id
        df = recluster_locus(locus_dir, locus_id, grid_key, max_gap_kb, max_span_kb, min_proteins)
        if df.empty:
            continue
        # Save per-locus
        df[['locus_id','genome','block_id','scaffold','strand','start','end','span_kb','num_target_proteins','num_query_matches']].to_csv(
            out02 / f'{locus_id}_synteny_blocks.tsv', sep='\t', index=False
        )
        all_blocks.append(df)

    if not all_blocks:
        # Return empty file path
        combined = out02 / 'all_synteny_blocks.tsv'
        if not combined.exists():
            pd.DataFrame(columns=['locus_id','genome','block_id','scaffold','strand','start','end','span_kb','num_target_proteins','num_query_matches']).to_csv(combined, sep='\t', index=False)
        return combined

    combined_df = pd.concat(all_blocks, ignore_index=True)
    combined = out02 / 'all_synteny_blocks.tsv'
    combined_df[['locus_id','genome','block_id','scaffold','strand','start','end','span_kb','num_target_proteins','num_query_matches']].to_csv(
        combined, sep='\t', index=False
    )

    # Rerank by target overlap to create filtered set
    out03 = fam_dir / f'03_filtered_blocks_tuned_{grid_key}'
    out03.mkdir(parents=True, exist_ok=True)
    # Temporarily symlink/copy combined into expected location? Instead, let rerank read 02 per-locus tsv already saved.
    # rerank_family reads from 02_synteny_blocks; we want it to see tuned 02 dir.
    # Easiest: temporarily rename directories in object? Simpler: run rerank on fam_dir after copying tuned files into a temp
    # directory under fam_dir that mimics 02_synteny_blocks, then move output to tuned 03.

    temp02 = fam_dir / f'_temp_02_for_rerank_{grid_key}'
    if temp02.exists():
        for p in temp02.glob('*'):
            p.unlink()
    else:
        temp02.mkdir()

    # Copy tuned per-locus files
    for p in out02.glob('*_synteny_blocks.tsv'):
        (temp02 / p.name).write_text(p.read_text())

    # Monkey-patch family structure by pointing rerank at temp02 via symlink from expected dir name
    real02 = fam_dir / '02_synteny_blocks'
    backup02 = fam_dir / f'_backup_02_{grid_key}'
    swapped = False
    try:
        if real02.exists():
            real02.rename(backup02)
            swapped = True
        temp02.rename(real02)
        # Run rerank
        rerank_family(fam_dir)
        # Move reranked output into tuned 03 directory
        src = fam_dir / '03_filtered_blocks_reranked' / 'synteny_blocks_filtered.tsv'
        dst = out03 / 'synteny_blocks_filtered.tsv'
        if src.exists():
            dst.write_text(src.read_text())
    finally:
        # Restore original structure
        if real02.exists():
            # Move temp (now real02) back to temp02 name
            real02.rename(temp02)
        if swapped and backup02.exists():
            backup02.rename(real02)
        # Cleanup temp
        for p in temp02.glob('*'):
            p.unlink()
        temp02.rmdir()

    return combined


def find_families(outputs_root: Path):
    families = []
    for p in outputs_root.iterdir():
        if not p.is_dir():
            continue
        if (p / '02_synteny_blocks').exists() and (p / '04_target_genes' / 'all_target_loci.tsv').exists():
            families.append(p)
    return sorted(families)


def main():
    ap = argparse.ArgumentParser(description='Recluster synteny blocks from existing BLAST XML (no BLAST), then rerank by targets.')
    ap.add_argument('--family-dir', type=Path, help='Run only for this family (outputs/<FAM>)')
    ap.add_argument('--outputs-root', type=Path, default=Path('outputs'))
    ap.add_argument('--max-gap-kb', type=int, default=200)
    ap.add_argument('--max-span-kb', type=int, default=500)
    ap.add_argument('--min-proteins', type=int, default=2)
    args = ap.parse_args()

    fams = [args.family_dir] if args.family_dir else find_families(args.outputs_root)
    print(f'Found {len(fams)} families')

    rows = []
    for fam in fams:
        print(f'=== {fam.name} ===')
        combined = recluster_family(fam, args.max_gap_kb, args.max_span_kb, args.min_proteins)
        rows.append({'family': fam.name, 'combined_blocks': str(combined)})

    if rows:
        if len(fams) == 1:
            out = fams[0] / 'recluster_summary.tsv'
            pd.DataFrame(rows).to_csv(out, sep='\t', index=False)
            print(f'Summary: {out}')
        else:
            pd.DataFrame(rows).to_csv('recluster_summary.tsv', sep='\t', index=False)
            print('Summary: recluster_summary.tsv')


if __name__ == '__main__':
    main()
