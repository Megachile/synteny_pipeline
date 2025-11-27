#!/usr/bin/env python3
"""
Phase 2 (redesign): Detect synteny blocks using deduplicated flanking proteins.

Key changes vs legacy:
- Reads locus_scale from Phase 1 to scale max-gap/max-span and merge thresholds
- Keeps all blocks (no best-per-genome reduction)
- Merge-adjacent enabled by default to reduce fragmentation

Usage:
    python 02_synteny_detection.py \
        --locus-defs outputs/<family>/phase1_v2/locus_definitions.tsv \
        --genome-db-dir data/ragtag_dbs \
        --output-dir outputs/<family>/phase2_synteny_v2 \
        [--evalue 1e-5] [--base-max-targets 1000] [--threads 16]
"""

from __future__ import annotations

from pathlib import Path
import argparse
import subprocess
import xml.etree.ElementTree as ET
from collections import defaultdict
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Reuse normalize_scaffold if available
try:
    from synteny_utils import normalize_scaffold
except Exception:  # fallback
    def normalize_scaffold(name: str) -> str:
        return str(name).split()[0]


def is_valid_blast_xml(xml_file: Path) -> bool:
    try:
        with open(xml_file, 'rb') as f:
            f.seek(0, 2)
            n = f.tell()
            if n == 0:
                return False
            f.seek(max(0, n - 2048))
            tail = f.read().decode('utf-8', errors='ignore')
            return '</BlastOutput>' in tail
    except Exception:
        return False


def run_tblastn(query_file: Path, genome_db: Path, output_xml: Path, *, evalue: str, max_targets: int, threads: int) -> bool:
    cmd = [
        'tblastn',
        '-query', str(query_file),
        '-db', str(genome_db),
        '-outfmt', '5',
        '-evalue', str(evalue),
        '-max_target_seqs', str(max_targets),
        '-num_threads', str(threads),
    ]
    with open(output_xml, 'w') as out:
        res = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True)
    return res.returncode == 0


def parse_blast_xml(xml_file: Path) -> list[dict]:
    hits: list[dict] = []
    try:
        tree = ET.parse(xml_file)
        root = tree.getroot()
        for iteration in root.findall('.//Iteration'):
            query_id = iteration.find('Iteration_query-def').text.split()[0]
            for hit in iteration.findall('.//Hit'):
                hit_def = hit.find('Hit_def').text
                hit_accession = hit.find('Hit_accession').text
                for hsp in hit.findall('.//Hsp'):
                    evalue = float(hsp.find('Hsp_evalue').text)
                    bitscore = float(hsp.find('Hsp_bit-score').text)
                    identity = int(hsp.find('Hsp_identity').text)
                    alen = int(hsp.find('Hsp_align-len').text)
                    sstart = int(hsp.find('Hsp_hit-from').text)
                    send = int(hsp.find('Hsp_hit-to').text)
                    h_hseq = hsp.find('Hsp_hseq').text
                    h_qseq = hsp.find('Hsp_qseq').text
                    if sstart <= send:
                        strand = '+'; start = sstart; end = send
                    else:
                        strand = '-'; start = send; end = sstart
                    hits.append({
                        'qseqid': query_id,
                        'sseqid': hit_accession,
                        'scaffold_desc': hit_def,
                        'strand': strand,
                        'coord_start': start,
                        'coord_end': end,
                        'evalue': evalue,
                        'bitscore': bitscore,
                        'pident': (identity / max(1, alen)) * 100.0,
                        'length': alen,
                        'hit_protein_seq': h_hseq,
                        'query_protein_seq': h_qseq,
                    })
    except Exception:
        pass
    return hits


def _build_query_order_map(flanking_faa: Path) -> dict[str, int]:
    order: dict[str, int] = {}
    try:
        for i, rec in enumerate(SeqIO.parse(str(flanking_faa), 'fasta')):
            order[rec.id.split()[0]] = i
    except Exception:
        pass
    return order


def _select_blocks_by_chain(scaffold_hits: list[dict], strand: str, q_order_map: dict[str, int], block_id_start: int,
                            max_gap_kb: int, max_span_kb: int, min_proteins: int) -> tuple[list[dict], int]:
    blocks: list[dict] = []
    block_id = block_id_start
    max_gap_bp = max_gap_kb * 1000
    max_span_bp = max_span_kb * 1000

    # choose best HSP per query id (lower evalue, tie by bitscore)
    best_by_q: dict[str, dict] = {}
    for h in scaffold_hits:
        q = h['qseqid']
        cur = best_by_q.get(q)
        if cur is None or (h['evalue'] < cur['evalue'] or (h['evalue'] == cur['evalue'] and h['bitscore'] > cur['bitscore'])):
            best_by_q[q] = h

    seq = []
    for q, h in best_by_q.items():
        qi = q_order_map.get(q)
        if qi is None:
            qb = q.split('|')[0]
            qi = q_order_map.get(qb)
        if qi is None:
            continue
        pos = h['coord_start'] if strand == '+' else -h['coord_start']
        seq.append((pos, qi, h))
    if not seq:
        return blocks, block_id

    seq.sort(key=lambda x: x[0])
    unused = seq[:]

    # increasing chains
    while unused:
        chain = []
        last_q = -1
        rest = []
        for pos, qi, h in unused:
            if qi > last_q:
                chain.append((pos, qi, h)); last_q = qi
            else:
                rest.append((pos, qi, h))
        if len(chain) < min_proteins:
            break

        seg_start = chain[0][2]['coord_start']
        seg_end = chain[0][2]['coord_end']
        current = [chain[0][2]]
        for _, _, h in chain[1:]:
            prev = current[-1]
            if h['coord_start'] > prev['coord_end']:
                gap_bp = h['coord_start'] - prev['coord_end']
            elif prev['coord_start'] > h['coord_end']:
                gap_bp = prev['coord_start'] - h['coord_end']
            else:
                gap_bp = 0
            prospective_start = min(seg_start, h['coord_start'])
            prospective_end = max(seg_end, h['coord_end'])
            if gap_bp > max_gap_bp or (prospective_end - prospective_start) > max_span_bp:
                if len(current) >= min_proteins:
                    blocks.append({
                        'block_id': f"block_{block_id:05d}",
                        'scaffold': normalize_scaffold(current[0]['sseqid']),
                        'strand': strand,
                        'start': min(x['coord_start'] for x in current),
                        'end': max(x['coord_end'] for x in current),
                        'num_target_proteins': len(set(x['qseqid'] for x in current)),
                        'num_query_matches': len(set(x['qseqid'] for x in current)),
                        'hits': current,
                    })
                    block_id += 1
                seg_start = h['coord_start']; seg_end = h['coord_end']; current = [h]
            else:
                current.append(h)
                seg_start = min(seg_start, h['coord_start'])
                seg_end = max(seg_end, h['coord_end'])
        if len(current) >= min_proteins:
            blocks.append({
                'block_id': f"block_{block_id:05d}",
                'scaffold': normalize_scaffold(current[0]['sseqid']),
                'strand': strand,
                'start': min(x['coord_start'] for x in current),
                'end': max(x['coord_end'] for x in current),
                'num_target_proteins': len(set(x['qseqid'] for x in current)),
                'num_query_matches': len(set(x['qseqid'] for x in current)),
                'hits': current,
            })
            block_id += 1
        used_q = set(h['qseqid'] for _, _, h in chain)
        unused = [t for t in rest if t[2]['qseqid'] not in used_q]

    # decreasing chains (inversions)
    while unused:
        chain = []
        last_q = float('inf')
        rest = []
        for pos, qi, h in unused:
            if qi < last_q:
                chain.append((pos, qi, h)); last_q = qi
            else:
                rest.append((pos, qi, h))
        if len(chain) < min_proteins:
            break
        seg_start = chain[0][2]['coord_start']
        seg_end = chain[0][2]['coord_end']
        current = [chain[0][2]]
        for _, _, h in chain[1:]:
            prev = current[-1]
            if h['coord_start'] > prev['coord_end']:
                gap_bp = h['coord_start'] - prev['coord_end']
            elif prev['coord_start'] > h['coord_end']:
                gap_bp = prev['coord_start'] - h['coord_end']
            else:
                gap_bp = 0
            prospective_start = min(seg_start, h['coord_start'])
            prospective_end = max(seg_end, h['coord_end'])
            if gap_bp > max_gap_bp or (prospective_end - prospective_start) > max_span_bp:
                if len(current) >= min_proteins:
                    blocks.append({
                        'block_id': f"block_{block_id:05d}",
                        'scaffold': normalize_scaffold(current[0]['sseqid']),
                        'strand': strand,
                        'start': min(x['coord_start'] for x in current),
                        'end': max(x['coord_end'] for x in current),
                        'num_target_proteins': len(set(x['qseqid'] for x in current)),
                        'num_query_matches': len(set(x['qseqid'] for x in current)),
                        'hits': current,
                    })
                    block_id += 1
                seg_start = h['coord_start']; seg_end = h['coord_end']; current = [h]
            else:
                current.append(h)
                seg_start = min(seg_start, h['coord_start'])
                seg_end = max(seg_end, h['coord_end'])
        if len(current) >= min_proteins:
            blocks.append({
                'block_id': f"block_{block_id:05d}",
                'scaffold': normalize_scaffold(current[0]['sseqid']),
                'strand': strand,
                'start': min(x['coord_start'] for x in current),
                'end': max(x['coord_end'] for x in current),
                'num_target_proteins': len(set(x['qseqid'] for x in current)),
                'num_query_matches': len(set(x['qseqid'] for x in current)),
                'hits': current,
            })
            block_id += 1
        used_q = set(h['qseqid'] for _, _, h in chain)
        unused = [t for t in rest if t[2]['qseqid'] not in used_q]

    return blocks, block_id


def cluster_into_blocks(hits: list[dict], q_order_map: dict[str, int], max_gap_kb: int, max_span_kb: int, min_proteins: int) -> list[dict]:
    grouped: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for h in hits:
        grouped[(h['sseqid'], h['strand'])].append(h)
    blocks: list[dict] = []
    block_id = 1
    for (scaf, strand), arr in grouped.items():
        b, block_id = _select_blocks_by_chain(arr, strand, q_order_map, block_id, max_gap_kb, max_span_kb, min_proteins)
        blocks.extend(b)
    return blocks


def merge_adjacent_blocks(blocks: list[dict], max_merge_gap_kb: int, max_merged_span_kb: int) -> list[dict]:
    if not blocks:
        return blocks
    grouped: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for b in blocks:
        grouped[(b['scaffold'], b['strand'])].append(b)
    merged_all: list[dict] = []
    for key, arr in grouped.items():
        arr.sort(key=lambda x: int(x['start']))
        cur = None
        for b in arr:
            if cur is None:
                cur = {**b}
                cur['hits'] = list(cur.get('hits', []))
                continue
            gap_bp = max(0, int(b['start']) - int(cur['end']))
            span_bp = max(int(cur['end']), int(b['end'])) - min(int(cur['start']), int(b['start']))
            if gap_bp <= max_merge_gap_kb * 1000 and span_bp <= max_merged_span_kb * 1000:
                cur['start'] = min(int(cur['start']), int(b['start']))
                cur['end'] = max(int(cur['end']), int(b['end']))
                cur['hits'] = list(cur.get('hits', [])) + list(b.get('hits', []))
                qset = set(h['qseqid'] for h in cur['hits'])
                cur['num_target_proteins'] = len(qset)
                cur['num_query_matches'] = len(qset)
            else:
                merged_all.append(cur)
                cur = {**b}; cur['hits'] = list(cur.get('hits', []))
        if cur is not None:
            merged_all.append(cur)
    return merged_all


def save_hit_sequences(hits: list[dict], genome_id: str, locus_output: Path) -> int:
    seqs_dir = locus_output / 'hit_sequences'
    seqs_dir.mkdir(exist_ok=True, parents=True)
    records: list[SeqRecord] = []
    for i, h in enumerate(hits, 1):
        protein_seq = (h.get('hit_protein_seq') or '').replace('-', '')
        if not protein_seq:
            continue
        seq_id = f"{genome_id}_{h['qseqid']}_{h['sseqid']}_{h['coord_start']}-{h['coord_end']}_hit{i}"
        rec = SeqRecord(Seq(protein_seq), id=seq_id,
                        description=f"{genome_id} | {h['qseqid']} hit | {h['sseqid']}:{h['coord_start']}-{h['coord_end']}({h['strand']}) | pident={h['pident']:.1f}% | len={len(protein_seq)}aa")
        records.append(rec)
    if records:
        out_fa = seqs_dir / f"{genome_id}_hit_proteins.fasta"
        SeqIO.write(records, out_fa, 'fasta')
    return len(records)


def _filter_flanking(flanking_file: Path, locus_id: str, locus_out: Path,
                     closest_each_side: int | None, closest_internal: int | None) -> Path:
    """Filter flanking FASTA to closest anchors and save as <locus>_flanking_filtered.faa.

    If closest_each_side is None, returns the original path (and still writes a copy for downstream).
    """
    locus_out.mkdir(parents=True, exist_ok=True)
    filtered = locus_out / f"{locus_id}_flanking_filtered.faa"
    try:
        records = list(SeqIO.parse(str(flanking_file), 'fasta'))
        if closest_each_side is None:
            # Just write a copy for downstream consistency
            with open(filtered, 'w') as w:
                SeqIO.write(records, w, 'fasta')
            return filtered

        u_sel, i_sel, d_sel = [], [], []
        for rec in records:
            desc = rec.description or ''
            parts = desc.split()
            label = parts[1] if len(parts) > 1 else ''
            if label.startswith('U') and label[1:].isdigit():
                idx = int(label[1:])
                if idx <= int(closest_each_side):
                    u_sel.append((idx, rec))
            elif label.startswith('D') and label[1:].isdigit():
                idx = int(label[1:])
                if idx <= int(closest_each_side):
                    d_sel.append((idx, rec))
            elif label.startswith('I') and label[1:].isdigit():
                if closest_internal and int(label[1:]) <= int(closest_internal):
                    i_sel.append((int(label[1:]), rec))

        u_sel.sort(key=lambda x: x[0])
        i_sel.sort(key=lambda x: x[0])
        d_sel.sort(key=lambda x: x[0])
        selected = [r for _, r in u_sel] + [r for _, r in i_sel] + [r for _, r in d_sel]
        if not selected:
            # Fallback to copy all
            selected = records
        with open(filtered, 'w') as w:
            SeqIO.write(selected, w, 'fasta')
        return filtered
    except Exception:
        # Last resort: copy
        try:
            with open(flanking_file, 'r') as rf, open(filtered, 'w') as wf:
                wf.write(rf.read())
            return filtered
        except Exception:
            return flanking_file


def process_locus(row: pd.Series, genome_dbs: dict[str, Path], out_dir: Path,
                  evalue: str, base_max_targets: int, threads: int,
                  base_max_gap_kb: int, base_max_span_kb: int,
                  base_merge_gap_kb: int, base_merge_span_kb: int,
                  min_proteins: int,
                  closest_each_side: int | None,
                  closest_internal: int | None) -> list[dict]:
    locus_id = str(row['locus_id'])
    flanking_file = Path(row['flanking_file'])
    expected_chr = str(row.get('expected_chromosome', '') or '')
    locus_scale = float(row.get('locus_scale', 1.0) or 1.0)
    flanking_span_kb = row.get('flanking_span_kb', None)

    # Use data-driven parameters from Phase 1 if available (Balanced formula)
    if flanking_span_kb is not None and not pd.isna(flanking_span_kb) and flanking_span_kb > 0:
        max_gap_kb = int(flanking_span_kb / 2.5)
        max_span_kb = int(flanking_span_kb * 1.1)
        merge_gap_kb = int(flanking_span_kb / 5)  # More conservative for merging
        merge_span_kb = int(flanking_span_kb * 1.3)
        param_source = "data-driven"
    else:
        # Fallback to scaled parameters (with caps) if flanking_span_kb not available
        def scale(v: int, lo: int, hi: int) -> int:
            return max(lo, min(hi, int(round(v * locus_scale))))

        max_gap_kb = scale(base_max_gap_kb, 100, 1000)
        max_span_kb = scale(base_max_span_kb, 300, 1500)
        merge_gap_kb = scale(base_merge_gap_kb, 10, 200)
        merge_span_kb = scale(base_merge_span_kb, 400, 2000)
        param_source = "scaled"

    max_targets = max(50, int(round(base_max_targets * locus_scale)))

    print(f"  {locus_id}: {param_source} gap={max_gap_kb}kb span={max_span_kb}kb merge_gap={merge_gap_kb}kb merge_span={merge_span_kb}kb targets={max_targets}")

    # Output dirs
    locus_out = out_dir / locus_id
    locus_out.mkdir(parents=True, exist_ok=True)
    blast_dir = locus_out / 'blast_xml'
    blast_dir.mkdir(exist_ok=True)
    # Filter flanking anchors and write filtered FASTA
    flanking_filtered = _filter_flanking(flanking_file, locus_id, locus_out, closest_each_side, closest_internal)

    # Build query order from filtered flanking file
    q_order = _build_query_order_map(flanking_filtered)

    all_blocks: list[dict] = []
    for genome_name, genome_db in genome_dbs.items():
        if genome_name.startswith('GCA_') or genome_name.startswith('GCF_'):
            genome_id = '_'.join(genome_name.split('_')[:2])
        else:
            genome_id = genome_name

        xml_path = blast_dir / f"{genome_name}.xml"
        if xml_path.exists() and is_valid_blast_xml(xml_path):
            print(f"    {genome_name}: Using existing BLAST")
        else:
            if xml_path.exists():
                try: xml_path.unlink()
                except Exception: pass
            print(f"    {genome_name}: tBLASTn...", end='')
            ok = run_tblastn(flanking_filtered, genome_db, xml_path, evalue=evalue, max_targets=max_targets, threads=threads)
            print(" done" if ok else " FAILED")
            if not ok:
                continue

        hits = parse_blast_xml(xml_path)
        if not hits:
            continue
        _ = save_hit_sequences(hits, genome_id, locus_out)

        blocks = cluster_into_blocks(hits, q_order, max_gap_kb, max_span_kb, min_proteins)
        blocks = merge_adjacent_blocks(blocks, merge_gap_kb, merge_span_kb)

        # Annotate and stash
        for b in blocks:
            b['genome'] = genome_id
            b['locus_id'] = locus_id
            b['span_kb'] = round((int(b['end']) - int(b['start'])) / 1000.0, 1)
            b['expected_chromosome'] = expected_chr
        all_blocks.extend(blocks)

    # Write hits TSV
    all_hits_rows = []
    for b in all_blocks:
        for h in b.get('hits', []):
            row = {k: v for k, v in h.items() if k not in ('hit_protein_seq', 'query_protein_seq')}
            row['genome'] = b['genome']; row['block_id'] = b['block_id']
            all_hits_rows.append(row)
    if all_hits_rows:
        pd.DataFrame(all_hits_rows).to_csv(locus_out / 'flanking_blast_all.tsv', sep='\t', index=False)

    # Write per-locus blocks file
    if all_blocks:
        out_rows = []
        for b in all_blocks:
            out_rows.append({
                'locus_id': b['locus_id'],
                'genome': b['genome'],
                'block_id': b['block_id'],
                'scaffold': b['scaffold'],
                'strand': b['strand'],
                'start': b['start'],
                'end': b['end'],
                'span_kb': b['span_kb'],
                'num_target_proteins': b['num_target_proteins'],
                'num_query_matches': b['num_query_matches'],
                'expected_chromosome': b.get('expected_chromosome', ''),
            })
        blocks_df = pd.DataFrame(out_rows)
        (out_dir / f"{locus_id}_synteny_blocks.tsv").write_text(blocks_df.to_csv(sep='\t', index=False))

    return all_blocks


def parse_args():
    p = argparse.ArgumentParser(description='Phase 2 (redesign): synteny detection with scaled parameters')
    p.add_argument('--locus-defs', type=Path, required=True, help='Path to phase1_v2 locus_definitions.tsv')
    p.add_argument('--genome-db-dir', type=Path, required=True, help='Directory with *.nhr BLAST DBs')
    p.add_argument('--output-dir', type=Path, required=True, help='Directory for per-locus synteny outputs')
    p.add_argument('--evalue', type=str, default='1e-5')
    p.add_argument('--base-max-targets', type=int, default=1000)
    p.add_argument('--threads', type=int, default=16)
    p.add_argument('--base-max-gap-kb', type=int, default=300)
    p.add_argument('--base-max-span-kb', type=int, default=800)
    p.add_argument('--base-merge-gap-kb', type=int, default=50)
    p.add_argument('--base-merge-span-kb', type=int, default=800)
    p.add_argument('--min-proteins', type=int, default=3)
    p.add_argument('--closest-each-side', type=int, help='Select closest U/D anchors (e.g., 12) to filter flanking set')
    p.add_argument('--closest-internal', type=int, default=0, help='Select closest internal anchors (e.g., 4)')
    p.add_argument('--locus', type=str, help='Optional: process only this locus_id')
    p.add_argument('--exclude-genomes', type=Path, help='File with genome IDs to exclude (one per line, # comments)')
    return p.parse_args()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    print('=' * 80)
    print('PHASE 2 (REDESIGN): SYNTENY DETECTION (SCALED)')
    print('=' * 80)
    print(f"Locus defs: {args.locus_defs}")
    print(f"DB dir:     {args.genome_db_dir}")
    print(f"Out dir:    {args.output_dir}")

    loci_df = pd.read_csv(args.locus_defs, sep='\t')
    if args.locus:
        loci_df = loci_df[loci_df['locus_id'] == args.locus]
        if loci_df.empty:
            print(f"ERROR: locus {args.locus} not found")
            return
        print(f"Processing single locus: {args.locus}")
    else:
        print(f"Loaded loci: {len(loci_df)}")

    # Load genome exclusion list if provided
    exclude_genomes: set[str] = set()
    if args.exclude_genomes and args.exclude_genomes.exists():
        with open(args.exclude_genomes) as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    exclude_genomes.add(line)
        print(f"Exclusion list: {len(exclude_genomes)} genomes")

    # Discover genome DBs
    genome_dbs: dict[str, Path] = {}
    for db in sorted(args.genome_db_dir.glob('*.nhr')):
        db_name = db.stem
        # Extract genome ID (GCA_XXXXXXXX.X or GCF_XXXXXXXX.X) for exclusion check
        if db_name.startswith('GCA_') or db_name.startswith('GCF_'):
            genome_id = '_'.join(db_name.split('_')[:2])
            if genome_id in exclude_genomes:
                continue
        genome_dbs[db_name] = args.genome_db_dir / db_name
    print(f"Genomes: {len(genome_dbs)}" + (f" ({len(exclude_genomes)} excluded)" if exclude_genomes else ""))

    all_blocks: list[dict] = []
    for _, row in loci_df.iterrows():
        blocks = process_locus(row, genome_dbs, args.output_dir,
                               args.evalue, args.base_max_targets, args.threads,
                               args.base_max_gap_kb, args.base_max_span_kb,
                               args.base_merge_gap_kb, args.base_merge_span_kb,
                               args.min_proteins,
                               args.closest_each_side, args.closest_internal)
        all_blocks.extend(blocks)

    # Combined file for convenience
    if all_blocks:
        rows = []
        for b in all_blocks:
            rows.append({
                'locus_id': b['locus_id'], 'genome': b['genome'], 'block_id': b['block_id'],
                'scaffold': b['scaffold'], 'strand': b['strand'],
                'start': b['start'], 'end': b['end'], 'span_kb': b['span_kb'],
                'num_target_proteins': b['num_target_proteins'], 'num_query_matches': b['num_query_matches'],
                'expected_chromosome': b.get('expected_chromosome', ''),
            })
        pd.DataFrame(rows).to_csv(args.output_dir / 'combined_synteny_blocks.tsv', sep='\t', index=False)

    print('\nDone.')


if __name__ == '__main__':
    main()
