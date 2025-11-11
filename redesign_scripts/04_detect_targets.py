#!/usr/bin/env python3
"""
Phase 4 (redesign): BLAST target proteins against genomes and cluster hits into loci.

Changes vs legacy:
- Build combined multi-query FASTA from Phase 1_v2 per-locus targets
- Run tBLASTn once per genome
- Cluster hits into loci per (locus_id, scaffold, strand) with split gaps scaled by locus_scale

Outputs:
- all_target_loci.tsv
- combined_targets.faa
- blast_xml/*.xml
"""

from __future__ import annotations

from pathlib import Path
import argparse
import subprocess
import xml.etree.ElementTree as ET
from collections import defaultdict
from typing import Dict, List, Tuple
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def run_tblastn(query_file: Path, genome_db: Path, output_xml: Path, *, evalue: str, max_targets: int, threads: int) -> bool:
    cmd = [
        'tblastn', '-query', str(query_file), '-db', str(genome_db),
        '-outfmt', '5', '-evalue', str(evalue), '-max_target_seqs', str(max_targets), '-num_threads', str(threads)
    ]
    with open(output_xml, 'w') as out:
        res = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True)
    return res.returncode == 0


def parse_blast_xml(xml_file: Path) -> List[Dict]:
    hits: List[Dict] = []
    try:
        tree = ET.parse(xml_file); root = tree.getroot()
        for iteration in root.findall('.//Iteration'):
            qid = iteration.find('Iteration_query-def').text.split()[0]
            qlen = int(iteration.find('Iteration_query-len').text)
            for hit in iteration.findall('.//Hit'):
                hit_def = hit.find('Hit_def').text; hit_acc = hit.find('Hit_accession').text
                for hsp in hit.findall('.//Hsp'):
                    evalue = float(hsp.find('Hsp_evalue').text)
                    bits = float(hsp.find('Hsp_bit-score').text)
                    ident = int(hsp.find('Hsp_identity').text)
                    alen = int(hsp.find('Hsp_align-len').text)
                    sstart = int(hsp.find('Hsp_hit-from').text)
                    send = int(hsp.find('Hsp_hit-to').text)
                    if sstart <= send:
                        strand = '+'; start = sstart; end = send
                    else:
                        strand = '-'; start = send; end = sstart
                    hits.append({
                        'query_id': qid,
                        'query_length': qlen,
                        'scaffold': hit_acc,
                        'scaffold_desc': hit_def,
                        'strand': strand,
                        'start': start,
                        'end': end,
                        'evalue': evalue,
                        'bitscore': bits,
                        'align_length': alen,
                        'pident': (ident / max(1, alen)) * 100.0,
                    })
    except Exception:
        pass
    return hits


def merge_intervals(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    if not intervals:
        return []
    arr = sorted(intervals, key=lambda x: x[0])
    merged = [arr[0]]
    for s, e in arr[1:]:
        ls, le = merged[-1]
        if s <= le:
            merged[-1] = (ls, max(le, e))
        else:
            merged.append((s, e))
    return merged


def cluster_hits_for_locus(hits: List[Dict], min_split_gap_kb: int) -> List[Dict]:
    if not hits:
        return []
    out: List[Dict] = []
    by_key: Dict[Tuple[str, str], List[Dict]] = defaultdict(list)
    for h in hits:
        by_key[(h['scaffold'], h['strand'])].append(h)
    gap_bp = min_split_gap_kb * 1000
    for (scf, strand), arr in by_key.items():
        intervals = [(h['start'], h['end']) for h in arr]
        cov = merge_intervals(intervals)
        # derive split points where gap >= gap_bp
        splits: List[Tuple[int, int]] = []
        for i in range(len(cov) - 1):
            a = cov[i]; b = cov[i+1]
            if b[0] - a[1] >= gap_bp:
                splits.append((a[1], b[0]))
        # segments between boundaries
        positions = sorted([p for it in cov for p in it])
        if not positions:
            continue
        segs: List[Tuple[int, int]] = []
        cur_s = positions[0]; cur_e = positions[-1]
        if splits:
            cur_s = positions[0]
            for (gap_s, gap_e) in splits:
                segs.append((cur_s, gap_s))
                cur_s = gap_e
            if cur_s < positions[-1]:
                segs.append((cur_s, positions[-1]))
        else:
            segs.append((positions[0], positions[-1]))

        for seg_s, seg_e in segs:
            seg_hits = [h for h in arr if not (h['end'] < seg_s or h['start'] > seg_e)]
            if not seg_hits:
                continue
            out.append({
                'scaffold': scf,
                'strand': strand,
                'start': min(h['start'] for h in seg_hits),
                'end': max(h['end'] for h in seg_hits),
                'span_kb': (max(h['end'] for h in seg_hits) - min(h['start'] for h in seg_hits)) / 1000.0,
                'num_hits': len(seg_hits),
                'best_evalue': min(h['evalue'] for h in seg_hits),
                'best_bitscore': max(h['bitscore'] for h in seg_hits),
            })
    return out


def parse_args():
    p = argparse.ArgumentParser(description='Phase 4 (redesign): detect targets (scaled split gap)')
    p.add_argument('--locus-defs', type=Path, required=True, help='Phase1_v2 locus_definitions.tsv (for locus_scale)')
    p.add_argument('--phase1-dir', type=Path, required=True, help='Phase1_v2 directory containing <locus>/*_targets.faa')
    p.add_argument('--genome-db-dir', type=Path, required=True, help='Directory with *.nhr BLAST DBs')
    p.add_argument('--output-dir', type=Path, required=True, help='Phase4_v2 output directory')
    p.add_argument('--evalue', type=str, default='1e-5')
    p.add_argument('--max-targets', type=int, default=10000)
    p.add_argument('--threads', type=int, default=16)
    p.add_argument('--base-min-split-kb', type=int, default=40, help='Base coverage gap (kb) to split loci')
    return p.parse_args()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    blast_dir = args.output_dir / 'blast_xml'; blast_dir.mkdir(exist_ok=True)

    loci_df = pd.read_csv(args.locus_defs, sep='\t')
    locus_scale_map = {str(r['locus_id']): float(r.get('locus_scale', 1.0) or 1.0) for _, r in loci_df.iterrows()}

    # Build combined targets FASTA and map query_id -> locus_id
    combined = args.output_dir / 'combined_targets.faa'
    mapping: Dict[str, str] = {}
    records: List[SeqRecord] = []
    for _, row in loci_df.iterrows():
        locus_id = str(row['locus_id'])
        target_fa = args.phase1_dir / locus_id / f"{locus_id}_targets.faa"
        if not target_fa.exists():
            continue
        for rec in SeqIO.parse(str(target_fa), 'fasta'):
            qid = rec.id.split()[0]
            mapping[qid] = locus_id
            records.append(rec)
    if not records:
        print('ERROR: no targets found in phase1_v2 locus directories')
        return
    SeqIO.write(records, str(combined), 'fasta')

    # Discover genome DBs
    genome_dbs = {db.stem: args.genome_db_dir / db.stem for db in sorted(args.genome_db_dir.glob('*.nhr'))}

    all_out_rows: List[Dict] = []
    for genome_name, genome_db in genome_dbs.items():
        xml_path = blast_dir / f"{genome_name}.xml"
        if xml_path.exists():
            print(f"{genome_name}: using existing BLAST")
        else:
            print(f"{genome_name}: tBLASTn combined targets...", end='')
            ok = run_tblastn(combined, genome_db, xml_path, evalue=args.evalue, max_targets=args.max_targets, threads=args.threads)
            print(' done' if ok else ' FAILED')
            if not ok:
                continue

        hits = parse_blast_xml(xml_path)
        # attach locus_id via query_id mapping
        for h in hits:
            locus_id = mapping.get(h['query_id']) or mapping.get(h['query_id'].split('|')[0])
            if not locus_id:
                continue
            h['locus_id'] = locus_id
            h['genome'] = genome_name

        # group hits by locus
        by_locus: Dict[str, List[Dict]] = defaultdict(list)
        for h in hits:
            lid = h.get('locus_id')
            if lid:
                by_locus[lid].append(h)

        # cluster per locus with scaled split gaps
        for locus_id, lh in by_locus.items():
            scale = float(locus_scale_map.get(locus_id, 1.0) or 1.0)
            min_split_kb = max(10, min(200, int(round(args.base_min_split_kb * scale))))
            clusters = cluster_hits_for_locus(lh, min_split_kb)
            for c in clusters:
                all_out_rows.append({
                    'locus_id': locus_id,
                    'parent_locus': locus_id,
                    'genome': genome_name,
                    'scaffold': c['scaffold'],
                    'strand': c['strand'],
                    'start': c['start'],
                    'end': c['end'],
                    'span_kb': c['span_kb'],
                    'num_hits': c['num_hits'],
                    'best_evalue': c['best_evalue'],
                    'best_bitscore': c['best_bitscore'],
                })

    if all_out_rows:
        out_tsv = args.output_dir / 'all_target_loci.tsv'
        pd.DataFrame(all_out_rows).to_csv(out_tsv, sep='\t', index=False)
        print(f'Wrote: {out_tsv} ({len(all_out_rows)} loci)')
    else:
        print('No target loci found')


if __name__ == '__main__':
    main()

