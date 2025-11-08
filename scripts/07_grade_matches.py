#!/usr/bin/env python3
"""
Phase 6 grading: assign match quality grades by combining BLAST strength and
Exonerate integrity into a concise label per syntenic target.

Inputs:
  --syntenic: 05_classified/*/syntenic_targets.tsv (Phase 5 output)
  --extracted-dir: 06_extracted_sequences (Phase 6 output dir)
  --query-proteins: 04_target_genes/combined_targets.faa (to get query lengths)

Output:
  A TSV written next to the syntenic file as <dir>/syntenic_targets_graded.tsv with:
    grade ∈ {intact, degraded_fragment, degraded_pseudogene, degraded_no_cds, bad_match}
    status_exonerate ∈ {intact, fragment, pseudogene, no_cds, unknown}
    coverage_pct (if available), plus BLAST fields for context.

Grading rules (default thresholds; adjust as needed):
  Intact: Exonerate status == intact
  Degraded pseudogene: Exonerate == pseudogene
  Degraded fragment: Exonerate == fragment AND (coverage >= 70% OR BLAST medium)
  Degraded no_cds: Exonerate == no_cds AND BLAST strong
  Bad match: otherwise

BLAST tiers:
  Strong: best_evalue <= 1e-30 OR best_bitscore >= 100
  Medium: best_evalue <= 1e-10 OR best_bitscore >= 70
"""

from pathlib import Path
import argparse
import pandas as pd
from Bio import SeqIO


def load_query_lengths(query_faa: Path) -> dict:
    qlen = {}
    try:
        for rec in SeqIO.parse(str(query_faa), 'fasta'):
            qlen[rec.id.split()[0]] = len(rec.seq)
    except Exception:
        pass
    return qlen


def scan_extraction_dir(extracted_root: Path, genome: str, locus_name: str) -> dict:
    """Inspect Phase 6 outputs for a specific target to recover Exonerate status and lengths.

    Returns dict with keys: status_exonerate, protein_len, cds_len, coverage_pct (optional)
    """
    info = {
        'status_exonerate': 'unknown',
        'protein_len': None,
        'cds_len': None,
        'coverage_pct': None,
    }

    target_dir = extracted_root / genome / locus_name
    if not target_dir.exists():
        return info

    # Prefer gene1_* files; fall back to any *_protein/_cds
    protein = None
    cds = None
    gene1_p = list(target_dir.glob('*_gene1_protein.fasta'))
    gene1_c = list(target_dir.glob('*_gene1_cds.fasta'))
    if gene1_p:
        protein = gene1_p[0]
    else:
        any_p = list(target_dir.glob('*_protein.fasta'))
        if any_p:
            protein = any_p[0]

    if gene1_c:
        cds = gene1_c[0]
    else:
        any_c = list(target_dir.glob('*_cds.fasta'))
        if any_c:
            cds = any_c[0]

    # Parse protein length from sequence (header may include length, but sequence is definitive)
    if protein and protein.exists():
        try:
            for rec in SeqIO.parse(str(protein), 'fasta'):
                info['protein_len'] = len(rec.seq)
                break
        except Exception:
            pass

    # Parse exonerate status from CDS header (status:... embedded), and get CDS length
    if cds and cds.exists():
        try:
            for rec in SeqIO.parse(str(cds), 'fasta'):
                header = rec.description
                info['cds_len'] = len(rec.seq)
                # Extract status: <word>
                # Header template: status:<functional_status> exons:<n>
                for token in header.split():
                    if token.startswith('status:'):
                        info['status_exonerate'] = token.split(':', 1)[1]
                        break
                break
        except Exception:
            pass

    return info


def grade_row(row, ex_info, qlen_map):
    # BLAST tiers
    e = row.get('best_evalue', 1.0)
    bs = row.get('best_bitscore', 0)
    try:
        e = float(e)
    except Exception:
        e = 1.0
    try:
        bs = float(bs)
    except Exception:
        bs = 0.0
    strong = (e <= 1e-30) or (bs >= 100)
    medium = (e <= 1e-10) or (bs >= 70)

    # Coverage if available
    cov = None
    qid = str(row.get('query_id', '')).split()[0]
    qlen = qlen_map.get(qid)
    if ex_info.get('protein_len') and qlen:
        cov = ex_info['protein_len'] / qlen * 100.0 if qlen > 0 else None

    status = ex_info.get('status_exonerate', 'unknown')

    # Grade assignment
    if status == 'intact':
        grade = 'intact'
    elif status == 'pseudogene':
        grade = 'degraded_pseudogene'
    elif status == 'fragment':
        if cov is not None and cov >= 70:
            grade = 'degraded_fragment'
        elif medium:
            grade = 'degraded_fragment'
        else:
            grade = 'bad_match'
    elif status == 'no_cds':
        grade = 'degraded_no_cds' if strong or medium else 'bad_match'
    else:
        # No exonerate info: accept only if strong; otherwise bad
        grade = 'degraded_no_cds' if strong else 'bad_match'

    return grade, status, cov


def main():
    ap = argparse.ArgumentParser(description='Assign quality grades to Phase 6 matches (syntenic targets)')
    ap.add_argument('--syntenic', type=Path, required=True, help='Path to syntenic_targets.tsv')
    ap.add_argument('--extracted-dir', type=Path, required=True, help='Path to 06_extracted_sequences directory')
    ap.add_argument('--query-proteins', type=Path, required=True, help='Path to combined_targets.faa')
    ap.add_argument('--output', type=Path, help='Optional path for graded TSV output')
    args = ap.parse_args()

    df = pd.read_csv(args.syntenic, sep='\t')
    qlen = load_query_lengths(args.query_proteins)
    extracted_root = args.extracted_dir

    grades = []
    for _, row in df.iterrows():
        genome = str(row.get('genome', ''))
        locus_name = str(row.get('locus_name', ''))
        ex_info = scan_extraction_dir(extracted_root, genome, locus_name)
        g, status, cov = grade_row(row, ex_info, qlen)
        grades.append({
            'locus_name': locus_name,
            'genome': genome,
            'assigned_to': row.get('assigned_to', row.get('parent_locus', '')),
            'scaffold': row.get('scaffold', ''),
            'start': row.get('start', ''),
            'end': row.get('end', ''),
            'query_id': row.get('query_id', ''),
            'best_evalue': row.get('best_evalue', ''),
            'best_bitscore': row.get('best_bitscore', ''),
            'status_exonerate': status,
            'coverage_pct': None if cov is None else round(cov, 1),
            'grade': g
        })

    out_df = pd.DataFrame(grades)
    if args.output:
        out_path = args.output
    else:
        out_path = args.syntenic.parent / 'syntenic_targets_graded.tsv'
    out_df.to_csv(out_path, sep='\t', index=False)
    print(f'Wrote: {out_path}')


if __name__ == '__main__':
    main()

