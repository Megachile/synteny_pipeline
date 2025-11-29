#!/usr/bin/env python3
"""
Phase 7 (Helixer): Validate synteny by SwissProt annotation of flanking genes.

For each flanking gene match (Helixer flanking → BK flanking), check if both
genes have the same SwissProt annotation. This validates whether the synteny
calls are real (same gene) or spurious (different genes that happen to match).

Inputs:
- Phase 4b flanking proteins (flanking_proteins.faa)
- Phase 5 v3 flanking matches (flanking_matches.tsv)
- Phase 1 BK flanking proteins
- SwissProt database

Outputs:
- flanking_swissprot.tsv: SwissProt annotations for all flanking proteins
- synteny_validation.tsv: Match pairs with agreement/disagreement status
- validation_summary.tsv: Aggregate validation stats
"""

from __future__ import annotations

import argparse
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Dict, Set, Tuple

import pandas as pd


def clean_fasta_for_diamond(input_faa: Path, output_faa: Path) -> int:
    """Clean FASTA file for DIAMOND by replacing invalid characters."""
    n_written = 0
    with open(input_faa) as f_in, open(output_faa, 'w') as f_out:
        current_header = None
        current_seq = []
        for line in f_in:
            line = line.strip()
            if line.startswith('>'):
                if current_header and current_seq:
                    seq = ''.join(current_seq).replace('.', 'X')
                    if len(seq) > 0:
                        f_out.write(f"{current_header}\n{seq}\n")
                        n_written += 1
                current_header = line
                current_seq = []
            else:
                current_seq.append(line)
        if current_header and current_seq:
            seq = ''.join(current_seq).replace('.', 'X')
            if len(seq) > 0:
                f_out.write(f"{current_header}\n{seq}\n")
                n_written += 1
    return n_written


def run_diamond_swissprot(
    query_faa: Path,
    swissprot_db: Path,
    output_tsv: Path,
    evalue: float = 1e-5,
    threads: int = 8,
) -> bool:
    """Run DIAMOND blastp against SwissProt."""
    cmd = [
        'diamond', 'blastp',
        '--query', str(query_faa),
        '--db', str(swissprot_db).replace('.dmnd', ''),
        '--out', str(output_tsv),
        '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length', 'evalue', 'bitscore', 'stitle',
        '--evalue', str(evalue),
        '--max-target-seqs', '1',
        '--threads', str(threads),
        '--quiet'
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.returncode == 0


def parse_swissprot_hit(sseqid: str, stitle: str) -> Dict:
    """Parse SwissProt hit ID and description."""
    # sseqid format could be:
    # - sp|P12345|GENE_ORGANISM
    # - P12345.1 (just accession with version)
    parts = sseqid.split('|')
    if len(parts) >= 3:
        accession = parts[1]
        gene_org = parts[2]
        gene = gene_org.split('_')[0] if '_' in gene_org else gene_org
    else:
        # Just accession, possibly with version
        accession = sseqid.split('.')[0] if '.' in sseqid else sseqid
        gene = ''

    # Extract protein name from description
    # Format: P12345.1 RecName: Full=Protein name; ...
    protein_name = ''
    if 'RecName: Full=' in stitle:
        name_start = stitle.find('RecName: Full=') + len('RecName: Full=')
        name_end = stitle.find(';', name_start)
        if name_end == -1:
            name_end = stitle.find(' [', name_start)
        if name_end == -1:
            name_end = len(stitle)
        protein_name = stitle[name_start:name_end].strip()

    # Extract description (before OS=)
    desc = stitle.split(' OS=')[0] if ' OS=' in stitle else stitle

    return {
        'swissprot_accession': accession,
        'swissprot_gene': gene,
        'swissprot_protein_name': protein_name,
        'swissprot_desc': desc,
    }


def collect_bk_flanking_proteins(phase1_dir: Path) -> Tuple[Dict[str, str], Path]:
    """
    Collect all BK flanking proteins from Phase 1 and write combined FASTA.
    Returns: (protein_id -> sequence dict, combined_faa_path)
    """
    proteins = {}

    locus_defs = phase1_dir / 'locus_definitions.tsv'
    if not locus_defs.exists():
        return proteins, None

    df = pd.read_csv(locus_defs, sep='\t')

    for _, row in df.iterrows():
        locus_id = row['locus_id']
        flanking_file = phase1_dir / f"{locus_id}_flanking.faa"

        if not flanking_file.exists():
            continue

        current_id = None
        current_seq = []

        with open(flanking_file) as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_id and current_seq:
                        proteins[current_id] = ''.join(current_seq)
                    # Extract XP ID
                    header = line[1:].split()[0]
                    current_id = header.split('|')[0] if '|' in header else header
                    current_seq = []
                else:
                    current_seq.append(line)

            if current_id and current_seq:
                proteins[current_id] = ''.join(current_seq)

    return proteins


def main():
    parser = argparse.ArgumentParser(
        description="Phase 7 (Helixer): Validate synteny via SwissProt annotation"
    )
    parser.add_argument("--family", required=True)
    parser.add_argument("--phase4b-dir", type=Path, required=True,
                        help="Phase 4b output with flanking_proteins.faa")
    parser.add_argument("--phase5-dir", type=Path, required=True,
                        help="Phase 5 v3 output with flanking_matches.tsv")
    parser.add_argument("--phase1-dir", type=Path, required=True,
                        help="Phase 1 directory with BK/LB flanking proteins")
    parser.add_argument("--swissprot-db", type=Path, required=True,
                        help="SwissProt DIAMOND database (.dmnd)")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--evalue", type=float, default=1e-5)
    parser.add_argument("--threads", type=int, default=8)

    args = parser.parse_args()

    print("=" * 80)
    print("PHASE 7 (HELIXER): VALIDATE SYNTENY BY SWISSPROT")
    print("=" * 80)
    print()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Load flanking matches from Phase 5 v3
    matches_file = args.phase5_dir / "flanking_matches.tsv"
    if not matches_file.exists():
        print(f"[ERROR] No flanking matches file: {matches_file}")
        return 1

    matches_df = pd.read_csv(matches_file, sep='\t')
    print(f"Loaded {len(matches_df)} flanking gene matches")

    # Get unique Helixer flanking IDs that had matches
    helixer_ids_with_matches = set(matches_df['helixer_flanking'].unique())
    bk_xp_ids_matched = set(matches_df['bk_xp'].unique())

    print(f"  Unique Helixer flanking: {len(helixer_ids_with_matches)}")
    print(f"  Unique BK XP IDs: {len(bk_xp_ids_matched)}")

    # Load and clean Helixer flanking proteins
    helixer_faa = args.phase4b_dir / "flanking_proteins.faa"
    if not helixer_faa.exists():
        print(f"[ERROR] No flanking proteins: {helixer_faa}")
        return 1

    helixer_cleaned = args.output_dir / "helixer_flanking_cleaned.faa"
    n_helixer = clean_fasta_for_diamond(helixer_faa, helixer_cleaned)
    print(f"Cleaned {n_helixer} Helixer flanking proteins")

    # Collect BK flanking proteins
    print("Collecting BK flanking proteins from Phase 1...")
    bk_proteins = collect_bk_flanking_proteins(args.phase1_dir)
    print(f"  Found {len(bk_proteins)} unique BK flanking proteins")

    # Write BK flanking proteins that were matched
    bk_matched_faa = args.output_dir / "bk_flanking_matched.faa"
    with open(bk_matched_faa, 'w') as f:
        for xp_id in bk_xp_ids_matched:
            if xp_id in bk_proteins:
                seq = bk_proteins[xp_id]
                f.write(f">{xp_id}\n{seq}\n")

    n_bk_written = sum(1 for xp in bk_xp_ids_matched if xp in bk_proteins)
    print(f"  Wrote {n_bk_written} matched BK proteins")

    # Run SwissProt on both sets
    print("\nRunning SwissProt annotation...")

    # Helixer flanking
    print("  DIAMOND: Helixer flanking vs SwissProt...")
    helixer_sp_out = args.output_dir / "helixer_swissprot.tsv"
    success = run_diamond_swissprot(
        helixer_cleaned, args.swissprot_db, helixer_sp_out,
        evalue=args.evalue, threads=args.threads
    )
    if not success:
        print("[ERROR] DIAMOND failed for Helixer flanking")
        return 1

    # BK flanking
    print("  DIAMOND: BK flanking vs SwissProt...")
    bk_sp_out = args.output_dir / "bk_swissprot.tsv"
    success = run_diamond_swissprot(
        bk_matched_faa, args.swissprot_db, bk_sp_out,
        evalue=args.evalue, threads=args.threads
    )
    if not success:
        print("[ERROR] DIAMOND failed for BK flanking")
        return 1

    # Parse SwissProt results
    def load_swissprot_results(tsv_path: Path) -> Dict[str, Dict]:
        """Load SwissProt results into dict by query ID."""
        results = {}
        try:
            df = pd.read_csv(tsv_path, sep='\t', header=None,
                           names=['qseqid', 'sseqid', 'pident', 'length', 'evalue', 'bitscore', 'stitle'])
            for _, row in df.iterrows():
                parsed = parse_swissprot_hit(row['sseqid'], row['stitle'])
                parsed['evalue'] = row['evalue']
                parsed['pident'] = row['pident']
                results[row['qseqid']] = parsed
        except pd.errors.EmptyDataError:
            pass
        return results

    helixer_sp = load_swissprot_results(helixer_sp_out)
    bk_sp = load_swissprot_results(bk_sp_out)

    print(f"\nSwissProt hits:")
    print(f"  Helixer flanking: {len(helixer_sp)}")
    print(f"  BK flanking: {len(bk_sp)}")

    # Validate each match
    validation_rows = []

    for _, match in matches_df.iterrows():
        helixer_id = match['helixer_flanking']
        bk_xp = match['bk_xp']
        target_id = match['target_gene_id']
        locus_id = match['locus_id']

        hel_hit = helixer_sp.get(helixer_id, {})
        bk_hit = bk_sp.get(bk_xp, {})

        hel_acc = hel_hit.get('swissprot_accession', '')
        bk_acc = bk_hit.get('swissprot_accession', '')
        hel_name = hel_hit.get('swissprot_protein_name', '')
        bk_name = bk_hit.get('swissprot_protein_name', '')

        # Determine agreement
        if not hel_acc and not bk_acc:
            status = 'both_no_hit'
        elif not hel_acc:
            status = 'helixer_no_hit'
        elif not bk_acc:
            status = 'bk_no_hit'
        elif hel_acc == bk_acc:
            status = 'same_accession'  # Strongest agreement - same protein
        elif hel_name and bk_name and hel_name == bk_name:
            status = 'same_protein_name'  # Same functional annotation
        else:
            status = 'different'  # Potential false positive

        validation_rows.append({
            'target_gene_id': target_id,
            'locus_id': locus_id,
            'helixer_flanking': helixer_id,
            'bk_xp': bk_xp,
            'helixer_swissprot_acc': hel_acc,
            'helixer_swissprot_name': hel_name,
            'bk_swissprot_acc': bk_acc,
            'bk_swissprot_name': bk_name,
            'agreement_status': status,
            'diamond_pident': match.get('pident', ''),
            'diamond_evalue': match.get('evalue', ''),
        })

    validation_df = pd.DataFrame(validation_rows)
    validation_df.to_csv(args.output_dir / "synteny_validation.tsv", sep='\t', index=False)

    # Summary statistics
    print("\n" + "=" * 60)
    print("VALIDATION SUMMARY")
    print("=" * 60)

    status_counts = validation_df['agreement_status'].value_counts()
    total = len(validation_df)

    print(f"\nTotal flanking gene matches: {total}")
    print("\nAgreement status:")
    for status, count in status_counts.items():
        pct = 100 * count / total
        print(f"  {status}: {count} ({pct:.1f}%)")

    # Calculate validation rate (same_accession + same_protein_name as "validated")
    n_validated = status_counts.get('same_accession', 0) + status_counts.get('same_protein_name', 0)
    n_different = status_counts.get('different', 0)
    n_with_both_hits = n_validated + n_different

    if n_with_both_hits > 0:
        validation_rate = 100 * n_validated / n_with_both_hits
        print(f"\nValidation rate (same gene when both have hits): {validation_rate:.1f}%")
        print(f"  ({n_validated} validated / {n_with_both_hits} with both hits)")

    # Summary by target
    target_summary = []
    for target_id in validation_df['target_gene_id'].unique():
        target_matches = validation_df[validation_df['target_gene_id'] == target_id]
        n_matches = len(target_matches)
        n_same = len(target_matches[target_matches['agreement_status'].isin(['same_accession', 'same_protein_name'])])
        n_diff = len(target_matches[target_matches['agreement_status'] == 'different'])

        target_summary.append({
            'target_gene_id': target_id,
            'n_flanking_matches': n_matches,
            'n_validated': n_same,
            'n_different': n_diff,
            'validation_pct': 100 * n_same / n_matches if n_matches > 0 else 0,
        })

    summary_df = pd.DataFrame(target_summary)
    summary_df.to_csv(args.output_dir / "validation_summary.tsv", sep='\t', index=False)

    print(f"\n[OUTPUT] {args.output_dir / 'synteny_validation.tsv'}")
    print(f"[OUTPUT] {args.output_dir / 'validation_summary.tsv'}")

    return 0


if __name__ == "__main__":
    exit(main())
