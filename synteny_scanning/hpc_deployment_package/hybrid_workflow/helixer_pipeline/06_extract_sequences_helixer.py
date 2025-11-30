#!/usr/bin/env python3
"""
Phase 6 (Helixer): Extract and organize target sequences.

Creates per-genome/per-locus directory structure matching the original tblastn pipeline:
    {output_dir}/{genome}/{locus}/gene1_protein.fasta

For Helixer, we extract complete gene predictions directly from the Helixer proteome
(no Exonerate needed since Helixer already provides complete gene models).

Inputs:
- Phase 4 Helixer targets (all_target_loci.tsv)
- Phase 5 Helixer classifications (syntenic_targets.tsv, unplaceable_targets.tsv)
- Helixer proteomes ({genome}_helixer_proteins.faa)
- Phase 1 query proteins (for length comparison)

Outputs:
- Per-genome/per-locus protein FASTA files
- all_extracted_genes.faa (aggregated)
- length_qc.tsv (QC report for unusual lengths)
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def compute_length_flag(query_length_aa: int, extracted_length_aa: int) -> str:
    """
    Classify extracted protein length relative to query.

    Returns: 'ok' (0.8-1.3x), 'short' (<0.8x), 'long' (>1.3x), or 'unknown'
    """
    if query_length_aa <= 0 or extracted_length_aa <= 0:
        return "unknown"
    ratio = extracted_length_aa / query_length_aa
    if 0.8 <= ratio <= 1.3:
        return "ok"
    if ratio < 0.8:
        return "short"
    return "long"


def load_query_protein_lengths(phase1_dir: Path) -> Dict[str, int]:
    """
    Load query protein lengths from Phase 1 target FASTA files.

    Returns: {locus_id: protein_length_aa}
    """
    lengths = {}

    # Find all target FASTA files
    for fasta_file in phase1_dir.glob("*_targets.faa"):
        locus_id = fasta_file.stem.replace("_targets", "")
        for record in SeqIO.parse(fasta_file, 'fasta'):
            # Use first protein in each locus as query length
            lengths[locus_id] = len(str(record.seq).replace('*', '').replace('.', ''))
            break

    return lengths


def load_helixer_proteome(helixer_dir: Path, genome: str) -> Dict[str, str]:
    """Load Helixer protein sequences for a genome."""
    proteome_file = helixer_dir / genome / f"{genome}_helixer_proteins.faa"

    if not proteome_file.exists():
        # Try alternate naming
        alt_file = helixer_dir / genome / "helixer_proteins.faa"
        if alt_file.exists():
            proteome_file = alt_file
        else:
            return {}

    sequences = {}
    for record in SeqIO.parse(proteome_file, 'fasta'):
        sequences[record.id] = str(record.seq)

    return sequences


def get_protein_coordinates(gff_file: Path, gene_id: str) -> Optional[Tuple[str, int, int, str]]:
    """
    Get coordinates for a gene from GFF3 file.

    Returns: (scaffold, start, end, strand) or None
    """
    if not gff_file.exists():
        return None

    with open(gff_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            feature_type = parts[2]
            if feature_type != 'gene':
                continue

            attrs = parts[8]
            # Check if this gene matches
            if f"ID={gene_id}" in attrs or f"Name={gene_id}" in attrs:
                scaffold = parts[0]
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                return scaffold, start, end, strand

    return None


def write_protein_fasta(
    output_file: Path,
    gene_id: str,
    protein_seq: str,
    scaffold: str,
    start: int,
    end: int,
    strand: str,
    query_length_aa: int = 0,
) -> Tuple[int, str]:
    """
    Write protein FASTA with metadata header matching original pipeline format.

    Returns: (extracted_length_aa, length_flag)
    """
    length_aa = len(protein_seq.replace('*', '').replace('.', ''))

    # Calculate coverage and length flag if query length provided
    if query_length_aa > 0:
        coverage = (length_aa / query_length_aa) * 100
        cov_str = f" cov:{coverage:.1f}%"
        length_flag = compute_length_flag(query_length_aa, length_aa)
        flag_str = f" len_flag:{length_flag}"
    else:
        cov_str = ""
        length_flag = "unknown"
        flag_str = ""

    header = (
        f">gene1 {scaffold}:{start}-{end} "
        f"strand:{strand} "
        f"length:{length_aa}aa"
        f"{cov_str}"
        f"{flag_str}"
    )

    # Clean sequence
    clean_seq = protein_seq.replace('.', 'X').replace('*', '')

    with open(output_file, 'w') as f:
        f.write(f"{header}\n{clean_seq}\n")

    return length_aa, length_flag


def main():
    parser = argparse.ArgumentParser(
        description="Phase 6 (Helixer): Extract and organize target sequences"
    )
    parser.add_argument("--family", required=True, help="Gene family name")
    parser.add_argument("--phase1-dir", type=Path,
                        help="Phase 1 directory with query proteins (for length QC)")
    parser.add_argument("--phase4-dir", type=Path, required=True,
                        help="Phase 4 Helixer output directory")
    parser.add_argument("--phase5-dir", type=Path, required=True,
                        help="Phase 5 Helixer output directory")
    parser.add_argument("--helixer-dir", type=Path, required=True,
                        help="Directory with Helixer proteomes")
    parser.add_argument("--output-dir", type=Path, required=True,
                        help="Output directory")
    parser.add_argument("--verbose", action="store_true")

    args = parser.parse_args()

    print("=" * 80)
    print("PHASE 6 (HELIXER): EXTRACT AND ORGANIZE SEQUENCES")
    print("=" * 80)
    print()
    print(f"Family: {args.family}")
    print(f"Phase 1 dir: {args.phase1_dir}")
    print(f"Phase 4 dir: {args.phase4_dir}")
    print(f"Phase 5 dir: {args.phase5_dir}")
    print(f"Output dir: {args.output_dir}")
    print()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Load query protein lengths from Phase 1 for length QC
    query_lengths = {}
    if args.phase1_dir and args.phase1_dir.exists():
        query_lengths = load_query_protein_lengths(args.phase1_dir)
        print(f"Loaded query lengths for {len(query_lengths)} loci")
    else:
        print("No Phase 1 dir provided - length QC will be skipped")

    # Load Phase 4 targets for coordinate info
    targets_file = args.phase4_dir / "all_target_loci.tsv"
    if not targets_file.exists():
        print(f"[ERROR] Missing Phase 4 targets: {targets_file}")
        return 1

    targets_df = pd.read_csv(targets_file, sep='\t')
    print(f"Loaded {len(targets_df)} targets from Phase 4")

    # Create lookup for target coordinates and query info
    target_info = {}
    for _, row in targets_df.iterrows():
        target_id = row['target_gene_id']
        target_info[target_id] = {
            'scaffold': row['scaffold'],
            'start': int(row['start']),
            'end': int(row['end']),
            'strand': row['strand'],
            'query_id': row.get('query_id', ''),
            'parent_locus': row.get('parent_locus', ''),
        }

    # Track QC data for length report
    qc_records: List[Dict] = []

    # Load Phase 5 classifications
    syntenic_file = args.phase5_dir / "syntenic_targets.tsv"
    unplaceable_file = args.phase5_dir / "unplaceable_targets.tsv"

    syntenic_df = pd.DataFrame()
    unplaceable_df = pd.DataFrame()

    if syntenic_file.exists() and syntenic_file.stat().st_size > 0:
        try:
            syntenic_df = pd.read_csv(syntenic_file, sep='\t')
            print(f"Loaded {len(syntenic_df)} syntenic targets")
        except pd.errors.EmptyDataError:
            print("  0 syntenic targets (empty file)")

    if unplaceable_file.exists() and unplaceable_file.stat().st_size > 0:
        try:
            unplaceable_df = pd.read_csv(unplaceable_file, sep='\t')
            print(f"Loaded {len(unplaceable_df)} unplaceable targets")
        except pd.errors.EmptyDataError:
            print("  0 unplaceable targets (empty file)")

    # Track extraction stats
    extracted_count = 0
    failed_count = 0
    genome_proteomes = {}  # Cache loaded proteomes

    # Process syntenic targets
    print("\n[1] Extracting syntenic targets...")

    for _, row in syntenic_df.iterrows():
        genome = row['genome']
        target_id = row['target_gene_id']
        assigned_locus = row.get('assigned_to', row.get('locus_id', row.get('parent_locus', 'unknown')))
        parent_locus = row.get('parent_locus', assigned_locus)

        # Get target info
        info = target_info.get(target_id, {})
        scaffold = info.get('scaffold', row.get('scaffold', 'unknown'))
        start = info.get('start', row.get('start', 0))
        end = info.get('end', row.get('end', 0))
        strand = info.get('strand', row.get('strand', '+'))

        # Get query length for this locus
        query_length = query_lengths.get(parent_locus, 0)

        # Create output directory: {output}/{genome}/{assigned_locus}/
        target_dir = args.output_dir / genome / assigned_locus
        target_dir.mkdir(parents=True, exist_ok=True)

        # Load proteome if not cached
        if genome not in genome_proteomes:
            genome_proteomes[genome] = load_helixer_proteome(args.helixer_dir, genome)

        proteome = genome_proteomes[genome]

        # Extract protein sequence
        # Target ID format: {genome}_{scaffold}_{gene_index}
        # Helixer proteome IDs may have .1 suffix: {genome}_{scaffold}_{gene_index}.1
        protein_seq = proteome.get(target_id)

        if not protein_seq:
            # Try with .1 suffix (common for Helixer proteins)
            protein_seq = proteome.get(f"{target_id}.1")

        if not protein_seq:
            # Try .2 suffix
            protein_seq = proteome.get(f"{target_id}.2")

        if not protein_seq:
            if args.verbose:
                print(f"  [WARNING] No protein for {target_id} in {genome}")
            failed_count += 1
            continue

        # Write protein FASTA
        output_file = target_dir / f"{assigned_locus}_gene1_protein.fasta"
        extracted_length, length_flag = write_protein_fasta(
            output_file, target_id, protein_seq,
            scaffold, start, end, strand,
            query_length_aa=query_length
        )
        extracted_count += 1

        # Track QC data
        qc_records.append({
            'genome': genome,
            'target_id': target_id,
            'locus': assigned_locus,
            'placement': 'syntenic',
            'query_length_aa': query_length,
            'extracted_length_aa': extracted_length,
            'length_flag': length_flag,
            'ratio': round(extracted_length / query_length, 2) if query_length > 0 else 0,
        })

        if args.verbose:
            print(f"  {genome}/{assigned_locus}: {extracted_length}aa ({length_flag})")

    # Process unplaceable targets
    print("\n[2] Extracting unplaceable targets...")

    for _, row in unplaceable_df.iterrows():
        genome = row['genome']
        target_id = row['target_gene_id']
        parent_locus = row.get('parent_locus', 'unknown')

        # Get target info
        info = target_info.get(target_id, {})
        scaffold = info.get('scaffold', row.get('scaffold', 'unknown'))
        start = info.get('start', row.get('start', 0))
        end = info.get('end', row.get('end', 0))
        strand = info.get('strand', row.get('strand', '+'))

        # Get query length for this locus
        query_length = query_lengths.get(parent_locus, 0)

        # Create unique tag for unplaceable
        unique_tag = f"{scaffold}_{start}_{end}"

        # Create output directory: {output}/{genome}/UNPLACED/{unique_tag}/
        target_dir = args.output_dir / genome / "UNPLACED" / unique_tag
        target_dir.mkdir(parents=True, exist_ok=True)

        # Load proteome if not cached
        if genome not in genome_proteomes:
            genome_proteomes[genome] = load_helixer_proteome(args.helixer_dir, genome)

        proteome = genome_proteomes[genome]

        # Extract protein sequence (try with .1 suffix)
        protein_seq = proteome.get(target_id)
        if not protein_seq:
            protein_seq = proteome.get(f"{target_id}.1")
        if not protein_seq:
            protein_seq = proteome.get(f"{target_id}.2")

        if not protein_seq:
            if args.verbose:
                print(f"  [WARNING] No protein for {target_id} in {genome}")
            failed_count += 1
            continue

        # Write protein FASTA
        output_file = target_dir / f"{unique_tag}_gene1_protein.fasta"
        extracted_length, length_flag = write_protein_fasta(
            output_file, target_id, protein_seq,
            scaffold, start, end, strand,
            query_length_aa=query_length
        )
        extracted_count += 1

        # Track QC data
        qc_records.append({
            'genome': genome,
            'target_id': target_id,
            'locus': f"UNPLACED/{unique_tag}",
            'placement': 'unplaceable',
            'query_length_aa': query_length,
            'extracted_length_aa': extracted_length,
            'length_flag': length_flag,
            'ratio': round(extracted_length / query_length, 2) if query_length > 0 else 0,
        })

        if args.verbose:
            print(f"  {genome}/UNPLACED/{unique_tag}: {extracted_length}aa ({length_flag})")

    # Aggregate all proteins
    print("\n[3] Aggregating extracted proteins...")

    all_proteins_file = args.output_dir / "all_extracted_genes.faa"
    protein_count = 0

    with open(all_proteins_file, 'w') as outfile:
        for genome_dir in args.output_dir.iterdir():
            if not genome_dir.is_dir() or genome_dir.name.startswith('.'):
                continue

            for target_dir in genome_dir.iterdir():
                if not target_dir.is_dir():
                    continue

                # Handle UNPLACED subdirectories
                if target_dir.name == 'UNPLACED':
                    for unplaced_dir in target_dir.iterdir():
                        if unplaced_dir.is_dir():
                            for protein_file in unplaced_dir.glob("*_protein.fasta"):
                                for record in SeqIO.parse(protein_file, 'fasta'):
                                    SeqIO.write(record, outfile, 'fasta')
                                    protein_count += 1
                else:
                    for protein_file in target_dir.glob("*_protein.fasta"):
                        for record in SeqIO.parse(protein_file, 'fasta'):
                            SeqIO.write(record, outfile, 'fasta')
                            protein_count += 1

    print(f"  Aggregated {protein_count} proteins to {all_proteins_file.name}")

    # Write QC report
    print("\n[4] Generating length QC report...")

    if qc_records:
        qc_df = pd.DataFrame(qc_records)
        qc_file = args.output_dir / "length_qc.tsv"
        qc_df.to_csv(qc_file, sep='\t', index=False)

        # Calculate QC statistics
        n_ok = len(qc_df[qc_df['length_flag'] == 'ok'])
        n_short = len(qc_df[qc_df['length_flag'] == 'short'])
        n_long = len(qc_df[qc_df['length_flag'] == 'long'])
        n_unknown = len(qc_df[qc_df['length_flag'] == 'unknown'])

        print(f"  Length QC summary:")
        print(f"    ok (0.8-1.3x):  {n_ok}")
        print(f"    short (<0.8x):  {n_short}")
        print(f"    long (>1.3x):   {n_long}")
        print(f"    unknown:        {n_unknown}")
        print(f"  Saved to: {qc_file.name}")

        # Flag unusual lengths for review
        if n_short > 0 or n_long > 0:
            unusual = qc_df[qc_df['length_flag'].isin(['short', 'long'])].sort_values('ratio')
            unusual_file = args.output_dir / "length_qc_flagged.tsv"
            unusual.to_csv(unusual_file, sep='\t', index=False)
            print(f"\n  [QC] {len(unusual)} targets flagged for unusual length:")
            print(f"       Saved to: {unusual_file.name}")

            # Show extreme examples
            if n_short > 0:
                shortest = unusual[unusual['length_flag'] == 'short'].head(3)
                for _, r in shortest.iterrows():
                    print(f"       SHORT: {r['genome']}/{r['locus']} - {r['extracted_length_aa']}aa vs {r['query_length_aa']}aa query ({r['ratio']}x)")

            if n_long > 0:
                longest = unusual[unusual['length_flag'] == 'long'].tail(3)
                for _, r in longest.iterrows():
                    print(f"       LONG:  {r['genome']}/{r['locus']} - {r['extracted_length_aa']}aa vs {r['query_length_aa']}aa query ({r['ratio']}x)")

    # Summary
    print("\n" + "=" * 80)
    print("PHASE 6 (HELIXER) COMPLETE")
    print("=" * 80)
    print(f"\nExtracted: {extracted_count} targets")
    print(f"Failed: {failed_count} targets")
    print(f"Aggregated proteins: {protein_count}")
    print(f"\nOutput directory: {args.output_dir}")

    return 0


if __name__ == "__main__":
    exit(main())
