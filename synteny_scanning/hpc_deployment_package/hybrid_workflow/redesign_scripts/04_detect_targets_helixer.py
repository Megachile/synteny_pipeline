#!/usr/bin/env python3
"""
Phase 4 (Helixer version): DIAMOND protein search against Helixer-predicted proteomes.

Key differences vs tblastn version:
- Uses DIAMOND blastp instead of tblastn (much faster)
- Gene coordinates come from Helixer GFF3 instead of HSP clustering
- No complex HSP merging logic needed

Outputs (per family):
- all_target_loci.tsv       # per-(query,gene)-level hits (same format as tblastn)
- diamond_results/*.tsv     # raw DIAMOND output per genome
"""

from __future__ import annotations

import argparse
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import pandas as pd
import re


def clean_fasta_for_diamond(input_faa: Path, output_faa: Path) -> int:
    """
    Clean FASTA file for DIAMOND by removing/replacing invalid characters.

    DIAMOND doesn't accept '.' in sequences. Some Helixer outputs have:
    - '.' (possibly stop codons or gaps)
    - Other invalid characters

    Returns number of sequences written.
    """
    valid_aa = set('ACDEFGHIKLMNPQRSTVWYX*')
    n_written = 0
    n_skipped = 0

    with open(input_faa) as f_in, open(output_faa, 'w') as f_out:
        current_header = None
        current_seq = []

        for line in f_in:
            line = line.strip()
            if line.startswith('>'):
                # Write previous sequence if valid
                if current_header:
                    seq = ''.join(current_seq)
                    # Replace '.' with 'X' (unknown amino acid)
                    seq = seq.replace('.', 'X')
                    # Check if sequence has valid characters
                    if len(seq) > 0 and all(c.upper() in valid_aa for c in seq):
                        f_out.write(f"{current_header}\n")
                        # Write sequence in 60-char lines
                        for i in range(0, len(seq), 60):
                            f_out.write(f"{seq[i:i+60]}\n")
                        n_written += 1
                    else:
                        n_skipped += 1

                current_header = line
                current_seq = []
            else:
                current_seq.append(line)

        # Write last sequence
        if current_header:
            seq = ''.join(current_seq)
            seq = seq.replace('.', 'X')
            if len(seq) > 0 and all(c.upper() in valid_aa for c in seq):
                f_out.write(f"{current_header}\n")
                for i in range(0, len(seq), 60):
                    f_out.write(f"{seq[i:i+60]}\n")
                n_written += 1
            else:
                n_skipped += 1

    if n_skipped > 0:
        print(f"    Cleaned FASTA: {n_written} sequences written, {n_skipped} skipped (invalid chars)")

    return n_written


def parse_helixer_gff3(gff3_path: Path) -> Dict[str, Dict]:
    """
    Parse Helixer GFF3 to extract gene coordinates.

    Returns dict: gene_id -> {scaffold, strand, start, end, mrna_id}
    """
    genes = {}

    with open(gff3_path) as f:
        for line in f:
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            scaffold, source, feature, start, end, score, strand, phase, attributes = parts

            # Only extract gene features (not mRNA, CDS, etc.)
            if feature != 'gene':
                continue

            # Parse attributes to get gene ID
            attr_dict = {}
            for attr in attributes.split(';'):
                if '=' in attr:
                    key, value = attr.split('=', 1)
                    attr_dict[key] = value

            gene_id = attr_dict.get('ID', '')
            if not gene_id:
                continue

            genes[gene_id] = {
                'scaffold': scaffold,
                'strand': strand,
                'start': int(start),
                'end': int(end),
            }

    return genes


def protein_to_gene_id(protein_id: str) -> str:
    """
    Convert protein ID (e.g., Aulacidea_tavakolii_CM021338.1_RagTag_000001.1)
    to gene ID (e.g., Aulacidea_tavakolii_CM021338.1_RagTag_000001).

    Helixer protein IDs have .1 (transcript number) suffix; gene IDs don't.
    """
    # Remove the last .N suffix (transcript number)
    if '.' in protein_id:
        parts = protein_id.rsplit('.', 1)
        if parts[1].isdigit():
            return parts[0]
    return protein_id


def build_diamond_db(proteome_faa: Path, db_path: Path, clean_dir: Path = None) -> bool:
    """
    Build DIAMOND database from Helixer proteome FASTA.

    Automatically cleans FASTA to handle invalid characters from Helixer.
    """
    db_path.parent.mkdir(parents=True, exist_ok=True)

    # Clean FASTA file first (handle '.' and other invalid chars)
    if clean_dir is None:
        clean_dir = db_path.parent / "cleaned_fasta"
    clean_dir.mkdir(parents=True, exist_ok=True)
    cleaned_faa = clean_dir / proteome_faa.name

    if not cleaned_faa.exists():
        print(f"    Cleaning FASTA for DIAMOND...")
        n_seqs = clean_fasta_for_diamond(proteome_faa, cleaned_faa)
        if n_seqs == 0:
            print(f"  [ERROR] No valid sequences after cleaning")
            return False
        print(f"    Cleaned: {n_seqs} sequences ready")

    cmd = [
        "diamond", "makedb",
        "--in", str(cleaned_faa),
        "--db", str(db_path),
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  [ERROR] DIAMOND makedb failed: {result.stderr}")
        return False
    return True


def run_diamond_blastp(
    query_faa: Path,
    db_path: Path,
    output_tsv: Path,
    *,
    evalue: float = 1e-5,
    max_targets: int = 50,
    threads: int = 4,
) -> bool:
    """Run DIAMOND blastp against Helixer proteome database."""

    cmd = [
        "diamond", "blastp",
        "--query", str(query_faa),
        "--db", str(db_path),
        "--out", str(output_tsv),
        "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "mismatch",
                        "gapopen", "qstart", "qend", "sstart", "send", "evalue",
                        "bitscore", "qlen", "slen",
        "--max-target-seqs", str(max_targets),
        "--evalue", str(evalue),
        "--threads", str(threads),
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  [ERROR] DIAMOND blastp failed: {result.stderr}")
        return False
    return True


def parse_diamond_results(tsv_path: Path) -> pd.DataFrame:
    """Parse DIAMOND output TSV into DataFrame."""
    columns = [
        'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
        'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen'
    ]

    if not tsv_path.exists() or tsv_path.stat().st_size == 0:
        return pd.DataFrame(columns=columns)

    return pd.read_csv(tsv_path, sep='\t', names=columns, header=None)


def aggregate_hits_per_gene(
    diamond_df: pd.DataFrame,
    gff3_data: Dict[str, Dict],
    genome_id: str,
    query_to_locus: Dict[str, str],
    gene_family: str,
) -> List[Dict]:
    """
    Aggregate DIAMOND hits per target gene.

    For each unique (query_id, target_gene) pair:
    - Get coordinates from GFF3
    - Calculate query coverage
    - Keep best e-value/bitscore
    """
    if diamond_df.empty:
        return []

    results = []

    # Group by query and target gene
    for (qseqid, sseqid), group in diamond_df.groupby(['qseqid', 'sseqid']):
        # Convert protein ID to gene ID
        gene_id = protein_to_gene_id(sseqid)

        # Look up gene coordinates from GFF3
        if gene_id not in gff3_data:
            print(f"    Warning: Gene {gene_id} not found in GFF3")
            continue

        gene_info = gff3_data[gene_id]

        # Get query length from first hit
        query_len = group['qlen'].iloc[0]

        # Calculate query coverage (union of all HSPs)
        query_spans = [(row['qstart'], row['qend']) for _, row in group.iterrows()]
        query_covered = calculate_union_coverage(query_spans)
        query_coverage = query_covered / query_len if query_len > 0 else 0

        # Get best hit stats
        best_idx = group['evalue'].idxmin()
        best_evalue = group.loc[best_idx, 'evalue']
        best_bitscore = group.loc[best_idx, 'bitscore']

        # Map query to locus
        # Query format: XP_033208249.1|LOC117167432
        parent_locus = query_to_locus.get(qseqid, '')
        if not parent_locus:
            # Try extracting LOC ID from after the pipe
            if '|' in qseqid:
                loc_id = qseqid.split('|')[1]
                parent_locus = query_to_locus.get(loc_id, '')
        if not parent_locus:
            # Try base query ID
            base_qid = qseqid.split('|')[0]
            parent_locus = query_to_locus.get(base_qid, 'unknown')

        # Calculate span
        span_kb = (gene_info['end'] - gene_info['start'] + 1) / 1000.0

        results.append({
            'genome': genome_id,
            'scaffold': gene_info['scaffold'],
            'strand': gene_info['strand'],
            'start': gene_info['start'],
            'end': gene_info['end'],
            'span_kb': span_kb,
            'num_hsps': len(group),
            'best_evalue': best_evalue,
            'best_bitscore': best_bitscore,
            'query_id': qseqid,
            'parent_locus': parent_locus,
            'query_start_min': group['qstart'].min(),
            'query_end_max': group['qend'].max(),
            'query_span_coverage': (group['qend'].max() - group['qstart'].min() + 1) / query_len if query_len > 0 else 0,
            'combined_query_coverage': query_coverage,
            'gene_family': gene_family,
            'target_gene_id': gene_id,  # Extra column for Helixer version
        })

    return results


def calculate_union_coverage(spans: List[Tuple[int, int]]) -> int:
    """Calculate union coverage of query spans (handles overlaps)."""
    if not spans:
        return 0

    # Sort by start position
    sorted_spans = sorted(spans)

    # Merge overlapping spans
    merged = [sorted_spans[0]]
    for start, end in sorted_spans[1:]:
        if start <= merged[-1][1] + 1:  # Overlapping or adjacent
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))

    # Sum covered positions
    return sum(end - start + 1 for start, end in merged)


def load_query_to_locus_map(phase1_dir: Path) -> Dict[str, str]:
    """
    Load mapping from query protein ID to locus ID.

    Query IDs have format: XP_033208249.1|LOC117167432
    locus_definitions.tsv has: locus_id, cluster_members (LOC IDs)

    Maps both full query ID and LOC ID to locus.
    """
    mapping = {}

    # Primary source: locus_definitions.tsv
    locus_defs = phase1_dir / "locus_definitions.tsv"
    if locus_defs.exists():
        df = pd.read_csv(locus_defs, sep='\t')
        if 'locus_id' in df.columns and 'cluster_members' in df.columns:
            for _, row in df.iterrows():
                locus_id = str(row['locus_id'])
                members = str(row['cluster_members']).split(';')
                for member in members:
                    member = member.strip()
                    if member:
                        # Map LOC ID to locus
                        mapping[member] = locus_id
                        # Also map with LOC prefix variations
                        if not member.startswith('LOC'):
                            mapping[f'LOC{member}'] = locus_id

    # Also try cluster_members.tsv
    members_file = phase1_dir / "cluster_members.tsv"
    if members_file.exists():
        df = pd.read_csv(members_file, sep='\t')
        if 'protein_id' in df.columns and 'locus_id' in df.columns:
            for _, row in df.iterrows():
                mapping[str(row['protein_id'])] = str(row['locus_id'])

    # Also try flanking_genes.tsv
    flanking_file = phase1_dir / "flanking_genes.tsv"
    if flanking_file.exists():
        df = pd.read_csv(flanking_file, sep='\t')
        if 'protein_id' in df.columns and 'locus_id' in df.columns:
            for _, row in df.iterrows():
                mapping[str(row['protein_id'])] = str(row['locus_id'])

    return mapping


def main():
    parser = argparse.ArgumentParser(
        description="Phase 4 (Helixer): DIAMOND target detection against Helixer proteomes"
    )
    parser.add_argument(
        "--family", required=True,
        help="Gene family name (e.g., ferritin_MC102)"
    )
    parser.add_argument(
        "--phase1-dir", type=Path, required=True,
        help="Path to phase1_v2 directory with query proteins"
    )
    parser.add_argument(
        "--helixer-dir", type=Path, required=True,
        help="Directory containing Helixer outputs (proteome FASTA and GFF3)"
    )
    parser.add_argument(
        "--genome-list", type=Path, required=True,
        help="TSV file with genome_id column listing genomes to process"
    )
    parser.add_argument(
        "--output-dir", type=Path, required=True,
        help="Output directory for Phase 4 Helixer results"
    )
    parser.add_argument(
        "--diamond-db-dir", type=Path, default=None,
        help="Directory for DIAMOND databases (default: output-dir/diamond_dbs)"
    )
    parser.add_argument(
        "--evalue", type=float, default=1e-5,
        help="E-value threshold for DIAMOND (default: 1e-5)"
    )
    parser.add_argument(
        "--max-targets", type=int, default=50,
        help="Max target sequences per query (default: 50)"
    )
    parser.add_argument(
        "--threads", type=int, default=4,
        help="Number of threads for DIAMOND (default: 4)"
    )

    args = parser.parse_args()

    print("=" * 80)
    print("PHASE 4 (HELIXER): DIAMOND TARGET DETECTION")
    print("=" * 80)
    print()

    # Setup directories
    args.output_dir.mkdir(parents=True, exist_ok=True)
    diamond_db_dir = args.diamond_db_dir or (args.output_dir / "diamond_dbs")
    diamond_db_dir.mkdir(parents=True, exist_ok=True)
    diamond_results_dir = args.output_dir / "diamond_results"
    diamond_results_dir.mkdir(parents=True, exist_ok=True)

    # Load query proteins
    query_faa = args.phase1_dir / "query_proteins.faa"
    if not query_faa.exists():
        # Fallback to combined_targets.faa
        query_faa = args.phase1_dir.parent / "phase4_v2" / "combined_targets.faa"

    if not query_faa.exists():
        print(f"[ERROR] Query proteins not found: {query_faa}")
        return 1

    print(f"Query proteins: {query_faa}")

    # Load query-to-locus mapping
    query_to_locus = load_query_to_locus_map(args.phase1_dir)
    print(f"Loaded {len(query_to_locus)} query-to-locus mappings")

    # Load genome list
    if args.genome_list.exists():
        genome_df = pd.read_csv(args.genome_list, sep='\t')
        if 'genome_id' in genome_df.columns:
            genomes = genome_df['genome_id'].tolist()
        else:
            genomes = genome_df.iloc[:, 0].tolist()
    else:
        # Auto-detect from helixer directory
        genomes = [d.name for d in args.helixer_dir.iterdir()
                   if d.is_dir() and (d / f"{d.name}_helixer_proteins.faa").exists()]

    print(f"Processing {len(genomes)} genomes")
    print()

    # Process each genome
    all_results = []

    for genome_id in genomes:
        print(f"\n[{genome_id}]")

        # Find Helixer files
        genome_dir = args.helixer_dir / genome_id
        proteome_faa = genome_dir / f"{genome_id}_helixer_proteins.faa"
        gff3_file = genome_dir / f"{genome_id}_helixer.gff3"

        if not proteome_faa.exists():
            print(f"  Skipping: no Helixer proteome found")
            continue

        if not gff3_file.exists():
            print(f"  Skipping: no Helixer GFF3 found")
            continue

        # Build DIAMOND database if needed
        db_path = diamond_db_dir / f"{genome_id}_helixer"
        if not (db_path.parent / f"{db_path.name}.dmnd").exists():
            print(f"  Building DIAMOND database...")
            if not build_diamond_db(proteome_faa, db_path):
                continue

        # Run DIAMOND blastp
        diamond_out = diamond_results_dir / f"{genome_id}.tsv"
        if not diamond_out.exists():
            print(f"  Running DIAMOND blastp...")
            if not run_diamond_blastp(
                query_faa, db_path, diamond_out,
                evalue=args.evalue,
                max_targets=args.max_targets,
                threads=args.threads,
            ):
                continue

        # Parse results
        diamond_df = parse_diamond_results(diamond_out)
        if diamond_df.empty:
            print(f"  No DIAMOND hits")
            continue

        print(f"  DIAMOND hits: {len(diamond_df)}")

        # Parse GFF3
        print(f"  Parsing GFF3...")
        gff3_data = parse_helixer_gff3(gff3_file)
        print(f"  Genes in GFF3: {len(gff3_data)}")

        # Aggregate hits per gene
        genome_results = aggregate_hits_per_gene(
            diamond_df, gff3_data, genome_id, query_to_locus, args.family
        )

        print(f"  Target genes: {len(genome_results)}")
        all_results.extend(genome_results)

    # Write output
    if all_results:
        output_df = pd.DataFrame(all_results)

        # Reorder columns to match tblastn format
        column_order = [
            'genome', 'scaffold', 'strand', 'start', 'end', 'span_kb',
            'num_hsps', 'best_evalue', 'best_bitscore', 'query_id', 'parent_locus',
            'query_start_min', 'query_end_max', 'query_span_coverage',
            'combined_query_coverage', 'gene_family', 'target_gene_id'
        ]
        output_df = output_df[[c for c in column_order if c in output_df.columns]]

        output_file = args.output_dir / "all_target_loci.tsv"
        output_df.to_csv(output_file, sep='\t', index=False)
        print(f"\n[OUTPUT] Wrote {len(output_df)} target genes to {output_file}")
    else:
        print("\n[WARNING] No target genes found!")

    print("\nPhase 4 (Helixer) complete.")
    return 0


if __name__ == "__main__":
    exit(main())
