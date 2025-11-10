#!/usr/bin/env python3
"""
Phase 1 REDESIGNED: tblastn-first locus discovery.

**Problem Fixed:**
Previous workflow used GFF annotation coordinates as starting point, but BRAKER3
annotations can be misplaced by 600kb+. This caused downstream failures when
tblastn found the ACTUAL target location far from the synteny block.

**New Approach:**
1. Run tblastn on genome FIRST → precise target location (ground truth)
2. Define flanking window around THAT location (±500kb)
3. Query GFF for genes within that window
4. Extract flanking proteins

**Key Principle:**
"Build locus around where target ACTUALLY IS, not where annotation THINKS it is"

Usage:
    python 01_discover_paralogs_tblastn_first.py \\
        --query-protein XP_015672091.1 \\
        --output-dir outputs/ferritin_phase1_tblastn \\
        --gene-family ferritin_MC102 \\
        --window-size 500 \\
        --evalue 1e-10 \\
        --min-coverage 0.5 \\
        --min-identity 40
"""

import sys
import re
import subprocess
import tempfile
from pathlib import Path
from collections import defaultdict
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import json
from datetime import datetime
import argparse

# Genome configuration
GENOMES = {
    'BK': {
        'name': 'Belonocnema kinseyi',
        'accession': 'GCF_010883055.1',
        'genome': 'data/ragtag_output/GCA_010883055.1/ragtag.scaffold.fasta',
        'proteome': 'data/proteomes/GCF_010883055.1_Bkinseyi_NCBI.faa',
        'gff': 'data/genomes/GCF_010883055.1_Bkinseyi.gff',
        'type': 'NCBI'
    },
    'LB': {
        'name': 'Leptopilina boulardi',
        'accession': 'GCF_019393585.1',
        'genome': 'data/ragtag_output/GCA_019393585.1/ragtag.scaffold.fasta',
        'proteome': 'data/proteomes/GCF_019393585.1_Lboulardi_NCBI.faa',
        'gff': 'data/genomes/GCF_019393585.1_Lboulardi.gff',
        'type': 'NCBI'
    },
    'TR': {
        'name': 'Torymus regularis',
        'accession': 'GCA_020615435.1',
        'genome': 'data/ragtag_output/GCA_020615435.1/ragtag.scaffold.fasta',
        'proteome': 'data/proteomes/GCA_020615435.1.faa',
        'gff': 'data/genomes/GCA_020615435.1.gff3',
        'type': 'BRAKER3'
    },
    'DR': {
        'name': 'Druon quercuslanigerum',
        'accession': 'GCA_030998225.1',
        'genome': 'data/ragtag_output/GCA_030998225.1/ragtag.scaffold.fasta',
        'proteome': 'data/proteomes/GCA_030998225.1.faa',
        'gff': 'data/genomes/GCA_030998225.1.gff3',
        'type': 'BRAKER3'
    }
}

# Default parameters (can be overridden by command line)
WINDOW_SIZE_KB = 500  # Flanking window: ±500kb from target
EVALUE_THRESHOLD = 1e-5  # E-value for target detection (relaxed for tblastn)
MIN_COVERAGE = 0.25  # Hit must cover 25% of query (tblastn is fragmented)
MIN_IDENTITY = 35.0  # Minimum % identity (relaxed for tblastn)


def run_tblastn(query_protein, genome_fasta, evalue=1e-10):
    """
    Run tblastn to find protein query in genome.

    Returns:
        pandas DataFrame with hits (columns: qseqid, sseqid, pident, length,
        mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qlen, slen)
    """

    print(f"  Running tblastn on {genome_fasta}...")

    # Create temp files
    with tempfile.NamedTemporaryFile(mode='w', suffix='.faa', delete=False) as qf:
        SeqIO.write([query_protein], qf, 'fasta')
        query_file = qf.name

    with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as of:
        output_file = of.name

    try:
        # Run tblastn
        cmd = [
            'tblastn',
            '-query', query_file,
            '-subject', genome_fasta,
            '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen',
            '-evalue', str(evalue),
            '-max_target_seqs', '100'
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, check=True)

        # Write output
        with open(output_file, 'w') as f:
            f.write(result.stdout)

        # Parse results
        if result.stdout.strip():
            hits = pd.read_csv(output_file, sep='\t', header=None,
                             names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch',
                                   'gapopen', 'qstart', 'qend', 'sstart', 'send',
                                   'evalue', 'bitscore', 'qlen', 'slen'])
            print(f"    Found {len(hits)} hits")
            return hits
        else:
            print(f"    No hits found")
            return pd.DataFrame()

    finally:
        # Cleanup
        Path(query_file).unlink(missing_ok=True)
        Path(output_file).unlink(missing_ok=True)


def filter_tblastn_hits(hits, min_coverage=0.5, min_identity=40.0):
    """
    Filter tblastn hits by coverage and identity.

    Returns:
        Filtered DataFrame with additional 'coverage' column
    """

    if hits.empty:
        return hits

    # Calculate query coverage
    hits['coverage'] = (hits['qend'] - hits['qstart'] + 1) / hits['qlen']

    # Filter
    filtered = hits[
        (hits['coverage'] >= min_coverage) &
        (hits['pident'] >= min_identity)
    ].copy()

    print(f"    After filtering (cov≥{min_coverage}, id≥{min_identity}%): {len(filtered)} hits")

    return filtered


def parse_hit_coordinates(hit):
    """
    Parse hit coordinates and determine strand.

    Returns:
        dict with scaffold, start, end, strand, score
    """

    scaffold = hit['sseqid']
    sstart = int(hit['sstart'])
    send = int(hit['send'])

    # Determine strand and fix coordinates
    if sstart <= send:
        strand = '+'
        start = sstart
        end = send
    else:
        strand = '-'
        start = send
        end = sstart

    return {
        'scaffold': scaffold,
        'start': start,
        'end': end,
        'strand': strand,
        'evalue': float(hit['evalue']),
        'bitscore': float(hit['bitscore']),
        'pident': float(hit['pident']),
        'coverage': float(hit['coverage']) if 'coverage' in hit else 0.0
    }


def define_flanking_window(target_coords, window_size_kb=500):
    """
    Define flanking window around target.

    Args:
        target_coords: dict with scaffold, start, end, strand
        window_size_kb: size of flanking window in kb

    Returns:
        dict with window_start, window_end
    """

    window_bp = window_size_kb * 1000

    # Window extends from both sides
    window_start = max(1, target_coords['start'] - window_bp)
    window_end = target_coords['end'] + window_bp

    return {
        'scaffold': target_coords['scaffold'],
        'window_start': window_start,
        'window_end': window_end,
        'target_start': target_coords['start'],
        'target_end': target_coords['end'],
        'target_strand': target_coords['strand'],
        'window_size_kb': window_size_kb
    }


def parse_gff_genes_in_window(gff_file, window, annotation_type='BRAKER3'):
    """
    Parse GFF to find genes within window.

    Args:
        gff_file: path to GFF file
        window: dict from define_flanking_window()
        annotation_type: 'NCBI' or 'BRAKER3'

    Returns:
        list of dicts: [{gene_id, protein_id, start, end, strand, distance_from_target}, ...]
    """

    print(f"  Parsing {annotation_type} GFF for genes in window...")
    print(f"    Window: {window['scaffold']}:{window['window_start']}-{window['window_end']}")
    print(f"    Target: {window['target_start']}-{window['target_end']} ({window['target_strand']})")

    genes = []
    protein_map = {}  # gene_id -> protein_id

    with open(gff_file) as f:
        for line in f:
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            scaffold = parts[0]
            feature_type = parts[2]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            attributes = parts[8]

            # Only process features on our target scaffold
            if scaffold != window['scaffold']:
                continue

            # Only process features within window
            if end < window['window_start'] or start > window['window_end']:
                continue

            # Extract gene features
            if feature_type == 'gene':
                if annotation_type == 'NCBI':
                    # NCBI: Name=LOC117167432
                    match = re.search(r'Name=(LOC\d+)', attributes)
                    if match:
                        gene_id = match.group(1)
                        gene_mid = (start + end) // 2
                        target_mid = (window['target_start'] + window['target_end']) // 2
                        distance = abs(gene_mid - target_mid)

                        genes.append({
                            'gene_id': gene_id,
                            'start': start,
                            'end': end,
                            'strand': strand,
                            'distance': distance
                        })

                elif annotation_type == 'BRAKER3':
                    # BRAKER3: ID=g10440;
                    match = re.search(r'ID=(g\d+);', attributes)
                    if match:
                        gene_id = match.group(1)
                        gene_mid = (start + end) // 2
                        target_mid = (window['target_start'] + window['target_end']) // 2
                        distance = abs(gene_mid - target_mid)

                        genes.append({
                            'gene_id': gene_id,
                            'start': start,
                            'end': end,
                            'strand': strand,
                            'distance': distance
                        })

            # Map proteins to genes for NCBI
            elif feature_type == 'CDS' and annotation_type == 'NCBI':
                protein_match = re.search(r'protein_id=([^;\s]+)', attributes)
                gene_match = re.search(r'gene=(LOC\d+)', attributes)
                if protein_match and gene_match:
                    protein_map[gene_match.group(1)] = protein_match.group(1)

            # Map transcripts to genes for BRAKER3
            elif feature_type == 'mRNA' and annotation_type == 'BRAKER3':
                transcript_match = re.search(r'ID=(g\d+\.t\d+);', attributes)
                parent_match = re.search(r'Parent=(g\d+);', attributes)
                if transcript_match and parent_match:
                    protein_map[parent_match.group(1)] = transcript_match.group(1)

    # Add protein IDs to gene records
    for gene in genes:
        gene['protein_id'] = protein_map.get(gene['gene_id'], gene['gene_id'] + '.t1')

    print(f"    Found {len(genes)} genes in window")

    return genes


def assign_upstream_downstream(genes, target_coords):
    """
    Assign genes as upstream or downstream relative to target.
    Handles strand correctly.

    Args:
        genes: list from parse_gff_genes_in_window()
        target_coords: dict with target start, end, strand

    Returns:
        tuple of (upstream_genes, downstream_genes) sorted by distance
    """

    target_mid = (target_coords['start'] + target_coords['end']) // 2
    target_strand = target_coords['strand']

    upstream = []
    downstream = []

    for gene in genes:
        gene_mid = (gene['start'] + gene['end']) // 2

        # Skip genes that overlap the target
        if not (gene['end'] < target_coords['start'] or gene['start'] > target_coords['end']):
            continue

        if target_strand == '+':
            # For + strand: upstream is lower coordinates, downstream is higher
            if gene_mid < target_mid:
                upstream.append(gene)
            else:
                downstream.append(gene)
        else:
            # For - strand: upstream is higher coordinates, downstream is lower
            if gene_mid > target_mid:
                upstream.append(gene)
            else:
                downstream.append(gene)

    # Sort by distance from target (closest first)
    upstream.sort(key=lambda g: g['distance'])
    downstream.sort(key=lambda g: g['distance'])

    return upstream, downstream


def extract_protein_sequences(genes, proteome_file, annotation_type='BRAKER3'):
    """
    Extract protein sequences for genes.

    Args:
        genes: list of gene dicts with 'protein_id'
        proteome_file: path to proteome FASTA
        annotation_type: 'NCBI' or 'BRAKER3'

    Returns:
        dict mapping gene_id -> SeqRecord
    """

    proteome = SeqIO.to_dict(SeqIO.parse(proteome_file, 'fasta'))
    proteins = {}

    for gene in genes:
        protein_id = gene['protein_id']
        gene_id = gene['gene_id']

        if protein_id in proteome:
            proteins[gene_id] = proteome[protein_id]
        elif f"{protein_id}.1" in proteome:
            # Try with version suffix
            proteins[gene_id] = proteome[f"{protein_id}.1"]
        elif annotation_type == 'BRAKER3' and gene_id + '.t1' in proteome:
            # BRAKER3 fallback
            proteins[gene_id] = proteome[gene_id + '.t1']

    return proteins


def save_flanking_proteins(upstream_genes, downstream_genes, proteins, output_file):
    """
    Save flanking proteins with U/D position labels.

    Format:
        >protein_id|gene_id U1 [description]
        >protein_id|gene_id D1 [description]
    """

    with open(output_file, 'w') as f:
        # Write upstream (closest to target = U1)
        for i, gene in enumerate(upstream_genes, 1):
            gene_id = gene['gene_id']
            if gene_id in proteins:
                seq = proteins[gene_id]
                # Create new record with modified ID and description
                new_record = SeqRecord(
                    seq.seq,
                    id=f"{seq.id}|{gene_id}",
                    description=f"U{i} {seq.description}"
                )
                SeqIO.write([new_record], f, 'fasta')

        # Write downstream (closest to target = D1)
        for i, gene in enumerate(downstream_genes, 1):
            gene_id = gene['gene_id']
            if gene_id in proteins:
                seq = proteins[gene_id]
                new_record = SeqRecord(
                    seq.seq,
                    id=f"{seq.id}|{gene_id}",
                    description=f"D{i} {seq.description}"
                )
                SeqIO.write([new_record], f, 'fasta')


def name_locus(genome_code, scaffold, locus_number):
    """Create systematic locus name."""

    # Shorten scaffold name
    if scaffold.startswith('CM'):
        # RagTag scaffolds like CM021340.1_RagTag
        scaffold_short = scaffold.split('_')[0].replace('CM', 'CM')[:10]
    elif scaffold.startswith('scaffold'):
        scaffold_short = 'scf' + scaffold.replace('scaffold_', '')[:7]
    elif scaffold.startswith('NC_'):
        scaffold_short = 'chr' + scaffold.split('.')[0][-2:]
    elif scaffold.startswith('NW_'):
        scaffold_short = 'scf' + scaffold.split('.')[0][-4:]
    else:
        scaffold_short = scaffold[:10]

    locus_letter = chr(ord('a') + locus_number)

    return f"{genome_code}_{scaffold_short}_{locus_letter}"


def process_genome(genome_code, genome_info, query_protein, output_dir, args):
    """
    Process one genome: run tblastn, find targets, extract flanking genes.

    Returns:
        list of locus dicts
    """

    print(f"\n{'='*80}")
    print(f"Processing {genome_code}: {genome_info['name']}")
    print(f"{'='*80}")

    # Check files exist
    if not Path(genome_info['genome']).exists():
        print(f"  ERROR: Genome file not found: {genome_info['genome']}")
        return []

    if not Path(genome_info['gff']).exists():
        print(f"  WARNING: GFF file not found: {genome_info['gff']}")
        print(f"  Will save target coordinates but no flanking genes")

    # Step 1: Run tblastn
    hits = run_tblastn(query_protein, genome_info['genome'], args.evalue)

    if hits.empty:
        print(f"  No tblastn hits found in {genome_code}")
        return []

    # Step 2: Filter hits
    filtered_hits = filter_tblastn_hits(hits, args.min_coverage, args.min_identity)

    if filtered_hits.empty:
        print(f"  No hits passed filtering in {genome_code}")
        return []

    # Step 3: Process each hit as a potential locus
    loci = []
    locus_counter = 0

    for idx, hit in filtered_hits.iterrows():
        # Parse hit coordinates
        target_coords = parse_hit_coordinates(hit)

        # Define flanking window
        window = define_flanking_window(target_coords, args.window_size)

        # Create locus name
        locus_name = name_locus(genome_code, target_coords['scaffold'], locus_counter)
        locus_counter += 1

        # Create locus directory
        locus_dir = output_dir / locus_name
        locus_dir.mkdir(exist_ok=True)

        print(f"\n  Locus: {locus_name}")
        print(f"    Target: {target_coords['scaffold']}:{target_coords['start']}-{target_coords['end']} ({target_coords['strand']})")
        print(f"    E-value: {target_coords['evalue']:.2e}, Identity: {target_coords['pident']:.1f}%")

        # Save target seed info
        target_seed_df = pd.DataFrame([{
            'locus_id': locus_name,
            'scaffold': target_coords['scaffold'],
            'start': target_coords['start'],
            'end': target_coords['end'],
            'strand': target_coords['strand'],
            'evalue': target_coords['evalue'],
            'bitscore': target_coords['bitscore'],
            'pident': target_coords['pident'],
            'coverage': target_coords['coverage']
        }])
        target_seed_df.to_csv(locus_dir / 'target_seed.tsv', sep='\t', index=False)

        # Save window definition
        window_df = pd.DataFrame([window])
        window_df.to_csv(locus_dir / 'window_definition.tsv', sep='\t', index=False)

        # Step 4: Parse GFF for genes in window (if GFF exists)
        upstream_genes = []
        downstream_genes = []
        flanking_count = 0

        if Path(genome_info['gff']).exists():
            genes_in_window = parse_gff_genes_in_window(
                genome_info['gff'],
                window,
                genome_info['type']
            )

            if genes_in_window:
                # Assign upstream/downstream
                upstream_genes, downstream_genes = assign_upstream_downstream(
                    genes_in_window, target_coords
                )

                flanking_count = len(upstream_genes) + len(downstream_genes)

                print(f"    Flanking: {len(upstream_genes)}U + {len(downstream_genes)}D = {flanking_count}")

                # Save flanking genes info
                flanking_df = pd.DataFrame(genes_in_window)
                if not flanking_df.empty:
                    flanking_df.to_csv(locus_dir / 'flanking_genes.tsv', sep='\t', index=False)

                # Step 5: Extract protein sequences
                if flanking_count > 0:
                    proteins = extract_protein_sequences(
                        genes_in_window,
                        genome_info['proteome'],
                        genome_info['type']
                    )

                    print(f"    Extracted {len(proteins)} protein sequences")

                    # Save flanking proteins with U/D labels
                    flanking_faa = locus_dir / 'flanking_proteins.faa'
                    save_flanking_proteins(upstream_genes, downstream_genes, proteins, flanking_faa)

        # Step 6: Create locus summary
        summary = f"""Locus: {locus_name}
Genome: {genome_code} ({genome_info['name']})
Gene Family: {args.gene_family}

TARGET LOCATION (from tblastn):
  Scaffold: {target_coords['scaffold']}
  Coordinates: {target_coords['start']}-{target_coords['end']} ({target_coords['strand']} strand)
  E-value: {target_coords['evalue']:.2e}
  % Identity: {target_coords['pident']:.1f}%
  Coverage: {target_coords['coverage']:.1%}

FLANKING WINDOW:
  Window size: ±{args.window_size} kb
  Window coordinates: {window['window_start']}-{window['window_end']}
  Genes found: {flanking_count}
    Upstream: {len(upstream_genes)}
    Downstream: {len(downstream_genes)}

FILES:
  target_seed.tsv - tblastn hit details
  window_definition.tsv - flanking window coordinates
  flanking_genes.tsv - genes found in window
  flanking_proteins.faa - protein sequences with U/D labels
"""

        with open(locus_dir / 'locus_summary.txt', 'w') as f:
            f.write(summary)

        # Add to loci list
        loci.append({
            'locus_id': locus_name,
            'gene_family': args.gene_family,
            'genome': genome_code,
            'scaffold': target_coords['scaffold'],
            'target_start': target_coords['start'],
            'target_end': target_coords['end'],
            'target_strand': target_coords['strand'],
            'target_evalue': target_coords['evalue'],
            'upstream_count': len(upstream_genes),
            'downstream_count': len(downstream_genes),
            'flanking_count': flanking_count,
            'flanking_file': str(locus_dir / 'flanking_proteins.faa')
        })

    return loci


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="tblastn-first locus discovery with precise target coordinates",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
  python 01_discover_paralogs_tblastn_first.py \\
      --query-protein XP_015672091.1 \\
      --output-dir outputs/ferritin_phase1_v2 \\
      --gene-family ferritin_MC102
        """
    )

    parser.add_argument('--query-protein', required=True, type=str,
                       help='Query protein ID (must exist in BK proteome)')
    parser.add_argument('--output-dir', required=True, type=Path,
                       help='Output directory for locus definitions')
    parser.add_argument('--gene-family', type=str, default='unknown',
                       help='Gene family name (default: unknown)')
    parser.add_argument('--window-size', type=int, default=WINDOW_SIZE_KB,
                       help=f'Flanking window size in kb (default: {WINDOW_SIZE_KB})')
    parser.add_argument('--evalue', type=float, default=EVALUE_THRESHOLD,
                       help=f'tblastn e-value threshold (default: {EVALUE_THRESHOLD})')
    parser.add_argument('--min-coverage', type=float, default=MIN_COVERAGE,
                       help=f'Minimum query coverage (default: {MIN_COVERAGE})')
    parser.add_argument('--min-identity', type=float, default=MIN_IDENTITY,
                       help=f'Minimum percent identity (default: {MIN_IDENTITY})')
    parser.add_argument('--genomes', type=str, default='BK,LB,TR,DR',
                       help='Comma-separated list of genomes to search (default: BK,LB,TR,DR)')

    return parser.parse_args()


def main():
    args = parse_args()

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    print("="*80)
    print("PHASE 1: tblastn-FIRST LOCUS DISCOVERY")
    print("="*80)
    print(f"\nQuery protein: {args.query_protein}")
    print(f"Gene family: {args.gene_family}")
    print(f"Output directory: {args.output_dir}")
    print(f"\nParameters:")
    print(f"  Flanking window: ±{args.window_size} kb")
    print(f"  E-value threshold: {args.evalue}")
    print(f"  Min coverage: {args.min_coverage}")
    print(f"  Min identity: {args.min_identity}%")

    # Load query protein from BK proteome
    bk_proteome = SeqIO.to_dict(SeqIO.parse(GENOMES['BK']['proteome'], 'fasta'))

    if args.query_protein not in bk_proteome:
        # Try with version suffix
        if f"{args.query_protein}.1" in bk_proteome:
            query_protein = bk_proteome[f"{args.query_protein}.1"]
        else:
            print(f"\nERROR: Query protein {args.query_protein} not found in BK proteome")
            print(f"Available proteins: {list(bk_proteome.keys())[:10]}...")
            sys.exit(1)
    else:
        query_protein = bk_proteome[args.query_protein]

    print(f"\nQuery protein loaded: {query_protein.id} ({len(query_protein.seq)} aa)")

    # Process selected genomes
    genome_codes = [g.strip() for g in args.genomes.split(',')]
    all_loci = []

    for genome_code in genome_codes:
        if genome_code not in GENOMES:
            print(f"\nWARNING: Unknown genome code: {genome_code}")
            continue

        loci = process_genome(
            genome_code,
            GENOMES[genome_code],
            query_protein,
            args.output_dir,
            args
        )

        all_loci.extend(loci)

    # Save locus definitions (compatible with downstream phases)
    print(f"\n{'='*80}")
    print("SAVING RESULTS")
    print(f"{'='*80}")

    if all_loci:
        locus_df = pd.DataFrame(all_loci)
        locus_definitions_file = args.output_dir / 'locus_definitions.tsv'
        locus_df.to_csv(locus_definitions_file, sep='\t', index=False)
        print(f"  locus_definitions.tsv: {len(all_loci)} loci")

        # Save NEW locus_coordinates.tsv for Phase 2/6 validation
        coords_df = locus_df[[
            'locus_id', 'genome', 'scaffold', 'target_start', 'target_end',
            'target_strand', 'target_evalue', 'upstream_count', 'downstream_count'
        ]].copy()
        coords_df.columns = [
            'locus_id', 'genome', 'scaffold', 'start', 'end',
            'strand', 'target_evalue', 'num_flanking_up', 'num_flanking_down'
        ]
        locus_coords_file = args.output_dir / 'locus_coordinates.tsv'
        coords_df.to_csv(locus_coords_file, sep='\t', index=False)
        print(f"  locus_coordinates.tsv: NEW - target locations for validation")
    else:
        print(f"  WARNING: No loci found - no output files created")
        print(f"  Try relaxing --evalue, --min-coverage, or --min-identity thresholds")

    # Summary
    print(f"\n{'='*80}")
    print("SUMMARY")
    print(f"{'-'*80}")

    if all_loci:
        for genome_code in genome_codes:
            genome_loci = [l for l in all_loci if l['genome'] == genome_code]
            if genome_loci:
                avg_flanking = sum(l['flanking_count'] for l in genome_loci) / len(genome_loci)
                print(f"{genome_code}: {len(genome_loci)} loci (avg {avg_flanking:.1f} flanking genes)")

        print(f"\nTotal loci discovered: {len(all_loci)}")
        print(f"Output saved to: {args.output_dir}/")
    else:
        print(f"No loci discovered across {len(genome_codes)} genomes")
        print(f"Consider adjusting thresholds or checking query protein sequence")


if __name__ == "__main__":
    main()
