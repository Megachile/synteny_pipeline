#!/usr/bin/env python3
"""
Step 04: BLAST target genes against all genomes (OPTIMIZED VERSION).

Creates a combined multi-query FASTA with all unique target proteins,
runs BLAST once per genome, then assigns hits to loci with best-match tracking.

Usage:
    python 04_blast_targets.py \\
        --blocks <path/to/synteny_blocks_filtered.tsv> \\
        --locus-defs <path/to/locus_definitions.tsv> \\
        --query-proteins <path/to/01_extracted_proteins> \\
        --genome-db-dir <path/to/ragtag_dbs> \\
        --output-dir <path/to/target_genes> \\
        [--evalue 1e-5] \\
        [--max-targets 10000] \\
        [--threads 16]
"""

from pathlib import Path
import subprocess
import pandas as pd
from collections import defaultdict
import xml.etree.ElementTree as ET
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
import sys
sys.path.insert(0, str(Path(__file__).parent))
import extract_with_exonerate

def run_tblastn(query_file, genome_db, output_xml, evalue, max_targets, threads):
    """Run tBLASTn search for target proteins."""
    cmd = [
        'tblastn',
        '-query', str(query_file),
        '-db', str(genome_db),
        '-outfmt', '5',  # XML format
        '-evalue', str(evalue),
        '-max_target_seqs', str(max_targets),
        '-num_threads', str(threads)
    ]

    with open(output_xml, 'w') as f:
        result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, text=True)

    return result.returncode == 0

def parse_blast_xml(xml_file):
    """Parse BLAST XML and extract hits with query info."""
    hits = []

    try:
        tree = ET.parse(xml_file)
        root = tree.getroot()

        for iteration in root.findall('.//Iteration'):
            query_id = iteration.find('Iteration_query-def').text.split()[0]
            query_len = int(iteration.find('Iteration_query-len').text)

            for hit in iteration.findall('.//Hit'):
                hit_def = hit.find('Hit_def').text
                hit_accession = hit.find('Hit_accession').text

                for hsp in hit.findall('.//Hsp'):
                    evalue = float(hsp.find('Hsp_evalue').text)
                    bitscore = float(hsp.find('Hsp_bit-score').text)
                    identity = int(hsp.find('Hsp_identity').text)
                    align_length = int(hsp.find('Hsp_align-len').text)

                    sstart = int(hsp.find('Hsp_hit-from').text)
                    send = int(hsp.find('Hsp_hit-to').text)

                    # Extract reading frame from tBLASTn results
                    frame = int(hsp.find('Hsp_hit-frame').text)

                    if sstart < send:
                        strand = '+'
                        start = sstart
                        end = send
                    else:
                        strand = '-'
                        start = send
                        end = sstart

                    hits.append({
                        'query_id': query_id,
                        'scaffold': hit_accession,
                        'scaffold_desc': hit_def,
                        'strand': strand,
                        'frame': frame,
                        'start': start,
                        'end': end,
                        'evalue': evalue,
                        'bitscore': bitscore,
                        'identity': identity,
                        'align_length': align_length,
                        'pident': identity / align_length * 100,
                        'query_length': query_len
                    })

    except Exception as e:
        print(f"    Parse error: {e}")

    return hits

def cluster_into_loci(hits, locus_id, gene_family, gap_kb=50, min_hits=1):
    """Cluster target hits into loci based on proximity.

    Args:
        gap_kb: Maximum gap in kb between hits to cluster together
        min_hits: Minimum number of hits to call a locus (default 1 for target genes)
    """
    if not hits:
        return []

    # Group by scaffold, strand, and reading frame
    # Reading frame is critical - HSPs in different frames cannot be the same gene
    scaffold_hits = defaultdict(list)
    for hit in hits:
        frame = hit.get('frame', 0)  # Default to 0 if frame not present
        key = (hit['scaffold'], hit['strand'], frame)
        scaffold_hits[key].append(hit)

    loci = []
    locus_num = 1

    for (scaffold, strand, frame), s_hits in scaffold_hits.items():
        # Sort by position
        s_hits.sort(key=lambda x: x['start'])

        # Cluster by gap
        current_cluster = [s_hits[0]]

        for hit in s_hits[1:]:
            prev_end = current_cluster[-1]['end']
            gap = (hit['start'] - prev_end) / 1000  # kb

            if gap <= gap_kb:
                current_cluster.append(hit)
            else:
                # Save current cluster as locus
                if len(current_cluster) >= min_hits:
                    locus_name = f"{locus_id}_{scaffold}_{locus_num:03d}"
                    start = min(h['start'] for h in current_cluster)
                    end = max(h['end'] for h in current_cluster)
                    best_evalue = min(h['evalue'] for h in current_cluster)

                    loci.append({
                        'locus_name': locus_name,
                        'scaffold': scaffold,
                        'strand': strand,
                        'frame': frame,
                        'start': start,
                        'end': end,
                        'span_kb': (end - start) / 1000,
                        'num_hits': len(current_cluster),
                        'gene_family': gene_family,
                        'best_evalue': best_evalue
                    })
                    locus_num += 1

                # Start new cluster
                current_cluster = [hit]

        # Save final cluster
        if len(current_cluster) >= min_hits:
            locus_name = f"{locus_id}_{scaffold}_{locus_num:03d}"
            start = min(h['start'] for h in current_cluster)
            end = max(h['end'] for h in current_cluster)
            best_evalue = min(h['evalue'] for h in current_cluster)

            loci.append({
                'locus_name': locus_name,
                'scaffold': scaffold,
                'strand': strand,
                'frame': frame,
                'start': start,
                'end': end,
                'span_kb': (end - start) / 1000,
                'num_hits': len(current_cluster),
                'gene_family': gene_family,
                'best_evalue': best_evalue
            })
            locus_num += 1

    return loci

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="BLAST target genes against all genomes",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('--blocks', required=True, type=Path,
                        help='Path to synteny_blocks_filtered.tsv')
    parser.add_argument('--locus-defs', required=True, type=Path,
                        help='Path to locus_definitions.tsv')
    parser.add_argument('--query-proteins', required=True, type=Path,
                        help='Directory containing query protein FASTA files')
    parser.add_argument('--genome-db-dir', required=True, type=Path,
                        help='Directory with BLAST databases for genomes')
    parser.add_argument('--output-dir', required=True, type=Path,
                        help='Output directory for target genes')
    parser.add_argument('--evalue', type=str, default="1e-5",
                        help='E-value threshold (default: 1e-5)')
    parser.add_argument('--max-targets', type=int, default=10000,
                        help='Maximum BLAST targets per query (default: 10000)')
    parser.add_argument('--threads', type=int, default=16,
                        help='Number of threads for BLAST (default: 16)')

    return parser.parse_args()

def main():
    """Main execution function."""
    args = parse_args()

    print("=" * 80, flush=True)
    print("STEP 04: BLAST TARGET GENES (OPTIMIZED)", flush=True)
    print("=" * 80, flush=True)
    print(f"\nInput files:")
    print(f"  Locus definitions: {args.locus_defs}")
    print(f"  Query proteins dir: {args.query_proteins}")
    print(f"  Genome DB dir: {args.genome_db_dir}")
    print(f"  Output directory: {args.output_dir}")
    print(f"\nParameters:")
    print(f"  E-value threshold: {args.evalue}")
    print(f"  Max BLAST targets: {args.max_targets}")
    print(f"  Threads: {args.threads}")

    # Create output directory
    args.output_dir.mkdir(exist_ok=True, parents=True)

    # Load locus definitions
    print("\n[1] Loading locus definitions...", flush=True)
    loci_df = pd.read_csv(args.locus_defs, sep='\t')
    print(f"  Loaded {len(loci_df)} loci", flush=True)

    # Find unique target proteins
    print("\n[2] Finding unique target proteins...", flush=True)

    target_proteins = {}  # seq -> (protein_id, list of loci using it)
    locus_to_protein = {}  # locus_id -> protein_id

    for _, locus_row in loci_df.iterrows():
        locus_id = locus_row['locus_id']
        targets_file = args.query_proteins / locus_id / f"{locus_id}_targets.faa"

        if not targets_file.exists():
            print(f"  WARNING: Target file not found for {locus_id}", flush=True)
            continue

        # Read target protein
        with open(targets_file) as f:
            records = list(SeqIO.parse(f, 'fasta'))
            if records:
                target_seq = str(records[0].seq)
                target_id = records[0].id

                if target_seq not in target_proteins:
                    target_proteins[target_seq] = {
                        'protein_id': target_id,
                        'loci': [],
                        'record': records[0]
                    }
                target_proteins[target_seq]['loci'].append(locus_id)
                locus_to_protein[locus_id] = target_id

    print(f"  Found {len(target_proteins)} unique target proteins for {len(loci_df)} loci", flush=True)
    for i, (seq, info) in enumerate(target_proteins.items(), 1):
        print(f"    Protein {i} ({info['protein_id']}): used by {len(info['loci'])} loci {info['loci']}", flush=True)

    # Create combined query file
    print("\n[3] Creating combined query file...", flush=True)
    combined_query = args.output_dir / "combined_targets.faa"

    with open(combined_query, 'w') as f:
        for info in target_proteins.values():
            SeqIO.write(info['record'], f, 'fasta')

    print(f"  Saved {len(target_proteins)} proteins to {combined_query.name}", flush=True)

    # Find genome databases
    print("\n[4] Finding genome databases...", flush=True)
    genome_dbs = {}
    for db_file in args.genome_db_dir.glob("*.nhr"):
        genome_name = db_file.stem
        genome_dbs[genome_name] = args.genome_db_dir / genome_name
    print(f"  Found {len(genome_dbs)} genome databases", flush=True)

    # Run BLAST once per genome
    print(f"\n[5] Running tBLASTn ({len(genome_dbs)} searches)...", flush=True)
    blast_dir = args.output_dir / "blast_xml"
    blast_dir.mkdir(exist_ok=True)

    genome_hits = {}  # genome -> list of hits

    for genome_name, genome_db in genome_dbs.items():
        output_xml = blast_dir / f"{genome_name}.xml"

        if output_xml.exists():
            print(f"  {genome_name}: Using existing results", flush=True)
        else:
            print(f"  {genome_name}: Running tBLASTn...", end='', flush=True)
            if run_tblastn(combined_query, genome_db, output_xml, args.evalue, args.max_targets, args.threads):
                print(" done", flush=True)
            else:
                print(" FAILED", flush=True)
                continue

        # Parse hits
        hits = parse_blast_xml(output_xml)
        genome_hits[genome_name] = hits

    print(f"  Completed {len(genome_hits)} genomes", flush=True)

    # Infer gap_kb from config or use default
    gap_kb = 50  # Default clustering gap

    # Assign hits to loci and track best matches
    print("\n[6] Assigning hits to loci...", flush=True)
    all_target_loci = []

    for _, locus_row in loci_df.iterrows():
        locus_id = locus_row['locus_id']
        gene_family = locus_row['gene_family']
        expected_protein = locus_to_protein.get(locus_id)

        print(f"\n  Processing {locus_id} (expects {expected_protein})...", flush=True)

        locus_output = args.output_dir / locus_id
        locus_output.mkdir(exist_ok=True, parents=True)

        # Collect all hits for this locus across genomes
        for genome_name, hits in genome_hits.items():
            # Filter to hits for this locus's expected protein
            # BUT also check for query mismatches (hits to other proteins)
            locus_hits = []
            query_mismatches = []

            for hit in hits:
                if hit['query_id'] == expected_protein:
                    locus_hits.append(hit)
                else:
                    # Track hits to other query proteins (potential ortho-paralogs)
                    query_mismatches.append(hit)

            if not locus_hits:
                continue

            # Cluster hits into target loci (min_hits=1 for single-protein targets)
            target_loci = cluster_into_loci(locus_hits, locus_id, gene_family, gap_kb, min_hits=1)

            for target_locus in target_loci:
                target_locus['genome'] = genome_name
                target_locus['parent_locus'] = locus_id
                target_locus['expected_query'] = expected_protein
                target_locus['query_matched'] = expected_protein  # Matches as expected

                # Check if there are strong hits to other queries in this region
                region_mismatches = [h for h in query_mismatches if
                                   h['scaffold'] == target_locus['scaffold'] and
                                   h['start'] >= target_locus['start'] - 50000 and
                                   h['end'] <= target_locus['end'] + 50000]

                if region_mismatches:
                    best_mismatch = max(region_mismatches, key=lambda x: x['bitscore'])
                    target_locus['also_matches'] = best_mismatch['query_id']
                    target_locus['mismatch_bitscore'] = best_mismatch['bitscore']
                else:
                    target_locus['also_matches'] = None
                    target_locus['mismatch_bitscore'] = None

                all_target_loci.append(target_locus)

    # Save results
    print("\n[7] Saving target loci...", flush=True)
    loci_file = args.output_dir / "all_target_loci.tsv"

    if all_target_loci:
        loci_df_out = pd.DataFrame(all_target_loci)
        loci_df_out.to_csv(loci_file, sep='\t', index=False)
        print(f"  Saved {len(loci_df_out)} target loci to {loci_file.name}", flush=True)

        # Summary
        print("\n[8] Summary by locus:", flush=True)
        for parent_locus in sorted(loci_df_out['parent_locus'].unique()):
            parent_loci = loci_df_out[loci_df_out['parent_locus'] == parent_locus]
            genomes_with_targets = parent_loci['genome'].nunique()
            total_loci = len(parent_loci)
            with_mismatches = parent_loci['also_matches'].notna().sum()
            print(f"  {parent_locus}: {total_loci} loci in {genomes_with_targets} genomes ({with_mismatches} with query mismatches)", flush=True)
    else:
        print("  No target loci found", flush=True)

    # Overall summary
    print("\n" + "=" * 80, flush=True)
    print("TARGET BLAST COMPLETE", flush=True)
    print("=" * 80, flush=True)

    if all_target_loci:
        print(f"\nTotal target loci found: {len(all_target_loci)}", flush=True)
        print(f"Genomes with targets: {loci_df_out['genome'].nunique()}", flush=True)
        print(f"Average loci per genome: {len(all_target_loci) / loci_df_out['genome'].nunique():.1f}", flush=True)

        with_mismatches = loci_df_out['also_matches'].notna().sum()
        print(f"Loci with query mismatches: {with_mismatches} ({with_mismatches/len(all_target_loci)*100:.1f}%)", flush=True)
    else:
        print("\nNo target loci found", flush=True)

    print(f"\nOutputs saved to: {args.output_dir}", flush=True)
    print("\nNext step: 05_classify_targets.py", flush=True)

if __name__ == "__main__":
    main()
