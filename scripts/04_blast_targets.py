#!/usr/bin/env python3
"""
Step 04: BLAST target genes against all genomes (OPTIMIZED VERSION).

Creates a combined multi-query FASTA with all unique target proteins,
runs BLAST once per genome, then assigns hits to loci with best-match tracking.

Usage:
    python 04_blast_targets.py \\
        --locus-defs <path/to/locus_definitions.tsv> \\
        --query-proteins <path/to/01_extracted_proteins> \\
        --genome-db-dir <path/to/ragtag_dbs> \\
        --output-dir <path/to/target_genes> \\
        [--evalue 1e-5] \\
        [--max-targets 10000] \\
        [--threads 16]
"""

from pathlib import Path
import re
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


def _normalize_scaffold(name: str) -> str:
    s = str(name).split()[0]
    # Remove RagTag suffix
    if s.endswith('_RagTag'):
        s = s[: -len('_RagTag')]
    # Strip version for common accessions
    for pat in (r'^(CM\d+)(?:\.\d+)?$', r'^(NC_\d+)(?:\.\d+)?$', r'^(NW_\d+)(?:\.\d+)?$', r'^(NT_\d+)(?:\.\d+)?$'):
        m = re.match(pat, s)
        if m:
            return m.group(1)
    return s

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
                        'scaffold': _normalize_scaffold(hit_accession),
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

def merge_intervals(intervals):
    """Merge overlapping intervals to create coverage track.

    Args:
        intervals: List of (start, end) tuples

    Returns:
        List of non-overlapping (start, end) tuples sorted by start position
    """
    if not intervals:
        return []

    # Sort by start position
    sorted_intervals = sorted(intervals, key=lambda x: x[0])
    merged = [sorted_intervals[0]]

    for current in sorted_intervals[1:]:
        last_merged = merged[-1]

        # If current overlaps with last merged, extend the last merged
        if current[0] <= last_merged[1]:
            merged[-1] = (last_merged[0], max(last_merged[1], current[1]))
        else:
            # No overlap, start new interval
            merged.append(current)

    return merged


def find_coverage_gaps(coverage_track, min_gap_kb=10):
    """Find zero-coverage gaps >= min_gap_kb between coverage segments.

    Args:
        coverage_track: List of (start, end) coverage intervals
        min_gap_kb: Minimum gap size in kb to be considered a breakpoint

    Returns:
        List of (gap_start, gap_end, gap_size_kb) tuples
    """
    gaps = []
    for i in range(len(coverage_track) - 1):
        gap_start = coverage_track[i][1]  # end of current segment
        gap_end = coverage_track[i+1][0]  # start of next segment
        gap_size_kb = (gap_end - gap_start) / 1000

        if gap_size_kb >= min_gap_kb:
            gaps.append((gap_start, gap_end, gap_size_kb))

    return gaps


def chain_hsps_simple(hits, segment_start, segment_end):
    """Chain HSPs within a coverage segment using simple colinear scoring.

    Args:
        hits: List of HSP hits within this segment
        segment_start: Start of coverage segment
        segment_end: End of coverage segment

    Returns:
        List of chains (each chain is a list of hits)
    """
    if not hits:
        return []

    # Sort by subject start position
    hits = sorted(hits, key=lambda h: h['start'])

    # For now, use simple greedy chaining:
    # All hits in the segment belong to one chain (they share coverage)
    # Future enhancement: DP-based chaining considering bitscore and gaps

    # Filter hits that actually fall within segment bounds
    segment_hits = [h for h in hits if h['start'] >= segment_start and h['end'] <= segment_end]

    if segment_hits:
        return [segment_hits]  # Single chain containing all hits
    else:
        return []


def cluster_into_loci(hits, locus_id, gene_family, min_split_gap_kb=10, min_hits=1):
    """Cluster target hits into loci using coverage-track-based splitting.

    NEW ALGORITHM:
    1. Group by scaffold and strand ONLY (ignore frame - it flips at exon boundaries)
    2. Build HSP coverage track (union of all HSP intervals, any frame)
    3. Split at zero-coverage gaps >= min_split_gap_kb (hard breakpoints)
    4. Within each segment, chain HSPs by colinear order
    5. Each chain becomes a locus

    Args:
        hits: List of BLAST HSP hits
        locus_id: Base locus ID for naming
        gene_family: Gene family name
        min_split_gap_kb: Minimum zero-coverage gap (kb) to split clusters (default 10kb)
        min_hits: Minimum number of hits to call a locus (default 1)

    Returns:
        List of locus dictionaries
    """
    if not hits:
        return []

    # Group by scaffold and strand ONLY (no frame)
    scaffold_hits = defaultdict(list)
    for hit in hits:
        key = (hit['scaffold'], hit['strand'])
        scaffold_hits[key].append(hit)

    loci = []
    locus_num = 1

    for (scaffold, strand), s_hits in scaffold_hits.items():
        # Build coverage track: union of all HSP intervals (any frame)
        intervals = [(h['start'], h['end']) for h in s_hits]
        coverage_track = merge_intervals(intervals)

        # Find coverage gaps >= min_split_gap_kb
        gaps = find_coverage_gaps(coverage_track, min_split_gap_kb)

        # Create split points from gaps
        split_points = []
        for gap_start, gap_end, gap_size_kb in gaps:
            split_points.append(gap_start)  # End of segment before gap
            split_points.append(gap_end)    # Start of segment after gap

        # If no splits, all hits belong to one locus
        if not split_points:
            # No large gaps - all hits form one locus
            if len(s_hits) >= min_hits:
                locus_name = f"{locus_id}_{scaffold}_{locus_num:03d}"
                start = min(h['start'] for h in s_hits)
                end = max(h['end'] for h in s_hits)
                best_evalue = min(h['evalue'] for h in s_hits)
                best_bitscore = max(h.get('bitscore', 0) for h in s_hits)
                frames = {h.get('frame', 0) for h in s_hits}

                loci.append({
                    'locus_name': locus_name,
                    'scaffold': scaffold,
                    'strand': strand,
                    'frame': list(frames)[0] if len(frames) == 1 else 0,
                    'frames_present': ','.join(map(str, sorted(frames))),
                    'start': start,
                    'end': end,
                    'span_kb': (end - start) / 1000,
                    'num_hits': len(s_hits),
                    'gene_family': gene_family,
                    'best_evalue': best_evalue,
                    'best_bitscore': best_bitscore
                })
                locus_num += 1
        else:
            # Split hits based on gap boundaries
            split_points.sort()

            # Add boundaries at start and end
            all_positions = sorted([h['start'] for h in s_hits] + [h['end'] for h in s_hits])
            region_start = all_positions[0]
            region_end = all_positions[-1]

            # Create segments between split points
            segments = []
            current_start = region_start

            for i in range(0, len(split_points), 2):
                segment_end = split_points[i]  # End before gap
                segments.append((current_start, segment_end))

                if i + 1 < len(split_points):
                    current_start = split_points[i + 1]  # Start after gap

            # Add final segment
            if current_start < region_end:
                segments.append((current_start, region_end))

            # Assign hits to segments and create loci
            for seg_start, seg_end in segments:
                segment_hits = [
                    h for h in s_hits
                    if not (h['end'] < seg_start or h['start'] > seg_end)
                ]

                if len(segment_hits) >= min_hits:
                    locus_name = f"{locus_id}_{scaffold}_{locus_num:03d}"
                    start = min(h['start'] for h in segment_hits)
                    end = max(h['end'] for h in segment_hits)
                    best_evalue = min(h['evalue'] for h in segment_hits)
                    best_bitscore = max(h.get('bitscore', 0) for h in segment_hits)
                    frames = {h.get('frame', 0) for h in segment_hits}

                    loci.append({
                        'locus_name': locus_name,
                        'scaffold': scaffold,
                        'strand': strand,
                        'frame': list(frames)[0] if len(frames) == 1 else 0,
                        'frames_present': ','.join(map(str, sorted(frames))),
                        'start': start,
                        'end': end,
                        'span_kb': (end - start) / 1000,
                        'num_hits': len(segment_hits),
                        'gene_family': gene_family,
                        'best_evalue': best_evalue,
                        'best_bitscore': best_bitscore
                    })
                    locus_num += 1

    return loci

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="BLAST target genes against all genomes",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

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
    parser.add_argument('--gene-family', type=str, default='unknown',
                        help='Gene family name for these targets (default: unknown)')

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
    def _normalize_genome_id(name: str) -> str:
        # Match 02_synteny_detection normalization: GCA_XXXXXX.Y -> keep prefix+version
        if name.startswith('GCA_') or name.startswith('GCF_'):
            parts = name.split('_')
            if len(parts) >= 2:
                return parts[0] + '_' + parts[1]
        return name

    for db_file in args.genome_db_dir.glob("*.nhr"):
        genome_name = db_file.stem
        genome_dbs[genome_name] = {
            'db': args.genome_db_dir / genome_name,
            'genome_id': _normalize_genome_id(genome_name)
        }
    print(f"  Found {len(genome_dbs)} genome databases", flush=True)

    # Run BLAST once per genome
    print(f"\n[5] Running tBLASTn ({len(genome_dbs)} searches)...", flush=True)
    blast_dir = args.output_dir / "blast_xml"
    blast_dir.mkdir(exist_ok=True)

    genome_hits = {}  # genome -> list of hits

    for genome_name, meta in genome_dbs.items():
        genome_db = meta['db']
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

    # Set minimum coverage gap for splitting clusters
    # Gaps >= this size (kb) with zero coverage will split loci
    min_split_gap_kb = 10  # Default: 10kb zero-coverage gap splits clusters

    # Process all hits - NO per-locus duplication
    print("\n[6] Clustering target hits (coverage-track-based)...", flush=True)
    print(f"    Algorithm: Split at zero-coverage gaps >= {min_split_gap_kb}kb", flush=True)
    all_target_loci = []

    for genome_name, hits in genome_hits.items():
        if not hits:
            continue

        print(f"  Processing {genome_name}: {len(hits)} hits", flush=True)

        # Filter multi-frame artifacts: keep only best frame per genomic location
        # NOTE: Frame is no longer used for clustering, but we still filter obvious artifacts
        hits_by_location = defaultdict(list)
        for hit in hits:
            # Group by scaffold + approximate location (5kb bins)
            location_key = (hit['scaffold'], hit['start'] // 5000)
            hits_by_location[location_key].append(hit)

        filtered_hits = []
        for location_key, location_hits in hits_by_location.items():
            if len(location_hits) == 1:
                filtered_hits.append(location_hits[0])
            else:
                # Multiple hits at same location - check if multi-frame artifacts
                frames = {h['frame'] for h in location_hits}
                if len(frames) > 1:
                    # Multi-frame artifact - keep only strongest
                    best_hit = min(location_hits, key=lambda h: h['evalue'])
                    filtered_hits.append(best_hit)
                    print(f"    Multi-frame artifact at {location_key}: kept frame {best_hit['frame']} (e-value {best_hit['evalue']:.2e}), filtered {len(location_hits)-1} others", flush=True)
                else:
                    # Same frame, might be legitimate multi-exon hits
                    filtered_hits.extend(location_hits)

        print(f"    After multi-frame filtering: {len(filtered_hits)} hits", flush=True)

        # Cluster hits by query protein
        for query_id in {h['query_id'] for h in filtered_hits}:
            query_hits = [h for h in filtered_hits if h['query_id'] == query_id]

            # Create a locus_id placeholder (will be assigned by Phase 5)
            locus_id_placeholder = f"{query_id.split('.')[0]}_targets"

            # Cluster into loci using coverage-track-based algorithm
            target_loci = cluster_into_loci(query_hits, locus_id_placeholder, args.gene_family, min_split_gap_kb, min_hits=1)

            # Normalize genome id to match Phase 2/3 outputs
            genome_id = genome_dbs[genome_name]['genome_id'] if genome_name in genome_dbs else _normalize_genome_id(genome_name)

            for target_locus in target_loci:
                target_locus['genome'] = genome_id
                target_locus['query_id'] = query_id
                # NO parent_locus assignment - Phase 5 will assign based on synteny proximity

                all_target_loci.append(target_locus)

    # Save results
    print("\n[7] Saving target loci...", flush=True)
    loci_file = args.output_dir / "all_target_loci.tsv"

    if all_target_loci:
        loci_df_out = pd.DataFrame(all_target_loci)
        loci_df_out.to_csv(loci_file, sep='\t', index=False)
        print(f"  Saved {len(loci_df_out)} target loci to {loci_file.name}", flush=True)

        # Summary
        print("\n[8] Summary by query protein:", flush=True)
        for query_id in sorted(loci_df_out['query_id'].unique()):
            query_loci = loci_df_out[loci_df_out['query_id'] == query_id]
            genomes_with_targets = query_loci['genome'].nunique()
            total_loci = len(query_loci)
            print(f"  {query_id}: {total_loci} loci in {genomes_with_targets} genomes", flush=True)
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

    else:
        print("\nNo target loci found", flush=True)

    print(f"\nOutputs saved to: {args.output_dir}", flush=True)
    print("\nNext step: 05_classify_targets.py", flush=True)

if __name__ == "__main__":
    main()
