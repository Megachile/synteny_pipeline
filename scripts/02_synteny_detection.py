#!/usr/bin/env python3
"""
Step 02: Detect synteny blocks using flanking proteins.

Uses tBLASTn to search flanking proteins against target genomes.
Clusters BLAST hits into synteny blocks based on proximity.

Usage:
    python 02_synteny_detection.py \\
        --locus-defs <path/to/locus_definitions.tsv> \\
        --genome-db-dir <path/to/blast_dbs> \\
        --output-dir <path/to/02_synteny_blocks> \\
        --evalue 1e-5 \\
        --max-targets 50 \\
        --threads 16 \\
        [--locus <locus_id>]  # Optional: for SLURM array jobs
"""

from pathlib import Path
import subprocess
import pandas as pd
from collections import defaultdict
import xml.etree.ElementTree as ET
import csv
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def run_tblastn(query_file, genome_db, output_xml, evalue="1e-5", max_targets=50, threads=16):
    """Run tBLASTn search."""
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
    """Parse BLAST XML and extract hits."""
    hits = []

    try:
        tree = ET.parse(xml_file)
        root = tree.getroot()

        for iteration in root.findall('.//Iteration'):
            query_id = iteration.find('Iteration_query-def').text.split()[0]  # e.g., "locus_U1"

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

                    # Extract HSP sequences (already translated by tBLASTn)
                    hsp_hseq = hsp.find('Hsp_hseq').text  # Translated hit sequence
                    hsp_qseq = hsp.find('Hsp_qseq').text  # Query sequence

                    # Determine strand
                    if sstart < send:
                        strand = '+'
                        start = sstart
                        end = send
                    else:
                        strand = '-'
                        start = send
                        end = sstart

                    hits.append({
                        'qseqid': query_id,
                        'sseqid': hit_accession,
                        'scaffold_desc': hit_def,
                        'strand': strand,
                        'coord_start': start,
                        'coord_end': end,
                        'evalue': evalue,
                        'bitscore': bitscore,
                        'pident': identity / align_length * 100,
                        'length': align_length,
                        'hit_protein_seq': hsp_hseq,
                        'query_protein_seq': hsp_qseq
                    })

    except Exception as e:
        print(f"    Parse error: {e}")

    return hits

def cluster_into_blocks(hits, max_gap_kb=500):
    """Cluster BLAST hits into synteny blocks."""

    # Group by scaffold and strand
    grouped = defaultdict(list)
    for hit in hits:
        key = (hit['sseqid'], hit['strand'])
        grouped[key].append(hit)

    blocks = []
    block_id = 1

    for (scaffold, strand), scaffold_hits in grouped.items():
        # Sort by position
        scaffold_hits = sorted(scaffold_hits, key=lambda x: x['coord_start'])

        # Cluster into blocks
        current_block_hits = [scaffold_hits[0]]

        for hit in scaffold_hits[1:]:
            # Check if within max_gap of current block
            block_end = max(h['coord_end'] for h in current_block_hits)
            gap = hit['coord_start'] - block_end

            if gap <= max_gap_kb * 1000:  # Convert kb to bp
                current_block_hits.append(hit)
            else:
                # Save current block and start new one
                # Count unique TARGET proteins (sseqid) in the genome, not query proteins
                unique_target_proteins = len(set((h['sseqid'], h['coord_start'], h['coord_end']) for h in current_block_hits))
                if unique_target_proteins >= config.MIN_PROTEINS_FOR_BLOCK:
                    blocks.append({
                        'block_id': f"block_{block_id:05d}",
                        'scaffold': scaffold,
                        'strand': strand,
                        'start': min(h['coord_start'] for h in current_block_hits),
                        'end': max(h['coord_end'] for h in current_block_hits),
                        'num_target_proteins': unique_target_proteins,
                        'num_query_matches': len(set(h['qseqid'] for h in current_block_hits)),
                        'hits': current_block_hits
                    })
                    block_id += 1

                current_block_hits = [hit]

        # Don't forget the last block
        # Count unique TARGET proteins (sseqid) in the genome, not query proteins
        unique_target_proteins = len(set((h['sseqid'], h['coord_start'], h['coord_end']) for h in current_block_hits))
        if unique_target_proteins >= config.MIN_PROTEINS_FOR_BLOCK:
            blocks.append({
                'block_id': f"block_{block_id:05d}",
                'scaffold': scaffold,
                'strand': strand,
                'start': min(h['coord_start'] for h in current_block_hits),
                'end': max(h['coord_end'] for h in current_block_hits),
                'num_target_proteins': unique_target_proteins,
                'num_query_matches': len(set(h['qseqid'] for h in current_block_hits)),
                'hits': current_block_hits
            })
            block_id += 1

    return blocks



def save_hit_sequences(hits, genome_id, locus_output):
    """
    Save protein sequences from BLAST HSPs.
    
    Args:
        hits: List of BLAST hit dictionaries (with hit_protein_seq field)
        genome_id: Standardized genome ID
        locus_output: Output directory for this locus
    
    Returns:
        Number of sequences saved
    """
    seqs_dir = locus_output / "hit_sequences"
    seqs_dir.mkdir(exist_ok=True)
    
    sequences = []
    
    for i, hit in enumerate(hits, 1):
        # Remove gaps from the aligned sequence
        protein_seq = hit['hit_protein_seq'].replace('-', '')
        
        # Create SeqRecord
        seq_id = f"{genome_id}_{hit['qseqid']}_{hit['sseqid']}_{hit['coord_start']}-{hit['coord_end']}_hit{i}"
        seq_record = SeqRecord(
            seq=Seq(protein_seq),
            id=seq_id,
            description=f"{genome_id} | {hit['qseqid']} hit | {hit['sseqid']}:{hit['coord_start']}-{hit['coord_end']}({hit['strand']}) | pident={hit['pident']:.1f}% | len={len(protein_seq)}aa"
        )
        
        sequences.append(seq_record)
    
    # Save sequences for this genome
    if sequences:
        output_fasta = seqs_dir / f"{genome_id}_hit_proteins.fasta"
        SeqIO.write(sequences, output_fasta, "fasta")
    
    return len(sequences)

def filter_flanking_proteins(flanking_file, locus_id):
    """
    Filter flanking proteins to top N for synteny detection.

    Reads only the first MAX_FLANKING_FOR_SYNTENY proteins from the file
    and saves to a filtered temporary file in the output directory.

    Args:
        flanking_file: Path to original flanking proteins file
        locus_id: Locus identifier for naming

    Returns:
        Path to filtered flanking file, or None on error
    """
    try:
        # Read first N proteins
        filtered_seqs = []
        with open(flanking_file, 'r') as f:
            for i, record in enumerate(SeqIO.parse(f, 'fasta')):
                if i >= config.MAX_FLANKING_FOR_SYNTENY:
                    break
                filtered_seqs.append(record)

        if not filtered_seqs:
            print(f"    WARNING: No sequences found in {flanking_file}")
            return None

        # Save filtered file in output directory
        filtered_path = config.STEP02_SYNTENY / locus_id / f"{locus_id}_flanking_filtered.faa"
        filtered_path.parent.mkdir(exist_ok=True, parents=True)

        with open(filtered_path, 'w') as f:
            SeqIO.write(filtered_seqs, f, 'fasta')

        print(f"    Using {len(filtered_seqs)} flanking proteins (filtered from {flanking_file.name})")

        return filtered_path

    except Exception as e:
        print(f"    ERROR filtering flanking proteins: {e}")
        return None

def process_locus(locus_id, flanking_file, genome_dbs, output_dir, evalue="1e-5", max_targets=50, threads=16, max_gap_kb=500):
    """Process synteny detection for one locus."""

    print(f"\n  Processing {locus_id}...")

    # Input files - flanking_file path now passed as argument from locus_definitions.tsv
    if not flanking_file.exists():
        print(f"    ERROR: Flanking proteins not found: {flanking_file}")
        return []

    # Filter flanking proteins to top N for efficiency
    filtered_flanking = filter_flanking_proteins(flanking_file, locus_id)
    if filtered_flanking is None:
        print(f"    ERROR: Failed to filter flanking proteins")
        return []

    # Output directory
    locus_output = output_dir / locus_id
    locus_output.mkdir(exist_ok=True, parents=True)
    blast_dir = locus_output / "blast_xml"
    blast_dir.mkdir(exist_ok=True)

    all_blocks = []

    # Run BLAST against each genome
    for genome_name, genome_db in genome_dbs.items():
        # Extract GCA accession if it's in the name, otherwise use full name
        # Handle formats like: GCA_020615435.1_ASM2061543v1_genomic or species_name
        if genome_name.startswith('GCA_') or genome_name.startswith('GCF_'):
            genome_id = genome_name.split('_')[0] + '_' + genome_name.split('_')[1]
        else:
            genome_id = genome_name

        output_xml = blast_dir / f"{genome_name}.xml"

        # Run tBLASTn
        if output_xml.exists():
            print(f"    {genome_name}: Using existing BLAST results")
        else:
            print(f"    {genome_name}: Running tBLASTn...", end='')
            if run_tblastn(filtered_flanking, genome_db, output_xml, evalue, max_targets, threads):
                print(" done")
            else:
                print(" FAILED")
                continue

        # Parse hits (includes HSP sequences)
        hits = parse_blast_xml(output_xml)

        if not hits:
            print("")
            continue

        # Save HSP sequences for flanking proteins (for SwissProt validation)
        num_seqs = save_hit_sequences(hits, genome_id, locus_output)

        # Cluster into blocks
        blocks = cluster_into_blocks(hits, max_gap_kb)
        print(f" - {len(blocks)} synteny blocks, {num_seqs} HSP sequences saved")

        # Add genome and locus info to blocks
        for block in blocks:
            block['genome'] = genome_id  # Use standardized genome ID
            block['locus_id'] = locus_id
            block['span_kb'] = round((block['end'] - block['start']) / 1000, 1)
            all_blocks.append(block)

    # Summary for this locus
    print(f"\n    Total for {locus_id}:")
    print(f"      Synteny blocks: {len(all_blocks)}")

    # Save detailed hits
    hits_file = locus_output / "flanking_blast_all.tsv"
    all_hits = []
    for block in all_blocks:
        for hit in block['hits']:
            # Don't save the full sequences to TSV (too large)
            hit_copy = {k: v for k, v in hit.items() if k not in ['hit_protein_seq', 'query_protein_seq']}
            hit_copy['genome'] = block['genome']
            hit_copy['block_id'] = block['block_id']
            all_hits.append(hit_copy)

    if all_hits:
        hits_df = pd.DataFrame(all_hits)
        hits_df.to_csv(hits_file, sep='\t', index=False)
        print(f"    Saved detailed hits: {hits_file.name}")

    return all_blocks

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Detect synteny blocks using flanking proteins",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('--locus-defs', required=True, type=Path,
                        help='Path to locus_definitions.tsv (with flanking_file column)')
    parser.add_argument('--genome-db-dir', required=True, type=Path,
                        help='Directory containing BLAST databases (*.nhr files)')
    parser.add_argument('--output-dir', required=True, type=Path,
                        help='Output directory for synteny blocks')
    parser.add_argument('--evalue', type=str, default="1e-5",
                        help='E-value threshold for BLAST (default: 1e-5)')
    parser.add_argument('--max-targets', type=int, default=50,
                        help='Maximum number of BLAST targets (default: 50)')
    parser.add_argument('--max-gap-kb', type=int, default=500,
                        help='Maximum gap (kb) between hits in same block (default: 500)')
    parser.add_argument('--threads', type=int, default=16,
                        help='Number of BLAST threads (default: 16)')
    parser.add_argument('--locus', type=str,
                        help='Process only this specific locus (for SLURM arrays)')

    return parser.parse_args()

def main():
    """Main execution function."""
    args = parse_args()

    print("=" * 80)
    print("STEP 02: SYNTENY DETECTION USING FLANKING PROTEINS")
    print("=" * 80)
    print(f"\nInput files:")
    print(f"  Locus definitions: {args.locus_defs}")
    print(f"  Genome databases: {args.genome_db_dir}")
    print(f"  Output directory: {args.output_dir}")
    print(f"\nParameters:")
    print(f"  E-value: {args.evalue}")
    print(f"  Max targets: {args.max_targets}")
    print(f"  Max gap: {args.max_gap_kb} kb")
    print(f"  Threads: {args.threads}")

    # Create output directory
    args.output_dir.mkdir(exist_ok=True, parents=True)

    # Load locus definitions
    print("\n[1] Loading locus definitions...")
    loci_df = pd.read_csv(args.locus_defs, sep='\t')

    # Filter to specific locus if requested
    if args.locus:
        loci_df = loci_df[loci_df['locus_id'] == args.locus]
        if loci_df.empty:
            print(f"  ERROR: Locus '{args.locus}' not found in {args.locus_defs}")
            return
        print(f"  Processing single locus: {args.locus}")
    else:
        print(f"  Loaded {len(loci_df)} loci")

    # Find genome databases
    print("\n[2] Finding genome databases...")
    genome_dbs = {}
    for db_file in args.genome_db_dir.glob("*.nhr"):
        genome_name = db_file.stem
        genome_dbs[genome_name] = args.genome_db_dir / genome_name
    print(f"  Found {len(genome_dbs)} genome databases")

    # Process each locus
    print("\n[3] Detecting synteny blocks...")
    all_blocks = []

    for _, row in loci_df.iterrows():
        locus_id = row['locus_id']
        flanking_file = Path(row['flanking_file'])
        locus_blocks = process_locus(
            locus_id, flanking_file, genome_dbs, args.output_dir,
            args.evalue, args.max_targets, args.threads, args.max_gap_kb
        )
        all_blocks.extend(locus_blocks)

    # Filter to best block per genome per locus
    print("\n[4] Filtering to best block per genome...")
    if all_blocks:
        # Group by locus and genome, keep block with most query matches
        blocks_df = pd.DataFrame(all_blocks)
        best_blocks = blocks_df.loc[blocks_df.groupby(['locus_id', 'genome'])['num_query_matches'].idxmax()]
        print(f"  Filtered {len(blocks_df)} blocks → {len(best_blocks)} best blocks (one per genome per locus)")
        all_blocks = best_blocks.to_dict('records')

    # Save blocks per locus (for array jobs)
    print("\n[5] Saving synteny blocks...")

    # Determine which locus to use for filename (should be only one if --locus specified)
    loci_list = loci_df['locus_id'].unique().tolist()
    if len(loci_list) == 1:
        locus_id = loci_list[0]
    else:
        # Multiple loci - shouldn't happen in array mode
        locus_id = "combined"

    blocks_file = args.output_dir / f"{locus_id}_synteny_blocks.tsv"

    if all_blocks:
        # Convert to simple format for saving
        rows = []
        for block in all_blocks:
            rows.append({
                'locus_id': block['locus_id'],
                'genome': block['genome'],
                'block_id': block['block_id'],
                'scaffold': block['scaffold'],
                'strand': block['strand'],
                'start': block['start'],
                'end': block['end'],
                'span_kb': block['span_kb'],
                'num_target_proteins': block['num_target_proteins'],
                'num_query_matches': block['num_query_matches']
            })

        blocks_df = pd.DataFrame(rows)
        blocks_df.to_csv(blocks_file, sep='\t', index=False)
        print(f"  Saved {len(blocks_df)} blocks to {blocks_file.name}")

    # Summary
    print("\n" + "=" * 80)
    print("SYNTENY DETECTION COMPLETE")
    print("=" * 80)
    print(f"\nTotal synteny blocks found: {len(all_blocks)}")

    if all_blocks:
        blocks_df = pd.DataFrame([{
            'locus_id': b['locus_id'],
            'genome': b['genome'],
            'num_query_matches': b.get('num_query_matches', 0)
        } for b in all_blocks])

        print(f"\nBlocks by locus:")
        for locus_id in loci_list:
            locus_blocks = blocks_df[blocks_df['locus_id'] == locus_id]
            print(f"  {locus_id}: {len(locus_blocks)} blocks")

    print(f"\nOutputs saved to: {args.output_dir}")
    print("\nNext step: 03_filter_blocks.py")

if __name__ == "__main__":
    main()