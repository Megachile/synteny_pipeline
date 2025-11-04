#!/usr/bin/env python3
"""
Step 04: BLAST target genes against all genomes.

Uses tBLASTn to search target proteins against nucleotide databases.
Clusters hits into loci based on proximity.
"""

from pathlib import Path
import subprocess
import pandas as pd
from collections import defaultdict
import xml.etree.ElementTree as ET
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import config
import sys
sys.path.insert(0, str(Path(__file__).parent))
import extract_with_exonerate

def run_tblastn(query_file, genome_db, output_xml):
    """Run tBLASTn search for target proteins."""
    cmd = [
        'tblastn',  # Protein query vs nucleotide database
        '-query', str(query_file),
        '-db', str(genome_db),
        '-outfmt', '5',  # XML format
        '-evalue', config.TARGET_BLAST_EVALUE,
        '-max_target_seqs', str(config.TARGET_BLAST_MAX_TARGETS),
        '-num_threads', str(config.BLAST_THREADS)
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
            query_id = iteration.find('Iteration_query-def').text.split()[0]

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

                    # Extract HSP sequences (already translated by tBLASTn)
                    hsp_hseq = hsp.find('Hsp_hseq').text  # Translated hit sequence
                    hsp_qseq = hsp.find('Hsp_qseq').text  # Query sequence
                    query_len = int(iteration.find('Iteration_query-len').text)  # Reference protein length

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
                        'hit_protein_seq': hsp_hseq,
                        'query_protein_seq': hsp_qseq,
                        'query_length': query_len
                    })

    except Exception as e:
        print(f"    Parse error: {e}")

    return hits


def save_target_sequences(hits, genome_id, locus_id, locus_output):
    """
    Save target protein sequences from BLAST HSPs.

    Args:
        hits: List of BLAST hit dictionaries (with hit_protein_seq field)
        genome_id: Standardized genome ID
        locus_id: Locus identifier
        locus_output: Output directory for this locus

    Returns:
        Number of sequences saved
    """
    seqs_dir = config.STEP04_TARGETS / "target_sequences"
    seqs_dir.mkdir(exist_ok=True, parents=True)

    sequences = []

    for i, hit in enumerate(hits, 1):
        # Remove gaps from the aligned sequence
        protein_seq = hit['hit_protein_seq'].replace('-', '')

        # Create SeqRecord
        seq_id = f"{genome_id}_{locus_id}_{hit['query_id']}_{hit['scaffold']}_{hit['start']}-{hit['end']}_hit{i}"
        seq_record = SeqRecord(
            seq=Seq(protein_seq),
            id=seq_id,
            description=f"{genome_id} | {locus_id} | {hit['query_id']} hit | {hit['scaffold']}:{hit['start']}-{hit['end']}({hit['strand']}) | pident={hit['pident']:.1f}% | len={len(protein_seq)}aa"
        )

        sequences.append(seq_record)

    # Save sequences for this genome+locus
    if sequences:
        output_fasta = seqs_dir / f"{genome_id}_{locus_id}_target_proteins.fasta"
        SeqIO.write(sequences, output_fasta, "fasta")

    return len(sequences)

def cluster_into_loci(hits, locus_id, gene_family, max_gap_kb=10):
    """Cluster BLAST hits into loci based on proximity."""

    # Group by scaffold
    grouped = defaultdict(list)
    for hit in hits:
        grouped[hit['scaffold']].append(hit)

    loci = []
    locus_num = 1

    for scaffold, scaffold_hits in grouped.items():
        # Sort by position
        scaffold_hits = sorted(scaffold_hits, key=lambda x: x['start'])

        # Cluster overlapping/nearby hits
        clusters = []
        current_cluster = [scaffold_hits[0]]

        for hit in scaffold_hits[1:]:
            # Check if overlaps or within max_gap of current cluster
            cluster_end = max(h['end'] for h in current_cluster)

            # Check if reading frame matches current cluster (for tBLASTn)
            # HSPs in different frames represent different reading frames
            cluster_frame = current_cluster[0].get('frame', None)
            hit_frame = hit.get('frame', None)
            frame_compatible = (cluster_frame is None or hit_frame is None or cluster_frame == hit_frame)

            if hit['start'] <= cluster_end + (max_gap_kb * 1000) and frame_compatible:
                # Add to current cluster
                current_cluster.append(hit)
            else:
                # Start new cluster
                clusters.append(current_cluster)
                current_cluster = [hit]

        clusters.append(current_cluster)

        # Create locus for each cluster
        for cluster in clusters:
            # Get cluster boundaries
            cluster_start = min(h['start'] for h in cluster)
            cluster_end = max(h['end'] for h in cluster)

            # Get best e-value
            best_evalue = min(h['evalue'] for h in cluster)

            # Get all unique queries
            unique_queries = set(h['query_id'] for h in cluster)

            # Determine strand (majority vote)
            plus_count = sum(1 for h in cluster if h['strand'] == '+')
            minus_count = sum(1 for h in cluster if h['strand'] == '-')
            strand = '+' if plus_count >= minus_count else '-'

            # Get reading frame (should be same for all in cluster now)
            frame = cluster[0].get('frame', None)

            loci.append({
                'locus_name': f"{locus_id}_target_{locus_num:03d}",
                'gene_family': gene_family,
                'scaffold': scaffold,
                'strand': strand,
                'frame': frame,
                'start': cluster_start,
                'end': cluster_end,
                'num_hits': len(cluster),
                'num_queries': len(unique_queries),
                'query_ids': ','.join(sorted(unique_queries)),
                'best_evalue': best_evalue
            })

            locus_num += 1

    return loci

def process_locus_targets(locus_id, gene_family, genome_dbs):
    """Process target gene BLAST for one locus."""

    print(f"\n  Processing {locus_id} targets...")

    # Input file
    targets_file = config.STEP01_PROTEINS / locus_id / f"{locus_id}_targets.faa"

    if not targets_file.exists():
        print(f"    ERROR: Target proteins not found: {targets_file}")
        return []

    # Output directory
    locus_output = config.STEP04_TARGETS / locus_id
    locus_output.mkdir(exist_ok=True, parents=True)
    blast_dir = locus_output / "blast_xml"
    blast_dir.mkdir(exist_ok=True)

    all_loci = []
    all_hits = []

    # Run BLAST against each genome
    for genome_name, genome_db in genome_dbs.items():
        output_xml = blast_dir / f"{genome_name}.xml"

        # Run tBLASTn
        if output_xml.exists():
            print(f"    {genome_name}: Using existing BLAST results")
        else:
            if config.VERBOSE:
                print(f"    {genome_name}: Running tBLASTn...", end='')

            if run_tblastn(targets_file, genome_db, output_xml):
                if config.VERBOSE:
                    print(" done")
            else:
                if config.VERBOSE:
                    print(" FAILED")
                continue

        # Parse hits
        hits = parse_blast_xml(output_xml)

        # Extract genome ID
        if genome_name.startswith('GCA_') or genome_name.startswith('GCF_'):
            genome_id = genome_name.split('_')[0] + '_' + genome_name.split('_')[1]
        else:
            genome_id = genome_name

        if not hits:
            continue

        # Add genome info
        for hit in hits:
            hit['genome'] = genome_name
            hit['locus_id'] = locus_id
            all_hits.append(hit)

        # Cluster into loci
        loci = cluster_into_loci(hits, locus_id, gene_family, config.TARGET_CLUSTER_GAP_KB)

        # Find genome FASTA file for Exonerate extraction
        genome_fasta = config.RAGTAG_FASTA_DIR / genome_name / "ragtag.scaffold.fasta"
        if not genome_fasta.exists():
            print(f"    Warning: Genome FASTA not found for {genome_name}, skipping Exonerate")
            # Add genome info without extraction
            for locus in loci:
                locus['genome'] = genome_name
                locus['parent_locus'] = locus_id
                all_loci.append(locus)
            continue

        # Extract gene structures for each target locus using Exonerate
        total_genes = 0
        for locus in loci:
            # Create a pseudo-block structure for extract_block_genes
            # Get all hits for this locus
            locus_hits = [h for h in hits if
                         h['scaffold'] == locus['scaffold'] and
                         h['start'] >= locus['start'] and
                         h['end'] <= locus['end']]

            block_structure = {
                'block_id': locus['locus_name'],
                'scaffold': locus['scaffold'],
                'start': locus['start'],
                'end': locus['end'],
                'strand': locus['strand'],
                'hits': locus_hits
            }

            # Output directory for this target locus
            locus_exonerate_output = locus_output / "target_genes_exonerate" / locus['locus_name']

            # Extract genes with Exonerate
            extracted_genes = extract_with_exonerate.extract_block_genes(
                block=block_structure,
                query_protein_file=targets_file,
                genome_fasta=genome_fasta,
                output_dir=locus_exonerate_output
            )

            # Store extracted genes in locus
            locus['extracted_genes'] = extracted_genes
            locus['num_genes_extracted'] = len(extracted_genes)
            total_genes += len(extracted_genes)

            # Add genome info
            locus['genome'] = genome_name
            locus['parent_locus'] = locus_id
            all_loci.append(locus)

        if config.VERBOSE and total_genes > 0:
            print(f" - Exonerate extracted {total_genes} target genes from {len(loci)} loci")

    print(f"    Found {len(all_loci)} target loci across {len(set(l['genome'] for l in all_loci))} genomes")

    # Save detailed hits
    if all_hits:
        hits_file = locus_output / "target_blast_hits.tsv"
        hits_df = pd.DataFrame(all_hits)
        hits_df.to_csv(hits_file, sep='\t', index=False)

    return all_loci

def main():
    """Main execution function."""
    print("=" * 80)
    print("STEP 04: BLAST TARGET GENES")
    print("=" * 80)

    # Create output directory
    config.STEP04_TARGETS.mkdir(exist_ok=True, parents=True)

    # Load locus definitions
    print("\n[1] Loading locus definitions...")
    loci_df = pd.read_csv(config.LOCI_DEFINITIONS_FILE, sep='\t')
    print(f"  Loaded {len(loci_df)} loci")

    # Find genome databases
    print("\n[2] Finding genome databases...")
    genome_dbs = {}
    for db_file in config.RAGTAG_DB_DIR.glob("*.nhr"):
        genome_name = db_file.stem
        genome_dbs[genome_name] = config.RAGTAG_DB_DIR / genome_name
    print(f"  Found {len(genome_dbs)} genome databases")

    # Group loci by target protein to avoid redundant BLASTs
    print("\n[3] Grouping loci by target protein...")
    from Bio import SeqIO

    target_to_loci = {}  # Maps target protein sequence to list of loci using it

    for _, locus_row in loci_df.iterrows():
        locus_id = locus_row['locus_id']
        targets_file = config.STEP01_PROTEINS / locus_id / f"{locus_id}_targets.faa"

        if not targets_file.exists():
            print(f"  WARNING: Target file not found for {locus_id}")
            continue

        # Read target protein sequence
        with open(targets_file) as f:
            records = list(SeqIO.parse(f, 'fasta'))
            if records:
                target_seq = str(records[0].seq)
                target_id = records[0].id

                if target_seq not in target_to_loci:
                    target_to_loci[target_seq] = {
                        'protein_id': target_id,
                        'loci': [],
                        'gene_family': locus_row['gene_family']
                    }
                target_to_loci[target_seq]['loci'].append(locus_id)

    print(f"  Found {len(target_to_loci)} unique target proteins for {len(loci_df)} loci")
    for i, (seq, info) in enumerate(target_to_loci.items(), 1):
        print(f"    Protein {i} ({info['protein_id']}): used by {len(info['loci'])} loci")

    # BLAST each unique target protein once
    print("\n[4] BLASTing unique target proteins...")
    all_loci = []

    for target_seq, info in target_to_loci.items():
        protein_id = info['protein_id']
        loci_using_protein = info['loci']
        gene_family = info['gene_family']

        print(f"\n  Processing {protein_id} (used by {len(loci_using_protein)} loci)...")

        # Use the first locus's target file as the query
        first_locus = loci_using_protein[0]
        targets_file = config.STEP01_PROTEINS / first_locus / f"{first_locus}_targets.faa"

        # Create shared BLAST output directory
        shared_blast_dir = config.STEP04_TARGETS / f"shared_blast_{protein_id}"
        shared_blast_dir.mkdir(exist_ok=True, parents=True)

        # Run BLAST once for this protein
        for genome_name, genome_db in genome_dbs.items():
            output_xml = shared_blast_dir / f"{genome_name}.xml"

            if not output_xml.exists():
                if config.VERBOSE:
                    print(f"    {genome_name}: Running tBLASTn...", end='')

                if run_tblastn(targets_file, genome_db, output_xml):
                    if config.VERBOSE:
                        print(" done")
                else:
                    if config.VERBOSE:
                        print(" FAILED")

        # Now process results for ALL loci using this protein
        for locus_id in loci_using_protein:
            locus_targets = process_locus_targets_from_blast(
                locus_id, gene_family, genome_dbs, shared_blast_dir
            )
            all_loci.extend(locus_targets)

    # Save all target loci
    print("\n[4] Saving target loci...")
    loci_file = config.STEP04_TARGETS / "all_target_loci.tsv"

    if all_loci:
        loci_df = pd.DataFrame(all_loci)
        loci_df.to_csv(loci_file, sep='\t', index=False)
        print(f"  Saved {len(loci_df)} target loci to {loci_file.name}")

        # Summary by parent locus
        print("\n[5] Summary by locus:")
        for parent_locus in loci_df['parent_locus'].unique():
            parent_loci = loci_df[loci_df['parent_locus'] == parent_locus]
            genomes_with_targets = parent_loci['genome'].nunique()
            total_loci = len(parent_loci)
            print(f"  {parent_locus}: {total_loci} loci in {genomes_with_targets} genomes")

    # Overall summary
    print("\n" + "=" * 80)
    print("TARGET BLAST COMPLETE")
    print("=" * 80)

    if all_loci:
        print(f"\nTotal target loci found: {len(all_loci)}")
        print(f"Genomes with targets: {loci_df['genome'].nunique()}")
        print(f"Average loci per genome: {len(all_loci) / loci_df['genome'].nunique():.1f}")
    else:
        print("\nNo target loci found")

    print(f"\nOutputs saved to: {config.STEP04_TARGETS}")
    print("\nNext step: 05_classify_targets.py")

if __name__ == "__main__":
    main()