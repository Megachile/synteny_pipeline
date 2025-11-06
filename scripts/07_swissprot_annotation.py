#!/usr/bin/env python3
"""
Step 07: Annotate flanking proteins with SwissProt - ONLY proteins in synteny blocks.

CRITICAL FIX: Only annotates proteins that are actually members of the filtered synteny blocks,
not all DIAMOND hits. Uses synteny_blocks_filtered.tsv to identify which block to use per genome,
then filters flanking_blast_all.tsv to get only proteins in those blocks.

Usage:
    python 07_swissprot_annotation.py \
        --synteny-dir <path/to/02_synteny_blocks> \
        --filtered-blocks <path/to/synteny_blocks_filtered.tsv> \
        --swissprot-db <path/to/swissprot.dmnd> \
        --output-dir <path/to/swissprot_annotations> \
        [--evalue 1e-5] \
        [--threads 16]
"""

from pathlib import Path
import subprocess
import pandas as pd
import xml.etree.ElementTree as ET
from Bio import SeqIO
import argparse
import re
from collections import defaultdict

def run_diamond(query_file, output_xml, swissprot_db, evalue, threads, use_diamond=True):
    """Run DIAMOND against SwissProt database (100-1000x faster than BLASTP)."""

    if use_diamond:
        cmd = [
            'diamond', 'blastp',
            '--query', str(query_file),
            '--db', str(swissprot_db),
            '--outfmt', '5',  # XML format (same as BLAST)
            '--evalue', str(evalue),
            '--max-target-seqs', '1',  # Only need top hit
            '--threads', str(threads),
            '--quiet'  # Suppress progress messages
        ]
    else:
        # Fall back to traditional BLASTP
        cmd = [
            'blastp',
            '-query', str(query_file),
            '-db', str(swissprot_db).replace('.dmnd', ''),  # Remove .dmnd extension for BLAST
            '-outfmt', '5',  # XML format
            '-evalue', str(evalue),
            '-max_target_seqs', '1',  # Only need top hit
            '-num_threads', str(threads)
        ]

    with open(output_xml, 'w') as f:
        result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, text=True)

    return result.returncode == 0

def parse_swissprot_xml(xml_file):
    """Parse SwissProt BLAST results."""
    annotations = {}

    try:
        tree = ET.parse(xml_file)
        root = tree.getroot()

        for iteration in root.findall('.//Iteration'):
            query_def = iteration.find('Iteration_query-def').text
            # Extract query ID (e.g., "GCA_XXXXXX_U1_scaffold_1000-2000")
            query_id = query_def.split()[0]

            # Check if hit found
            message = iteration.find('Iteration_message')
            if message is not None and message.text == 'No hits found':
                annotations[query_id] = {
                    'query_id': query_id,
                    'swissprot_id': 'No match',
                    'swissprot_description': 'No SwissProt match',
                    'percent_identity': 0,
                    'evalue': 1.0,
                    'bitscore': 0
                }
                continue

            # Get best hit
            hit = iteration.find('.//Hit')
            if hit is not None:
                hit_def = hit.find('Hit_def').text
                hit_id = hit.find('Hit_id').text

                # Get HSP details
                hsp = hit.find('.//Hsp')
                if hsp is not None:
                    identity = int(hsp.find('Hsp_identity').text)
                    align_len = int(hsp.find('Hsp_align-len').text)
                    evalue = float(hsp.find('Hsp_evalue').text)
                    bitscore = float(hsp.find('Hsp_bit-score').text)

                    annotations[query_id] = {
                        'query_id': query_id,
                        'swissprot_id': hit_id.split('|')[1] if '|' in hit_id else hit_id,
                        'swissprot_description': hit_def,
                        'percent_identity': identity / align_len * 100,
                        'evalue': evalue,
                        'bitscore': bitscore
                    }

    except Exception as e:
        print(f"    Parse error: {e}")

    return annotations

def load_synteny_block_hits(locus_id, synteny_dir, filtered_blocks_df):
    """
    Load the specific BLAST hits that are in the filtered synteny blocks.

    Returns dict: genome_id -> set of (qseqid, sseqid, coord_start, coord_end) tuples
    """
    # Get the filtered blocks for this locus
    locus_blocks = filtered_blocks_df[filtered_blocks_df['locus_id'] == locus_id]

    if locus_blocks.empty:
        print(f"  WARNING: No filtered blocks found for {locus_id}")
        return {}

    # Load flanking_blast_all.tsv which has block assignments
    blast_file = synteny_dir / locus_id / "flanking_blast_all.tsv"
    if not blast_file.exists():
        print(f"  WARNING: Flanking BLAST file not found: {blast_file}")
        return {}

    blast_df = pd.read_csv(blast_file, sep='\t')

    # For each genome, get the SPECIFIC HITS in its best block
    genome_to_hits = {}

    for _, block_row in locus_blocks.iterrows():
        genome = block_row['genome']
        block_id = block_row['block_id']

        # Get all BLAST hits in this specific block
        block_hits = blast_df[
            (blast_df['genome'] == genome) &
            (blast_df['block_id'] == block_id)
        ]

        # Create set of (qseqid, sseqid, coord_start, coord_end) tuples
        hit_set = set()
        for _, hit in block_hits.iterrows():
            hit_set.add((
                hit['qseqid'],
                hit['sseqid'],
                int(hit['coord_start']),
                int(hit['coord_end'])
            ))

        genome_to_hits[genome] = hit_set

    return genome_to_hits

def filter_fasta_by_hits(input_fasta, hit_set, output_fasta):
    """
    Filter a FASTA file to only include sequences matching specific BLAST hits.

    Matches by (qseqid, sseqid, coords) from BLAST results against FASTA headers.

    Header format: genome_qseqid_sseqid_coords_hitN
    Example: GCA_011634705.1_XP_033209112.1|LOC117167954_JAABKJ010027599_782-1942_hit1
    Description: GCA_011634705.1 | XP_033209112.1|LOC117167954 hit | JAABKJ010027599:782-1942(+) | ...

    CRITICAL: qseqid contains a pipe (e.g., XP_033209112.1|LOC117167954), so when split by '|':
      desc_parts[0] = "genome_id ..."
      desc_parts[1] = " genome_id again "
      desc_parts[2] = " XP_033209112.1"
      desc_parts[3] = "LOC117167954 hit "
      desc_parts[4] = " scaffold:coords(strand) "
    """
    kept_seqs = []

    for record in SeqIO.parse(input_fasta, "fasta"):
        # Parse FASTA header to extract qseqid, sseqid, coordinates
        desc_parts = record.description.split('|')

        if len(desc_parts) < 5:  # Need at least 5 parts now
            continue

        try:
            # Extract qseqid by reconstructing from desc_parts[2] and [3]
            # desc_parts[2] = " XP_033209112.1"
            # desc_parts[3] = "LOC117167954 hit "
            qseqid = desc_parts[2].strip() + '|' + desc_parts[3].replace(' hit', '').strip()

            # Extract sseqid and coords from desc_parts[4] (e.g., " JAABKJ010027599:782-1942(+) ")
            scaffold_info = desc_parts[4].strip().split(':')
            if len(scaffold_info) < 2:
                continue

            sseqid = scaffold_info[0].strip()
            coords_strand = scaffold_info[1].strip()  # e.g., "782-1942(+)"

            # Extract coordinates (remove strand info)
            coords = coords_strand.split('(')[0]  # Remove "(+)" or "(-)"
            coord_parts = coords.split('-')
            if len(coord_parts) != 2:
                continue

            coord_start = int(coord_parts[0])
            coord_end = int(coord_parts[1])

            # Check if this hit is in our set
            hit_tuple = (qseqid, sseqid, coord_start, coord_end)
            if hit_tuple in hit_set:
                kept_seqs.append(record)

        except (ValueError, IndexError) as e:
            # Skip malformed headers
            continue

    if kept_seqs:
        SeqIO.write(kept_seqs, output_fasta, "fasta")

    return len(kept_seqs)

def load_flanking_protein_mapping(locus_id, synteny_dir):
    """Load flanking proteins and create position mapping for this locus."""
    from Bio import SeqIO

    # Construct path to flanking protein file
    # File is in synteny_dir/locus_id/locus_id_flanking_filtered.faa
    flanking_file = synteny_dir / locus_id / f"{locus_id}_flanking_filtered.faa"

    if not flanking_file.exists():
        print(f"    Warning: Flanking file not found: {flanking_file}")
        return {}

    # Read proteins in order
    protein_ids = []
    for record in SeqIO.parse(flanking_file, "fasta"):
        # Extract protein ID (first part before pipe or space)
        prot_id = record.id.split('|')[0].split()[0]
        protein_ids.append(prot_id)

    if len(protein_ids) == 0:
        return {}

    # Create mapping: U1, U2, ... and D1, D2, ...
    # Proteins are in order: upstream (furthest to closest), then downstream (closest to furthest)
    # So we need to split them. Assume equal up/down for now
    n_proteins = len(protein_ids)
    n_up = n_proteins // 2
    n_down = n_proteins - n_up

    position_to_protein = {}

    # Upstream proteins (in reverse order for labeling)
    for i in range(n_up):
        pos_label = f"U{n_up - i}"
        position_to_protein[pos_label] = protein_ids[i]

    # Downstream proteins
    for i in range(n_down):
        pos_label = f"D{i + 1}"
        position_to_protein[pos_label] = protein_ids[n_up + i]

    return position_to_protein

def annotate_locus_hits(locus_id, position_to_protein, synteny_dir, output_dir, filtered_blocks_df,
                        swissprot_db, evalue, threads, use_diamond):
    """Annotate hit sequences for one locus - ONLY hits in synteny blocks."""

    print(f"\n  Processing {locus_id}...", flush=True)

    # Load which BLAST hits are in the filtered synteny blocks
    print(f"    Loading synteny block hit assignments...", flush=True)
    genome_to_hits = load_synteny_block_hits(locus_id, synteny_dir, filtered_blocks_df)

    if not genome_to_hits:
        print(f"    No synteny blocks found for {locus_id}", flush=True)
        return []

    total_hits = sum(len(hits) for hits in genome_to_hits.values())
    print(f"    Found {total_hits} BLAST hits in synteny blocks across {len(genome_to_hits)} genomes", flush=True)

    # Find hit sequences directory
    hit_seqs_dir = synteny_dir / locus_id / "hit_sequences"

    if not hit_seqs_dir.exists():
        print(f"    ERROR: Hit sequences not found: {hit_seqs_dir}", flush=True)
        return []

    # Get all genome hit sequence files
    hit_files = list(hit_seqs_dir.glob("*_hit_proteins.fasta"))

    if not hit_files:
        print(f"    No hit sequence files found", flush=True)
        return []

    # Output directory
    locus_output = output_dir / locus_id
    locus_output.mkdir(exist_ok=True, parents=True)

    all_annotations = []

    # Process each genome's hits
    for hit_file in hit_files:
        genome_id = hit_file.stem.replace('_hit_proteins', '')

        if genome_id not in genome_to_hits:
            # This genome has no synteny block for this locus
            continue

        block_hit_set = genome_to_hits[genome_id]

        print(f"    {genome_id} ({len(block_hit_set)} hits in block)...", end='', flush=True)

        # Filter FASTA to only hits in the synteny block
        filtered_fasta = locus_output / f"{genome_id}_block_hits.fasta"
        n_kept = filter_fasta_by_hits(hit_file, block_hit_set, filtered_fasta)

        if n_kept == 0:
            print(f" no proteins kept after filtering - SKIP", flush=True)
            continue

        print(f" {n_kept} seqs filtered...", end='', flush=True)

        # Run SwissProt search (DIAMOND or BLAST) on filtered set
        swissprot_xml = locus_output / f"{genome_id}_swissprot.xml"

        if swissprot_xml.exists():
            print(" using existing result", end='', flush=True)
        else:
            if run_diamond(filtered_fasta, swissprot_xml, swissprot_db, evalue, threads, use_diamond):
                print(" search done", end='', flush=True)
            else:
                print(" search FAILED", flush=True)
                continue

        # Parse annotations
        annotations = parse_swissprot_xml(swissprot_xml)

        # Add genome and locus info
        for query_id, annot in annotations.items():
            # Extract BK protein ID from query_id for position mapping
            # Format: GCA_037103525.1_XP_033209112.1|LOC117167954_CM021339.1_RagTag_...
            bk_protein_id = "unknown"
            flanking_pos = "unknown"

            try:
                # Extract XP_/NP_ accession number for position lookup
                match = re.search(r'([XN]P_\d+\.\d+)', query_id)
                if match:
                    bk_protein_id = match.group(1)

                # Try to find flanking position from position_to_protein mapping
                # by reverse lookup using the protein accession
                for pos, prot_id in position_to_protein.items():
                    if prot_id == bk_protein_id:
                        flanking_pos = pos
                        break
            except Exception as e:
                print(f"    Warning: Could not parse query_id {query_id}: {e}")

            # Use SwissProt functional annotation as bk_protein
            # This gives us readable names like "Soma ferritin" instead of accessions
            bk_protein_annotation = annot.get('swissprot_description', 'No annotation')

            annot_record = {
                'locus': locus_id,
                'genome': genome_id,
                'bk_protein_id': bk_protein_id,  # Landmark protein accession (e.g., XP_033209112.1) - KEY FOR PHASE 8
                'bk_protein_annotation': bk_protein_annotation,  # SwissProt functional name (e.g., "Ferritin")
                'query_id': query_id,
                'flanking_position': flanking_pos,
                **annot
            }
            all_annotations.append(annot_record)

        print(f" - {len(annotations)} proteins annotated", flush=True)

    return all_annotations

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Annotate flanking proteins with SwissProt (ONLY proteins in synteny blocks)",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('--synteny-dir', required=True, type=Path,
                        help='Path to 02_synteny_blocks directory')
    parser.add_argument('--filtered-blocks', required=True, type=Path,
                        help='Path to synteny_blocks_filtered.tsv (from Phase 3)')
    parser.add_argument('--swissprot-db', required=True, type=Path,
                        help='Path to SwissProt DIAMOND database (.dmnd)')
    parser.add_argument('--output-dir', required=True, type=Path,
                        help='Output directory for annotations')
    parser.add_argument('--evalue', type=str, default="1e-5",
                        help='E-value threshold (default: 1e-5)')
    parser.add_argument('--threads', type=int, default=16,
                        help='Number of threads for DIAMOND (default: 16)')
    parser.add_argument('--locus', type=str,
                        help='Process only this specific locus (for SLURM arrays)')

    return parser.parse_args()

def main():
    """Main execution function."""
    args = parse_args()

    print("=" * 80, flush=True)
    print("STEP 07: SWISSPROT ANNOTATION (SYNTENY BLOCK PROTEINS ONLY)", flush=True)
    print("=" * 80, flush=True)
    print(f"\nInput files:")
    print(f"  Synteny directory: {args.synteny_dir}")
    print(f"  Filtered blocks: {args.filtered_blocks}")
    print(f"  SwissProt database: {args.swissprot_db}")
    print(f"  Output directory: {args.output_dir}")
    print(f"\nParameters:")
    print(f"  E-value threshold: {args.evalue}")
    print(f"  Threads: {args.threads}")

    # Create output directory
    args.output_dir.mkdir(exist_ok=True, parents=True)

    # Determine if using DIAMOND based on file extension
    use_diamond = str(args.swissprot_db).endswith('.dmnd')
    db_type = "DIAMOND (.dmnd)" if use_diamond else "BLAST"

    # Check if SwissProt database exists
    print("\n[1] Checking SwissProt database...", flush=True)
    if use_diamond:
        db_exists = args.swissprot_db.exists()
    else:
        db_exists = (Path(str(args.swissprot_db) + ".phr").exists() or
                     args.swissprot_db.exists())

    print(f"  Using {db_type} database", flush=True)

    if not db_exists:
        print(f"  WARNING: SwissProt database not found: {args.swissprot_db}", flush=True)
        print("  Skipping SwissProt annotation...", flush=True)
        return

    # Load filtered synteny blocks
    print("\n[2] Loading filtered synteny blocks...", flush=True)
    if not args.filtered_blocks.exists():
        print(f"  ERROR: Filtered blocks file not found: {args.filtered_blocks}", flush=True)
        return

    filtered_blocks_df = pd.read_csv(args.filtered_blocks, sep='\t')
    print(f"  Loaded {len(filtered_blocks_df)} synteny blocks", flush=True)

    # Get unique loci
    loci = filtered_blocks_df['locus_id'].unique().tolist()

    # Filter to specific locus if requested
    if args.locus:
        if args.locus in loci:
            loci = [args.locus]
            print(f"  Processing single locus: {args.locus}", flush=True)
        else:
            print(f"  ERROR: Locus '{args.locus}' not found in filtered blocks", flush=True)
            return
    else:
        print(f"  Found {len(loci)} loci", flush=True)

    # Process each locus
    print("\n[3] Annotating synteny block proteins...", flush=True)
    all_annotations = []

    for locus_id in loci:
        # Load position → protein mapping for this locus
        position_to_protein = load_flanking_protein_mapping(locus_id, args.synteny_dir)
        print(f"  Loaded {len(position_to_protein)} protein position mappings for {locus_id}", flush=True)

        locus_annotations = annotate_locus_hits(
            locus_id, position_to_protein, args.synteny_dir, args.output_dir,
            filtered_blocks_df, args.swissprot_db, args.evalue, args.threads, use_diamond
        )
        all_annotations.extend(locus_annotations)

    # Save all annotations
    print("\n[4] Saving annotations...", flush=True)

    # Use locus-specific filename if processing single locus
    if args.locus:
        output_file = args.output_dir / f"{args.locus}_swissprot_annotations.tsv"
    else:
        output_file = args.output_dir / "genome_specific_swissprot_annotations.tsv"

    if all_annotations:
        annot_df = pd.DataFrame(all_annotations)
        annot_df.to_csv(output_file, sep='\t', index=False)
        print(f"  Saved {len(annot_df)} annotations to {output_file.name}")

        # Summary
        print("\n[5] Summary:")
        print(f"  Total annotations: {len(annot_df)}")
        print(f"  Unique genomes: {annot_df['genome'].nunique()}")
        print(f"  Unique loci: {annot_df['locus'].nunique()}")
        print(f"  With SwissProt match: {len(annot_df[annot_df['swissprot_id'] != 'No match'])}")
        print(f"  No match: {len(annot_df[annot_df['swissprot_id'] == 'No match'])}")

        # Per-locus summary
        print("\n[6] Matches per locus:")
        for locus_id in loci:
            locus_annot = annot_df[annot_df['locus'] == locus_id]
            matches = len(locus_annot[locus_annot['swissprot_id'] != 'No match'])
            total = len(locus_annot)
            if total > 0:
                print(f"  {locus_id}: {matches}/{total} ({matches/total*100:.1f}%) with SwissProt match")

    else:
        # Create empty file
        pd.DataFrame(columns=['locus', 'genome', 'bk_protein', 'query_id', 'flanking_position',
                             'swissprot_id', 'swissprot_description',
                             'percent_identity', 'evalue', 'bitscore']).to_csv(output_file, sep='\t', index=False)
        print(f"  No annotations found, created empty file: {output_file.name}")

    # Overall summary
    print("\n" + "=" * 80, flush=True)
    print("SWISSPROT ANNOTATION COMPLETE", flush=True)
    print("=" * 80, flush=True)
    print(f"\nAnnotations saved to: {args.output_dir}", flush=True)
    print(f"Output file: {output_file.name}", flush=True)
    print("\nNext step: 08a_generate_locus_matrices.py", flush=True)

if __name__ == "__main__":
    main()
