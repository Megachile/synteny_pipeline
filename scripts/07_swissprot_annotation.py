#!/usr/bin/env python3
"""
Step 07: Annotate flanking proteins with SwissProt.

Searches extracted hit sequences from each genome against SwissProt
to get functional annotations for comparative analysis.

Usage:
    python 07_swissprot_annotation.py \\
        --synteny-dir <path/to/02_synteny_blocks> \\
        --swissprot-db <path/to/swissprot.dmnd> \\
        --output-dir <path/to/swissprot_annotations> \\
        [--evalue 1e-5] \\
        [--threads 16]
"""

from pathlib import Path
import subprocess
import pandas as pd
import xml.etree.ElementTree as ET
from Bio import SeqIO
import argparse
import re

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

def annotate_locus_hits(locus_id, position_to_protein, synteny_dir, output_dir, swissprot_db, evalue, threads, use_diamond):
    """Annotate hit sequences for one locus."""

    print(f"\n  Processing {locus_id}...", flush=True)

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
        print(f"    {genome_id}...", end='', flush=True)

        # Run SwissProt search (DIAMOND or BLAST)
        swissprot_xml = locus_output / f"{genome_id}_swissprot.xml"

        if swissprot_xml.exists():
            print(" using existing result", end='', flush=True)
        else:
            if run_diamond(hit_file, swissprot_xml, swissprot_db, evalue, threads, use_diamond):
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
                'bk_protein': bk_protein_annotation,  # SwissProt functional name
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
        description="Annotate flanking proteins with SwissProt",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('--synteny-dir', required=True, type=Path,
                        help='Path to 02_synteny_blocks directory')
    parser.add_argument('--swissprot-db', required=True, type=Path,
                        help='Path to SwissProt DIAMOND database (.dmnd)')
    parser.add_argument('--output-dir', required=True, type=Path,
                        help='Output directory for annotations')
    parser.add_argument('--evalue', type=str, default="1e-5",
                        help='E-value threshold (default: 1e-5)')
    parser.add_argument('--threads', type=int, default=16,
                        help='Number of threads for DIAMOND (default: 16)')

    return parser.parse_args()

def main():
    """Main execution function."""
    args = parse_args()

    print("=" * 80, flush=True)
    print("STEP 07: SWISSPROT ANNOTATION OF HIT SEQUENCES", flush=True)
    print("=" * 80, flush=True)
    print(f"\nInput files:")
    print(f"  Synteny directory: {args.synteny_dir}")
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
        if use_diamond:
            print("\n  NOTE: Install DIAMOND SwissProt database:")
            print("    1. Download: ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz")
            print("    2. Extract: gunzip uniprot_sprot.fasta.gz")
            print("    3. Format: diamond makedb --in uniprot_sprot.fasta --db swissprot.dmnd")
        else:
            print("\n  NOTE: Install SwissProt for functional annotations:")
            print("    1. Download: ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz")
            print("    2. Extract and format: makeblastdb -in uniprot_sprot.fasta -dbtype prot")
        return

    # Load locus definitions from synteny directory
    print("\n[2] Loading locus definitions...", flush=True)
    loci_def_file = args.synteny_dir.parent / "locus_definitions.tsv"

    if not loci_def_file.exists():
        print(f"  ERROR: Locus definitions not found: {loci_def_file}", flush=True)
        print("  Expected at parent of synteny-dir", flush=True)
        return

    loci_df = pd.read_csv(loci_def_file, sep='\t')
    loci = loci_df['locus_id'].tolist()
    print(f"  Loaded {len(loci)} loci", flush=True)

    # Process each locus
    print("\n[3] Annotating hit sequences...", flush=True)
    all_annotations = []

    for locus_id in loci:
        # Load position → protein mapping for this locus
        position_to_protein = load_flanking_protein_mapping(locus_id, args.synteny_dir)
        print(f"  Loaded {len(position_to_protein)} protein position mappings for {locus_id}", flush=True)

        locus_annotations = annotate_locus_hits(
            locus_id, position_to_protein, args.synteny_dir, args.output_dir,
            args.swissprot_db, args.evalue, args.threads, use_diamond
        )
        all_annotations.extend(locus_annotations)

    # Save all annotations
    print("\n[4] Saving annotations...", flush=True)
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
