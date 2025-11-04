#!/usr/bin/env python3
"""
Step 07: Annotate flanking proteins with SwissProt.

Searches extracted hit sequences from each genome against SwissProt
to get functional annotations for comparative analysis.
"""

print("DEBUG: Starting imports...", flush=True)
from pathlib import Path
print("DEBUG: Imported pathlib", flush=True)
import subprocess
print("DEBUG: Imported subprocess", flush=True)
import pandas as pd
print("DEBUG: Imported pandas", flush=True)
import xml.etree.ElementTree as ET
print("DEBUG: Imported xml", flush=True)
from Bio import SeqIO
print("DEBUG: Imported Bio.SeqIO", flush=True)
import config
print("DEBUG: Imported config", flush=True)

def run_diamond(query_file, output_xml):
    """Run DIAMOND against SwissProt database (100-1000x faster than BLASTP)."""

    if config.USE_DIAMOND:
        cmd = [
            'diamond', 'blastp',
            '--query', str(query_file),
            '--db', str(config.SWISSPROT_DB),
            '--outfmt', '5',  # XML format (same as BLAST)
            '--evalue', config.SWISSPROT_BLAST_EVALUE,
            '--max-target-seqs', '1',  # Only need top hit
            '--threads', str(config.BLAST_THREADS),
            '--quiet'  # Suppress progress messages
        ]
    else:
        # Fall back to traditional BLASTP
        cmd = [
            'blastp',
            '-query', str(query_file),
            '-db', str(config.SWISSPROT_DB).replace('.dmnd', ''),  # Remove .dmnd extension for BLAST
            '-outfmt', '5',  # XML format
            '-evalue', config.SWISSPROT_BLAST_EVALUE,
            '-max_target_seqs', '1',  # Only need top hit
            '-num_threads', str(config.BLAST_THREADS)
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

def load_flanking_protein_mapping(locus_id, loci_df):
    """Load flanking proteins and create position mapping for this locus."""
    from Bio import SeqIO

    # Find this locus in loci_df
    locus_row = loci_df[loci_df['locus_id'] == locus_id]
    if locus_row.empty:
        return {}

    flanking_file = locus_row.iloc[0]['flanking_file']
    if not Path(flanking_file).exists():
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

def annotate_locus_hits(locus_id, position_to_protein):
    """Annotate hit sequences for one locus."""

    print(f"\n  Processing {locus_id}...")

    # Find hit sequences directory
    hit_seqs_dir = config.STEP02_SYNTENY / locus_id / "hit_sequences"

    if not hit_seqs_dir.exists():
        print(f"    ERROR: Hit sequences not found: {hit_seqs_dir}")
        return []

    # Get all genome hit sequence files
    hit_files = list(hit_seqs_dir.glob("*_hit_proteins.fasta"))

    if not hit_files:
        print(f"    No hit sequence files found")
        return []

    # Output directory
    locus_output = config.STEP07_SWISSPROT / locus_id
    locus_output.mkdir(exist_ok=True, parents=True)

    all_annotations = []

    # Process each genome's hits
    for hit_file in hit_files:
        genome_id = hit_file.stem.replace('_hit_proteins', '')
        print(f"    {genome_id}...", end='')

        # Run SwissProt search (DIAMOND or BLAST)
        swissprot_xml = locus_output / f"{genome_id}_swissprot.xml"

        if swissprot_xml.exists():
            print(" using existing result", end='')
        else:
            if run_diamond(hit_file, swissprot_xml):
                print(" search done", end='')
            else:
                print(" search FAILED")
                continue

        # Parse annotations
        annotations = parse_swissprot_xml(swissprot_xml)

        # Add genome and locus info
        for query_id, annot in annotations.items():
            # Parse query ID to extract components
            # Format: GCA_XXXXXX_U1_scaffold_1000-2000
            parts = query_id.split('_')
            if len(parts) >= 3:
                flanking_pos = parts[2]  # U1, D1, etc.
            else:
                flanking_pos = "unknown"

            # Map position to BK protein ID
            bk_protein = position_to_protein.get(flanking_pos, "unknown")

            annot_record = {
                'locus': locus_id,
                'genome': genome_id,
                'bk_protein': bk_protein,
                'query_id': query_id,
                'flanking_position': flanking_pos,
                **annot
            }
            all_annotations.append(annot_record)

        print(f" - {len(annotations)} proteins annotated")

    return all_annotations

def main():
    """Main execution function."""
    print("DEBUG: Entered main() function", flush=True)
    print("DEBUG: About to print header", flush=True)
    print("=" * 80, flush=True)
    print("STEP 07: SWISSPROT ANNOTATION OF HIT SEQUENCES", flush=True)
    print("=" * 80, flush=True)

    # Create output directory
    print("DEBUG: About to create output directory", flush=True)
    config.STEP07_SWISSPROT.mkdir(exist_ok=True, parents=True)
    print("DEBUG: Created output directory", flush=True)

    # Load locus definitions
    print("\n[1] Loading locus definitions...", flush=True)
    loci_df = pd.read_csv(config.LOCI_DEFINITIONS_FILE, sep='\t')
    loci = loci_df['locus_id'].tolist()
    print(f"  Loaded {len(loci)} loci", flush=True)

    # Check if SwissProt database exists
    print("\n[2] Checking SwissProt database...", flush=True)
    if config.USE_DIAMOND:
        db_exists = config.SWISSPROT_DB.exists()
        db_type = "DIAMOND (.dmnd)"
    else:
        db_exists = (Path(str(config.SWISSPROT_DB).replace('.dmnd', '') + ".phr").exists() or
                     Path(str(config.SWISSPROT_DB).replace('.dmnd', '') + ".fasta").exists())
        db_type = "BLAST"

    print(f"  Using {db_type} database", flush=True)

    if not db_exists:
        print(f"  WARNING: SwissProt database not found: {config.SWISSPROT_DB}")
        print("  Skipping SwissProt annotation...")
        if config.USE_DIAMOND:
            print("\n  NOTE: Install DIAMOND SwissProt database:")
            print("    1. Download: ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz")
            print("    2. Extract: gunzip uniprot_sprot.fasta.gz")
            print("    3. Format: diamond makedb --in uniprot_sprot.fasta --db swissprot.dmnd")
        else:
            print("\n  NOTE: Install SwissProt for functional annotations:")
            print("    1. Download: ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz")
            print("    2. Extract and format: makeblastdb -in uniprot_sprot.fasta -dbtype prot")
        print("    3. Update config.py with database path")
        return

    # Process each locus
    print("\n[3] Annotating hit sequences...")
    all_annotations = []

    for locus_id in loci:
        # Load position → protein mapping for this locus
        position_to_protein = load_flanking_protein_mapping(locus_id, loci_df)
        print(f"  Loaded {len(position_to_protein)} protein position mappings for {locus_id}")

        locus_annotations = annotate_locus_hits(locus_id, position_to_protein)
        all_annotations.extend(locus_annotations)

    # Save all annotations
    print("\n[4] Saving annotations...")
    output_file = config.STEP07_SWISSPROT / "genome_specific_swissprot_annotations.tsv"

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
    print("\n" + "=" * 80)
    print("SWISSPROT ANNOTATION COMPLETE")
    print("=" * 80)
    print(f"\nAnnotations saved to: {config.STEP07_SWISSPROT}")
    print(f"Output file: {output_file.name}")
    print("\nNext step: 08a_generate_locus_matrices.py")

if __name__ == "__main__":
    print("DEBUG: In __main__ block, about to call main()", flush=True)
    main()
    print("DEBUG: main() returned", flush=True)
