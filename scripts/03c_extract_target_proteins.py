#!/usr/bin/env python3
"""
Phase 3c: Extract target proteins from proteomes.

Reads locus_definitions.tsv and extracts the target protein for each locus
from the appropriate proteome FASTA file, using GFF annotations to map
LOC IDs to protein IDs.

Usage:
    python 03c_extract_target_proteins.py \\
        --locus-defs <path/to/locus_definitions.tsv> \\
        --proteome-dir <path/to/proteomes> \\
        --genome-dir <path/to/genomes> \\
        --output-dir <path/to/extracted_proteins>
"""

from pathlib import Path
import pandas as pd
from Bio import SeqIO
import argparse

# Genome metadata
GENOME_INFO = {
    'BK': {
        'proteome': 'GCF_010883055.1_Bkinseyi_NCBI.faa',
        'gff': 'GCF_010883055.1_Bkinseyi.gff'
    },
    'LB': {
        'proteome': 'GCF_019393585.1_Lboulardi_NCBI.faa',
        'gff': 'GCF_019393585.1_Lboulardi.gff'
    },
    'TR': {
        'proteome': 'GCA_020615435.1.faa',
        'gff': 'GCA_020615435.1.gff3'
    },
    'DR': {
        'proteome': 'GCA_030998225.1.faa',
        'gff': 'GCA_030998225.1.gff3'
    }
}

def parse_gff_for_protein_id(gff_file, loc_id):
    """Extract protein ID (XP_* or similar) for a given LOC ID from GFF."""
    protein_id = None

    with open(gff_file) as f:
        for line in f:
            if line.startswith('#'):
                continue

            if loc_id in line and 'CDS' in line:
                # Extract protein_id from attributes
                fields = line.strip().split('\t')
                if len(fields) >= 9:
                    attributes = fields[8]

                    # Look for protein_id in attributes
                    for attr in attributes.split(';'):
                        if 'protein_id=' in attr:
                            protein_id = attr.split('=')[1].strip()
                            return protein_id
                        elif 'Name=' in attr and (attr.startswith('Name=XP_') or attr.startswith('Name=NP_')):
                            protein_id = attr.split('=')[1].strip()
                            return protein_id

    return protein_id

def extract_protein_from_proteome(proteome_file, protein_id):
    """Extract a protein sequence by ID from proteome FASTA."""
    with open(proteome_file) as f:
        for record in SeqIO.parse(f, 'fasta'):
            # Check if protein_id matches record ID or is in description
            if protein_id in record.id or protein_id in record.description:
                return record

    return None

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Extract target proteins from proteomes",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('--locus-defs', required=True, type=Path,
                        help='Path to locus_definitions.tsv')
    parser.add_argument('--proteome-dir', required=True, type=Path,
                        help='Directory containing proteome FASTA files')
    parser.add_argument('--genome-dir', required=True, type=Path,
                        help='Directory containing genome GFF annotation files')
    parser.add_argument('--output-dir', required=True, type=Path,
                        help='Output directory for extracted target proteins')

    return parser.parse_args()

def main():
    """Main execution function."""
    args = parse_args()

    print("=" * 80)
    print("PHASE 3c: EXTRACT TARGET PROTEINS")
    print("=" * 80)
    print(f"\nInput files:")
    print(f"  Locus definitions: {args.locus_defs}")
    print(f"  Proteome directory: {args.proteome_dir}")
    print(f"  Genome directory: {args.genome_dir}")
    print(f"  Output directory: {args.output_dir}")

    # Create output directory
    args.output_dir.mkdir(exist_ok=True, parents=True)

    # Load locus definitions
    print("\n[1] Loading locus definitions...")
    loci_df = pd.read_csv(args.locus_defs, sep='\t')
    print(f"  Loaded {len(loci_df)} loci")

    # Extract target protein for each locus
    print("\n[2] Extracting target proteins...")

    extracted_count = 0
    failed_count = 0

    for _, locus in loci_df.iterrows():
        locus_id = locus['locus_id']
        genome = locus['genome']
        target_gene = locus['target_gene']  # LOC ID

        print(f"\n  Processing {locus_id}...")
        print(f"    Genome: {genome}")
        print(f"    Target gene: {target_gene}")

        # Check if genome is supported
        if genome not in GENOME_INFO:
            print(f"    WARNING: Genome {genome} not in GENOME_INFO, skipping")
            failed_count += 1
            continue

        # Get genome files
        genome_info = GENOME_INFO[genome]
        proteome_file = args.proteome_dir / genome_info['proteome']
        gff_file = args.genome_dir / genome_info['gff']

        # Check files exist
        if not proteome_file.exists():
            print(f"    ERROR: Proteome file not found: {proteome_file}")
            failed_count += 1
            continue

        if not gff_file.exists():
            print(f"    ERROR: GFF file not found: {gff_file}")
            failed_count += 1
            continue

        # Step 1: Map LOC ID to protein ID using GFF
        print(f"    Parsing GFF to find protein ID...")
        protein_id = parse_gff_for_protein_id(gff_file, target_gene)

        if not protein_id:
            print(f"    WARNING: Could not find protein ID for {target_gene} in GFF")
            failed_count += 1
            continue

        print(f"    Found protein ID: {protein_id}")

        # Step 2: Extract protein from proteome
        print(f"    Extracting from proteome...")
        protein_record = extract_protein_from_proteome(proteome_file, protein_id)

        if not protein_record:
            print(f"    WARNING: Could not find {protein_id} in proteome")
            failed_count += 1
            continue

        # Step 3: Save to locus-specific file
        locus_dir = args.output_dir / locus_id
        locus_dir.mkdir(exist_ok=True, parents=True)

        output_file = locus_dir / f"{locus_id}_targets.faa"
        with open(output_file, 'w') as f:
            SeqIO.write(protein_record, f, 'fasta')

        print(f"    ✓ Saved to {output_file}")
        extracted_count += 1

    # Summary
    print("\n" + "=" * 80)
    print("TARGET PROTEIN EXTRACTION COMPLETE")
    print("=" * 80)
    print(f"\nSuccessfully extracted: {extracted_count}/{len(loci_df)}")
    if failed_count > 0:
        print(f"Failed: {failed_count}/{len(loci_df)}")

    print(f"\nOutput directory: {args.output_dir}")
    print("\nNext step: 04_blast_targets.py")

if __name__ == "__main__":
    main()
