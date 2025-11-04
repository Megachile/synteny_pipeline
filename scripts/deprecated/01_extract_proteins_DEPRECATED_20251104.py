#!/usr/bin/env python3
"""
Step 01: Extract flanking and target proteins from reference genome.

This script extracts:
1. Flanking proteins (upstream and downstream) for synteny detection
2. Target proteins for homolog search

Outputs both individual files (for debugging) and concatenated files (for pipeline).
"""

from pathlib import Path
import pandas as pd
from Bio import SeqIO
import csv
import config

def extract_proteins(loci_df, proteins_file, output_dir):
    """Extract proteins from reference genome based on locus definitions."""

    # Load all proteins into memory
    print("  Loading reference proteins...")
    proteins_dict = {}
    for record in SeqIO.parse(proteins_file, "fasta"):
        # Extract protein ID (e.g., XP_033228573.1 from full header)
        protein_id = record.id.split()[0]
        proteins_dict[protein_id] = record
    print(f"  Loaded {len(proteins_dict)} proteins")

    # Process each locus
    for _, locus_row in loci_df.iterrows():
        locus_id = locus_row['locus_id']
        print(f"\n  Processing {locus_id}...")

        # Create locus output directory
        locus_dir = output_dir / locus_id
        locus_dir.mkdir(exist_ok=True, parents=True)

        # Extract upstream proteins
        upstream_proteins = locus_row['true_upstream_proteins'].split(',')
        upstream_names = locus_row['true_upstream_names'].split(',')

        all_flanking_records = []  # For concatenated file

        for i, (protein_id, protein_name) in enumerate(zip(upstream_proteins, upstream_names), 1):
            if protein_id in proteins_dict:
                record = proteins_dict[protein_id]
                # Rename for clarity
                record.id = f"{locus_id}_U{i}"
                record.description = f"Upstream {i} - {protein_name} - {protein_id}"

                # Save individual file
                individual_file = locus_dir / f"{locus_id}_U{i}.fasta"
                SeqIO.write([record], individual_file, "fasta")

                # Add to concatenated list
                all_flanking_records.append(record)

                print(f"    Extracted U{i}: {protein_name} ({protein_id})")
            else:
                print(f"    WARNING: U{i} not found: {protein_name} ({protein_id})")

        # Extract downstream proteins
        downstream_proteins = locus_row['true_downstream_proteins'].split(',')
        downstream_names = locus_row['true_downstream_names'].split(',')

        for i, (protein_id, protein_name) in enumerate(zip(downstream_proteins, downstream_names), 1):
            if protein_id in proteins_dict:
                record = proteins_dict[protein_id]
                # Rename for clarity
                record.id = f"{locus_id}_D{i}"
                record.description = f"Downstream {i} - {protein_name} - {protein_id}"

                # Save individual file
                individual_file = locus_dir / f"{locus_id}_D{i}.fasta"
                SeqIO.write([record], individual_file, "fasta")

                # Add to concatenated list
                all_flanking_records.append(record)

                print(f"    Extracted D{i}: {protein_name} ({protein_id})")
            else:
                print(f"    WARNING: D{i} not found: {protein_name} ({protein_id})")

        # Save concatenated flanking proteins file
        concatenated_flanking = locus_dir / f"{locus_id}_flanking.faa"
        SeqIO.write(all_flanking_records, concatenated_flanking, "fasta")
        print(f"    Saved concatenated flanking: {concatenated_flanking.name}")

        # Extract target proteins
        target_proteins = locus_row['target_proteins'].split(',')
        target_names = locus_row['target_names'].split(',')

        all_target_records = []  # For concatenated file

        for i, (protein_id, protein_name) in enumerate(zip(target_proteins, target_names), 1):
            if protein_id in proteins_dict:
                record = proteins_dict[protein_id]
                # Rename for clarity
                record.id = f"{locus_id}_target_{i}"
                record.description = f"Target {i} - {protein_name} - {protein_id}"

                # Save individual file
                individual_file = locus_dir / f"{locus_id}_target_{i}.fasta"
                SeqIO.write([record], individual_file, "fasta")

                # Add to concatenated list
                all_target_records.append(record)

                print(f"    Extracted Target {i}: {protein_name} ({protein_id})")
            else:
                print(f"    WARNING: Target {i} not found: {protein_name} ({protein_id})")

        # Save concatenated target proteins file
        concatenated_targets = locus_dir / f"{locus_id}_targets.faa"
        SeqIO.write(all_target_records, concatenated_targets, "fasta")
        print(f"    Saved concatenated targets: {concatenated_targets.name}")

        print(f"    Total: {len(all_flanking_records)} flanking, {len(all_target_records)} target proteins")

def main():
    """Main execution function."""
    print("=" * 80)
    print("STEP 01: EXTRACT PROTEINS FROM REFERENCE GENOME")
    print("=" * 80)

    # Create output directory
    config.STEP01_PROTEINS.mkdir(exist_ok=True, parents=True)

    # Load locus definitions
    print("\n[1] Loading locus definitions...")
    loci_df = pd.read_csv(config.LOCI_DEFINITIONS_FILE, sep='\t')
    print(f"  Loaded {len(loci_df)} loci")

    # Extract proteins
    print("\n[2] Extracting proteins...")
    extract_proteins(loci_df, config.BK_PROTEINS_FILE, config.STEP01_PROTEINS)

    # Summary
    print("\n" + "=" * 80)
    print("EXTRACTION COMPLETE")
    print("=" * 80)
    print(f"\nOutputs saved to: {config.STEP01_PROTEINS}")
    print("\nFor each locus, created:")
    print("  - Individual protein files (U1.fasta, D1.fasta, etc.)")
    print("  - Concatenated flanking proteins (*_flanking.faa)")
    print("  - Concatenated target proteins (*_targets.faa)")
    print("\nNext step: 02_synteny_detection.py")

if __name__ == "__main__":
    main()