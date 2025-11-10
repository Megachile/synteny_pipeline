#!/usr/bin/env python3
"""
One-off script to fix flanking_position labels in existing SwissProt annotation files.

The bug: Phase 7 was assuming equal upstream/downstream splits instead of reading
actual position labels (U1, U2, D1, D2, etc.) from FASTA headers.

This script reads the correct positions from flanking_filtered.faa files and updates
the genome_specific_swissprot_annotations.tsv files.
"""

from pathlib import Path
import pandas as pd
from Bio import SeqIO
import sys

def load_correct_positions(flanking_file):
    """Load correct protein->position mapping from flanking FASTA file."""
    protein_to_position = {}

    if not flanking_file.exists():
        print(f"    Warning: Flanking file not found: {flanking_file}")
        return protein_to_position

    for record in SeqIO.parse(flanking_file, "fasta"):
        # Extract protein ID (first part before pipe)
        prot_id = record.id.split('|')[0].split()[0]

        # Extract position label from description
        # Format: "XP_033225730.1|LOC117178413 D1 XP_033225730.1 lazarillo protein-like [species]"
        description = record.description
        parts = description.split()

        pos_label = None
        for part in parts[1:4]:  # Check first few fields after ID
            if part.startswith('U') or part.startswith('D'):
                if len(part) > 1 and part[1:].isdigit():
                    pos_label = part
                    break

        if pos_label:
            protein_to_position[prot_id] = pos_label
        else:
            print(f"    Warning: Could not parse position for {prot_id}")

    return protein_to_position


def fix_annotation_file(annot_file, synteny_dir):
    """Fix flanking_position labels in one annotation file."""

    print(f"\nProcessing: {annot_file}")

    # Read annotation file
    df = pd.read_csv(annot_file, sep='\t')

    if df.empty:
        print("  Empty file, skipping")
        return

    # Group by locus to load correct positions
    fixes_made = 0
    position_errors = 0

    for locus_id in df['locus'].unique():
        # Load correct positions from flanking file
        flanking_file = synteny_dir / locus_id / f"{locus_id}_flanking_filtered.faa"
        protein_to_position = load_correct_positions(flanking_file)

        if not protein_to_position:
            print(f"  Warning: No position mapping for {locus_id}")
            continue

        # Fix positions for this locus
        locus_mask = df['locus'] == locus_id

        for idx in df[locus_mask].index:
            bk_protein_id = df.at[idx, 'bk_protein_id']
            current_pos = df.at[idx, 'flanking_position']

            # Look up correct position
            correct_pos = protein_to_position.get(bk_protein_id, 'unknown')

            if correct_pos != current_pos:
                if current_pos != 'unknown' and correct_pos != 'unknown':
                    print(f"    {locus_id} {bk_protein_id}: {current_pos} -> {correct_pos}")
                df.at[idx, 'flanking_position'] = correct_pos
                fixes_made += 1

            if correct_pos == 'unknown':
                position_errors += 1

    # Write corrected file
    df.to_csv(annot_file, sep='\t', index=False)

    print(f"  Fixed {fixes_made} positions, {position_errors} still unknown")

    return fixes_made, position_errors


def main():
    """Main execution."""

    # Find all annotation files
    base_dir = Path("outputs")

    if not base_dir.exists():
        print(f"Error: {base_dir} not found")
        sys.exit(1)

    # Process all families
    total_fixes = 0
    total_errors = 0

    for family_dir in sorted(base_dir.glob("*_BK_LB_batch")):
        annot_file = family_dir / "07_swissprot_annotations" / "genome_specific_swissprot_annotations.tsv"
        synteny_dir = family_dir / "02_synteny_blocks"

        if not annot_file.exists():
            continue

        fixes, errors = fix_annotation_file(annot_file, synteny_dir)
        total_fixes += fixes
        total_errors += errors

    print("\n" + "="*80)
    print(f"COMPLETE: Fixed {total_fixes} position labels across all families")
    print(f"Remaining unknown positions: {total_errors}")
    print("="*80)


if __name__ == "__main__":
    main()
