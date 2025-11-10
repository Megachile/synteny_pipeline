#!/usr/bin/env python3
"""
Summarize Phase 2 (Synteny Detection) results across all gene families.

Phase 2 uses tBLASTn to search flanking proteins against all genomes,
then clusters hits into synteny blocks based on genomic proximity.

Generates a TSV report with:
- input_loci: Number of loci from Phase 1
- total_blocks: Total synteny blocks found across all genomes
- genomes_with_blocks: Number of genomes with at least one block
- avg_blocks_per_locus: Average blocks per locus
- avg_proteins_per_block: Average proteins (hits) per block

Usage:
    python summarize_phase2.py
"""

import pandas as pd
from pathlib import Path

def summarize_phase2():
    """Generate Phase 2 summary for all families."""

    results = []

    # Find all family outputs
    outputs_dir = Path("outputs")
    family_dirs = sorted(outputs_dir.glob("*_BK_LB_batch"))

    for family_dir in family_dirs:
        family_name = family_dir.name.replace("_BK_LB_batch", "")

        # Check if Phase 2 completed
        phase2_dir = family_dir / "02_synteny_blocks"
        blocks_file = phase2_dir / "all_synteny_blocks.tsv"

        if not blocks_file.exists():
            results.append({
                'family': family_name,
                'status': 'INCOMPLETE',
                'input_loci': 0,
                'total_blocks': 0,
                'genomes_with_blocks': 0,
                'avg_blocks_per_locus': 0.0,
                'avg_proteins_per_block': 0.0
            })
            continue

        # Read synteny blocks
        blocks_df = pd.read_csv(blocks_file, sep='\t')

        # Count input loci (from Phase 1)
        locus_defs = family_dir / "phase12_landmark" / "locus_definitions.tsv"
        if locus_defs.exists():
            locus_df = pd.read_csv(locus_defs, sep='\t')
            input_loci = len(locus_df)
        else:
            input_loci = 0

        # Total blocks
        total_blocks = len(blocks_df)

        # Genomes with blocks
        genomes_with_blocks = blocks_df['genome'].nunique() if len(blocks_df) > 0 else 0

        # Average blocks per locus
        if input_loci > 0:
            avg_blocks_per_locus = total_blocks / input_loci
        else:
            avg_blocks_per_locus = 0.0

        # Average proteins per block (use num_query_matches column)
        if len(blocks_df) > 0 and 'num_query_matches' in blocks_df.columns:
            avg_proteins_per_block = blocks_df['num_query_matches'].mean()
        else:
            avg_proteins_per_block = 0.0

        results.append({
            'family': family_name,
            'status': 'COMPLETE',
            'input_loci': input_loci,
            'total_blocks': total_blocks,
            'genomes_with_blocks': genomes_with_blocks,
            'avg_blocks_per_locus': round(avg_blocks_per_locus, 1),
            'avg_proteins_per_block': round(avg_proteins_per_block, 1)
        })

    # Create DataFrame
    summary_df = pd.DataFrame(results)

    # Sort by family name
    summary_df = summary_df.sort_values('family')

    # Save to file
    output_file = "outputs/phase2_summary.tsv"
    summary_df.to_csv(output_file, sep='\t', index=False)

    # Print summary
    print("=" * 110)
    print("PHASE 2 SUMMARY: SYNTENY DETECTION")
    print("=" * 110)
    print()
    print(summary_df.to_string(index=False))
    print()
    print("=" * 110)
    print(f"Total families: {len(summary_df)}")
    print(f"Completed: {len(summary_df[summary_df['status'] == 'COMPLETE'])}")
    print(f"Incomplete: {len(summary_df[summary_df['status'] == 'INCOMPLETE'])}")
    print()
    if len(summary_df[summary_df['status'] == 'COMPLETE']) > 0:
        complete_df = summary_df[summary_df['status'] == 'COMPLETE']
        print(f"Total synteny blocks: {complete_df['total_blocks'].sum()}")
        print(f"Average blocks per family: {complete_df['total_blocks'].mean():.1f}")
    print()
    print(f"Summary saved to: {output_file}")
    print("=" * 110)

if __name__ == "__main__":
    summarize_phase2()
