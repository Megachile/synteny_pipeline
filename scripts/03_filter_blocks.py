#!/usr/bin/env python3
"""
Step 03: Filter synteny blocks to select best block per genome per locus.

Uses chromosome-aware scoring to select the best synteny block when multiple
blocks exist for the same genome/locus combination.

Usage:
    python 03_filter_blocks.py \\
        --input <path/to/all_synteny_blocks.tsv> \\
        --output-dir <path/to/03_filtered_blocks> \\
        [--verbose]
"""

from pathlib import Path
import pandas as pd
import argparse

def calculate_block_score(block, expected_chromosome=None):
    """
    Calculate score for a synteny block.

    Scoring:
    - Base score = number of proteins found
    - Bonus +20 if on expected chromosome (when known)
    - Penalty -10 if on fragmented scaffold (e.g., NW_*, JAACYO*)
    """
    score = block['num_proteins']

    # Chromosome bonus
    if expected_chromosome and block['scaffold'] == expected_chromosome:
        score += 20

    # Fragmented scaffold penalty
    if block['scaffold'].startswith(('NW_', 'JAACYO')):
        score -= 10

    return score

def filter_best_blocks(blocks_df, verbose=False):
    """Select best block per genome per locus."""

    filtered_blocks = []

    # Group by locus and genome
    grouped = blocks_df.groupby(['locus_id', 'genome'])

    for (locus_id, genome), group in grouped:
        if len(group) == 1:
            # Only one block for this genome/locus
            filtered_blocks.append(group.iloc[0].to_dict())
        else:
            # Multiple blocks - select best
            # Calculate scores
            group['score'] = group.apply(lambda x: calculate_block_score(x), axis=1)

            # Select block with highest score
            best_idx = group['score'].idxmax()
            best_block = group.loc[best_idx].to_dict()

            # Remove temporary score column from output
            del best_block['score']

            filtered_blocks.append(best_block)

            if verbose:
                print(f"    {genome}/{locus_id}: {len(group)} blocks -> selected block on {best_block['scaffold']}")

    return pd.DataFrame(filtered_blocks)

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Filter synteny blocks to select best per genome per locus",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('--input', required=True, type=Path,
                        help='Path to all_synteny_blocks.tsv from Phase 2')
    parser.add_argument('--output-dir', required=True, type=Path,
                        help='Output directory for filtered blocks')
    parser.add_argument('--verbose', action='store_true',
                        help='Print verbose output for block selection')

    return parser.parse_args()

def main():
    """Main execution function."""
    args = parse_args()

    print("=" * 80)
    print("STEP 03: FILTER SYNTENY BLOCKS (BEST PER GENOME)")
    print("=" * 80)
    print(f"\nInput: {args.input}")
    print(f"Output: {args.output_dir}")

    # Create output directory
    args.output_dir.mkdir(exist_ok=True, parents=True)

    # Load synteny blocks
    print("\n[1] Loading synteny blocks...")

    if not args.input.exists():
        print(f"  ERROR: Synteny blocks file not found: {args.input}")
        print("  Run Step 02 first.")
        return

    blocks_df = pd.read_csv(args.input, sep='\t')
    print(f"  Loaded {len(blocks_df)} synteny blocks")

    # Get unique loci
    loci = blocks_df['locus_id'].unique()
    print(f"  Loci: {', '.join(loci)}")

    # Filter blocks
    print("\n[2] Filtering to best block per genome...")
    filtered_df = filter_best_blocks(blocks_df, args.verbose)

    reduction = (1 - len(filtered_df) / len(blocks_df)) * 100
    print(f"  Filtered: {len(blocks_df)} -> {len(filtered_df)} blocks ({reduction:.1f}% reduction)")

    # Save filtered blocks
    print("\n[3] Saving filtered blocks...")

    # Save combined file
    all_filtered_file = args.output_dir / "synteny_blocks_filtered.tsv"
    filtered_df.to_csv(all_filtered_file, sep='\t', index=False)
    print(f"  Saved combined: {all_filtered_file.name}")

    # Also save per-locus files for easier access
    for locus_id in loci:
        locus_blocks = filtered_df[filtered_df['locus_id'] == locus_id]

        if len(locus_blocks) > 0:
            locus_file = args.output_dir / f"{locus_id}_synteny_blocks_filtered.tsv"
            locus_blocks.to_csv(locus_file, sep='\t', index=False)
            print(f"  Saved {locus_id}: {len(locus_blocks)} blocks")

    # Summary statistics
    print("\n[4] Summary statistics:")

    for locus_id in loci:
        locus_blocks = filtered_df[filtered_df['locus_id'] == locus_id]
        print(f"\n  {locus_id}:")
        print(f"    Genomes with blocks: {locus_blocks['genome'].nunique()}")

        if len(locus_blocks) > 0:
            # Check column name (might be num_proteins or num_query_matches)
            num_col = 'num_proteins' if 'num_proteins' in locus_blocks.columns else 'num_query_matches'
            print(f"    Proteins found: {locus_blocks[num_col].min()}-{locus_blocks[num_col].max()} (median: {locus_blocks[num_col].median():.0f})")

            # Chromosome distribution
            chrom_dist = locus_blocks['scaffold'].value_counts().head(5)
            print(f"    Top scaffolds:")
            for scaffold, count in chrom_dist.items():
                pct = count / len(locus_blocks) * 100
                print(f"      {scaffold}: {count} ({pct:.1f}%)")

    # Overall summary
    print("\n" + "=" * 80)
    print("FILTERING COMPLETE")
    print("=" * 80)
    print(f"\nFiltered blocks saved to: {args.output_dir}")
    print(f"  Total blocks: {len(filtered_df)}")
    print(f"  Unique genomes: {filtered_df['genome'].nunique()}")
    print(f"  Unique loci: {filtered_df['locus_id'].nunique()}")
    print("\nNext step: 04_blast_targets.py")

if __name__ == "__main__":
    main()