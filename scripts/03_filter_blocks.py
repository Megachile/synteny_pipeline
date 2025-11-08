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
    - Base score = number of matched queries (preferred) or unique hit segments
    - Bonus +20 if on expected chromosome (when known)
    - Penalty -10 if on fragmented scaffold (common contig prefixes)
    """
    # Prefer num_query_matches; fall back to num_target_proteins or num_proteins
    base_count = (
        block['num_query_matches']
        if 'num_query_matches' in block
        else (block['num_target_proteins'] if 'num_target_proteins' in block else block.get('num_proteins', 0))
    )

    score = base_count

    # Chromosome bonus
    if expected_chromosome and block['scaffold'] == expected_chromosome:
        score += 20

    # Fragmented scaffold penalty (broaden beyond a single prefix)
    scaf = str(block['scaffold'])
    if scaf.startswith(('NW_', 'NODE_', 'JA')):
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
            # Multiple blocks - select best with tie-breakers:
            #  1) Highest score (num_query_matches with penalties/bonuses)
            #  2) Smaller span_kb preferred
            #  3) Higher density (matches per kb) preferred
            group = group.copy()
            group['score'] = group.apply(lambda x: calculate_block_score(x), axis=1)

            # Compute span_kb if missing
            if 'span_kb' not in group.columns:
                group['span_kb'] = (group['end'] - group['start']) / 1000.0

            # Base count for density calculation
            def _base_count(x):
                if 'num_query_matches' in x:
                    return x['num_query_matches']
                elif 'num_target_proteins' in x:
                    return x['num_target_proteins']
                else:
                    return x.get('num_proteins', 0)

            group['base_count'] = group.apply(lambda r: _base_count(r), axis=1)
            group['density'] = group.apply(lambda r: (r['base_count'] / max(r['span_kb'], 1.0)), axis=1)

            # Sort by score desc, span asc, density desc
            sorted_group = group.sort_values(by=['score', 'span_kb', 'density'], ascending=[False, True, False])
            best_idx = sorted_group.index[0]
            best_block = sorted_group.loc[best_idx].to_dict()

            # Remove temporary columns
            for tmp_col in ('score', 'base_count', 'density'):
                if tmp_col in best_block:
                    del best_block[tmp_col]

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

    # Check if we have any blocks
    if len(blocks_df) == 0:
        print("\n[2] No synteny blocks found - creating empty output")
        filtered_df = blocks_df.copy()
    else:
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

    # Sync locus_definitions.tsv with actual filtered loci
    print("\n" + "=" * 80)
    print("SYNCING LOCUS DEFINITIONS")
    print("=" * 80)
    sync_locus_definitions(args.output_dir, filtered_df['locus_id'].unique())

    print("\nNext step: 04_blast_targets.py")


def sync_locus_definitions(output_dir, actual_loci):
    """
    Sync locus_definitions.tsv with actual filtered loci.

    Phase 1/12 creates initial locus_definitions.tsv, but Phase 2-3 may change
    locus IDs due to tandem clustering, re-ranking, or filtering. This ensures
    locus_definitions.tsv matches reality before downstream phases use it.
    """
    import shutil

    # Find locus_definitions files
    locus_defs_file = output_dir / "locus_definitions.tsv"
    phase12_locus_defs = output_dir / "phase12_landmark" / "locus_definitions.tsv"
    phase12_dir = output_dir / "phase12_landmark"
    synteny_dir = output_dir / "02_synteny_blocks"

    if not phase12_locus_defs.exists():
        print("  ⚠ No phase12 locus_definitions.tsv found - skipping sync")
        return

    # Load phase12 definitions
    try:
        df = pd.read_csv(phase12_locus_defs, sep='\t')
        locus_ids_defs = set(df['locus_id'].values)
    except Exception as e:
        print(f"  ⚠ Could not read locus_definitions.tsv: {e}")
        return

    locus_ids_actual = set(actual_loci)

    # Check if sync needed
    if locus_ids_defs == locus_ids_actual:
        print(f"  ✓ Already synced ({len(locus_ids_actual)} loci)")
        # Copy to root if needed
        if not locus_defs_file.exists():
            df.to_csv(locus_defs_file, sep='\t', index=False)
            print(f"  ✓ Copied to root directory")
        return

    print(f"  Syncing {len(locus_ids_defs)} → {len(locus_ids_actual)} loci...")

    # Map old → new locus IDs by matching chromosome base
    old_to_new = {}
    unmatched_actual = set(locus_ids_actual)

    for old_id in locus_ids_defs:
        parts = old_id.split('_')
        if len(parts) < 3:
            continue
        old_base = f"{parts[0]}_{parts[1]}"  # e.g., "BK_chr1"

        for new_id in locus_ids_actual:
            new_parts = new_id.split('_')
            if len(new_parts) < 3:
                continue
            new_base = f"{new_parts[0]}_{new_parts[1]}"

            if old_base == new_base and new_id in unmatched_actual:
                old_to_new[old_id] = new_id
                unmatched_actual.remove(new_id)
                break

    # Update locus IDs in dataframe
    for old_id, new_id in old_to_new.items():
        if old_id != new_id:
            df.loc[df['locus_id'] == old_id, 'locus_id'] = new_id
            # Update flanking_file path
            old_path = str(df.loc[df['locus_id'] == new_id, 'flanking_file'].iloc[0])
            new_path = old_path.replace(f"{old_id}_flanking.faa", f"{new_id}_flanking.faa")
            df.loc[df['locus_id'] == new_id, 'flanking_file'] = new_path
            print(f"    {old_id} → {new_id}")

            # Copy flanking file with new name
            if phase12_dir.exists():
                old_file = phase12_dir / f"{old_id}_flanking.faa"
                new_file = phase12_dir / f"{new_id}_flanking.faa"
                if old_file.exists() and not new_file.exists():
                    shutil.copy(old_file, new_file)

    # Remove loci that no longer exist
    df = df[df['locus_id'].isin(locus_ids_actual)]

    # Save synced definitions to root directory
    df.to_csv(locus_defs_file, sep='\t', index=False)
    print(f"  ✓ Synced {len(df)} loci to {locus_defs_file}")


if __name__ == "__main__":
    main()
