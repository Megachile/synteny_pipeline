#!/usr/bin/env python3
"""
Step 02b: Aggregate per-locus synteny blocks into combined file.

Run after array job completes to combine all locus-specific outputs.
"""

from pathlib import Path
import pandas as pd
import config

def main():
    """Aggregate all locus-specific synteny block files."""
    print("=" * 80)
    print("STEP 02b: AGGREGATE SYNTENY BLOCKS")
    print("=" * 80)

    # Find all locus-specific files
    print("\n[1] Finding locus-specific synteny block files...")
    locus_files = list(config.STEP02_SYNTENY.glob("*_synteny_blocks.tsv"))

    # Filter out the combined file if it exists
    locus_files = [f for f in locus_files if f.name != "all_synteny_blocks.tsv"]

    if not locus_files:
        print("  ERROR: No locus-specific files found!")
        print(f"  Looking in: {config.STEP02_SYNTENY}")
        return

    print(f"  Found {len(locus_files)} locus files:")
    for f in sorted(locus_files):
        print(f"    - {f.name}")

    # Load and combine
    print("\n[2] Combining blocks...")
    all_blocks = []

    for locus_file in sorted(locus_files):
        df = pd.read_csv(locus_file, sep='\t')
        all_blocks.append(df)
        print(f"  {locus_file.stem}: {len(df)} blocks")

    combined_df = pd.concat(all_blocks, ignore_index=True)
    print(f"\n  Total: {len(combined_df)} blocks")

    # Save combined file
    print("\n[3] Saving combined file...")
    output_file = config.STEP02_SYNTENY / "all_synteny_blocks.tsv"
    combined_df.to_csv(output_file, sep='\t', index=False)
    print(f"  Saved to: {output_file.name}")

    # Summary
    print("\n[4] Summary:")
    print(f"  Total blocks: {len(combined_df)}")
    print(f"  Unique loci: {combined_df['locus_id'].nunique()}")
    print(f"  Unique genomes: {combined_df['genome'].nunique()}")

    print("\n  Blocks per locus:")
    for locus_id in sorted(combined_df['locus_id'].unique()):
        count = len(combined_df[combined_df['locus_id'] == locus_id])
        print(f"    {locus_id}: {count}")

    print("\n" + "=" * 80)
    print("AGGREGATION COMPLETE")
    print("=" * 80)

if __name__ == "__main__":
    main()
