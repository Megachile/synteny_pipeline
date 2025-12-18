#!/usr/bin/env python3
"""
Collect BUSCO scores from Helixer proteomes and apply exclusion logic.

This script:
1. Parses all busco_v6_hym short_summary files
2. Applies exclusion rules (hard floor + redundancy check)
3. Updates gca_to_species.tsv with BUSCO columns
4. Optionally moves excluded genomes to ragtag_output_excluded/

Usage:
    python collect_busco_and_exclude.py --dry-run    # Preview only
    python collect_busco_and_exclude.py --execute    # Actually move excluded genomes
"""

import argparse
import glob
import os
import re
import shutil
import sys
from pathlib import Path

import pandas as pd

# =============================================================================
# EXCLUSION PARAMETERS - adjust these as needed
# =============================================================================
HARD_FLOOR = 40.0           # Exclude genomes below this (always)
REDUNDANCY_GAP = 14.0       # Exclude if this many pts below group median
MIN_GOOD_IN_GROUP = 2       # Need this many good genomes to trigger redundancy exclusion
GOOD_THRESHOLD = 75.0       # What counts as "good" alternative

# =============================================================================
# PATHS
# =============================================================================
BASE_DIR = Path("/carc/scratch/projects/emartins/2016456/adam/synteny_scanning/hpc_deployment_package")
RAGTAG_DIR = BASE_DIR / "hybrid_workflow" / "data" / "ragtag_output"
EXCLUDED_DIR = BASE_DIR / "hybrid_workflow" / "data" / "ragtag_output_excluded"
GCA_FILE = BASE_DIR / "data" / "gca_to_species.tsv"


def parse_busco_summary(summary_file: Path) -> dict:
    """Parse a BUSCO short_summary.txt file and return metrics."""
    with open(summary_file) as f:
        for line in f:
            if line.strip().startswith("C:"):
                # Parse C:74.0%[S:73.0%,D:1.0%],F:5.7%,M:20.4%
                match = re.search(
                    r'C:(\d+\.\d+)%.*S:(\d+\.\d+)%.*D:(\d+\.\d+)%.*F:(\d+\.\d+)%.*M:(\d+\.\d+)%',
                    line
                )
                if match:
                    complete, single, dup, frag, missing = map(float, match.groups())
                    return {
                        'busco_complete': complete,
                        'busco_single': single,
                        'busco_dup': dup,
                        'busco_frag': frag,
                        'busco_missing': missing
                    }
    return None


def collect_all_busco_scores(ragtag_dir: Path) -> pd.DataFrame:
    """Collect BUSCO scores from all genomes."""
    results = []

    for summary_file in ragtag_dir.glob("*/busco_v6_hym/short_summary*.txt"):
        folder_name = summary_file.parent.parent.name
        metrics = parse_busco_summary(summary_file)

        if metrics:
            metrics['folder_name'] = folder_name
            results.append(metrics)

    if not results:
        print("WARNING: No BUSCO summaries found!")
        return pd.DataFrame()

    return pd.DataFrame(results)


def load_species_mapping(gca_file: Path) -> pd.DataFrame:
    """Load gca_to_species.tsv and create lookup."""
    df = pd.read_csv(gca_file, sep='\t')
    # Normalize column names
    df.columns = ['accession', 'species', 'family', 'phylo_order', 'status']
    return df


def match_folder_to_species(folder_name: str, species_df: pd.DataFrame) -> dict:
    """Match a folder name to species info."""
    # Try direct accession match
    match = species_df[species_df['accession'] == folder_name]
    if len(match) == 1:
        return match.iloc[0].to_dict()

    # Try species name match (folder uses underscores)
    species_name = folder_name.replace('_', ' ')
    match = species_df[species_df['species'] == species_name]
    if len(match) == 1:
        return match.iloc[0].to_dict()

    # Try partial match on species column
    match = species_df[species_df['species'].str.replace(' ', '_') == folder_name]
    if len(match) == 1:
        return match.iloc[0].to_dict()

    return None


def compute_exclusion_decision(row: pd.Series, all_data: pd.DataFrame) -> tuple:
    """
    Determine if genome should be excluded and why.

    Returns: (decision, reason) where decision is 'keep' or 'exclude'
    """
    busco = row['busco_complete']
    phylo = row['phylo_order']

    # Rule 1: Hard floor
    if busco < HARD_FLOOR:
        return 'exclude', f'below {HARD_FLOOR}% floor'

    # Get group stats (only for genomes with BUSCO scores)
    group = all_data[all_data['phylo_order'] == phylo]
    n_in_group = len(group)

    # Rule 2: Only genome in group - keep regardless
    if n_in_group == 1:
        return 'keep', 'sole representative'

    # Rule 3: Redundancy check
    group_median = group['busco_complete'].median()
    n_good = len(group[group['busco_complete'] >= GOOD_THRESHOLD])
    gap_from_median = group_median - busco

    if n_good >= MIN_GOOD_IN_GROUP and gap_from_median > REDUNDANCY_GAP:
        return 'exclude', f'redundant ({busco:.1f}% vs {group_median:.1f}% median, {n_good} better alternatives)'

    return 'keep', 'passes filters'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dry-run', action='store_true',
                        help='Preview exclusions without moving files')
    parser.add_argument('--execute', action='store_true',
                        help='Actually move excluded genomes')
    parser.add_argument('--output-tsv', type=str, default=None,
                        help='Write results to this TSV (default: update gca_to_species.tsv)')
    args = parser.parse_args()

    if not args.dry_run and not args.execute:
        print("ERROR: Must specify --dry-run or --execute")
        sys.exit(1)

    # Collect BUSCO scores
    print(f"Scanning BUSCO outputs in {RAGTAG_DIR}...")
    busco_df = collect_all_busco_scores(RAGTAG_DIR)
    print(f"  Found {len(busco_df)} genomes with BUSCO scores")

    if busco_df.empty:
        print("No BUSCO data found. Exiting.")
        sys.exit(1)

    # Load species mapping
    print(f"\nLoading species mapping from {GCA_FILE}...")
    species_df = load_species_mapping(GCA_FILE)
    print(f"  Loaded {len(species_df)} species entries")

    # Match folders to species info
    print("\nMatching folders to species...")
    matched_info = []
    unmatched = []

    for _, row in busco_df.iterrows():
        folder = row['folder_name']
        info = match_folder_to_species(folder, species_df)

        if info:
            matched_info.append({
                'folder_name': folder,
                'accession': info['accession'],
                'species': info['species'],
                'family': info['family'],
                'phylo_order': info['phylo_order'],
                **{k: row[k] for k in ['busco_complete', 'busco_single', 'busco_dup', 'busco_frag', 'busco_missing']}
            })
        else:
            unmatched.append(folder)

    if unmatched:
        print(f"  WARNING: Could not match {len(unmatched)} folders: {unmatched}")

    matched_df = pd.DataFrame(matched_info)
    print(f"  Matched {len(matched_df)} genomes to species info")

    # Apply exclusion logic
    print("\nApplying exclusion logic...")
    print(f"  Parameters: floor={HARD_FLOOR}%, gap={REDUNDANCY_GAP}pts, min_good={MIN_GOOD_IN_GROUP}, good_threshold={GOOD_THRESHOLD}%")

    decisions = matched_df.apply(
        lambda r: compute_exclusion_decision(r, matched_df),
        axis=1
    )
    matched_df['busco_status'] = [d[0] for d in decisions]
    matched_df['busco_reason'] = [d[1] for d in decisions]

    # Sort by phylo order then BUSCO
    matched_df = matched_df.sort_values(['phylo_order', 'busco_complete'], ascending=[True, False])

    # Display results
    print("\n" + "=" * 100)
    print(f"{'Phylo':<5} {'Species':<35} {'BUSCO':>7} {'Status':<8} Reason")
    print("=" * 100)

    current_phylo = None
    for _, row in matched_df.iterrows():
        if current_phylo != row['phylo_order']:
            if current_phylo is not None:
                print("-" * 100)
            current_phylo = row['phylo_order']

        marker = "X" if row['busco_status'] == 'exclude' else " "
        print(f"{int(row['phylo_order']):<5} {row['species']:<35} {row['busco_complete']:>6.1f}% [{marker}] {row['busco_reason']}")

    print("=" * 100)

    # Summary
    excluded = matched_df[matched_df['busco_status'] == 'exclude']
    kept = matched_df[matched_df['busco_status'] == 'keep']

    print(f"\nSUMMARY: {len(kept)} keep, {len(excluded)} exclude")

    if len(excluded) > 0:
        print("\nGenomes to exclude:")
        for _, row in excluded.iterrows():
            print(f"  - {row['folder_name']}: {row['busco_reason']}")

    # Handle execution
    if args.execute:
        print("\n" + "=" * 50)
        print("EXECUTING EXCLUSIONS")
        print("=" * 50)

        # Create excluded directory
        EXCLUDED_DIR.mkdir(parents=True, exist_ok=True)

        # Move excluded genomes
        for _, row in excluded.iterrows():
            src = RAGTAG_DIR / row['folder_name']
            dst = EXCLUDED_DIR / row['folder_name']

            if src.exists():
                print(f"  Moving {row['folder_name']} -> ragtag_output_excluded/")
                shutil.move(str(src), str(dst))
            else:
                print(f"  WARNING: {src} not found, skipping")

        # Update gca_to_species.tsv
        output_file = args.output_tsv if args.output_tsv else GCA_FILE
        print(f"\nUpdating {output_file}...")

        # Merge BUSCO data back into species file
        species_df_updated = species_df.copy()
        busco_cols = ['busco_complete', 'busco_single', 'busco_dup', 'busco_frag', 'busco_missing', 'busco_status', 'busco_reason']

        for col in busco_cols:
            species_df_updated[col] = None

        for _, row in matched_df.iterrows():
            mask = species_df_updated['accession'] == row['accession']
            if not mask.any():
                # Try folder name match
                mask = species_df_updated['accession'] == row['folder_name']

            if mask.any():
                for col in busco_cols:
                    species_df_updated.loc[mask, col] = row[col]

        species_df_updated.to_csv(output_file, sep='\t', index=False)
        print(f"  Written updated species table with {len(busco_cols)} new columns")

    else:
        print("\n[DRY RUN] No files modified. Use --execute to apply changes.")


if __name__ == '__main__':
    main()
