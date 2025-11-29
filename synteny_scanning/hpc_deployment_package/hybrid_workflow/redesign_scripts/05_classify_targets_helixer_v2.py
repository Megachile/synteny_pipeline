#!/usr/bin/env python3
"""
Phase 5 (Helixer v2): Classify targets using existing tblastn synteny blocks.

Instead of re-computing synteny from Helixer flanking genes, use the
synteny blocks already computed by the tblastn pipeline (Phase 5_v2).
Check if Helixer target coordinates fall within known synteny blocks.

This is simpler and more consistent with the existing workflow.

Outputs:
- syntenic_targets.tsv      # Targets inside synteny blocks
- unplaceable_targets.tsv   # Targets outside synteny blocks
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd


def load_synteny_blocks(phase5_dir: Path) -> pd.DataFrame:
    """Load synteny blocks from tblastn pipeline Phase 5."""
    blocks_file = phase5_dir / "selected_synteny_blocks.tsv"
    if not blocks_file.exists():
        return pd.DataFrame()

    df = pd.read_csv(blocks_file, sep='\t')
    return df


def check_overlap(
    target_scaffold: str,
    target_start: int,
    target_end: int,
    block_scaffold: str,
    block_start: int,
    block_end: int,
) -> bool:
    """Check if target overlaps with synteny block."""
    if target_scaffold != block_scaffold:
        return False

    # Check coordinate overlap
    return not (target_end < block_start or target_start > block_end)


def classify_target(
    target: pd.Series,
    blocks_df: pd.DataFrame,
) -> Tuple[str, str, float]:
    """
    Classify a target as syntenic or unplaceable.

    Returns (classification, matched_locus, synteny_score)
    """
    if blocks_df.empty:
        return 'unplaceable', '', 0.0

    genome = target['genome']
    scaffold = target['scaffold']
    start = target['start']
    end = target['end']

    # Check each synteny block for this genome
    genome_blocks = blocks_df[blocks_df['genome'] == genome]

    for _, block in genome_blocks.iterrows():
        if check_overlap(scaffold, start, end,
                        block['scaffold'], block['start'], block['end']):
            # Get synteny score from block
            synteny_score = block.get('synteny_score', block.get('flanking_matches', 0))
            if isinstance(synteny_score, str):
                try:
                    synteny_score = float(synteny_score.rstrip('%')) / 100
                except:
                    synteny_score = 0.5
            elif synteny_score > 1:
                synteny_score = synteny_score / 100  # Convert percentage

            return 'syntenic', block['locus_id'], synteny_score

    return 'unplaceable', '', 0.0


def main():
    parser = argparse.ArgumentParser(
        description="Phase 5 (Helixer v2): Classify targets using tblastn synteny blocks"
    )
    parser.add_argument(
        "--family", required=True,
        help="Gene family name"
    )
    parser.add_argument(
        "--phase4-targets", type=Path, required=True,
        help="Phase 4 Helixer all_target_loci.tsv"
    )
    parser.add_argument(
        "--phase5-tblastn", type=Path, required=True,
        help="Phase 5 tblastn directory (phase5_v2)"
    )
    parser.add_argument(
        "--output-dir", type=Path, required=True,
        help="Output directory"
    )

    args = parser.parse_args()

    print("=" * 80)
    print("PHASE 5 (HELIXER V2): CLASSIFY TARGETS BY SYNTENY BLOCK OVERLAP")
    print("=" * 80)
    print()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Load Phase 4 Helixer targets
    targets_df = pd.read_csv(args.phase4_targets, sep='\t')
    print(f"Loaded {len(targets_df)} targets from Phase 4 Helixer")

    # Load tblastn synteny blocks
    blocks_df = load_synteny_blocks(args.phase5_tblastn)
    print(f"Loaded {len(blocks_df)} synteny blocks from tblastn Phase 5")

    if blocks_df.empty:
        print("[WARNING] No synteny blocks found - all targets will be unplaceable")

    # Classify each target
    results = []
    for _, target in targets_df.iterrows():
        classification, matched_locus, synteny_score = classify_target(target, blocks_df)

        result = target.to_dict()
        result['classification'] = classification
        result['matched_block_locus'] = matched_locus
        result['synteny_score'] = synteny_score
        results.append(result)

    results_df = pd.DataFrame(results)

    # Split into syntenic and unplaceable
    syntenic = results_df[results_df['classification'] == 'syntenic']
    unplaceable = results_df[results_df['classification'] == 'unplaceable']

    # Add locus_id column for compatibility with downstream scripts
    if 'matched_block_locus' in syntenic.columns:
        syntenic = syntenic.copy()
        syntenic['locus_id'] = syntenic['matched_block_locus']

    # Write outputs
    syntenic_file = args.output_dir / "syntenic_targets.tsv"
    syntenic.to_csv(syntenic_file, sep='\t', index=False)
    print(f"\n[OUTPUT] {len(syntenic)} syntenic targets → {syntenic_file}")

    unplaceable_file = args.output_dir / "unplaceable_targets.tsv"
    unplaceable.to_csv(unplaceable_file, sep='\t', index=False)
    print(f"[OUTPUT] {len(unplaceable)} unplaceable targets → {unplaceable_file}")

    # Summary by genome
    print("\nSummary by genome:")
    for genome in targets_df['genome'].unique():
        genome_targets = results_df[results_df['genome'] == genome]
        n_syntenic = len(genome_targets[genome_targets['classification'] == 'syntenic'])
        n_unplaceable = len(genome_targets[genome_targets['classification'] == 'unplaceable'])
        print(f"  {genome}: {n_syntenic} syntenic, {n_unplaceable} unplaceable")

    print("\nPhase 5 (Helixer v2) complete.")
    return 0


if __name__ == "__main__":
    exit(main())
