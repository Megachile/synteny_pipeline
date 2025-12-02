#!/usr/bin/env python3
"""
Phase 5b (Helixer): Refine unplaceable target classification.

Takes Phase 5 outputs and further classifies unplaceable targets:
- tandem: Within proximity of a syntenic target on same scaffold
- weak_match: Has partial synteny signal but below threshold
- orphan: No synteny signal at all

Inputs:
- Phase 5 syntenic_targets.tsv
- Phase 5 unplaceable_targets.tsv

Outputs:
- unplaceable_refined.tsv: Unplaceables with refined classification
- tandem_candidates.tsv: Targets classified as potential tandems
- classification_summary.tsv: Summary of classification counts
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd


def find_nearby_syntenic_targets(
    unplaceable_row: pd.Series,
    syntenic_df: pd.DataFrame,
    max_distance_kb: float = 50.0,
) -> List[Dict]:
    """
    Find syntenic targets on the same scaffold within max_distance_kb.

    Returns list of nearby syntenic targets with distance info.
    """
    # Handle empty syntenic DataFrame
    if syntenic_df.empty:
        return []

    genome = unplaceable_row['genome']
    scaffold = unplaceable_row['scaffold']
    target_start = unplaceable_row['start']
    target_end = unplaceable_row['end']
    target_mid = (target_start + target_end) / 2

    max_distance_bp = max_distance_kb * 1000

    nearby = []

    # Filter to same genome and scaffold
    same_location = syntenic_df[
        (syntenic_df['genome'] == genome) &
        (syntenic_df['scaffold'] == scaffold)
    ]

    for _, syn_row in same_location.iterrows():
        syn_start = syn_row['start']
        syn_end = syn_row['end']
        syn_mid = (syn_start + syn_end) / 2

        # Calculate distance between midpoints
        distance_bp = abs(target_mid - syn_mid)

        if distance_bp <= max_distance_bp:
            nearby.append({
                'syntenic_target_id': syn_row['target_gene_id'],
                'syntenic_locus': syn_row['locus_id'],
                'distance_kb': round(distance_bp / 1000, 2),
                'syntenic_start': syn_start,
                'syntenic_end': syn_end,
            })

    # Sort by distance
    nearby.sort(key=lambda x: x['distance_kb'])

    return nearby


def classify_unplaceable(
    row: pd.Series,
    nearby_syntenic: List[Dict],
    tandem_distance_kb: float = 50.0,
    weak_match_threshold: float = 0.05,
) -> Tuple[str, str, str]:
    """
    Classify an unplaceable target into refined categories.

    Returns: (refined_class, tandem_locus, notes)
    """
    unplaceable_reason = row.get('unplaceable_reason', '')
    best_locus_score = row.get('best_locus_score', 0)
    best_locus_match = row.get('best_locus_match', '')

    # Check for tandem (near a syntenic target)
    if nearby_syntenic:
        nearest = nearby_syntenic[0]
        if nearest['distance_kb'] <= tandem_distance_kb:
            return (
                'tandem',
                nearest['syntenic_locus'],
                f"within {nearest['distance_kb']}kb of {nearest['syntenic_target_id']}"
            )

    # Check for weak match (has partial synteny but below threshold)
    if best_locus_score > weak_match_threshold:
        return (
            'weak_match',
            best_locus_match,
            f"score {best_locus_score:.3f} for {best_locus_match}"
        )

    # Check reason-based classification
    if 'below_threshold' in unplaceable_reason:
        return (
            'weak_match',
            best_locus_match,
            f"below threshold: {unplaceable_reason}"
        )

    # Otherwise orphan
    if unplaceable_reason == 'no_flanking_genes':
        return ('orphan', '', 'no flanking genes extracted')
    elif unplaceable_reason == 'no_bk_hits':
        return ('orphan', '', 'flanking genes have no BK orthologs')
    elif unplaceable_reason == 'no_locus_matches':
        return ('orphan', '', 'BK orthologs not in any reference locus')
    else:
        return ('orphan', '', unplaceable_reason or 'unknown')


def main():
    parser = argparse.ArgumentParser(
        description="Phase 5b (Helixer): Refine unplaceable classification"
    )
    parser.add_argument("--family", required=True, help="Gene family name")
    parser.add_argument(
        "--phase5-dir", type=Path, required=True,
        help="Phase 5 output directory"
    )
    parser.add_argument(
        "--output-dir", type=Path, required=True,
        help="Output directory for refined classifications"
    )
    parser.add_argument(
        "--tandem-distance-kb", type=float, default=50.0,
        help="Max distance (kb) from syntenic target to classify as tandem (default: 50)"
    )
    parser.add_argument(
        "--weak-match-threshold", type=float, default=0.05,
        help="Minimum synteny score to classify as weak_match (default: 0.05)"
    )

    args = parser.parse_args()

    print("=" * 80)
    print("PHASE 5B (HELIXER): REFINE UNPLACEABLE CLASSIFICATION")
    print("=" * 80)
    print()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Load Phase 5 outputs
    syntenic_file = args.phase5_dir / "syntenic_targets.tsv"
    unplaceable_file = args.phase5_dir / "unplaceable_targets.tsv"

    if not syntenic_file.exists():
        print(f"[WARNING] No syntenic targets file: {syntenic_file}")
        syntenic_df = pd.DataFrame()
    else:
        try:
            syntenic_df = pd.read_csv(syntenic_file, sep='\t')
        except pd.errors.EmptyDataError:
            print(f"[WARNING] Syntenic targets file is empty: {syntenic_file}")
            syntenic_df = pd.DataFrame()

    if not unplaceable_file.exists():
        print(f"[ERROR] No unplaceable targets file: {unplaceable_file}")
        return 1

    unplaceable_df = pd.read_csv(unplaceable_file, sep='\t')

    print(f"Loaded {len(syntenic_df)} syntenic targets")
    print(f"Loaded {len(unplaceable_df)} unplaceable targets")

    if unplaceable_df.empty:
        print("\nNo unplaceables to refine. Exiting.")
        return 0

    # Classify each unplaceable
    refined_rows = []
    tandem_rows = []

    for _, row in unplaceable_df.iterrows():
        # Find nearby syntenic targets
        nearby = find_nearby_syntenic_targets(
            row, syntenic_df, args.tandem_distance_kb
        )

        # Classify
        refined_class, tandem_locus, notes = classify_unplaceable(
            row, nearby, args.tandem_distance_kb, args.weak_match_threshold
        )

        # Build refined row
        refined_row = row.to_dict()
        refined_row['refined_class'] = refined_class
        refined_row['tandem_of_locus'] = tandem_locus
        refined_row['classification_notes'] = notes

        # Add nearby syntenic info
        if nearby:
            refined_row['nearest_syntenic_target'] = nearby[0]['syntenic_target_id']
            refined_row['nearest_syntenic_locus'] = nearby[0]['syntenic_locus']
            refined_row['nearest_distance_kb'] = nearby[0]['distance_kb']
        else:
            refined_row['nearest_syntenic_target'] = ''
            refined_row['nearest_syntenic_locus'] = ''
            refined_row['nearest_distance_kb'] = -1

        refined_rows.append(refined_row)

        # Track tandems separately
        if refined_class == 'tandem':
            tandem_rows.append({
                'target_gene_id': row['target_gene_id'],
                'genome': row['genome'],
                'scaffold': row['scaffold'],
                'start': row['start'],
                'end': row['end'],
                'tandem_of_locus': tandem_locus,
                'nearest_syntenic_target': nearby[0]['syntenic_target_id'],
                'distance_kb': nearby[0]['distance_kb'],
                'gene_family': args.family,
            })

    # Write outputs
    refined_df = pd.DataFrame(refined_rows)
    refined_df.to_csv(
        args.output_dir / "unplaceable_refined.tsv", sep='\t', index=False
    )

    if tandem_rows:
        pd.DataFrame(tandem_rows).to_csv(
            args.output_dir / "tandem_candidates.tsv", sep='\t', index=False
        )

    # Summary
    class_counts = refined_df['refined_class'].value_counts()
    summary_rows = [
        {'classification': 'syntenic', 'count': len(syntenic_df)},
    ]
    for cls, count in class_counts.items():
        summary_rows.append({'classification': cls, 'count': count})
    summary_rows.append({'classification': 'TOTAL', 'count': len(syntenic_df) + len(unplaceable_df)})

    pd.DataFrame(summary_rows).to_csv(
        args.output_dir / "classification_summary.tsv", sep='\t', index=False
    )

    # Print summary
    print(f"\n[RESULTS]")
    print(f"  Syntenic: {len(syntenic_df)}")
    for cls in ['tandem', 'weak_match', 'orphan']:
        count = class_counts.get(cls, 0)
        print(f"  {cls}: {count}")
    print(f"  TOTAL: {len(syntenic_df) + len(unplaceable_df)}")

    # Show tandems
    if tandem_rows:
        print(f"\nTandem candidates ({len(tandem_rows)}):")
        for t in tandem_rows[:10]:
            print(f"  {t['target_gene_id']}: tandem of {t['tandem_of_locus']} "
                  f"({t['distance_kb']}kb from {t['nearest_syntenic_target']})")
        if len(tandem_rows) > 10:
            print(f"  ... and {len(tandem_rows) - 10} more")

    # Show orphan breakdown
    orphans = refined_df[refined_df['refined_class'] == 'orphan']
    if not orphans.empty:
        print(f"\nOrphan breakdown:")
        for reason, count in orphans['classification_notes'].value_counts().items():
            print(f"  {reason}: {count}")

    print(f"\n[OUTPUT] {args.output_dir}/unplaceable_refined.tsv")
    if tandem_rows:
        print(f"[OUTPUT] {args.output_dir}/tandem_candidates.tsv")
    print(f"[OUTPUT] {args.output_dir}/classification_summary.tsv")

    return 0


if __name__ == "__main__":
    exit(main())
