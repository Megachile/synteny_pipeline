#!/usr/bin/env python3
"""
Phase 8 (Helixer): Generate summary matrices for Helixer pipeline.

Creates summary matrices showing target counts and classifications per genome
for comparison with tblastn pipeline results.

Inputs:
- Phase 4 Helixer targets (all_target_loci.tsv)
- Phase 5 Helixer classifications (syntenic_targets.tsv, unplaceable_targets.tsv)
- Phase 6-7 Helixer annotations (annotated_targets.tsv) - optional

Outputs:
- {family}_helixer_summary.tsv: Per-genome target counts and classifications
- {family}_helixer_vs_tblastn.tsv: Comparison with tblastn pipeline
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd


def load_helixer_targets(
    phase4_dir: Path,
    phase5_dir: Path,
    phase67_dir: Optional[Path] = None,
) -> pd.DataFrame:
    """Load and merge all Helixer pipeline outputs."""
    # Load Phase 4 targets
    targets_file = phase4_dir / "all_target_loci.tsv"
    if not targets_file.exists():
        return pd.DataFrame()

    targets_df = pd.read_csv(targets_file, sep='\t')

    # Load Phase 5 classifications
    syntenic_file = phase5_dir / "syntenic_targets.tsv"
    unplaceable_file = phase5_dir / "unplaceable_targets.tsv"

    classification_map = {}
    locus_map = {}

    if syntenic_file.exists():
        syntenic_df = pd.read_csv(syntenic_file, sep='\t')
        for _, row in syntenic_df.iterrows():
            classification_map[row['target_gene_id']] = 'syntenic'
            if 'locus_id' in row:
                locus_map[row['target_gene_id']] = row['locus_id']

    if unplaceable_file.exists():
        unplaceable_df = pd.read_csv(unplaceable_file, sep='\t')
        for _, row in unplaceable_df.iterrows():
            if row['target_gene_id'] not in classification_map:
                classification_map[row['target_gene_id']] = 'unplaceable'

    targets_df['classification'] = targets_df['target_gene_id'].map(classification_map).fillna('unclassified')
    targets_df['assigned_locus'] = targets_df['target_gene_id'].map(locus_map).fillna('')

    # Load Phase 6-7 annotations if available
    if phase67_dir and phase67_dir.exists():
        annotations_file = phase67_dir / "annotated_targets.tsv"
        if annotations_file.exists():
            annotations_df = pd.read_csv(annotations_file, sep='\t')
            # Only keep annotation columns we need
            ann_cols = ['target_gene_id', 'swissprot_id', 'swissprot_gene', 'swissprot_description']
            ann_cols = [c for c in ann_cols if c in annotations_df.columns]
            if ann_cols:
                targets_df = targets_df.merge(
                    annotations_df[ann_cols],
                    on='target_gene_id',
                    how='left'
                )

    return targets_df


def load_tblastn_results(phase8_dir: Path, family: str) -> pd.DataFrame:
    """Load tblastn Phase 8 summary matrix for comparison."""
    summary_file = phase8_dir / f"{family}_summary_matrix.tsv"
    if not summary_file.exists():
        return pd.DataFrame()

    return pd.read_csv(summary_file, sep='\t')


def generate_helixer_summary(
    targets_df: pd.DataFrame,
    family: str,
) -> pd.DataFrame:
    """Generate summary matrix for Helixer targets."""
    if targets_df.empty:
        return pd.DataFrame()

    # Get unique loci for column headers
    loci = sorted(targets_df['parent_locus'].dropna().unique())

    # Build summary per genome
    rows = []
    for genome in sorted(targets_df['genome'].unique()):
        genome_df = targets_df[targets_df['genome'] == genome]

        row = {
            'genome': genome,
            'total_targets': len(genome_df),
            'syntenic': len(genome_df[genome_df['classification'] == 'syntenic']),
            'unplaceable': len(genome_df[genome_df['classification'] == 'unplaceable']),
        }

        # Count targets by parent locus
        for locus in loci:
            locus_df = genome_df[genome_df['parent_locus'] == locus]
            row[f'{locus}_count'] = len(locus_df)
            row[f'{locus}_syntenic'] = len(locus_df[locus_df['classification'] == 'syntenic'])

        # Get top SwissProt annotation if available
        if 'swissprot_gene' in genome_df.columns:
            top_gene = genome_df['swissprot_gene'].dropna()
            if len(top_gene) > 0:
                row['top_swissprot'] = top_gene.mode().iloc[0] if len(top_gene.mode()) > 0 else top_gene.iloc[0]
            else:
                row['top_swissprot'] = ''

        rows.append(row)

    return pd.DataFrame(rows)


def compare_with_tblastn(
    helixer_summary: pd.DataFrame,
    tblastn_summary: pd.DataFrame,
    family: str,
) -> pd.DataFrame:
    """Compare Helixer vs tblastn target counts."""
    if helixer_summary.empty:
        return pd.DataFrame()

    # Get genomes from Helixer results
    comparison_rows = []

    for _, row in helixer_summary.iterrows():
        genome = row['genome']
        helixer_total = row['total_targets']
        helixer_syntenic = row['syntenic']
        helixer_unplaceable = row['unplaceable']

        # Try to find matching tblastn results
        tblastn_total = 0
        if not tblastn_summary.empty:
            # Match by genome name (may need fuzzy matching)
            tblastn_match = tblastn_summary[
                tblastn_summary['genome_id'].str.contains(genome.split('_')[0], case=False, na=False) |
                tblastn_summary['species'].str.contains(genome.replace('_', ' '), case=False, na=False)
            ]
            if not tblastn_match.empty:
                tblastn_total = tblastn_match.iloc[0].get('total', 0)

        comparison_rows.append({
            'genome': genome,
            'helixer_total': helixer_total,
            'helixer_syntenic': helixer_syntenic,
            'helixer_unplaceable': helixer_unplaceable,
            'tblastn_total': tblastn_total,
            'difference': helixer_total - tblastn_total,
            'ratio': round(helixer_total / tblastn_total, 2) if tblastn_total > 0 else float('inf'),
        })

    return pd.DataFrame(comparison_rows)


def main():
    parser = argparse.ArgumentParser(
        description="Phase 8 (Helixer): Generate summary matrices"
    )
    parser.add_argument(
        "--family", required=True,
        help="Gene family name"
    )
    parser.add_argument(
        "--phase4-dir", type=Path, required=True,
        help="Phase 4 Helixer output directory"
    )
    parser.add_argument(
        "--phase5-dir", type=Path, required=True,
        help="Phase 5 Helixer output directory"
    )
    parser.add_argument(
        "--phase67-dir", type=Path,
        help="Phase 6-7 Helixer output directory (optional)"
    )
    parser.add_argument(
        "--tblastn-phase8", type=Path,
        help="tblastn Phase 8 directory for comparison (optional)"
    )
    parser.add_argument(
        "--output-dir", type=Path, required=True,
        help="Output directory"
    )

    args = parser.parse_args()

    print("=" * 80)
    print("PHASE 8 (HELIXER): GENERATE SUMMARY MATRICES")
    print("=" * 80)
    print()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Load Helixer targets
    targets_df = load_helixer_targets(
        args.phase4_dir,
        args.phase5_dir,
        args.phase67_dir,
    )

    if targets_df.empty:
        print(f"[WARNING] No targets found for {args.family}")
        return 0

    print(f"Loaded {len(targets_df)} targets for {args.family}")
    print(f"  Genomes: {targets_df['genome'].nunique()}")
    print(f"  Syntenic: {len(targets_df[targets_df['classification'] == 'syntenic'])}")
    print(f"  Unplaceable: {len(targets_df[targets_df['classification'] == 'unplaceable'])}")

    # Generate Helixer summary
    helixer_summary = generate_helixer_summary(targets_df, args.family)

    summary_file = args.output_dir / f"{args.family}_helixer_summary.tsv"
    helixer_summary.to_csv(summary_file, sep='\t', index=False)
    print(f"\n[OUTPUT] Helixer summary: {summary_file}")

    # Compare with tblastn if available
    if args.tblastn_phase8 and args.tblastn_phase8.exists():
        tblastn_summary = load_tblastn_results(args.tblastn_phase8, args.family)
        if not tblastn_summary.empty:
            comparison = compare_with_tblastn(helixer_summary, tblastn_summary, args.family)

            comparison_file = args.output_dir / f"{args.family}_helixer_vs_tblastn.tsv"
            comparison.to_csv(comparison_file, sep='\t', index=False)
            print(f"[OUTPUT] Comparison: {comparison_file}")

            # Print comparison summary
            print("\nHelixer vs tblastn comparison:")
            for _, row in comparison.iterrows():
                ratio_str = f"{row['ratio']:.1f}x" if row['ratio'] != float('inf') else "inf"
                print(f"  {row['genome']}: {row['helixer_total']} vs {row['tblastn_total']} ({ratio_str})")

    # Save detailed targets with all annotations
    detailed_file = args.output_dir / f"{args.family}_all_targets_detailed.tsv"
    targets_df.to_csv(detailed_file, sep='\t', index=False)
    print(f"[OUTPUT] Detailed targets: {detailed_file}")

    print("\nPhase 8 (Helixer) complete.")
    return 0


if __name__ == "__main__":
    exit(main())
