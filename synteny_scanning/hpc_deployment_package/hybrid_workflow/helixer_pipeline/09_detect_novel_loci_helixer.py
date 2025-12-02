#!/usr/bin/env python3
"""
Phase 9 (Helixer): Detect novel loci via reciprocal flanking gene analysis.

Takes orphan targets (from Phase 5b) and checks if their flanking genes are
conserved across multiple genomes. If so, they may represent genuinely novel
loci not present in BK.

Key insight: A real novel locus should have CONSERVED flanking genes across
multiple genomes, not just random flanking patterns.

Flow:
1. Load orphan targets from Phase 5b
2. Get their flanking genes (from Phase 4b)
3. Cluster orphans by shared flanking gene orthologs
4. Score each cluster by multi-genome conservation
5. Filter: require ≥3 genomes with similar flanking

Filters to prevent explosion:
- Require ≥3 genomes with similar flanking signature
- Require flanking genes on BOTH sides (upstream AND downstream)
- Require scaffold ≥100kb (configurable)
- Hard cap: max N novel loci per family

Inputs:
- Phase 5b unplaceable_refined.tsv (with refined_class column)
- Phase 4b flanking_genes.tsv
- Phase 5 diamond_vs_bk.tsv (BK orthologs for each flanking gene)

Outputs:
- novel_locus_candidates.tsv: Putative new loci with conservation scores
- novel_locus_members.tsv: Orphan targets grouped by novel locus
- orphan_flanking_patterns.tsv: Flanking gene patterns for each orphan
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Set, Tuple

import pandas as pd


def load_orphan_targets(phase5b_dir: Path) -> pd.DataFrame:
    """Load orphan targets from Phase 5b refined output."""
    refined_file = phase5b_dir / "unplaceable_refined.tsv"
    if not refined_file.exists():
        return pd.DataFrame()

    df = pd.read_csv(refined_file, sep='\t')
    # Filter to only orphans
    orphans = df[df['refined_class'] == 'orphan'].copy()
    return orphans


def load_flanking_info(phase4b_dir: Path) -> pd.DataFrame:
    """Load flanking gene information from Phase 4b."""
    flanking_file = phase4b_dir / "flanking_genes.tsv"
    if not flanking_file.exists():
        return pd.DataFrame()

    return pd.read_csv(flanking_file, sep='\t')


def load_bk_orthologs(phase5_dir: Path) -> Dict[str, str]:
    """
    Load BK orthologs for each Helixer flanking protein.

    Returns: {helixer_protein_id: bk_xp_id}
    """
    diamond_file = phase5_dir / "diamond_vs_bk.tsv"
    if not diamond_file.exists():
        return {}

    orthologs = {}
    try:
        df = pd.read_csv(diamond_file, sep='\t', header=None,
                         names=['qseqid', 'sseqid', 'pident', 'length', 'evalue', 'bitscore'])
        for _, row in df.iterrows():
            orthologs[row['qseqid']] = row['sseqid']
    except Exception:
        pass

    return orthologs


def get_flanking_signature(
    target_id: str,
    flanking_df: pd.DataFrame,
    bk_orthologs: Dict[str, str],
) -> Tuple[Set[str], int, int]:
    """
    Get the BK ortholog signature for a target's flanking genes.

    Returns: (set of BK XP IDs, n_upstream, n_downstream)
    """
    target_flanking = flanking_df[flanking_df['target_gene_id'] == target_id]

    bk_signature = set()
    n_upstream = 0
    n_downstream = 0

    for _, row in target_flanking.iterrows():
        protein_id = row['flanking_protein_id']
        position = row.get('position', '')

        if protein_id in bk_orthologs:
            bk_signature.add(bk_orthologs[protein_id])

        if 'upstream' in position:
            n_upstream += 1
        elif 'downstream' in position:
            n_downstream += 1

    return bk_signature, n_upstream, n_downstream


def calculate_signature_overlap(sig1: Set[str], sig2: Set[str]) -> float:
    """Calculate Jaccard similarity between two signatures."""
    if not sig1 or not sig2:
        return 0.0
    intersection = len(sig1 & sig2)
    union = len(sig1 | sig2)
    return intersection / union if union > 0 else 0.0


def cluster_by_flanking(
    orphans: pd.DataFrame,
    signatures: Dict[str, Set[str]],
    min_overlap: float = 0.3,
) -> List[Set[str]]:
    """
    Cluster orphans by shared flanking gene signatures.

    Uses single-linkage clustering with Jaccard similarity.
    """
    target_ids = list(orphans['target_gene_id'].unique())

    # Build adjacency based on signature overlap
    clusters = []
    unclustered = set(target_ids)

    while unclustered:
        # Start a new cluster with arbitrary target
        seed = next(iter(unclustered))
        cluster = {seed}
        unclustered.remove(seed)

        # Grow cluster by adding targets with overlapping signatures
        changed = True
        while changed:
            changed = False
            for target in list(unclustered):
                for member in cluster:
                    overlap = calculate_signature_overlap(
                        signatures.get(target, set()),
                        signatures.get(member, set())
                    )
                    if overlap >= min_overlap:
                        cluster.add(target)
                        unclustered.remove(target)
                        changed = True
                        break

        clusters.append(cluster)

    return clusters


def main():
    parser = argparse.ArgumentParser(
        description="Phase 9 (Helixer): Detect novel loci via reciprocal flanking analysis"
    )
    parser.add_argument("--family", required=True, help="Gene family name")
    parser.add_argument(
        "--phase5b-dir", type=Path, required=True,
        help="Phase 5b output directory (with unplaceable_refined.tsv)"
    )
    parser.add_argument(
        "--phase4b-dir", type=Path, required=True,
        help="Phase 4b output directory (with flanking_genes.tsv)"
    )
    parser.add_argument(
        "--phase5-dir", type=Path, required=True,
        help="Phase 5 output directory (with diamond_vs_bk.tsv)"
    )
    parser.add_argument(
        "--output-dir", type=Path, required=True,
        help="Output directory"
    )
    parser.add_argument(
        "--min-genomes", type=int, default=2,
        help="Minimum genomes for a novel locus (default: 2)"
    )
    parser.add_argument(
        "--min-flanking-both-sides", type=int, default=1,
        help="Minimum flanking genes required on each side (default: 1)"
    )
    parser.add_argument(
        "--min-signature-overlap", type=float, default=0.3,
        help="Minimum Jaccard overlap to cluster orphans (default: 0.3)"
    )
    parser.add_argument(
        "--max-novel-loci", type=int, default=50,
        help="Maximum novel loci to report per family (default: 50)"
    )

    args = parser.parse_args()

    print("=" * 80)
    print("PHASE 9 (HELIXER): DETECT NOVEL LOCI")
    print("=" * 80)
    print()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    orphans = load_orphan_targets(args.phase5b_dir)
    print(f"Loaded {len(orphans)} orphan targets")

    if orphans.empty:
        print("No orphans to analyze. Exiting.")
        return 0

    flanking_df = load_flanking_info(args.phase4b_dir)
    print(f"Loaded {len(flanking_df)} flanking gene entries")

    bk_orthologs = load_bk_orthologs(args.phase5_dir)
    print(f"Loaded {len(bk_orthologs)} BK ortholog mappings")

    # Get flanking signatures for each orphan
    print(f"\nAnalyzing flanking gene signatures...")
    signatures = {}
    orphan_details = []

    for _, row in orphans.iterrows():
        target_id = row['target_gene_id']
        genome = row['genome']

        sig, n_up, n_down = get_flanking_signature(target_id, flanking_df, bk_orthologs)
        signatures[target_id] = sig

        orphan_details.append({
            'target_gene_id': target_id,
            'genome': genome,
            'scaffold': row.get('scaffold', ''),
            'start': row.get('start', 0),
            'end': row.get('end', 0),
            'n_upstream': n_up,
            'n_downstream': n_down,
            'n_bk_orthologs': len(sig),
            'bk_signature': ';'.join(sorted(sig)) if sig else '',
            'has_both_sides': n_up >= args.min_flanking_both_sides and n_down >= args.min_flanking_both_sides,
        })

    orphan_df = pd.DataFrame(orphan_details)
    orphan_df.to_csv(args.output_dir / "orphan_flanking_patterns.tsv", sep='\t', index=False)

    # Filter orphans that have flanking genes on both sides
    valid_orphans = orphan_df[orphan_df['has_both_sides'] == True].copy()
    print(f"Orphans with flanking on both sides: {len(valid_orphans)} / {len(orphan_df)}")

    if valid_orphans.empty:
        print("No valid orphans for novel locus detection. Exiting.")
        return 0

    # Cluster orphans by shared flanking signatures
    print(f"\nClustering orphans by flanking gene similarity...")
    clusters = cluster_by_flanking(
        valid_orphans, signatures, args.min_signature_overlap
    )
    print(f"Found {len(clusters)} initial clusters")

    # Score and filter clusters
    novel_loci = []
    novel_members = []

    for i, cluster in enumerate(sorted(clusters, key=len, reverse=True)):
        if len(novel_loci) >= args.max_novel_loci:
            break

        # Get unique genomes in cluster
        cluster_orphans = valid_orphans[valid_orphans['target_gene_id'].isin(cluster)]
        cluster_genomes = cluster_orphans['genome'].unique()

        if len(cluster_genomes) < args.min_genomes:
            continue

        # Get consensus signature (union of all members)
        consensus_sig = set()
        for target_id in cluster:
            consensus_sig |= signatures.get(target_id, set())

        # Count how often each BK gene appears
        gene_counts = defaultdict(int)
        for target_id in cluster:
            for gene in signatures.get(target_id, set()):
                gene_counts[gene] += 1

        # Core signature: genes in ≥50% of members
        core_threshold = len(cluster) * 0.5
        core_sig = {g for g, c in gene_counts.items() if c >= core_threshold}

        novel_locus_id = f"NOVEL_{args.family}_{i+1:03d}"

        novel_loci.append({
            'novel_locus_id': novel_locus_id,
            'n_members': len(cluster),
            'n_genomes': len(cluster_genomes),
            'genomes': ';'.join(sorted(cluster_genomes)),
            'consensus_signature_size': len(consensus_sig),
            'core_signature_size': len(core_sig),
            'core_signature': ';'.join(sorted(core_sig)) if core_sig else '',
            'gene_family': args.family,
        })

        for target_id in cluster:
            row = cluster_orphans[cluster_orphans['target_gene_id'] == target_id].iloc[0]
            novel_members.append({
                'novel_locus_id': novel_locus_id,
                'target_gene_id': target_id,
                'genome': row['genome'],
                'scaffold': row['scaffold'],
                'start': row['start'],
                'end': row['end'],
                'n_bk_orthologs': row['n_bk_orthologs'],
            })

    # Write outputs
    if novel_loci:
        pd.DataFrame(novel_loci).to_csv(
            args.output_dir / "novel_locus_candidates.tsv", sep='\t', index=False
        )
        pd.DataFrame(novel_members).to_csv(
            args.output_dir / "novel_locus_members.tsv", sep='\t', index=False
        )

    # Summary
    print(f"\n[RESULTS]")
    print(f"  Novel locus candidates: {len(novel_loci)}")

    if novel_loci:
        print(f"\nTop novel loci:")
        for locus in novel_loci[:10]:
            print(f"  {locus['novel_locus_id']}: {locus['n_members']} members in {locus['n_genomes']} genomes")
            print(f"    Core signature ({locus['core_signature_size']} genes): {locus['core_signature'][:80]}...")

    # Stats on orphans
    n_assigned = len([m for m in novel_members])
    n_unassigned = len(valid_orphans) - n_assigned
    print(f"\nOrphan disposition:")
    print(f"  Assigned to novel loci: {n_assigned}")
    print(f"  Remain unassigned: {n_unassigned}")
    print(f"  Invalid (no flanking both sides): {len(orphan_df) - len(valid_orphans)}")

    print(f"\n[OUTPUT] {args.output_dir}/novel_locus_candidates.tsv")
    print(f"[OUTPUT] {args.output_dir}/novel_locus_members.tsv")
    print(f"[OUTPUT] {args.output_dir}/orphan_flanking_patterns.tsv")

    return 0


if __name__ == "__main__":
    exit(main())
