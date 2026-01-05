#!/usr/bin/env python3
"""
Phase 4b (Helixer): Extract flanking genes around target hits from GFF3.

For each TANDEM CLUSTER of targets (not individual targets), extract N genes
upstream and downstream from the Helixer GFF3 annotations. This prevents
family paralogs from being counted as flanking genes.

Key behavior:
- Targets within 50kb are grouped into tandem clusters
- Flanking genes are extracted for the CLUSTER boundaries
- Family members (other targets) are EXCLUDED from flanking counts
- Each target in a cluster shares the same flanking gene set

Outputs:
- flanking_genes.tsv       # Cluster + flanking gene info
- flanking_proteins.faa    # FASTA of flanking protein sequences
- tandem_clusters.tsv      # Tandem cluster assignments for targets
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import pandas as pd


def parse_helixer_gff3_genes(gff3_path: Path) -> List[Dict]:
    """
    Parse Helixer GFF3 to extract all gene features with coordinates.
    Returns list sorted by scaffold and position.
    """
    genes = []

    with open(gff3_path) as f:
        for line in f:
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            scaffold, source, feature, start, end, score, strand, phase, attributes = parts

            if feature != 'gene':
                continue

            # Parse attributes
            attr_dict = {}
            for attr in attributes.split(';'):
                if '=' in attr:
                    key, value = attr.split('=', 1)
                    attr_dict[key] = value

            gene_id = attr_dict.get('ID', '')
            if not gene_id:
                continue

            genes.append({
                'gene_id': gene_id,
                'scaffold': scaffold,
                'strand': strand,
                'start': int(start),
                'end': int(end),
            })

    # Sort by scaffold, then by start position
    genes.sort(key=lambda x: (x['scaffold'], x['start']))

    return genes


def build_gene_index(genes: List[Dict]) -> Dict[str, List[Dict]]:
    """Build index of genes by scaffold for efficient lookup."""
    index = defaultdict(list)
    for gene in genes:
        index[gene['scaffold']].append(gene)
    return index


# Tandem clustering thresholds (both must be satisfied)
MAX_INTERVENING_GENES = 1  # Max non-target genes between tandem targets
MAX_TANDEM_DISTANCE_KB = 200  # Max physical distance (sanity check)


def cluster_targets_by_gene_adjacency(
    targets_df: pd.DataFrame,
    genome_id: str,
    gene_index: Dict[str, List[Dict]],
) -> List[Dict]:
    """
    Cluster targets that are adjacent in gene order AND physically close.

    Two targets are "tandem" if BOTH conditions are met:
    1. <= MAX_INTERVENING_GENES non-target genes between them
    2. <= MAX_TANDEM_DISTANCE_KB physical distance

    This prevents calling distant "adjacent" genes tandem (e.g., if there's
    a large intergenic region or assembly gap).

    Returns list of cluster dicts with members, boundaries, etc.
    """
    genome_targets = targets_df[targets_df['genome'] == genome_id].copy()
    if genome_targets.empty:
        return []

    # Get all target gene IDs for this genome
    target_gene_ids = set(genome_targets['target_gene_id'].tolist())

    # Sort targets by scaffold and position
    genome_targets = genome_targets.sort_values(['scaffold', 'start'])

    clusters = []
    current_cluster = []

    for _, row in genome_targets.iterrows():
        if not current_cluster:
            current_cluster.append(row)
        else:
            prev_row = current_cluster[-1]

            # Check if same scaffold
            if row['scaffold'] != prev_row['scaffold']:
                clusters.append(_make_cluster_dict(current_cluster, genome_id, len(clusters)))
                current_cluster = [row]
                continue

            # Check physical distance
            physical_dist_kb = abs(row['start'] - prev_row['end']) / 1000.0

            # Count non-target genes between prev and current target
            n_intervening = count_intervening_genes(
                scaffold=row['scaffold'],
                start_pos=prev_row['end'],
                end_pos=row['start'],
                gene_index=gene_index,
                exclude_gene_ids=target_gene_ids,
            )

            # Both conditions must be met for tandem clustering
            is_tandem = (n_intervening <= MAX_INTERVENING_GENES and
                         physical_dist_kb <= MAX_TANDEM_DISTANCE_KB)

            if is_tandem:
                current_cluster.append(row)
            else:
                clusters.append(_make_cluster_dict(current_cluster, genome_id, len(clusters)))
                current_cluster = [row]

    # Don't forget the last cluster
    if current_cluster:
        clusters.append(_make_cluster_dict(current_cluster, genome_id, len(clusters)))

    return clusters


def count_intervening_genes(
    scaffold: str,
    start_pos: int,
    end_pos: int,
    gene_index: Dict[str, List[Dict]],
    exclude_gene_ids: set,
) -> int:
    """Count non-target genes between two positions."""
    scaffold_genes = gene_index.get(scaffold, [])
    count = 0
    for gene in scaffold_genes:
        # Gene is "between" if it starts after start_pos and ends before end_pos
        if gene['start'] > start_pos and gene['end'] < end_pos:
            if gene['gene_id'] not in exclude_gene_ids:
                count += 1
    return count


def cluster_targets_into_tandem_groups(
    targets_df: pd.DataFrame,
    genome_id: str,
    gene_index: Dict[str, List[Dict]] = None,
) -> List[Dict]:
    """
    Cluster targets into tandem groups.

    If gene_index is provided, uses gene adjacency (more accurate).
    Otherwise falls back to 50kb physical distance.
    """
    if gene_index is not None:
        return cluster_targets_by_gene_adjacency(targets_df, genome_id, gene_index)

    # Fallback: use physical distance (50kb threshold)
    MAX_TANDEM_DISTANCE_KB = 50
    genome_targets = targets_df[targets_df['genome'] == genome_id].copy()
    if genome_targets.empty:
        return []

    genome_targets = genome_targets.sort_values(['scaffold', 'start'])
    clusters = []
    current_cluster = []

    for _, row in genome_targets.iterrows():
        if not current_cluster:
            current_cluster.append(row)
        else:
            prev_row = current_cluster[-1]
            if (row['scaffold'] == prev_row['scaffold'] and
                abs(row['start'] - prev_row['end']) < MAX_TANDEM_DISTANCE_KB * 1000):
                current_cluster.append(row)
            else:
                clusters.append(_make_cluster_dict(current_cluster, genome_id, len(clusters)))
                current_cluster = [row]

    if current_cluster:
        clusters.append(_make_cluster_dict(current_cluster, genome_id, len(clusters)))

    return clusters


def _make_cluster_dict(members: List, genome_id: str, cluster_idx: int) -> Dict:
    """Create cluster dictionary from list of target rows."""
    member_df = pd.DataFrame(members)
    scaffold = member_df['scaffold'].iloc[0]

    # Cluster boundaries
    cluster_start = int(member_df['start'].min())
    cluster_end = int(member_df['end'].max())

    # Member gene IDs
    member_gene_ids = set(member_df['target_gene_id'].tolist())

    return {
        'cluster_id': f"{genome_id}_cluster_{cluster_idx}",
        'genome': genome_id,
        'scaffold': scaffold,
        'start': cluster_start,
        'end': cluster_end,
        'span_kb': (cluster_end - cluster_start) / 1000.0,
        'n_targets': len(member_df),
        'member_gene_ids': member_gene_ids,
        'member_targets': member_df,
    }


def get_flanking_genes_for_cluster(
    scaffold: str,
    cluster_start: int,
    cluster_end: int,
    gene_index: Dict[str, List[Dict]],
    exclude_gene_ids: set,
    n_upstream: int = 10,
    n_downstream: int = 10,
) -> Tuple[List[Dict], List[Dict]]:
    """
    Get N upstream and N downstream genes from a tandem cluster.

    Unlike per-target extraction, this:
    1. Uses cluster boundaries (not single gene position)
    2. EXCLUDES family members (other targets) from flanking counts
    3. Keeps extracting until we have N non-family flanking genes

    Returns (upstream_genes, downstream_genes) where genes are sorted
    by distance from cluster (closest first).
    """
    scaffold_genes = gene_index.get(scaffold, [])
    if not scaffold_genes:
        return [], []

    # Find genes that overlap or are within the cluster
    # These define the cluster boundaries in the gene index
    cluster_gene_indices = []
    for i, gene in enumerate(scaffold_genes):
        # Gene overlaps cluster if it doesn't end before or start after
        if gene['end'] >= cluster_start and gene['start'] <= cluster_end:
            cluster_gene_indices.append(i)

    if not cluster_gene_indices:
        # No genes found in cluster region - try to find insertion point
        for i, gene in enumerate(scaffold_genes):
            if gene['start'] > cluster_end:
                cluster_gene_indices = [i]
                break
        if not cluster_gene_indices:
            return [], []

    first_cluster_idx = min(cluster_gene_indices)
    last_cluster_idx = max(cluster_gene_indices)

    # Get upstream genes (before cluster), excluding family members
    upstream = []
    idx = first_cluster_idx - 1
    while len(upstream) < n_upstream and idx >= 0:
        gene = scaffold_genes[idx]
        if gene['gene_id'] not in exclude_gene_ids:
            upstream.append(gene)
        idx -= 1

    # Get downstream genes (after cluster), excluding family members
    downstream = []
    idx = last_cluster_idx + 1
    while len(downstream) < n_downstream and idx < len(scaffold_genes):
        gene = scaffold_genes[idx]
        if gene['gene_id'] not in exclude_gene_ids:
            downstream.append(gene)
        idx += 1

    return upstream, downstream


def get_flanking_genes(
    target_gene_id: str,
    scaffold: str,
    gene_index: Dict[str, List[Dict]],
    n_upstream: int = 10,
    n_downstream: int = 10,
) -> Tuple[List[Dict], List[Dict]]:
    """
    DEPRECATED: Use get_flanking_genes_for_cluster() instead.

    Get N upstream and N downstream genes from a target gene.
    """
    scaffold_genes = gene_index.get(scaffold, [])
    if not scaffold_genes:
        return [], []

    target_idx = None
    for i, gene in enumerate(scaffold_genes):
        if gene['gene_id'] == target_gene_id:
            target_idx = i
            break

    if target_idx is None:
        return [], []

    upstream = scaffold_genes[max(0, target_idx - n_upstream):target_idx]
    upstream = list(reversed(upstream))

    downstream = scaffold_genes[target_idx + 1:target_idx + 1 + n_downstream]

    return upstream, downstream


def load_helixer_proteins(faa_path: Path) -> Dict[str, str]:
    """Load Helixer protein sequences into dict."""
    proteins = {}
    current_id = None
    current_seq = []

    with open(faa_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    proteins[current_id] = ''.join(current_seq)
                # Protein ID is the header without '>'
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)

        if current_id:
            proteins[current_id] = ''.join(current_seq)

    return proteins


def gene_id_to_protein_id(gene_id: str) -> str:
    """Convert gene ID to protein ID (add .1 suffix)."""
    return f"{gene_id}.1"


def main():
    parser = argparse.ArgumentParser(
        description="Phase 4b (Helixer): Extract flanking genes from GFF3"
    )
    parser.add_argument(
        "--family", required=True,
        help="Gene family name"
    )
    parser.add_argument(
        "--phase4-output", type=Path, required=True,
        help="Path to Phase 4 Helixer output (all_target_loci.tsv)"
    )
    parser.add_argument(
        "--helixer-dir", type=Path, required=True,
        help="Directory containing Helixer outputs"
    )
    parser.add_argument(
        "--output-dir", type=Path, required=True,
        help="Output directory"
    )
    parser.add_argument(
        "--n-flanking", type=int, default=25,
        help="Number of flanking genes on each side (default: 25)"
    )

    args = parser.parse_args()

    print("=" * 80)
    print("PHASE 4B (HELIXER): EXTRACT FLANKING GENES")
    print("=" * 80)
    print()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Load Phase 4 targets
    targets_file = args.phase4_output
    if not targets_file.exists():
        print(f"[ERROR] Phase 4 output not found: {targets_file}")
        return 1

    targets_df = pd.read_csv(targets_file, sep='\t')
    print(f"Loaded {len(targets_df)} target genes from Phase 4")

    # Build global set of all target gene IDs (to exclude from flanking)
    all_target_gene_ids = set(targets_df['target_gene_id'].tolist())
    print(f"Will exclude {len(all_target_gene_ids)} family members from flanking counts")

    # Group targets by genome
    genomes = targets_df['genome'].unique()
    print(f"Processing {len(genomes)} genomes")
    print()

    all_flanking = []
    all_proteins = {}
    all_cluster_records = []  # Track tandem cluster assignments

    for genome_id in genomes:
        print(f"\n[{genome_id}]")

        # Find Helixer files
        genome_dir = args.helixer_dir / genome_id
        gff3_file = genome_dir / f"{genome_id}_helixer.gff3"
        faa_file = genome_dir / f"{genome_id}_helixer_proteins.faa"

        if not gff3_file.exists():
            print(f"  Skipping: no GFF3 found")
            continue

        # Parse GFF3
        print(f"  Parsing GFF3...")
        genes = parse_helixer_gff3_genes(gff3_file)
        gene_index = build_gene_index(genes)
        print(f"  Found {len(genes)} genes on {len(gene_index)} scaffolds")

        # Load proteins if available
        proteins = {}
        if faa_file.exists():
            proteins = load_helixer_proteins(faa_file)
            print(f"  Loaded {len(proteins)} protein sequences")

        # Cluster targets into tandem groups using gene adjacency
        clusters = cluster_targets_into_tandem_groups(targets_df, genome_id, gene_index)
        n_targets = len(targets_df[targets_df['genome'] == genome_id])
        n_singletons = sum(1 for c in clusters if c['n_targets'] == 1)
        n_tandem = len(clusters) - n_singletons
        print(f"  Clustered {n_targets} targets into {len(clusters)} groups ({n_tandem} tandem, {n_singletons} singletons)")

        # Process each tandem cluster
        for cluster in clusters:
            cluster_id = cluster['cluster_id']
            scaffold = cluster['scaffold']
            member_gene_ids = cluster['member_gene_ids']

            # Record cluster membership for each target
            for target_gene_id in member_gene_ids:
                all_cluster_records.append({
                    'genome': genome_id,
                    'target_gene_id': target_gene_id,
                    'tandem_cluster_id': cluster_id,
                    'cluster_scaffold': scaffold,
                    'cluster_start': cluster['start'],
                    'cluster_end': cluster['end'],
                    'cluster_span_kb': cluster['span_kb'],
                    'n_targets_in_cluster': cluster['n_targets'],
                    'gene_family': args.family,
                })

            # Get flanking genes for the CLUSTER (not individual targets)
            # Exclude ALL target gene IDs from this genome
            genome_target_ids = set(targets_df[targets_df['genome'] == genome_id]['target_gene_id'])

            upstream, downstream = get_flanking_genes_for_cluster(
                scaffold=scaffold,
                cluster_start=cluster['start'],
                cluster_end=cluster['end'],
                gene_index=gene_index,
                exclude_gene_ids=genome_target_ids,
                n_upstream=args.n_flanking,
                n_downstream=args.n_flanking,
            )

            # Record flanking info - associate with ALL targets in cluster
            for target_gene_id in member_gene_ids:
                # Get target info
                target_row = targets_df[
                    (targets_df['genome'] == genome_id) &
                    (targets_df['target_gene_id'] == target_gene_id)
                ].iloc[0]

                for i, gene in enumerate(upstream):
                    protein_id = gene_id_to_protein_id(gene['gene_id'])
                    flanking_record = {
                        'genome': genome_id,
                        'target_gene_id': target_gene_id,
                        'target_scaffold': scaffold,
                        'target_parent_locus': target_row['parent_locus'],
                        'tandem_cluster_id': cluster_id,
                        'flanking_gene_id': gene['gene_id'],
                        'flanking_protein_id': protein_id,
                        'flanking_scaffold': gene['scaffold'],
                        'flanking_start': gene['start'],
                        'flanking_end': gene['end'],
                        'flanking_strand': gene['strand'],
                        'position': f"upstream_{i+1}",
                        'distance_from_target': i + 1,
                        'gene_family': args.family,
                    }
                    all_flanking.append(flanking_record)

                    if protein_id in proteins:
                        all_proteins[protein_id] = proteins[protein_id]

                for i, gene in enumerate(downstream):
                    protein_id = gene_id_to_protein_id(gene['gene_id'])
                    flanking_record = {
                        'genome': genome_id,
                        'target_gene_id': target_gene_id,
                        'target_scaffold': scaffold,
                        'target_parent_locus': target_row['parent_locus'],
                        'tandem_cluster_id': cluster_id,
                        'flanking_gene_id': gene['gene_id'],
                        'flanking_protein_id': protein_id,
                        'flanking_scaffold': gene['scaffold'],
                        'flanking_start': gene['start'],
                        'flanking_end': gene['end'],
                        'flanking_strand': gene['strand'],
                        'position': f"downstream_{i+1}",
                        'distance_from_target': i + 1,
                        'gene_family': args.family,
                    }
                    all_flanking.append(flanking_record)

                    if protein_id in proteins:
                        all_proteins[protein_id] = proteins[protein_id]

        print(f"  Collected {len([f for f in all_flanking if f['genome'] == genome_id])} flanking records")

    # Write outputs
    if all_flanking:
        flanking_df = pd.DataFrame(all_flanking)
        flanking_file = args.output_dir / "flanking_genes.tsv"
        flanking_df.to_csv(flanking_file, sep='\t', index=False)
        print(f"\n[OUTPUT] Wrote {len(flanking_df)} flanking records to {flanking_file}")

        # Write protein FASTA
        faa_file = args.output_dir / "flanking_proteins.faa"
        with open(faa_file, 'w') as f:
            for protein_id, seq in all_proteins.items():
                f.write(f">{protein_id}\n")
                for i in range(0, len(seq), 60):
                    f.write(f"{seq[i:i+60]}\n")
        print(f"[OUTPUT] Wrote {len(all_proteins)} protein sequences to {faa_file}")
    else:
        print("\n[WARNING] No flanking genes found!")

    # Write tandem cluster assignments
    if all_cluster_records:
        cluster_df = pd.DataFrame(all_cluster_records)
        cluster_file = args.output_dir / "tandem_clusters.tsv"
        cluster_df.to_csv(cluster_file, sep='\t', index=False)

        # Summary stats
        n_tandem_targets = cluster_df[cluster_df['n_targets_in_cluster'] > 1]['target_gene_id'].nunique()
        n_singleton_targets = cluster_df[cluster_df['n_targets_in_cluster'] == 1]['target_gene_id'].nunique()
        max_cluster_size = cluster_df['n_targets_in_cluster'].max()

        print(f"[OUTPUT] Wrote {len(cluster_df)} cluster assignments to {cluster_file}")
        print(f"         Tandem targets: {n_tandem_targets}, Singletons: {n_singleton_targets}, Max cluster: {max_cluster_size}")

    print("\nPhase 4b (Helixer) complete.")
    return 0


if __name__ == "__main__":
    exit(main())
