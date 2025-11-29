#!/usr/bin/env python3
"""
Phase 4b (Helixer): Extract flanking genes around target hits from GFF3.

For each target gene found in Phase 4, extract N genes upstream and downstream
from the Helixer GFF3 annotations. This is much simpler than tblastn-based
flanking extraction since we already have gene predictions.

Outputs:
- flanking_genes.tsv       # Target + flanking gene info
- flanking_proteins.faa    # FASTA of flanking protein sequences
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


def get_flanking_genes(
    target_gene_id: str,
    scaffold: str,
    gene_index: Dict[str, List[Dict]],
    n_upstream: int = 10,
    n_downstream: int = 10,
) -> Tuple[List[Dict], List[Dict]]:
    """
    Get N upstream and N downstream genes from a target gene.

    Returns (upstream_genes, downstream_genes) where genes are sorted
    by distance from target (closest first).
    """
    scaffold_genes = gene_index.get(scaffold, [])
    if not scaffold_genes:
        return [], []

    # Find target gene index
    target_idx = None
    for i, gene in enumerate(scaffold_genes):
        if gene['gene_id'] == target_gene_id:
            target_idx = i
            break

    if target_idx is None:
        return [], []

    # Get upstream (genes before target)
    upstream = scaffold_genes[max(0, target_idx - n_upstream):target_idx]
    upstream = list(reversed(upstream))  # Closest first

    # Get downstream (genes after target)
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
        "--n-flanking", type=int, default=10,
        help="Number of flanking genes on each side (default: 10)"
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

    # Group targets by genome
    genomes = targets_df['genome'].unique()
    print(f"Processing {len(genomes)} genomes")
    print()

    all_flanking = []
    all_proteins = {}

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

        # Get targets for this genome
        genome_targets = targets_df[targets_df['genome'] == genome_id]
        print(f"  Processing {len(genome_targets)} targets...")

        for _, target in genome_targets.iterrows():
            target_gene_id = target['target_gene_id']
            scaffold = target['scaffold']

            # Get flanking genes
            upstream, downstream = get_flanking_genes(
                target_gene_id, scaffold, gene_index,
                n_upstream=args.n_flanking,
                n_downstream=args.n_flanking,
            )

            # Record flanking info
            for i, gene in enumerate(upstream):
                protein_id = gene_id_to_protein_id(gene['gene_id'])
                flanking_record = {
                    'genome': genome_id,
                    'target_gene_id': target_gene_id,
                    'target_scaffold': scaffold,
                    'target_parent_locus': target['parent_locus'],
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

                # Collect protein sequence
                if protein_id in proteins:
                    all_proteins[protein_id] = proteins[protein_id]

            for i, gene in enumerate(downstream):
                protein_id = gene_id_to_protein_id(gene['gene_id'])
                flanking_record = {
                    'genome': genome_id,
                    'target_gene_id': target_gene_id,
                    'target_scaffold': scaffold,
                    'target_parent_locus': target['parent_locus'],
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

    print("\nPhase 4b (Helixer) complete.")
    return 0


if __name__ == "__main__":
    exit(main())
