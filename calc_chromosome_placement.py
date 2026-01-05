#!/usr/bin/env python3
"""
Calculate chromosome placement statistics for Helixer annotations.

For each species, count how many genes are on the 10 RagTag-assigned chromosomes
(CM021338-CM021347) vs unplaced scaffolds.

This provides a quality metric showing how well each genome aligned to the
Belonocnema kinseyi reference during RagTag scaffolding.
"""

import sys
from pathlib import Path
from collections import defaultdict

# Add pipeline dir to path for imports
sys.path.insert(0, str(Path(__file__).parent))

from data_paths import get_data_dir, list_genomes, get_gff_path, get_species_mapping_file
import pandas as pd

# The 10 Belonocnema kinseyi reference chromosomes (base names)
BK_CHROMOSOME_BASES = [
    'CM021338.1', 'CM021339.1', 'CM021340.1', 'CM021341.1', 'CM021342.1',
    'CM021343.1', 'CM021344.1', 'CM021345.1', 'CM021346.1', 'CM021347.1'
]


def is_on_chromosome(seqid: str) -> bool:
    """Check if a sequence ID corresponds to a BK reference chromosome.

    Handles both:
    - Original BK names: CM021338.1
    - RagTag-scaffolded names: CM021338.1_RagTag
    """
    for chrom in BK_CHROMOSOME_BASES:
        if seqid == chrom or seqid == f"{chrom}_RagTag":
            return True
    return False


def count_genes_by_placement(gff_path: Path) -> dict:
    """
    Count genes on chromosomes vs unplaced scaffolds.

    Returns dict with:
        - on_chromosome: genes on CM* chromosomes
        - unplaced: genes on other scaffolds
        - total: total gene count
        - by_chrom: dict of counts per chromosome
    """
    on_chromosome = 0
    unplaced = 0
    by_chrom = defaultdict(int)

    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            seqid = parts[0]
            feature_type = parts[2]

            # Only count genes (not mRNA, CDS, etc.)
            if feature_type != 'gene':
                continue

            # Check if on a BK chromosome (handles _RagTag suffix)
            if is_on_chromosome(seqid):
                on_chromosome += 1
                by_chrom[seqid] += 1
            else:
                unplaced += 1
                # Track unplaced scaffold names too
                by_chrom[f"unplaced:{seqid[:20]}"] += 1

    total = on_chromosome + unplaced
    return {
        'on_chromosome': on_chromosome,
        'unplaced': unplaced,
        'total': total,
        'pct_on_chromosome': (on_chromosome / total * 100) if total > 0 else 0,
        'by_chrom': dict(by_chrom)
    }


def main():
    # Load species mapping for nice names
    species_file = get_species_mapping_file()
    species_df = pd.read_csv(species_file, sep='\t')

    # Create lookup from accession/folder to species name
    accession_to_species = dict(zip(species_df['accession'], species_df['species']))

    results = []

    print("Calculating chromosome placement for all genomes...")
    print("=" * 90)

    for genome_id in list_genomes():
        gff_path = get_gff_path(genome_id)
        if gff_path is None:
            continue

        stats = count_genes_by_placement(gff_path)

        # Get species name
        species_name = accession_to_species.get(genome_id, genome_id.replace('_', ' '))

        results.append({
            'genome_id': genome_id,
            'species': species_name,
            'total_genes': stats['total'],
            'on_chromosome': stats['on_chromosome'],
            'unplaced': stats['unplaced'],
            'pct_on_chromosome': stats['pct_on_chromosome']
        })

    # Sort by placement percentage
    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values('pct_on_chromosome', ascending=False)

    # Print summary
    print(f"{'Species':<40} {'Total':>8} {'On Chr':>8} {'Unplaced':>8} {'% Chr':>7}")
    print("-" * 90)

    for _, row in results_df.iterrows():
        print(f"{row['species'][:40]:<40} {row['total_genes']:>8} {row['on_chromosome']:>8} "
              f"{row['unplaced']:>8} {row['pct_on_chromosome']:>6.1f}%")

    print("=" * 90)

    # Overall stats
    total_genes = results_df['total_genes'].sum()
    total_on_chr = results_df['on_chromosome'].sum()
    overall_pct = total_on_chr / total_genes * 100 if total_genes > 0 else 0

    print(f"\nOverall: {total_on_chr:,} / {total_genes:,} genes on chromosomes ({overall_pct:.1f}%)")
    print(f"Mean per species: {results_df['pct_on_chromosome'].mean():.1f}%")
    print(f"Median per species: {results_df['pct_on_chromosome'].median():.1f}%")
    print(f"Range: {results_df['pct_on_chromosome'].min():.1f}% - {results_df['pct_on_chromosome'].max():.1f}%")

    # Save to TSV
    output_file = get_data_dir() / 'chromosome_placement_stats.tsv'
    results_df.to_csv(output_file, sep='\t', index=False)
    print(f"\nSaved to: {output_file}")


if __name__ == '__main__':
    main()
