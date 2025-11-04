#!/usr/bin/env python3
"""
Step 08b: Generate gene-type summary matrices.

Creates summary matrices by gene family showing:
- All genomes (rows) sorted by phylogenetic order
- All loci of that gene family (columns)
- Includes both syntenic and unplaceable targets
"""

from pathlib import Path
import pandas as pd
from collections import OrderedDict, defaultdict
import config

def load_species_and_phylo():
    """Load species mapping and phylogenetic order."""
    species_map = {}
    phylo_order_map = {}

    with open(config.SPECIES_MAP_FILE) as f:
        header = f.readline()  # Skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                genome_id = parts[0]
                species_map[genome_id] = parts[1]
                try:
                    phylo_order_map[genome_id] = int(parts[3])
                except (ValueError, IndexError):
                    phylo_order_map[genome_id] = 999

    return species_map, phylo_order_map

def create_gene_family_matrix(gene_family, loci_list, all_targets_df, synteny_blocks_df, species_map, phylo_order_map):
    """Create summary matrix for one gene family."""

    print(f"\n  Processing {gene_family}...")

    # Get all targets for this gene family
    family_targets = all_targets_df[all_targets_df['gene_family'] == gene_family]

    # Group targets by genome
    targets_by_genome = defaultdict(list)
    for _, target in family_targets.iterrows():
        targets_by_genome[target['genome']].append(target)

    # Get synteny blocks for this gene family's loci
    synteny_by_genome = {}
    for locus in loci_list:
        locus_blocks = synteny_blocks_df[synteny_blocks_df['locus_id'] == locus]
        for _, block in locus_blocks.iterrows():
            genome = block['genome']
            if genome not in synteny_by_genome:
                synteny_by_genome[genome] = []
            synteny_by_genome[genome].append(locus)

    # Get all unique locus categories
    # Syntenic loci (from loci_list) + unplaceable
    locus_categories = loci_list + [f"{gene_family}_unplaceable"]

    # Build matrix rows
    matrix_rows = []

    # Process all genomes
    all_genomes = set(species_map.keys())
    genomes_with_targets = set(targets_by_genome.keys())
    all_genomes.update(genomes_with_targets)

    for genome in all_genomes:
        row = OrderedDict()

        # Metadata
        row['genome_id'] = genome
        row['species'] = species_map.get(genome, genome)
        row['phylo_order'] = phylo_order_map.get(genome, 999)

        # Add synteny block presence
        has_synteny_block = genome in synteny_by_genome
        row[f'{gene_family}_synteny_block'] = 'Y' if has_synteny_block else 'N'

        # Count targets by category
        genome_targets = targets_by_genome.get(genome, [])

        # Initialize all columns
        for category in locus_categories:
            row[category] = ""

        # Fill in target counts
        category_counts = defaultdict(list)
        for target in genome_targets:
            assigned_to = target.get('assigned_to', target.get('description', 'unknown'))

            # Build target string with length and status
            length = int(target.get('orf_length_aa', 0))
            status = target.get('status', 'unknown')
            status_letter = 'F' if status == 'functional' else 'P'
            target_str = f"{length}{status_letter}"

            category_counts[assigned_to].append(target_str)

        # Format counts for each category
        for category, target_list in category_counts.items():
            if category in row:  # Only if it's a known category
                if target_list:
                    row[category] = f"[{'; '.join(target_list)}]"

        # Add total count
        row['total'] = len(genome_targets)

        matrix_rows.append(row)

    # Create DataFrame and sort
    if matrix_rows:
        matrix_df = pd.DataFrame(matrix_rows)
        matrix_df = matrix_df.sort_values('phylo_order', ascending=True)
        return matrix_df

    return pd.DataFrame()

def main():
    """Main execution function."""
    print("=" * 80)
    print("STEP 07b: GENERATE GENE-TYPE SUMMARY MATRICES")
    print("=" * 80)

    # Create output directory
    config.SUMMARY_MATRICES_DIR.mkdir(exist_ok=True, parents=True)

    # Load all required data
    print("\n[1] Loading data...")

    # Locus definitions
    loci_df = pd.read_csv(config.LOCI_DEFINITIONS_FILE, sep='\t')
    print(f"  Loaded {len(loci_df)} locus definitions")

    # Load synteny blocks to track presence/absence
    synteny_blocks_file = config.STEP03_FILTERED / "synteny_blocks_filtered.tsv"
    if synteny_blocks_file.exists():
        synteny_blocks_df = pd.read_csv(synteny_blocks_file, sep='\t')
        print(f"  Loaded {len(synteny_blocks_df)} synteny blocks")
    else:
        synteny_blocks_df = pd.DataFrame()
        print("  No synteny blocks found")

    # All classified targets (including unplaceable)
    targets_file = config.STEP05_CLASSIFIED / "all_targets_classified.tsv"
    if targets_file.exists():
        all_targets_df = pd.read_csv(targets_file, sep='\t')
        print(f"  Loaded {len(all_targets_df)} classified targets")
    else:
        all_targets_df = pd.DataFrame()
        print("  No classified targets found")

    # Species mapping
    species_map, phylo_order_map = load_species_and_phylo()
    print(f"  Loaded species mapping for {len(species_map)} genomes")

    # Get unique gene families
    gene_families = loci_df['gene_family'].unique()
    print(f"  Gene families: {', '.join(gene_families)}")

    # Process each gene family
    print("\n[2] Generating summary matrices...")

    for gene_family in gene_families:
        # Get loci for this gene family
        family_loci = loci_df[loci_df['gene_family'] == gene_family]['locus_id'].tolist()

        # Create matrix
        matrix_df = create_gene_family_matrix(
            gene_family, family_loci, all_targets_df, synteny_blocks_df,
            species_map, phylo_order_map
        )

        if not matrix_df.empty:
            # Save matrix
            output_file = config.SUMMARY_MATRICES_DIR / f"{gene_family}_summary_matrix.tsv"
            matrix_df.to_csv(output_file, sep='\t', index=False)

            print(f"    Saved: {output_file.name}")
            print(f"    Rows: {len(matrix_df)}")

            # Summary statistics
            with_targets = matrix_df[matrix_df['total'] > 0]
            print(f"    Genomes with targets: {len(with_targets)} ({len(with_targets)/len(matrix_df)*100:.1f}%)")

            # Count syntenic vs unplaceable
            syntenic_count = 0
            unplaceable_count = 0
            for col in family_loci:
                syntenic_count += len(matrix_df[matrix_df[col] != ""])

            unplaceable_col = f"{gene_family}_unplaceable"
            if unplaceable_col in matrix_df.columns:
                unplaceable_count = len(matrix_df[matrix_df[unplaceable_col] != ""])

            if syntenic_count + unplaceable_count > 0:
                print(f"    Syntenic: {syntenic_count} ({syntenic_count/(syntenic_count+unplaceable_count)*100:.1f}%)")
                print(f"    Unplaceable: {unplaceable_count} ({unplaceable_count/(syntenic_count+unplaceable_count)*100:.1f}%)")

    # NOTE: "All genes combined summary" removed as it's redundant for single gene families
    # and confusing for users. Each gene family has its own summary matrix.

    # Overall summary
    print("\n" + "=" * 80)
    print("SUMMARY MATRICES COMPLETE")
    print("=" * 80)
    print(f"\nMatrices saved to: {config.SUMMARY_MATRICES_DIR}")

    if not all_targets_df.empty:
        print(f"\nOverall statistics:")
        print(f"  Total targets: {len(all_targets_df)}")
        print(f"  Syntenic: {len(all_targets_df[all_targets_df['placement'] == 'synteny'])}")
        print(f"  Unplaceable: {len(all_targets_df[all_targets_df['placement'] == 'unplaceable'])}")
        print(f"  Genomes with targets: {all_targets_df['genome'].nunique()}")

    print("\nPipeline complete! Review outputs in:")
    print(f"  - Locus matrices: {config.LOCUS_MATRICES_DIR}")
    print(f"  - Summary matrices: {config.SUMMARY_MATRICES_DIR}")

if __name__ == "__main__":
    main()