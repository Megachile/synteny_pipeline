#!/usr/bin/env python3
"""
Step 08b: Generate gene-type summary matrices.

Creates summary matrices by gene family showing:
- All genomes (rows) sorted by phylogenetic order
- All loci of that gene family (columns)
- Includes both syntenic and unplaceable targets

Usage:
    python 08b_generate_summary_matrices.py \\
        --locus-defs <path/to/locus_definitions.tsv> \\
        --blocks <path/to/synteny_blocks_filtered.tsv> \\
        --targets <path/to/all_targets_classified.tsv> \\
        --species-map <path/to/gca_to_species.tsv> \\
        --extracted-seqs <path/to/06_extracted_sequences> \\
        --output-dir <path/to/gene_type_summaries>
"""

from pathlib import Path
import pandas as pd
from collections import OrderedDict, defaultdict
import argparse
import re

def load_species_and_phylo(species_map_file):
    """Load species mapping and phylogenetic order."""
    species_map = {}
    phylo_order_map = {}

    with open(species_map_file) as f:
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

def load_extracted_seq_metadata(extraction_dir):
    """
    Load metadata from Phase 6 extracted sequences.

    Returns dict: {locus_name: {'length_aa': int, 'status': str}}
    """
    metadata = {}

    if not extraction_dir or not Path(extraction_dir).exists():
        return metadata

    extraction_path = Path(extraction_dir)

    # Traverse genome directories
    for genome_dir in extraction_path.iterdir():
        if not genome_dir.is_dir():
            continue

        # Traverse locus directories
        for locus_dir in genome_dir.iterdir():
            if not locus_dir.is_dir():
                continue

            locus_name = locus_dir.name

            # Find CDS FASTA for status
            cds_files = list(locus_dir.glob("*_gene1_cds.fasta"))
            protein_files = list(locus_dir.glob("*_gene1_protein.fasta"))

            if cds_files and protein_files:
                # Parse status from CDS header
                with open(cds_files[0]) as f:
                    header = f.readline().strip()
                    status_match = re.search(r'status:(\w+)', header)
                    status = status_match.group(1) if status_match else 'unknown'

                # Parse length from protein header
                with open(protein_files[0]) as f:
                    header = f.readline().strip()
                    length_match = re.search(r'length:(\d+)aa', header)
                    length_aa = int(length_match.group(1)) if length_match else 0

                metadata[locus_name] = {
                    'length_aa': length_aa,
                    'status': status
                }

    return metadata

def create_gene_family_matrix(gene_family, loci_list, all_targets_df, synteny_blocks_df, species_map, phylo_order_map, seq_metadata, output_dir):
    """Create summary matrix for one gene family."""

    print(f"\n  Processing {gene_family}...")

    # Get all targets for this gene family
    family_targets = all_targets_df[all_targets_df['gene_family'] == gene_family]

    # Group targets by genome
    targets_by_genome = defaultdict(list)
    for _, target in family_targets.iterrows():
        targets_by_genome[target['genome']].append(target)

    # Read synteny percentages from Phase 8a locus matrices
    # Phase 8a already calculated correct synteny_pct based on flanking columns
    synteny_by_genome_locus = {}  # {genome: {locus: synteny_pct}}

    for locus in loci_list:
        # Load the locus matrix created by Phase 8a
        locus_matrix_file = output_dir / f"{locus}_genome_swissprot_matrix.tsv"

        if locus_matrix_file.exists():
            # Read synteny_pct from the locus matrix
            locus_matrix = pd.read_csv(locus_matrix_file, sep='\t')

            for _, row in locus_matrix.iterrows():
                genome = row['genome_id']
                synteny_pct = row.get('synteny_pct', 0.0)

                if genome not in synteny_by_genome_locus:
                    synteny_by_genome_locus[genome] = {}
                synteny_by_genome_locus[genome][locus] = synteny_pct
        else:
            print(f"    Warning: Locus matrix not found: {locus_matrix_file.name}")

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

        # Count targets by category
        genome_targets = targets_by_genome.get(genome, [])

        # Initialize all columns with "0%" (no synteny)
        for category in locus_categories:
            if category.endswith('_unplaceable'):
                row[category] = ""  # Unplaceable is special
            else:
                row[category] = "0%"

        # Fill in target counts
        category_counts = defaultdict(list)
        for target in genome_targets:
            assigned_to = target.get('assigned_to', target.get('description', 'unknown'))
            locus_name = target.get('locus_name', 'unknown')

            # Get metadata from extracted sequences
            meta = seq_metadata.get(locus_name, {})
            length = meta.get('length_aa', 0)
            status = meta.get('status', 'unknown')

            # Skip targets that failed extraction (Phase 6 deduplication removes these)
            if length == 0:
                continue

            # Status letters: I=intact, P=pseudogene, F=fragment
            if status == 'intact':
                status_letter = 'I'
            elif status == 'pseudogene':
                status_letter = 'P'
            elif status == 'fragment':
                status_letter = 'F'
            else:
                status_letter = '?'

            target_str = f"{length}{status_letter}"

            category_counts[assigned_to].append(target_str)

        # Format counts for each category with synteny percent
        for category, target_list in category_counts.items():
            if category in row:  # Only if it's a known category
                if target_list:
                    # Get synteny percent if available
                    if genome in synteny_by_genome_locus and category in synteny_by_genome_locus[genome]:
                        synteny_pct = synteny_by_genome_locus[genome][category]
                        row[category] = f"{synteny_pct}% [{'; '.join(target_list)}]"
                    else:
                        row[category] = f"[{'; '.join(target_list)}]"

        # Check for synteny blocks without targets
        # For loci where this genome has a synteny block but no extracted target
        if genome in synteny_by_genome_locus:
            for locus, synteny_pct in synteny_by_genome_locus[genome].items():
                if locus in row and row[locus] == "0%":
                    # Synteny block exists but no target found/extracted
                    row[locus] = f"{synteny_pct}% [empty]"

        # Add total count (only successfully extracted targets)
        row['total'] = sum(len(targets) for targets in category_counts.values())

        matrix_rows.append(row)

    # Create DataFrame and sort
    if matrix_rows:
        matrix_df = pd.DataFrame(matrix_rows)
        matrix_df = matrix_df.sort_values('phylo_order', ascending=True)
        return matrix_df

    return pd.DataFrame()

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate gene-type summary matrices",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('--locus-defs', required=True, type=Path,
                        help='Path to locus_definitions.tsv with gene_family column')
    parser.add_argument('--blocks', type=Path,
                        help='Path to synteny_blocks_filtered.tsv (optional)')
    parser.add_argument('--targets', type=Path,
                        help='Path to all_targets_classified.tsv (optional)')
    parser.add_argument('--species-map', required=True, type=Path,
                        help='Path to gca_to_species.tsv')
    parser.add_argument('--extracted-seqs', type=Path,
                        help='Path to 06_extracted_sequences directory (for length/status metadata)')
    parser.add_argument('--output-dir', required=True, type=Path,
                        help='Output directory for summary matrices')

    return parser.parse_args()

def main():
    """Main execution function."""
    args = parse_args()

    print("=" * 80)
    print("STEP 08b: GENERATE GENE-TYPE SUMMARY MATRICES")
    print("=" * 80)
    print(f"\nInput files:")
    print(f"  Locus definitions: {args.locus_defs}")
    print(f"  Synteny blocks: {args.blocks if args.blocks else 'None'}")
    print(f"  Targets: {args.targets if args.targets else 'None'}")
    print(f"  Species map: {args.species_map}")
    print(f"  Output directory: {args.output_dir}")

    # Create output directory
    args.output_dir.mkdir(exist_ok=True, parents=True)

    # Load all required data
    print("\n[1] Loading data...")

    # Locus definitions
    loci_df = pd.read_csv(args.locus_defs, sep='\t')
    print(f"  Loaded {len(loci_df)} locus definitions")

    # Load synteny blocks to track presence/absence
    if args.blocks and args.blocks.exists():
        synteny_blocks_df = pd.read_csv(args.blocks, sep='\t')
        print(f"  Loaded {len(synteny_blocks_df)} synteny blocks")
    else:
        synteny_blocks_df = pd.DataFrame()
        print("  No synteny blocks found")

    # All classified targets (including unplaceable)
    if args.targets and args.targets.exists():
        all_targets_df = pd.read_csv(args.targets, sep='\t')
        print(f"  Loaded {len(all_targets_df)} classified targets")
    else:
        all_targets_df = pd.DataFrame()
        print("  No classified targets found")

    # Species mapping
    species_map, phylo_order_map = load_species_and_phylo(args.species_map)
    print(f"  Loaded species mapping for {len(species_map)} genomes")

    # Load sequence metadata from extracted sequences
    if args.extracted_seqs and args.extracted_seqs.exists():
        seq_metadata = load_extracted_seq_metadata(args.extracted_seqs)
        print(f"  Loaded metadata for {len(seq_metadata)} extracted sequences")
    else:
        seq_metadata = {}
        print("  No extracted sequences found - lengths/status will be missing")

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
            species_map, phylo_order_map, seq_metadata, args.output_dir
        )

        if not matrix_df.empty:
            # Save matrix
            output_file = args.output_dir / f"{gene_family}_summary_matrix.tsv"
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

    # Overall summary
    print("\n" + "=" * 80)
    print("SUMMARY MATRICES COMPLETE")
    print("=" * 80)
    print(f"\nMatrices saved to: {args.output_dir}")

    if not all_targets_df.empty:
        print(f"\nOverall statistics:")
        print(f"  Total targets: {len(all_targets_df)}")
        print(f"  Syntenic: {len(all_targets_df[all_targets_df['placement'] == 'synteny'])}")
        print(f"  Unplaceable: {len(all_targets_df[all_targets_df['placement'] == 'unplaceable'])}")
        print(f"  Genomes with targets: {all_targets_df['genome'].nunique()}")

    print("\nPipeline complete!")

if __name__ == "__main__":
    main()