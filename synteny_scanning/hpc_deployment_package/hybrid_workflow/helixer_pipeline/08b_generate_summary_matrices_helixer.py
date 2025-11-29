#!/usr/bin/env python3
"""
Phase 8b (Helixer): Generate gene-family summary matrices.

Creates summary matrix by gene family showing:
- All genomes (rows) sorted by phylogenetic order
- All loci of that gene family (columns)
- Includes both syntenic and unplaceable targets

This produces the SAME output format as the original tblastn Phase 8b.

Inputs:
- Phase 1 locus definitions
- Phase 5 Helixer classifications
- Phase 6 Helixer extracted sequences
- Phase 8a locus matrices (for synteny percentages)
- Species/phylo mapping

Outputs:
- {gene_family}_summary_matrix.tsv
"""

from __future__ import annotations

import argparse
import re
from collections import OrderedDict, defaultdict
from pathlib import Path
from typing import Dict, Tuple

import pandas as pd


def load_species_and_phylo(species_map_file: Path) -> Tuple[Dict[str, str], Dict[str, int]]:
    """Load species mapping and phylogenetic order."""
    species_map = {}
    phylo_order_map = {}

    with open(species_map_file) as f:
        header = f.readline()
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


def load_extracted_seq_metadata(extraction_dir: Path) -> Dict[Tuple[str, str], Dict]:
    """Load metadata from Phase 6 extracted sequences."""
    metadata = {}

    if not extraction_dir or not extraction_dir.exists():
        return metadata

    for genome_dir in extraction_dir.iterdir():
        if not genome_dir.is_dir():
            continue

        genome_id = genome_dir.name

        for locus_dir in genome_dir.iterdir():
            if not locus_dir.is_dir():
                continue

            # Handle UNPLACED
            if locus_dir.name == 'UNPLACED':
                for unplaced_dir in locus_dir.iterdir():
                    if unplaced_dir.is_dir():
                        protein_files = list(unplaced_dir.glob("*_protein.fasta"))
                        if protein_files:
                            with open(protein_files[0]) as f:
                                header = f.readline().strip()
                                length_match = re.search(r'length:(\d+)aa', header)
                                length_aa = int(length_match.group(1)) if length_match else 0

                            metadata[(genome_id, unplaced_dir.name)] = {'length_aa': length_aa}
                continue

            locus_name = locus_dir.name

            protein_files = list(locus_dir.glob("*_protein.fasta"))
            if protein_files:
                with open(protein_files[0]) as f:
                    header = f.readline().strip()
                    length_match = re.search(r'length:(\d+)aa', header)
                    length_aa = int(length_match.group(1)) if length_match else 0

                metadata[(genome_id, locus_name)] = {'length_aa': length_aa}

    return metadata


def create_summary_matrix(
    gene_family: str,
    loci_list: list,
    targets_df: pd.DataFrame,
    species_map: Dict[str, str],
    phylo_order_map: Dict[str, int],
    seq_metadata: Dict,
    locus_matrices_dir: Path,
) -> pd.DataFrame:
    """Create summary matrix for gene family."""

    # Load synteny percentages from Phase 8a locus matrices
    synteny_by_genome_locus = {}
    for locus in loci_list:
        matrix_file = locus_matrices_dir / f"{locus}_genome_swissprot_matrix.tsv"
        if matrix_file.exists():
            matrix_df = pd.read_csv(matrix_file, sep='\t')
            for _, row in matrix_df.iterrows():
                genome = row['genome_id']
                synteny_pct = row.get('synteny_pct', 0)

                if genome not in synteny_by_genome_locus:
                    synteny_by_genome_locus[genome] = {}
                synteny_by_genome_locus[genome][locus] = synteny_pct

    # Group targets by genome
    targets_by_genome = defaultdict(list)
    if not targets_df.empty:
        for _, target in targets_df.iterrows():
            targets_by_genome[target['genome']].append(target)

    # Define columns: loci + unplaceable
    unplaceable_col = f"{gene_family}_unplaceable"
    columns = loci_list + [unplaceable_col]

    # Build matrix rows
    matrix_rows = []

    for genome in sorted(species_map.keys()):
        row = OrderedDict()

        row['genome_id'] = genome
        row['species'] = species_map.get(genome, genome)
        row['phylo_order'] = phylo_order_map.get(genome, 999)

        # Initialize all columns
        for col in columns:
            row[col] = ""

        # Get targets for this genome
        genome_targets = targets_by_genome.get(genome, [])

        # Track counts per category
        category_counts = defaultdict(list)

        for target in genome_targets:
            placement = target.get('placement', '')
            assigned_to = target.get('assigned_to', '')

            if placement == 'synteny' and assigned_to in loci_list:
                category = assigned_to
            else:
                category = unplaceable_col

            # Get length metadata
            if placement == 'synteny':
                lookup_key = (genome, assigned_to)
            else:
                # Unplaceable: use unique tag
                scaffold = str(target.get('scaffold', ''))
                start = str(target.get('start', ''))
                end = str(target.get('end', ''))
                unique_tag = f"{scaffold}_{start}_{end}"
                lookup_key = (genome, unique_tag)

            meta = seq_metadata.get(lookup_key, {})
            length = meta.get('length_aa', 0)

            if length > 0:
                category_counts[category].append(str(length))
            else:
                category_counts[category].append("hit")

        # Format cells
        for category, target_list in category_counts.items():
            if category in row:
                if category in loci_list:
                    # Syntenic: include synteny percentage
                    synteny_pct = synteny_by_genome_locus.get(genome, {}).get(category, 0)
                    row[category] = f"{synteny_pct}% [{'; '.join(target_list)}]"
                else:
                    # Unplaceable
                    row[category] = f"[{'; '.join(target_list)}]"

        # Fill in loci with synteny but no targets
        if genome in synteny_by_genome_locus:
            for locus, synteny_pct in synteny_by_genome_locus[genome].items():
                if locus in row and row[locus] == "":
                    if synteny_pct > 0:
                        row[locus] = f"{synteny_pct}% [empty]"
                    else:
                        row[locus] = "[not found]"

        # Mark remaining empty loci
        for locus in loci_list:
            if locus in row and row[locus] == "":
                row[locus] = "[not found]"

        # Total count
        row['total'] = sum(len(targets) for targets in category_counts.values())

        matrix_rows.append(row)

    # Create DataFrame
    if matrix_rows:
        matrix_df = pd.DataFrame(matrix_rows)
        matrix_df = matrix_df.sort_values('phylo_order', ascending=True)
        return matrix_df

    return pd.DataFrame()


def main():
    parser = argparse.ArgumentParser(
        description="Phase 8b (Helixer): Generate gene-family summary matrices"
    )
    parser.add_argument("--family", required=True, help="Gene family name")
    parser.add_argument("--phase1-dir", type=Path, required=True,
                        help="Phase 1 directory with locus definitions")
    parser.add_argument("--phase5-dir", type=Path, required=True,
                        help="Phase 5 Helixer output directory")
    parser.add_argument("--phase6-dir", type=Path, required=True,
                        help="Phase 6 Helixer extracted sequences")
    parser.add_argument("--phase8a-dir", type=Path, required=True,
                        help="Phase 8a locus matrices directory")
    parser.add_argument("--species-map", type=Path, required=True,
                        help="Path to gca_to_species_order.tsv")
    parser.add_argument("--output-dir", type=Path, required=True,
                        help="Output directory")

    args = parser.parse_args()

    print("=" * 80)
    print("PHASE 8b (HELIXER): GENERATE SUMMARY MATRICES")
    print("=" * 80)
    print()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    print("[1] Loading data...")

    # Species mapping
    species_map, phylo_order_map = load_species_and_phylo(args.species_map)
    print(f"  Species: {len(species_map)} genomes")

    # Locus definitions
    locus_defs_file = args.phase1_dir / "locus_definitions.tsv"
    if locus_defs_file.exists():
        locus_defs = pd.read_csv(locus_defs_file, sep='\t')
        gene_families = locus_defs['gene_family'].unique()
        loci_by_family = {
            gf: locus_defs[locus_defs['gene_family'] == gf]['locus_id'].tolist()
            for gf in gene_families
        }
        print(f"  Gene families: {len(gene_families)}")
    else:
        print("  [ERROR] No locus definitions found")
        return 1

    # Phase 5 targets
    targets_file = args.phase5_dir / "all_targets_classified.tsv"
    if targets_file.exists():
        targets_df = pd.read_csv(targets_file, sep='\t')
        print(f"  Targets: {len(targets_df)}")
    else:
        targets_df = pd.DataFrame()
        print("  No targets found")

    # Phase 6 metadata
    seq_metadata = load_extracted_seq_metadata(args.phase6_dir)
    print(f"  Extracted sequences: {len(seq_metadata)}")

    # Generate summary matrices
    print("\n[2] Generating summary matrices...")

    for gene_family, loci_list in loci_by_family.items():
        matrix_df = create_summary_matrix(
            gene_family, loci_list, targets_df,
            species_map, phylo_order_map, seq_metadata,
            args.phase8a_dir
        )

        if not matrix_df.empty:
            output_file = args.output_dir / f"{gene_family}_summary_matrix.tsv"
            matrix_df.to_csv(output_file, sep='\t', index=False)

            n_with_targets = len(matrix_df[matrix_df['total'] > 0])
            print(f"  {gene_family}: {len(matrix_df)} rows, {n_with_targets} with targets")

    # Overall statistics
    if not targets_df.empty:
        print("\n" + "=" * 80)
        print("SUMMARY STATISTICS")
        print("=" * 80)
        print(f"\nTotal targets: {len(targets_df)}")
        n_syntenic = len(targets_df[targets_df['placement'] == 'synteny'])
        n_unplaceable = len(targets_df[targets_df['placement'] == 'unplaceable'])
        print(f"  Syntenic: {n_syntenic} ({100*n_syntenic/len(targets_df):.1f}%)")
        print(f"  Unplaceable: {n_unplaceable} ({100*n_unplaceable/len(targets_df):.1f}%)")
        print(f"  Genomes with targets: {targets_df['genome'].nunique()}")

    print("\n" + "=" * 80)
    print("PHASE 8b (HELIXER) COMPLETE")
    print("=" * 80)
    print(f"\nSummary matrices saved to: {args.output_dir}")

    return 0


if __name__ == "__main__":
    exit(main())
