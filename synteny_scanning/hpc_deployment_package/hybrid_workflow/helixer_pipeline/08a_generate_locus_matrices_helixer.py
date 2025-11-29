#!/usr/bin/env python3
"""
Phase 8a (Helixer): Generate locus-specific matrices.

Creates one matrix per locus showing:
- All genomes (rows) sorted by phylogenetic order
- Flanking protein presence/absence (from Helixer Phase 4b)
- Target gene status (length)

This produces the SAME output format as the original tblastn Phase 8a.

Inputs:
- Phase 1 locus definitions
- Phase 4b flanking gene info
- Phase 5 Helixer classifications
- Phase 6 Helixer extracted sequences (for length metadata)
- Species/phylo mapping

Outputs:
- {locus_id}_genome_swissprot_matrix.tsv (one per locus)
"""

from __future__ import annotations

import argparse
import re
from collections import OrderedDict
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
    """
    Load metadata from Phase 6 extracted sequences.

    Returns dict keyed by (genome_id, locus_name):
      {(genome, locus): {'length_aa': int}}
    """
    metadata = {}

    if not extraction_dir or not extraction_dir.exists():
        return metadata

    for genome_dir in extraction_dir.iterdir():
        if not genome_dir.is_dir():
            continue

        genome_id = genome_dir.name

        for locus_dir in genome_dir.iterdir():
            if not locus_dir.is_dir() or locus_dir.name == 'UNPLACED':
                continue

            locus_name = locus_dir.name

            # Find protein FASTA
            protein_files = list(locus_dir.glob("*_protein.fasta"))
            if protein_files:
                with open(protein_files[0]) as f:
                    header = f.readline().strip()
                    length_match = re.search(r'length:(\d+)aa', header)
                    length_aa = int(length_match.group(1)) if length_match else 0

                metadata[(genome_id, locus_name)] = {'length_aa': length_aa}

    return metadata


def load_phase1_flanking_info(phase1_dir: Path) -> Dict[str, Dict]:
    """
    Load flanking gene information from Phase 1.

    Returns: {locus_id: {'upstream': [...], 'downstream': [...], 'gene_family': str}}
    """
    locus_info = {}

    locus_defs = phase1_dir / 'locus_definitions.tsv'
    if not locus_defs.exists():
        return {}

    df = pd.read_csv(locus_defs, sep='\t')

    for _, row in df.iterrows():
        locus_id = row['locus_id']
        gene_family = row.get('gene_family', 'unknown')

        # Load flanking proteins from Phase 1 file
        flanking_file = phase1_dir / f"{locus_id}_flanking.faa"
        if not flanking_file.exists():
            flanking_file = phase1_dir / f"{locus_id}_flanking_dedup.faa"

        upstream = []
        downstream = []

        if flanking_file.exists():
            with open(flanking_file) as f:
                for line in f:
                    if line.startswith('>'):
                        # Format: >XP_033209112.1|LOC117167954 U1 ...
                        parts = line[1:].strip().split()
                        if len(parts) >= 2:
                            protein_id = parts[0].split('|')[0]
                            position = parts[1] if len(parts) > 1 else ''

                            if position.startswith('U'):
                                upstream.append(protein_id)
                            elif position.startswith('D'):
                                downstream.append(protein_id)

        locus_info[locus_id] = {
            'upstream': upstream,
            'downstream': downstream,
            'gene_family': gene_family,
        }

    return locus_info


def create_locus_matrix(
    locus_id: str,
    locus_info: Dict,
    targets_df: pd.DataFrame,
    flanking_matches_df: pd.DataFrame,
    species_map: Dict[str, str],
    phylo_order_map: Dict[str, int],
    seq_metadata: Dict,
) -> pd.DataFrame:
    """Create matrix for one locus."""

    gene_family = locus_info.get('gene_family', locus_id)
    upstream_proteins = locus_info.get('upstream', [])
    downstream_proteins = locus_info.get('downstream', [])

    # Index targets by genome
    locus_targets = targets_df[
        (targets_df['placement'] == 'synteny') &
        (targets_df['assigned_to'] == locus_id)
    ]
    targets_by_genome = {row['genome']: row for _, row in locus_targets.iterrows()}

    # Build target_id -> genome mapping from targets_df
    target_to_genome = {}
    if not targets_df.empty:
        for _, row in targets_df.iterrows():
            target_to_genome[row['target_gene_id']] = row['genome']

    # Index flanking matches by genome with SwissProt annotations
    # flanking_by_genome[genome] = {bk_xp: annotation_string, ...}
    flanking_by_genome = {}
    if not flanking_matches_df.empty:
        for _, row in flanking_matches_df.iterrows():
            target_id = row['target_gene_id']
            bk_xp = row['bk_xp']
            genome = target_to_genome.get(target_id, 'unknown')

            if genome not in flanking_by_genome:
                flanking_by_genome[genome] = {}

            # Build annotation string with SwissProt info if available
            if 'swissprot_name' in row and row['swissprot_name'] and pd.notna(row['swissprot_name']):
                sp_name = str(row['swissprot_name'])
                sp_pident = row.get('swissprot_pident', 0)
                if sp_pident and sp_pident > 0:
                    annotation = f"{sp_name} ({int(sp_pident)}%)"
                else:
                    annotation = sp_name
            else:
                annotation = "match"

            flanking_by_genome[genome][bk_xp] = annotation

    # Build matrix rows
    matrix_rows = []

    for genome in sorted(species_map.keys()):
        row = OrderedDict()

        row['genome_id'] = genome
        row['species'] = species_map.get(genome, genome)
        row['phylo_order'] = phylo_order_map.get(genome, 999)

        # Check for synteny (flanking matches)
        genome_flanking = flanking_by_genome.get(genome, {})
        total_flanking = len(upstream_proteins) + len(downstream_proteins)
        all_flanking_proteins = set(upstream_proteins + downstream_proteins)
        n_matches = len(set(genome_flanking.keys()) & all_flanking_proteins)

        if total_flanking > 0:
            synteny_pct = round((n_matches / total_flanking) * 100, 1)
        else:
            synteny_pct = 0

        row['synteny_pct'] = synteny_pct
        row['num_proteins_found'] = n_matches

        # Add scaffold info if target exists
        if genome in targets_by_genome:
            target = targets_by_genome[genome]
            row['scaffold'] = target.get('scaffold', '')
            row['strand'] = target.get('strand', '')
            row['start'] = target.get('start', '')
            row['end'] = target.get('end', '')
        else:
            row['scaffold'] = ''
            row['strand'] = ''
            row['start'] = ''
            row['end'] = ''

        # Upstream proteins (U14, U13, ..., U1)
        for i, protein_id in enumerate(reversed(upstream_proteins), 1):
            col_name = f"U{len(upstream_proteins) - i + 1}_{protein_id[:15]}"
            if protein_id in genome_flanking:
                row[col_name] = genome_flanking[protein_id]  # SwissProt annotation or "match"
            else:
                row[col_name] = ""

        # TARGET column
        if genome in targets_by_genome:
            meta = seq_metadata.get((genome, locus_id), {})
            length = meta.get('length_aa', 0)
            if length > 0:
                row['TARGET'] = f"{gene_family} [{length}]"
            else:
                row['TARGET'] = f"{gene_family} [hit]"
        elif synteny_pct > 0:
            row['TARGET'] = f"{gene_family} [empty]"
        else:
            row['TARGET'] = f"{gene_family} [0%]"

        # Downstream proteins (D1, D2, ...)
        for i, protein_id in enumerate(downstream_proteins, 1):
            col_name = f"D{i}_{protein_id[:15]}"
            if protein_id in genome_flanking:
                row[col_name] = genome_flanking[protein_id]  # SwissProt annotation or "match"
            else:
                row[col_name] = ""

        matrix_rows.append(row)

    # Create DataFrame and sort
    if matrix_rows:
        matrix_df = pd.DataFrame(matrix_rows)
        matrix_df = matrix_df.sort_values('phylo_order', ascending=True)
        return matrix_df

    return pd.DataFrame()


def main():
    parser = argparse.ArgumentParser(
        description="Phase 8a (Helixer): Generate locus-specific matrices"
    )
    parser.add_argument("--family", required=True, help="Gene family name")
    parser.add_argument("--phase1-dir", type=Path, required=True,
                        help="Phase 1 directory with locus definitions")
    parser.add_argument("--phase5-dir", type=Path, required=True,
                        help="Phase 5 Helixer output directory")
    parser.add_argument("--phase6-dir", type=Path, required=True,
                        help="Phase 6 Helixer extracted sequences")
    parser.add_argument("--phase7-dir", type=Path, default=None,
                        help="Phase 7 SwissProt annotations (optional)")
    parser.add_argument("--species-map", type=Path, required=True,
                        help="Path to gca_to_species_order.tsv")
    parser.add_argument("--output-dir", type=Path, required=True,
                        help="Output directory")

    args = parser.parse_args()

    print("=" * 80)
    print("PHASE 8a (HELIXER): GENERATE LOCUS-SPECIFIC MATRICES")
    print("=" * 80)
    print()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    print("[1] Loading data...")

    # Species mapping
    species_map, phylo_order_map = load_species_and_phylo(args.species_map)
    print(f"  Species: {len(species_map)} genomes")

    # Phase 1 flanking info
    locus_info = load_phase1_flanking_info(args.phase1_dir)
    print(f"  Loci: {len(locus_info)}")

    # Phase 5 classifications
    targets_file = args.phase5_dir / "all_targets_classified.tsv"
    if targets_file.exists():
        targets_df = pd.read_csv(targets_file, sep='\t')
        print(f"  Targets: {len(targets_df)}")
    else:
        targets_df = pd.DataFrame()
        print("  No classified targets found")

    # Phase 5/7 flanking matches - prefer annotated version from Phase 7
    flanking_df = pd.DataFrame()
    swissprot_available = False

    if args.phase7_dir and (args.phase7_dir / "flanking_matches_annotated.tsv").exists():
        flanking_file = args.phase7_dir / "flanking_matches_annotated.tsv"
        flanking_df = pd.read_csv(flanking_file, sep='\t')
        swissprot_available = 'swissprot_name' in flanking_df.columns
        print(f"  Flanking matches (with SwissProt): {len(flanking_df)}")
    elif (args.phase5_dir / "flanking_matches.tsv").exists():
        flanking_file = args.phase5_dir / "flanking_matches.tsv"
        flanking_df = pd.read_csv(flanking_file, sep='\t')
        print(f"  Flanking matches: {len(flanking_df)}")
    else:
        print("  No flanking matches found")

    # Phase 6 extracted sequences metadata
    seq_metadata = load_extracted_seq_metadata(args.phase6_dir)
    print(f"  Extracted sequences: {len(seq_metadata)}")

    # Generate matrices
    print("\n[2] Generating matrices...")

    for locus_id, info in locus_info.items():
        matrix_df = create_locus_matrix(
            locus_id, info, targets_df, flanking_df,
            species_map, phylo_order_map, seq_metadata
        )

        if not matrix_df.empty:
            output_file = args.output_dir / f"{locus_id}_genome_swissprot_matrix.tsv"
            matrix_df.to_csv(output_file, sep='\t', index=False)

            n_with_synteny = len(matrix_df[matrix_df['synteny_pct'] > 0])
            print(f"  {locus_id}: {len(matrix_df)} rows, {n_with_synteny} with synteny")

    print("\n" + "=" * 80)
    print("PHASE 8a (HELIXER) COMPLETE")
    print("=" * 80)
    print(f"\nMatrices saved to: {args.output_dir}")

    return 0


if __name__ == "__main__":
    exit(main())
