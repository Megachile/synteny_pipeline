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
from typing import Dict, List, Tuple

import pandas as pd

from constants import (
    normalize_columns,
    validate_required_columns,
    COL_PLACEMENT, COL_ASSIGNED_TO,
    PLACEMENT_SYNTENY,
    LOWQUALITY_BITSCORE_THRESHOLD,
    FILE_PHASE2B_SYNTENY_BLOCKS, FILE_EMPTY_SYNTENY_BLOCKS_LEGACY,
)

# Additional threshold not in constants (specific to unplaceable classification)
MIN_FLANKING_FOR_DISTINCT = 3  # Need at least this many flanking genes to be "distinct"

# Reference genome definitions (hardcoded since they're not in gca_to_species)
# These are added inline with Helixer genomes at their phylogenetic positions
REFERENCE_GENOMES = {
    'Belonocnema_kinseyi_GCF': {
        'species': 'Belonocnema kinseyi (reference)',
        'phylo_order': 17,
        'genome_type': 'reference',
        'validation_key': 'BK',  # Key used in novel_locus_validation.tsv
    },
    'GCA_019393585.1': {
        'species': 'Leptopilina boulardi (outgroup)',
        'phylo_order': 6,
        'genome_type': 'reference',
        'validation_key': 'LB',  # Key used in novel_locus_validation.tsv
    },
}


def classify_unplaceable(target: dict) -> str:
    """
    Classify an unplaceable target into subcategories.

    Returns one of:
    - 'fragment': On a short scaffold with no/few flanking genes
    - 'lowquality': Hit bitscore below threshold
    - 'distinct': High-quality hit with flanking genes but no locus match
    """
    bitscore = target.get('best_bitscore', 0)
    n_flanking = target.get('n_flanking_genes', 0)
    reason = target.get('unplaceable_reason', '')

    # Fragment: no flanking genes available (short scaffold)
    if n_flanking == 0 or reason == 'no_flanking_genes':
        return 'fragment'

    # Low quality: weak hit regardless of flanking
    if bitscore < LOWQUALITY_BITSCORE_THRESHOLD:
        return 'lowquality'

    # Distinct: good hit with flanking but doesn't match known loci
    # This includes 'no_locus_matches' and 'no_bk_hits' reasons
    return 'distinct'


def load_species_and_phylo(species_map_file: Path) -> Tuple[Dict[str, str], Dict[str, int]]:
    """Load species mapping and phylogenetic order, excluding genomes marked as excluded."""
    species_map = {}
    phylo_order_map = {}

    with open(species_map_file) as f:
        header = f.readline().strip().split('\t')
        # Find status column index
        status_idx = header.index('status') if 'status' in header else -1

        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                # Skip excluded genomes
                if status_idx >= 0 and len(parts) > status_idx:
                    if parts[status_idx] == 'excluded':
                        continue

                genome_id = parts[0]
                species_map[genome_id] = parts[1]
                try:
                    phylo_order_map[genome_id] = int(parts[3])
                except (ValueError, IndexError):
                    phylo_order_map[genome_id] = 999

    return species_map, phylo_order_map


def load_novel_loci_clusters(phase5b_dir: Path) -> pd.DataFrame:
    """Load novel loci clustering results from Phase 5b."""
    if not phase5b_dir:
        return pd.DataFrame()

    cluster_file = phase5b_dir / "novel_loci_clustered.tsv"
    if not cluster_file.exists():
        return pd.DataFrame()

    try:
        df = pd.read_csv(cluster_file, sep='\t')
        return df
    except Exception as e:
        print(f"  [WARN] Could not load novel loci clusters: {e}")
        return pd.DataFrame()


def load_reference_validation(phase5b_dir: Path) -> Dict[str, Dict[str, Dict]]:
    """
    Load reference genome validation results from Phase 5b.

    Returns: {validation_key: {candidate_id: {n_hits, n_unique_flanking, synteny_score}}}
             e.g., {'BK': {'NOVEL_xxx': {...}, ...}, 'LB': {...}}
    """
    if not phase5b_dir:
        return {}

    validation_file = phase5b_dir / "novel_locus_validation.tsv"
    if not validation_file.exists():
        return {}

    try:
        df = pd.read_csv(validation_file, sep='\t')
        # Filter to reference genomes only
        ref_df = df[df['genome_type'] == 'reference']

        results = defaultdict(dict)
        for _, row in ref_df.iterrows():
            target_genome = row['target_genome']  # 'BK' or 'LB'
            candidate_id = row['candidate_id']
            results[target_genome][candidate_id] = {
                'n_hits': row.get('n_hits', 0),
                'n_unique_flanking': row.get('n_unique_flanking', 0),
                'synteny_score': row.get('synteny_score', 0.0),
            }

        return dict(results)
    except Exception as e:
        print(f"  [WARN] Could not load reference validation: {e}")
        return {}


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
    phase2b_blocks_df: pd.DataFrame,
    novel_loci_df: pd.DataFrame = None,
    ref_validation: Dict[str, Dict[str, Dict]] = None,
    novel_loci_id_to_name: Dict[str, str] = None,
) -> pd.DataFrame:
    """
    Create summary matrix for gene family.

    Shows all synteny blocks per cell with format:
    - Concordant (both methods): "75%* [318]"
    - Empty block (flanking only): "35% [empty]"
    - Target without flanking block: "[318]"
    """

    # Index Phase 2b blocks by (genome, locus)
    # blocks_by_genome_locus[genome][locus] = [{synteny_score, has_target_overlap, overlapping_target_id}, ...]
    blocks_by_genome_locus = defaultdict(lambda: defaultdict(list))
    if not phase2b_blocks_df.empty:
        for _, block in phase2b_blocks_df.iterrows():
            genome = block['genome']
            locus = block['locus_id']
            blocks_by_genome_locus[genome][locus].append({
                'synteny_score': block['synteny_score'],
                'has_target_overlap': block.get('has_target_overlap', False),
                'overlapping_target_id': block.get('overlapping_target_id', ''),
                'scaffold': block.get('scaffold', ''),
                'start': block.get('start', 0),
                'end': block.get('end', 0),
            })

    # Index targets by genome and target_id
    targets_by_genome = defaultdict(list)
    targets_by_id = {}
    if not targets_df.empty:
        for _, target in targets_df.iterrows():
            targets_by_genome[target['genome']].append(target)
            targets_by_id[target['target_gene_id']] = target

    # Define columns: loci + unplaceable + novel loci
    unplaceable_col = f"{gene_family}_unplaceable"
    columns = loci_list + [unplaceable_col]

    # Add novel loci columns if present
    novel_loci_list = []  # List of display names (locus_name)
    novel_loci_by_genome = {}  # {genome: {locus_name: 'member'|'empty'}}
    if novel_loci_df is not None and not novel_loci_df.empty:
        for _, row_nl in novel_loci_df.iterrows():
            # Use locus_name for column headers (e.g., "Synergus_coniferae_chr2")
            # Fall back to locus_id if locus_name not present
            locus_name = row_nl.get('locus_name', row_nl['locus_id'])
            novel_loci_list.append(locus_name)

            # Parse member genomes
            member_genomes = row_nl.get('member_genomes', '')
            if isinstance(member_genomes, str) and member_genomes:
                for g in member_genomes.split(';'):
                    g = g.strip()
                    if g:
                        if g not in novel_loci_by_genome:
                            novel_loci_by_genome[g] = {}
                        novel_loci_by_genome[g][locus_name] = 'member'

            # Parse empty genomes (flanking context but no target)
            empty_genomes = row_nl.get('empty_genomes', '')
            if isinstance(empty_genomes, str) and empty_genomes:
                for g in empty_genomes.split(';'):
                    g = g.strip()
                    if g:
                        if g not in novel_loci_by_genome:
                            novel_loci_by_genome[g] = {}
                        novel_loci_by_genome[g][locus_name] = 'empty'

        columns = columns + novel_loci_list

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

        genome_targets = targets_by_genome.get(genome, [])

        # Track which targets were matched to blocks
        matched_target_ids = set()

        # Process each locus
        for locus in loci_list:
            cell_entries = []

            # Get Phase 2b blocks for this genome/locus
            locus_blocks = blocks_by_genome_locus[genome].get(locus, [])

            for block in locus_blocks:
                synteny_pct = round(block['synteny_score'] * 100, 1)

                if block['has_target_overlap'] and block['overlapping_target_id']:
                    # Concordant block - both methods found it
                    target_id = block['overlapping_target_id']
                    matched_target_ids.add(target_id)

                    # Get AA length from metadata
                    meta = seq_metadata.get(target_id, {})
                    length = meta.get('length_aa', 0)

                    if length > 0:
                        cell_entries.append(f"{synteny_pct}%* [{length}]")
                    else:
                        cell_entries.append(f"{synteny_pct}%* [hit]")
                else:
                    # Empty block - flanking genes found, no target
                    cell_entries.append(f"{synteny_pct}% [empty]")

            # Check for targets assigned to this locus that weren't matched to blocks
            for target in genome_targets:
                if target.get('assigned_to') == locus and target['target_gene_id'] not in matched_target_ids:
                    # Target found by target-first but no matching Phase 2b block
                    target_id = target['target_gene_id']
                    meta = seq_metadata.get(target_id, {})
                    length = meta.get('length_aa', 0)

                    if length > 0:
                        cell_entries.append(f"[{length}]")
                    else:
                        cell_entries.append("[hit]")

            # Format cell
            if cell_entries:
                row[locus] = ", ".join(cell_entries)
            else:
                row[locus] = "-"

        # Handle unplaceable targets with subcategory counts
        unplaceable_entries = []
        n_fragment = 0
        n_lowquality = 0
        n_distinct = 0

        for target in genome_targets:
            if target.get(COL_PLACEMENT) != PLACEMENT_SYNTENY or target.get(COL_ASSIGNED_TO) not in loci_list:
                target_id = target['target_gene_id']
                if target_id not in matched_target_ids:
                    # Classify this unplaceable
                    target_dict = target.to_dict() if hasattr(target, 'to_dict') else dict(target)
                    category = classify_unplaceable(target_dict)

                    if category == 'fragment':
                        n_fragment += 1
                    elif category == 'lowquality':
                        n_lowquality += 1
                    else:  # distinct
                        n_distinct += 1

                    meta = seq_metadata.get(target_id, {})
                    length = meta.get('length_aa', 0)
                    if length > 0:
                        unplaceable_entries.append(f"[{length}]")
                    else:
                        unplaceable_entries.append("[hit]")

        if unplaceable_entries:
            row[unplaceable_col] = ", ".join(unplaceable_entries)
        else:
            row[unplaceable_col] = "-"

        # Add unplaceable count for easier analysis
        row['unplaceable_count'] = len(unplaceable_entries)

        # Add subcategory counts
        row['n_fragment'] = n_fragment
        row['n_lowquality'] = n_lowquality
        row['n_distinct'] = n_distinct

        # Total target count
        row['total'] = len(genome_targets)

        # Syntenic target count (total - unplaceable)
        row['syntenic_count'] = len(genome_targets) - len(unplaceable_entries)

        # Fill novel loci columns
        n_novel_hits = 0
        n_novel_empty = 0
        for novel_locus in novel_loci_list:
            genome_novel_status = novel_loci_by_genome.get(genome, {})
            status = genome_novel_status.get(novel_locus)
            if status == 'member':
                row[novel_locus] = "[hit]"
                n_novel_hits += 1
            elif status == 'empty':
                row[novel_locus] = "[empty]"
                n_novel_empty += 1
            else:
                row[novel_locus] = "-"

        # Add novel loci counts
        if novel_loci_list:
            row['n_novel_hits'] = n_novel_hits
            row['n_novel_empty'] = n_novel_empty

        matrix_rows.append(row)

    # Add reference genome validation data to novel loci columns
    # If ref genome already exists as a row (in species_map), update it; otherwise add new row
    if ref_validation and novel_loci_list and novel_loci_id_to_name:
        # Build lookup: validation_key -> ref_data for quick access
        for genome_id, ref_info in REFERENCE_GENOMES.items():
            validation_key = ref_info['validation_key']  # 'BK' or 'LB'

            # Check if we have validation data for this reference
            if validation_key not in ref_validation:
                continue

            ref_data = ref_validation[validation_key]

            # Check if this genome already exists in species_map (already has a row)
            if genome_id in species_map:
                # Find and update the existing row
                for existing_row in matrix_rows:
                    if existing_row['genome_id'] == genome_id:
                        # Update novel loci columns with reference validation data
                        n_novel_hits = 0
                        n_novel_empty = 0

                        for novel_locus in novel_loci_list:
                            matching_candidate = novel_loci_id_to_name.get(novel_locus)
                            if matching_candidate and matching_candidate in ref_data:
                                val_info = ref_data[matching_candidate]
                                synteny_score = val_info.get('synteny_score', 0.0)
                                n_hits = val_info.get('n_hits', 0)

                                if synteny_score > 0:
                                    pct = round(synteny_score * 100, 1)
                                    existing_row[novel_locus] = f"{pct}% [flanking]"
                                    n_novel_hits += 1
                                elif n_hits > 0:
                                    existing_row[novel_locus] = f"[{n_hits} hits]"
                                    n_novel_empty += 1

                        # Update counts
                        existing_row['n_novel_hits'] = n_novel_hits
                        existing_row['n_novel_empty'] = n_novel_empty
                        break
            else:
                # Add new row for reference genome not in species_map
                row = OrderedDict()
                row['genome_id'] = genome_id
                row['species'] = ref_info['species']
                row['phylo_order'] = ref_info['phylo_order']

                # Initialize all columns
                for col in columns:
                    row[col] = ""

                # Reference genomes don't have target hits through our pipeline
                for locus in loci_list:
                    row[locus] = "-"

                row[unplaceable_col] = "-"
                row['unplaceable_count'] = 0
                row['n_fragment'] = 0
                row['n_lowquality'] = 0
                row['n_distinct'] = 0
                row['total'] = 0
                row['syntenic_count'] = 0

                # Fill novel loci columns from validation results
                n_novel_hits = 0
                n_novel_empty = 0

                for novel_locus in novel_loci_list:
                    matching_candidate = novel_loci_id_to_name.get(novel_locus)
                    if matching_candidate and matching_candidate in ref_data:
                        val_info = ref_data[matching_candidate]
                        synteny_score = val_info.get('synteny_score', 0.0)
                        n_hits = val_info.get('n_hits', 0)

                        if synteny_score > 0:
                            pct = round(synteny_score * 100, 1)
                            row[novel_locus] = f"{pct}% [flanking]"
                            n_novel_hits += 1
                        elif n_hits > 0:
                            row[novel_locus] = f"[{n_hits} hits]"
                            n_novel_empty += 1
                        else:
                            row[novel_locus] = "-"
                    else:
                        row[novel_locus] = "-"

                if novel_loci_list:
                    row['n_novel_hits'] = n_novel_hits
                    row['n_novel_empty'] = n_novel_empty

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
    parser.add_argument("--phase2b-dir", type=Path, default=None,
                        help="Phase 2b synteny blocks directory (optional)")
    parser.add_argument("--phase5-dir", type=Path, required=True,
                        help="Phase 5 Helixer output directory")
    parser.add_argument("--phase6-dir", type=Path, required=True,
                        help="Phase 6 Helixer extracted sequences")
    parser.add_argument("--species-map", type=Path, required=True,
                        help="Path to gca_to_species_order.tsv")
    parser.add_argument("--output-dir", type=Path, required=True,
                        help="Output directory")
    parser.add_argument("--phase5b-dir", type=Path, default=None,
                        help="Phase 5b novel loci validation directory (optional)")

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

    # Phase 2b synteny blocks (optional)
    phase2b_blocks_df = pd.DataFrame()
    if args.phase2b_dir and args.phase2b_dir.exists():
        phase2b_file = args.phase2b_dir / FILE_PHASE2B_SYNTENY_BLOCKS
        if not phase2b_file.exists():
            phase2b_file = args.phase2b_dir / FILE_EMPTY_SYNTENY_BLOCKS_LEGACY
        if phase2b_file.exists():
            phase2b_blocks_df = pd.read_csv(phase2b_file, sep='\t')
            if 'has_target_overlap' in phase2b_blocks_df.columns:
                n_concordant = phase2b_blocks_df['has_target_overlap'].sum()
                n_empty = len(phase2b_blocks_df) - n_concordant
            else:
                n_concordant = 0
                n_empty = len(phase2b_blocks_df)
            print(f"  Phase 2b blocks: {len(phase2b_blocks_df)} ({n_concordant} concordant, {n_empty} empty)")
        else:
            print("  Phase 2b directory exists but no blocks file found")
    else:
        print("  Phase 2b not specified (target-only mode)")

    # Phase 5 targets - load and combine syntenic + unplaceable
    targets_df = pd.DataFrame()
    syntenic_file = args.phase5_dir / "syntenic_targets.tsv"
    unplaceable_file = args.phase5_dir / "unplaceable_targets.tsv"

    dfs_to_concat = []
    if syntenic_file.exists():
        try:
            syn_df = pd.read_csv(syntenic_file, sep='\t')
            if not syn_df.empty:
                dfs_to_concat.append(syn_df)
        except pd.errors.EmptyDataError:
            print("  [NOTE] Syntenic targets file is empty")
    if unplaceable_file.exists():
        try:
            unpl_df = pd.read_csv(unplaceable_file, sep='\t')
            if not unpl_df.empty:
                dfs_to_concat.append(unpl_df)
        except pd.errors.EmptyDataError:
            print("  [NOTE] Unplaceable targets file is empty")

    if dfs_to_concat:
        targets_df = pd.concat(dfs_to_concat, ignore_index=True)
        normalize_columns(targets_df)  # Standardize column names
        # Validate required columns are present
        required_cols = ['genome', 'target_gene_id', COL_PLACEMENT]
        try:
            validate_required_columns(targets_df, required_cols, "targets (syntenic + unplaceable)")
        except ValueError as e:
            print(f"  [WARNING] {e}")
        print(f"  Targets: {len(targets_df)}")
    else:
        print("  No targets found")

    # Phase 6 metadata - load from length_qc.tsv (keyed by target_id)
    seq_metadata = {}
    length_qc_file = args.phase6_dir / "length_qc.tsv"
    if length_qc_file.exists():
        length_df = pd.read_csv(length_qc_file, sep='\t')
        for _, row in length_df.iterrows():
            target_id = row['target_id']
            length_aa = int(row['extracted_length_aa']) if pd.notna(row['extracted_length_aa']) else 0
            seq_metadata[target_id] = {'length_aa': length_aa}
        print(f"  Extracted sequences: {len(seq_metadata)}")
    else:
        # Fallback to all_extracted_genes.faa
        all_genes_file = args.phase6_dir / "all_extracted_genes.faa"
        if all_genes_file.exists():
            current_id = None
            current_seq = []
            with open(all_genes_file) as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        if current_id and current_seq:
                            seq_metadata[current_id] = {'length_aa': len(''.join(current_seq))}
                        current_id = line[1:].split()[0]
                        current_seq = []
                    else:
                        current_seq.append(line)
                if current_id and current_seq:
                    seq_metadata[current_id] = {'length_aa': len(''.join(current_seq))}
            print(f"  Extracted sequences: {len(seq_metadata)}")

    # Phase 5b novel loci (optional)
    novel_loci_df = pd.DataFrame()
    ref_validation = {}
    novel_loci_id_to_name = {}

    if args.phase5b_dir:
        novel_loci_df = load_novel_loci_clusters(args.phase5b_dir)
        if not novel_loci_df.empty:
            print(f"  Novel loci: {len(novel_loci_df)} clusters")

            # Build locus_name -> representative mapping for reference genome lookup
            # The representative is the candidate_id used in validation file
            for _, row in novel_loci_df.iterrows():
                locus_name = row.get('locus_name', row['locus_id'])
                representative = row.get('representative', '')
                if representative:
                    novel_loci_id_to_name[locus_name] = representative
        else:
            print("  No novel loci clusters found")

        # Load reference genome validation results (BK, LB)
        ref_validation = load_reference_validation(args.phase5b_dir)
        if ref_validation:
            print(f"  Reference validation: {list(ref_validation.keys())}")
        else:
            print("  No reference validation data found")
    else:
        print("  Phase 5b not specified (no novel loci)")

    # Generate summary matrices
    print("\n[2] Generating summary matrices...")

    for gene_family, loci_list in loci_by_family.items():
        matrix_df = create_summary_matrix(
            gene_family, loci_list, targets_df,
            species_map, phylo_order_map, seq_metadata,
            phase2b_blocks_df, novel_loci_df,
            ref_validation, novel_loci_id_to_name
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
        n_syntenic = len(targets_df[targets_df[COL_PLACEMENT] == PLACEMENT_SYNTENY])
        n_unplaceable = len(targets_df[targets_df[COL_PLACEMENT] != PLACEMENT_SYNTENY])
        print(f"  Syntenic: {n_syntenic} ({100*n_syntenic/len(targets_df):.1f}%)")
        print(f"  Unplaceable: {n_unplaceable} ({100*n_unplaceable/len(targets_df):.1f}%)")
        print(f"  Genomes with targets: {targets_df['genome'].nunique()}")

        # Per-genome breakdown
        genomes_with_unplaceables = targets_df[targets_df[COL_PLACEMENT] != PLACEMENT_SYNTENY]['genome'].nunique()
        print(f"  Genomes with unplaceables: {genomes_with_unplaceables}")

        # Unplaceable subcategory breakdown
        unplaceable_df = targets_df[targets_df[COL_PLACEMENT] != PLACEMENT_SYNTENY].copy()
        if not unplaceable_df.empty:
            # Classify each unplaceable
            subcounts = {'fragment': 0, 'lowquality': 0, 'distinct': 0}
            for _, row in unplaceable_df.iterrows():
                cat = classify_unplaceable(row.to_dict())
                subcounts[cat] += 1

            print(f"\n  Unplaceable breakdown:")
            print(f"    Fragment (no flanking): {subcounts['fragment']}")
            print(f"    Low quality (bitscore<{LOWQUALITY_BITSCORE_THRESHOLD}): {subcounts['lowquality']}")
            print(f"    Distinct (good hit, no locus match): {subcounts['distinct']}")

    # Novel loci statistics
    if not novel_loci_df.empty:
        print("\n" + "-" * 40)
        print("NOVEL LOCI")
        print("-" * 40)
        print(f"  Total novel loci clusters: {len(novel_loci_df)}")
        total_members = novel_loci_df['n_members'].sum()
        total_empty = novel_loci_df['n_empty'].sum()
        print(f"  Total member genomes: {total_members}")
        print(f"  Total empty blocks: {total_empty}")
        for _, row in novel_loci_df.iterrows():
            locus_id = row['locus_id']
            locus_name = row.get('locus_name', locus_id)
            n_mem = row['n_members']
            n_emp = row['n_empty']
            print(f"    {locus_id} ({locus_name}): {n_mem} hits, {n_emp} empty")

    print("\n" + "=" * 80)
    print("PHASE 8b (HELIXER) COMPLETE")
    print("=" * 80)
    print(f"\nSummary matrices saved to: {args.output_dir}")

    return 0


if __name__ == "__main__":
    exit(main())
