#!/usr/bin/env python3
"""
Phase 8a (Helixer): Generate locus-specific matrices.

Creates one matrix per locus showing:
- All genomes (rows) sorted by phylogenetic order
- Flanking protein presence/absence with SwissProt annotations
- Target gene status (length)

This produces the SAME output format as the original tblastn Phase 8a.

Inputs:
- Phase 1 locus definitions (BK flanking proteins)
- Phase 5 Helixer classifications
- Phase 6 Helixer extracted sequences (for length metadata)
- Phase 7 SwissProt annotations (flanking_matches_annotated.tsv, phase2b_flanking_annotated.tsv)
- Phase 2b synteny blocks (phase2b_synteny_blocks.tsv)
- Species/phylo mapping

Key insight: Phase 7 already unions Phase 5 and Phase 2b flanking proteins and annotates
them with SwissProt. This script just displays the results in matrix format.

Outputs:
- {locus_id}_genome_swissprot_matrix.tsv (one per locus)
"""

from __future__ import annotations

import argparse
import re
from collections import OrderedDict
from pathlib import Path
from typing import Dict, Tuple, List

import pandas as pd


def shorten_swissprot_name(stitle: str, max_len: int = 25) -> str:
    """
    Shorten SwissProt stitle for column headers.

    Input format: "sp|P29981|CP4C1_BLADI Cytochrome P450 4C1 OS=Blaberus discoidalis..."
    Output: "Cytochrome P450 4C1"
    """
    if not stitle:
        return ''

    name = stitle

    # Handle sp|ID|GENE format - extract description after the ID parts
    if stitle.startswith('sp|') or stitle.startswith('tr|'):
        parts = stitle.split(' ', 1)
        if len(parts) > 1:
            name = parts[1]  # Everything after "sp|XXX|YYY "

    # Remove organism info "OS=Species name..."
    name = re.sub(r'\s+OS=.*$', '', name)

    # Remove common prefixes
    name = re.sub(r'^(Probable|Putative|Uncharacterized)\s+', '', name, flags=re.IGNORECASE)
    # Remove "RecName: Full=" prefix
    name = re.sub(r'^RecName:\s*Full=', '', name)
    # Remove AltName parts
    name = re.sub(r';\s*AltName:.*$', '', name)

    # Truncate if too long
    if len(name) > max_len:
        name = name[:max_len-2] + '..'

    return name.strip()


def load_bk_swissprot_annotations(bk_swissprot_file: Path) -> Dict[str, str]:
    """
    Load SwissProt annotations for BK proteins from DIAMOND output.

    Expects TSV format: qseqid sseqid pident length evalue bitscore stitle

    Returns: {bk_protein_id: short_swissprot_name}
    """
    bk_annotations = {}

    if not bk_swissprot_file or not bk_swissprot_file.exists():
        return bk_annotations

    try:
        df = pd.read_csv(bk_swissprot_file, sep='\t', header=None,
                         names=['qseqid', 'sseqid', 'pident', 'length', 'evalue', 'bitscore', 'stitle'])

        for _, row in df.iterrows():
            # qseqid format: XP_033229649.1|LOC117181192
            qseqid = str(row['qseqid'])
            protein_id = qseqid.split('|')[0]  # Just the XP ID

            if protein_id not in bk_annotations and pd.notna(row['stitle']):
                short_name = shorten_swissprot_name(str(row['stitle']))
                if short_name:
                    bk_annotations[protein_id] = short_name
    except Exception as e:
        print(f"  Warning: Could not load BK annotations from {bk_swissprot_file}: {e}")

    return bk_annotations


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


def load_extracted_seq_metadata(extraction_dir: Path) -> Dict[str, Dict]:
    """
    Load metadata from Phase 6 length_qc.tsv.

    Returns dict keyed by target_gene_id:
      {target_id: {'length_aa': int, 'length_flag': str}}
    """
    metadata = {}

    if not extraction_dir or not extraction_dir.exists():
        return metadata

    # Load from length_qc.tsv (preferred - has target_id directly)
    length_qc_file = extraction_dir / "length_qc.tsv"
    if length_qc_file.exists():
        df = pd.read_csv(length_qc_file, sep='\t')
        for _, row in df.iterrows():
            target_id = row.get('target_id', '')
            length_aa = int(row.get('extracted_length_aa', 0))
            length_flag = row.get('length_flag', '')
            if target_id:
                metadata[target_id] = {
                    'length_aa': length_aa,
                    'length_flag': length_flag
                }
        return metadata

    return metadata


def shorten_bk_description(desc: str, max_len: int = 30) -> str:
    """Shorten BK gene description for column headers."""
    if not desc:
        return ''

    # Remove species suffix like [Belonocnema kinseyi]
    if '[' in desc:
        desc = desc.split('[')[0].strip()

    # Remove common prefixes/suffixes
    desc = desc.replace(' isoform X1', '').replace(' isoform X2', '')
    desc = desc.replace('-like', '')

    # Truncate
    if len(desc) > max_len:
        desc = desc[:max_len-2] + '..'

    return desc


def load_phase1_flanking_info(phase1_dir: Path) -> Tuple[Dict[str, Dict], Dict[str, str]]:
    """
    Load flanking gene information from Phase 1, including BK annotations.

    Returns:
        locus_info: {locus_id: {'upstream': [...], 'downstream': [...], 'gene_family': str}}
        bk_annotations: {protein_id: short_description}
    """
    locus_info = {}
    bk_annotations = {}

    locus_defs = phase1_dir / 'locus_definitions.tsv'
    if not locus_defs.exists():
        return {}, {}

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
                        # Format: >XP_033232031.1|LOC117183006 U1 XP_033232031.1 description [species]
                        parts = line[1:].strip().split()
                        if len(parts) >= 2:
                            protein_id = parts[0].split('|')[0]
                            position = parts[1] if len(parts) > 1 else ''

                            # Extract description (everything after position and repeated XP ID)
                            if len(parts) >= 4:
                                # parts[2] is usually the XP ID again, parts[3:] is description
                                description = ' '.join(parts[3:])
                                short_desc = shorten_bk_description(description)
                                if short_desc and protein_id not in bk_annotations:
                                    bk_annotations[protein_id] = short_desc

                            if position.startswith('U'):
                                upstream.append(protein_id)
                            elif position.startswith('D'):
                                downstream.append(protein_id)

        locus_info[locus_id] = {
            'upstream': upstream,
            'downstream': downstream,
            'gene_family': gene_family,
        }

    return locus_info, bk_annotations


def create_locus_matrix(
    locus_id: str,
    locus_info: Dict,
    targets_df: pd.DataFrame,
    flanking_matches_df: pd.DataFrame,
    species_map: Dict[str, str],
    phylo_order_map: Dict[str, int],
    seq_metadata: Dict,
    bk_annotations: Dict[str, str] = None,
    empty_blocks_df: pd.DataFrame = None,
    phase2b_flanking_annotations: Dict[str, str] = None,
) -> pd.DataFrame:
    """Create matrix for one locus.

    Args:
        phase2b_flanking_annotations: Dict mapping BK protein ID -> Helixer SwissProt annotation.
            Built from Phase 7's phase2b_flanking_annotated.tsv.
    """

    if bk_annotations is None:
        bk_annotations = {}
    if empty_blocks_df is None:
        empty_blocks_df = pd.DataFrame()
    if phase2b_flanking_annotations is None:
        phase2b_flanking_annotations = {}

    gene_family = locus_info.get('gene_family', locus_id)
    upstream_proteins = locus_info.get('upstream', [])
    downstream_proteins = locus_info.get('downstream', [])

    # Index targets by genome - group for tandem detection
    locus_targets = targets_df[
        (targets_df['placement'] == 'synteny') &
        (targets_df['assigned_to'] == locus_id)
    ]

    # Group all targets per genome for tandem display
    targets_by_genome = {}
    for _, row in locus_targets.iterrows():
        genome = row['genome']
        if genome not in targets_by_genome:
            targets_by_genome[genome] = []
        targets_by_genome[genome].append(row.to_dict())

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

    # Index empty blocks by genome (for this locus)
    # empty_blocks_df has: locus_id, genome, synteny_score, flanking_genes_matched
    # Now also includes: selected, selection_reason, selection_confidence
    empty_blocks_by_genome = {}
    if not empty_blocks_df.empty:
        locus_empty = empty_blocks_df[empty_blocks_df['locus_id'] == locus_id]

        # If 'selected' column exists, filter to selected blocks only
        if 'selected' in locus_empty.columns:
            locus_empty = locus_empty[locus_empty['selected'] == True]

        for _, eb_row in locus_empty.iterrows():
            genome = eb_row['genome']
            if genome not in empty_blocks_by_genome:
                empty_blocks_by_genome[genome] = {
                    'synteny_score': eb_row['synteny_score'],
                    'flanking_genes': set(),
                    'scaffold': eb_row.get('scaffold', ''),
                    'start': eb_row.get('start', ''),
                    'end': eb_row.get('end', ''),
                    'selection_confidence': eb_row.get('selection_confidence', ''),
                    'selection_reason': eb_row.get('selection_reason', ''),
                }
            # Parse flanking genes matched (format: "XP_033229646.1;XP_033229883.1;...")
            if pd.notna(eb_row.get('flanking_genes_matched', '')):
                genes = str(eb_row['flanking_genes_matched']).split(';')
                empty_blocks_by_genome[genome]['flanking_genes'].update(genes)
            # Take best synteny score if multiple blocks
            if eb_row['synteny_score'] > empty_blocks_by_genome[genome]['synteny_score']:
                empty_blocks_by_genome[genome]['synteny_score'] = eb_row['synteny_score']
                empty_blocks_by_genome[genome]['selection_confidence'] = eb_row.get('selection_confidence', '')
                empty_blocks_by_genome[genome]['selection_reason'] = eb_row.get('selection_reason', '')

    # Build matrix rows
    matrix_rows = []

    for genome in sorted(species_map.keys()):
        row = OrderedDict()

        row['genome_id'] = genome
        row['species'] = species_map.get(genome, genome)
        row['phylo_order'] = phylo_order_map.get(genome, 999)

        # Check for synteny (flanking matches)
        # For concordant cases (both Phase 5 target AND Phase 2b block), union the flanking genes
        genome_flanking = flanking_by_genome.get(genome, {})
        total_flanking = len(upstream_proteins) + len(downstream_proteins)
        all_flanking_proteins = set(upstream_proteins + downstream_proteins)

        # Check if we have empty block data for this genome/locus
        has_empty_block = genome in empty_blocks_by_genome
        has_target = genome in targets_by_genome  # Check for actual Phase 5 target

        # Position discordance detection
        position_discordance = False
        if has_target and has_empty_block:
            target_scaffold = targets_by_genome[genome][0].get('scaffold', '')
            block_scaffold = empty_blocks_by_genome[genome].get('scaffold', '')
            if target_scaffold and block_scaffold and target_scaffold != block_scaffold:
                position_discordance = True
            elif target_scaffold == block_scaffold:
                # Check position overlap (targets should be within or near the block)
                target_start = int(targets_by_genome[genome][0].get('start', 0) or 0)
                target_end = int(targets_by_genome[genome][0].get('end', 0) or 0)
                block_start = int(empty_blocks_by_genome[genome].get('start', 0) or 0)
                block_end = int(empty_blocks_by_genome[genome].get('end', 0) or 0)
                # Check if target is within 100kb of block (generous overlap check)
                buffer = 100000  # 100kb
                if target_end < block_start - buffer or target_start > block_end + buffer:
                    position_discordance = True

        # Calculate Phase 5 synteny if we have target
        p5_synteny_pct = 0
        if has_target:
            n_p5_matches = len(set(genome_flanking.keys()) & all_flanking_proteins)
            if total_flanking > 0:
                p5_synteny_pct = round((n_p5_matches / total_flanking) * 100, 1)

        # Get Phase 2b data if available
        p2b_synteny_pct = 0
        if has_empty_block:
            eb_data = empty_blocks_by_genome[genome]
            p2b_synteny_pct = round(eb_data['synteny_score'] * 100, 1)
            # UNION: Add Phase 2b flanking genes to genome_flanking
            # Use Phase 7's Helixer SwissProt annotations (preferred) or BK annotations (fallback)
            for fg in eb_data['flanking_genes']:
                if fg in all_flanking_proteins and fg not in genome_flanking:
                    # phase2b_flanking_annotations maps BK protein -> Helixer SwissProt annotation
                    if fg in phase2b_flanking_annotations and phase2b_flanking_annotations[fg]:
                        genome_flanking[fg] = phase2b_flanking_annotations[fg]
                    elif fg in bk_annotations:
                        genome_flanking[fg] = bk_annotations[fg]
                    else:
                        genome_flanking[fg] = "match"

        # Use max synteny score (for concordant, this gives best evidence)
        if has_target or has_empty_block:
            synteny_pct = max(p5_synteny_pct, p2b_synteny_pct)
            n_matches = len(set(genome_flanking.keys()) & all_flanking_proteins)
        else:
            n_matches = 0
            synteny_pct = 0

        row['synteny_pct'] = synteny_pct
        row['num_proteins_found'] = n_matches

        # Add synteny confidence (from tiered block selection)
        if has_empty_block:
            eb_data = empty_blocks_by_genome[genome]
            row['synteny_confidence'] = eb_data.get('selection_confidence', '')
        else:
            row['synteny_confidence'] = ''

        # Add scaffold info - prefer target (use first/best), fallback to empty block
        if genome in targets_by_genome:
            targets_list = targets_by_genome[genome]
            target = targets_list[0]  # First target (highest score from Phase 5)
            row['scaffold'] = target.get('scaffold', '')
            row['strand'] = target.get('strand', '')
            row['start'] = target.get('start', '')
            row['end'] = target.get('end', '')
        elif has_empty_block:
            eb_data = empty_blocks_by_genome[genome]
            row['scaffold'] = eb_data.get('scaffold', '')
            row['strand'] = ''
            row['start'] = eb_data.get('start', '')
            row['end'] = eb_data.get('end', '')
        else:
            row['scaffold'] = ''
            row['strand'] = ''
            row['start'] = ''
            row['end'] = ''

        # Upstream proteins (U14, U13, ..., U1)
        for i, protein_id in enumerate(reversed(upstream_proteins), 1):
            pos_num = len(upstream_proteins) - i + 1
            # Use BK SwissProt annotation if available, otherwise protein ID
            if protein_id in bk_annotations:
                col_name = f"U{pos_num}_{bk_annotations[protein_id]}"
            else:
                col_name = f"U{pos_num}_{protein_id}"
            if protein_id in genome_flanking:
                row[col_name] = genome_flanking[protein_id]  # SwissProt annotation or "match"
            else:
                row[col_name] = ""

        # TARGET column
        # Include confidence indicator for LOW confidence synteny blocks
        confidence_suffix = ''
        if has_empty_block and row['synteny_confidence'] == 'LOW':
            confidence_suffix = ' *LOW*'
        # Add position discordance flag
        if position_discordance:
            confidence_suffix += ' *DISC*'

        if genome in targets_by_genome:
            targets_list = targets_by_genome[genome]
            n_targets = len(targets_list)

            # Get lengths for all targets
            lengths = []
            for t in targets_list:
                target_id = t.get('target_gene_id', '')
                meta = seq_metadata.get(target_id, {})
                length = meta.get('length_aa', 0)
                if length > 0:
                    lengths.append(int(length))

            # Sort lengths descending
            lengths.sort(reverse=True)

            # Format: synteny% [length] or synteny% [len1, len2] for tandems
            syn_pct = int(synteny_pct)
            if lengths:
                lengths_str = ', '.join(str(l) for l in lengths)
                row['TARGET'] = f"{gene_family} {syn_pct}% [{lengths_str}]{confidence_suffix}"
            else:
                row['TARGET'] = f"{gene_family} {syn_pct}% [{n_targets} hit{'s' if n_targets > 1 else ''}]{confidence_suffix}"
        elif synteny_pct > 0:
            row['TARGET'] = f"{gene_family} {int(synteny_pct)}% [empty]{confidence_suffix}"
        else:
            row['TARGET'] = f"{gene_family} [0%]"

        # Downstream proteins (D1, D2, ...)
        for i, protein_id in enumerate(downstream_proteins, 1):
            # Use BK SwissProt annotation if available, otherwise protein ID
            if protein_id in bk_annotations:
                col_name = f"D{i}_{bk_annotations[protein_id]}"
            else:
                col_name = f"D{i}_{protein_id}"
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
    parser.add_argument("--phase5b-dir", type=Path, default=None,
                        help="Phase 5b refined classifications (for tandems)")
    parser.add_argument("--phase2b-dir", type=Path, default=None,
                        help="Phase 2b empty blocks directory")
    parser.add_argument("--phase6-dir", type=Path, required=True,
                        help="Phase 6 Helixer extracted sequences")
    parser.add_argument("--phase7-dir", type=Path, default=None,
                        help="Phase 7 Helixer SwissProt annotations (optional)")
    parser.add_argument("--bk-swissprot", type=Path, default=None,
                        help="BK flanking proteins SwissProt annotations TSV (for column headers)")
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

    # Phase 1 flanking info (includes BK annotations from FASTA headers)
    locus_info, bk_annotations = load_phase1_flanking_info(args.phase1_dir)
    print(f"  Loci: {len(locus_info)}")
    print(f"  BK annotations (from Phase 1): {len(bk_annotations)}")

    # Phase 5 classifications - load syntenic targets
    targets_df = pd.DataFrame()
    syntenic_file = args.phase5_dir / "syntenic_targets.tsv"
    if syntenic_file.exists():
        targets_df = pd.read_csv(syntenic_file, sep='\t')
        # Normalize column names for compatibility
        # Only copy if destination column doesn't exist
        if 'classification' in targets_df.columns and 'placement' not in targets_df.columns:
            targets_df['placement'] = targets_df['classification']
        if 'locus_id' in targets_df.columns and 'assigned_to' not in targets_df.columns:
            targets_df['assigned_to'] = targets_df['locus_id']
        print(f"  Syntenic targets: {len(targets_df)}")
    else:
        print("  No syntenic targets found")

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

    # Phase 2b synteny blocks (try new filename first, then old for backward compatibility)
    empty_blocks_df = pd.DataFrame()
    if args.phase2b_dir:
        phase2b_blocks_file = args.phase2b_dir / "phase2b_synteny_blocks.tsv"
        if not phase2b_blocks_file.exists():
            phase2b_blocks_file = args.phase2b_dir / "empty_synteny_blocks.tsv"

        if phase2b_blocks_file.exists():
            empty_blocks_df = pd.read_csv(phase2b_blocks_file, sep='\t')
            print(f"  Phase 2b synteny blocks: {len(empty_blocks_df)}")
        else:
            print("  No Phase 2b synteny blocks data")
    else:
        print("  No Phase 2b directory provided")

    # Phase 7 Phase 2b flanking annotations (Helixer SwissProt annotations)
    # Build lookup: BK protein ID -> Helixer SwissProt annotation
    phase2b_flanking_annotations = {}
    if args.phase7_dir:
        phase2b_annot_file = args.phase7_dir / "phase2b_flanking_annotated.tsv"
        if phase2b_annot_file.exists():
            p2b_annot_df = pd.read_csv(phase2b_annot_file, sep='\t')
            for _, row in p2b_annot_df.iterrows():
                bk_id = str(row.get('bk_flanking_id', '')).split('|')[0]  # Get just XP ID
                sp_name = row.get('swissprot_name', '')
                sp_pident = row.get('swissprot_pident', 0)

                if bk_id and sp_name and pd.notna(sp_name):
                    if sp_pident and sp_pident > 0:
                        annotation = f"{sp_name} ({int(sp_pident)}%)"
                    else:
                        annotation = sp_name
                    # Keep best annotation per BK protein
                    if bk_id not in phase2b_flanking_annotations:
                        phase2b_flanking_annotations[bk_id] = annotation
            print(f"  Phase 2b flanking annotations (from Phase 7): {len(phase2b_flanking_annotations)}")

    # Optional: Override BK annotations with SwissProt file (usually not needed)
    if args.bk_swissprot:
        swissprot_annots = load_bk_swissprot_annotations(args.bk_swissprot)
        if swissprot_annots:
            bk_annotations.update(swissprot_annots)
            print(f"  BK SwissProt annotations (override): {len(swissprot_annots)}")

    # Generate matrices
    print("\n[2] Generating matrices...")

    for locus_id, info in locus_info.items():
        matrix_df = create_locus_matrix(
            locus_id, info, targets_df, flanking_df,
            species_map, phylo_order_map, seq_metadata, bk_annotations,
            empty_blocks_df, phase2b_flanking_annotations
        )

        if not matrix_df.empty:
            output_file = args.output_dir / f"{args.family}_{locus_id}_genome_swissprot_matrix.tsv"
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
