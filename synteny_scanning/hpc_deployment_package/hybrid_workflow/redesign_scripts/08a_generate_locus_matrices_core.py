#!/usr/bin/env python3
"""
Step 08a: Generate locus-specific matrices.

Creates one matrix per locus showing:
- All genomes (rows) sorted by phylogenetic order
- Flanking protein presence/absence with SwissProt annotations
- Target gene status (length + functional/pseudogene)

Usage:
    python 08a_generate_locus_matrices.py \\
        --locus-defs <path/to/locus_definitions.tsv> \\
        --synteny-dir <path/to/02_synteny_blocks> \\
        --blocks <path/to/synteny_blocks_filtered.tsv> \\
        --targets <path/to/all_targets_classified.tsv> \\
        [--graded-syntenic <path/to/syntenic_targets_graded.tsv>] \\
        [--keep-grades intact,degraded_fragment] \\
        --swissprot <path/to/genome_specific_swissprot_annotations.tsv> \\
        --reference-proteins <path/to/protein.faa> \\
        --species-map <path/to/gca_to_species.tsv> \\
        --output-dir <path/to/output/matrices>
"""

from pathlib import Path
import pandas as pd
from collections import OrderedDict
import argparse
import sys
import re

def load_extracted_seq_metadata(extraction_dir):
    """
    Load metadata from Phase 6 extracted sequences.

    Returns dict keyed by (genome_id, locus_name):
      {(genome, locus): {'length_aa': int, 'status': str}}
    """
    metadata = {}

    if not extraction_dir or not Path(extraction_dir).exists():
        return metadata

    extraction_path = Path(extraction_dir)

    # Traverse genome directories
    for genome_dir in extraction_path.iterdir():
        if not genome_dir.is_dir():
            continue

        genome_id = genome_dir.name

        # Traverse locus directories (skip UNPLACED - only load syntenic targets)
        for locus_dir in genome_dir.iterdir():
            if not locus_dir.is_dir() or locus_dir.name == 'UNPLACED':
                continue

            locus_name = locus_dir.name

            # Find CDS FASTA for status and protein for length
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

                metadata[(genome_id, locus_name)] = {
                    'length_aa': length_aa,
                    'status': status
                }

    return metadata

def get_protein_names_from_faa(proteins_file):
    """Extract protein names from protein.faa for better column headers."""
    protein_names = {}

    from Bio import SeqIO

    for record in SeqIO.parse(proteins_file, "fasta"):
        # Parse header: >XP_033228127.1 low density lipoprotein receptor adapter protein 1-like [Belonocnema kinseyi]
        protein_id = record.id
        description = record.description

        # Extract protein name from description
        if ' ' in description:
            # Remove protein ID and organism
            parts = description.split(' ', 1)[1]  # Skip the ID
            if '[' in parts:
                protein_name = parts.split('[')[0].strip()
            else:
                protein_name = parts.strip()

            # Simplify name - remove "uncharacterized protein " prefix
            if protein_name.startswith('uncharacterized protein '):
                protein_name = protein_name.replace('uncharacterized protein ', '')

            # Shorten long names
            if len(protein_name) > 30:
                protein_name = protein_name[:27] + '...'

            protein_names[protein_id] = protein_name
        else:
            protein_names[protein_id] = protein_id

    return protein_names

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
                    phylo_order_map[genome_id] = 999  # Default for missing

    return species_map, phylo_order_map

def parse_flanking_proteins_from_faa(flanking_file):
    """
    Parse flanking protein file to extract protein IDs, descriptions, and U/D positions.

    Returns tuple of:
    - upstream_proteins: list of (protein_id, protein_name) for upstream genes
    - downstream_proteins: list of (protein_id, protein_name) for downstream genes

    Expects Phase 1 format with U/D labels in descriptions:
      >XP_033209112.1|LOC117167954 U1 XP_033209112.1 description [species]
      >XP_033209139.1|LOC117167973 D1 XP_033209139.1 description [species]
    """
    from Bio import SeqIO

    upstream_proteins = []
    downstream_proteins = []

    if not flanking_file.exists():
        print(f"    Warning: Flanking file not found: {flanking_file}")
        return [], []

    for record in SeqIO.parse(flanking_file, "fasta"):
        # Header format: >XP_033209112.1|LOC117167954 U1 XP_033209112.1 description [species]
        protein_id = record.id.split('|')[0]  # Get XP_ ID

        # Extract description and position label
        description = record.description
        parts = description.split()

        if len(parts) < 2:
            # Malformed header, skip
            continue

        # Check for U/D label (should be 2nd element)
        position_label = parts[1] if len(parts) > 1 else None

        # Extract description (skip ID, position label, and protein ID repeat)
        if len(parts) >= 4:
            desc_parts = parts[3:]  # Everything after "U1 XP_..."
            desc = ' '.join(desc_parts)
            # Remove species name in brackets
            if '[' in desc:
                desc = desc.split('[')[0].strip()
        else:
            desc = protein_id

        # Add to appropriate list based on U/D label
        if position_label and position_label.startswith('U'):
            upstream_proteins.append((protein_id, desc))
        elif position_label and position_label.startswith('D'):
            downstream_proteins.append((protein_id, desc))
        else:
            # No label - legacy format, treat as upstream
            upstream_proteins.append((protein_id, desc))

    return upstream_proteins, downstream_proteins


def calculate_reference_locus_span(flanking_df, locus_id, synteny_dir):
    """
    Calculate the expected locus span from the reference genome's flanking hits.

    The reference genome is BK (Belonocnema_kinseyi_GCF), where loci are defined.
    Returns (expected_span_kb, ref_scaffold, ref_min_coord, ref_max_coord) or Nones if unavailable.
    """
    if flanking_df.empty:
        return None, None, None, None

    # Reference genome ID
    ref_genome = 'Belonocnema_kinseyi_GCF'

    # Filter to reference genome hits
    ref_hits = flanking_df[flanking_df['genome'] == ref_genome]

    if ref_hits.empty:
        return None, None, None, None

    # Get scaffold with most hits (in case of multiple scaffolds)
    scaffold_counts = ref_hits['sseqid'].value_counts()
    main_scaffold = scaffold_counts.index[0]

    # Filter to main scaffold
    scaffold_hits = ref_hits[ref_hits['sseqid'] == main_scaffold]

    min_coord = scaffold_hits['coord_start'].min()
    max_coord = scaffold_hits['coord_end'].max()
    span_kb = (max_coord - min_coord) / 1000.0

    return span_kb, main_scaffold, min_coord, max_coord


def create_locus_matrix(locus_id, locus_info, blocks_df, targets_df, swissprot_df, flanking_df, species_map, phylo_order_map, protein_names, synteny_dir, seq_metadata):
    """Create matrix for one locus."""

    print(f"\n  Processing {locus_id}...")

    # Derive gene family from locus_id (e.g., "BK_chr2_a" -> "ferritin" or use locus_id)
    # If gene_family column exists in locus_info, use it; otherwise use locus_id
    gene_family = locus_info.get('gene_family', locus_id.split('_')[0] if '_' in locus_id else locus_id)

    # Get flanking file path from locus_info, with robust fallbacks
    if 'flanking_file' in locus_info and pd.notna(locus_info['flanking_file']):
        flanking_file = Path(locus_info['flanking_file'])
    else:
        flanking_file = synteny_dir / locus_id / f"{locus_id}_flanking_filtered.faa"
    # If not found, try Phase 1 deduplicated, then non-dedup
    if not flanking_file.exists():
        try:
            import argparse as _ap, sys as _sys
            _p = _ap.ArgumentParser(add_help=False)
            _p.add_argument('--locus-defs')
            _known, _ = _p.parse_known_args(_sys.argv[1:])
            if _known.locus_defs:
                phase1_dir = Path(_known.locus_defs).parent
                cand = phase1_dir / f"{locus_id}_flanking_dedup.faa"
                if not cand.exists():
                    cand = phase1_dir / f"{locus_id}_flanking.faa"
                if cand.exists():
                    flanking_file = cand
        except Exception:
            pass

    # Parse flanking proteins from .faa file (now returns (upstream, downstream) tuples)
    upstream_data, downstream_data = parse_flanking_proteins_from_faa(flanking_file)

    # Unpack tuples
    upstream_proteins = [protein_id for protein_id, _ in upstream_data]
    upstream_names = [name for _, name in upstream_data]
    downstream_proteins = [protein_id for protein_id, _ in downstream_data]
    downstream_names = [name for _, name in downstream_data]

    total_expected = len(upstream_proteins) + len(downstream_proteins)

    # Calculate expected locus span from reference genome (BK)
    expected_span_kb, ref_scaffold, ref_min, ref_max = calculate_reference_locus_span(
        flanking_df, locus_id, synteny_dir
    )
    if expected_span_kb:
        print(f"  DEBUG: Reference locus span for {locus_id}: {expected_span_kb:.1f} kb ({ref_scaffold}: {ref_min}-{ref_max})")

    # Get blocks for this locus
    if not blocks_df.empty and 'locus_id' in blocks_df.columns:
        locus_blocks = blocks_df[blocks_df['locus_id'] == locus_id]
    else:
        locus_blocks = pd.DataFrame()

    # Index SwissProt annotations
    # KEY: Use bk_protein_id (the landmark protein accession) for matching
    swissprot_map = {}
    if not swissprot_df.empty:
        for _, row in swissprot_df.iterrows():
            if row['locus'] == locus_id:
                # Use bk_protein_id (e.g., XP_033209112.1) for the key, not bk_protein_annotation
                protein_id_key = row.get('bk_protein_id', row.get('bk_protein', 'unknown'))
                key = (row['genome'], row['locus'], protein_id_key)

                # Format annotation
                if row['swissprot_description'] == 'No SwissProt match' or pd.isna(row['swissprot_description']):
                    annotation = "no match"
                elif row['swissprot_description'] == 'SwissProt database not available':
                    annotation = "no match"  # Mark as present even if SwissProt not available
                else:
                    # Extract protein name
                    desc = row['swissprot_description']
                    if 'RecName: Full=' in desc:
                        start = desc.find('Full=') + 5
                        end_bracket = desc.find('[', start)
                        end_semi = desc.find(';', start)

                        if end_bracket > 0 and end_semi > 0:
                            end = min(end_bracket, end_semi)
                        elif end_bracket > 0:
                            end = end_bracket
                        elif end_semi > 0:
                            end = end_semi
                        else:
                            end = len(desc)

                        protein_name = desc[start:end].strip()
                    else:
                        protein_name = desc.split('[')[0].strip()[:30]  # Limit length

                    pident = round(float(row['percent_identity']))
                    annotation = f"{protein_name} ({pident}%)"

                swissprot_map[key] = annotation

    # Index target genes: synteny targets assigned to this locus + unplaceable targets from this parent locus
    targets_index = {}
    debug_unplaceable_count = 0
    debug_synteny_count = 0
    if not targets_df.empty:
        for _, target in targets_df.iterrows():
            # Include syntenic targets assigned to this locus
            if target.get('placement') == 'synteny' and target.get('assigned_to') == locus_id:
                key = target['genome']
                targets_index.setdefault(key, []).append(target)
                debug_synteny_count += 1
            # Also include unplaceable targets whose parent_locus matches this locus
            elif target.get('placement') == 'unplaceable' and target.get('parent_locus') == locus_id:
                key = target['genome']
                targets_index.setdefault(key, []).append(target)
                debug_unplaceable_count += 1

    print(f"  DEBUG: {locus_id} targets_df has {len(targets_df)} rows, indexed {debug_synteny_count} synteny + {debug_unplaceable_count} unplaceable = {len(targets_index)} genomes")

    # Build matrix rows
    matrix_rows = []

    # Get all other genomes (from blocks or from master list)
    all_genomes = set()

    # Add genomes with blocks
    all_genomes.update(locus_blocks['genome'].unique())

    # Add all genomes from species map to ensure complete coverage
    all_genomes.update(species_map.keys())

    # Exclude deprecated/renamed genomes
    for g in ('GCA_010883055.1',):
        if g in all_genomes:
            all_genomes.remove(g)

    for genome in all_genomes:
        # Skip GCA_010883055.1 - old Belonocnema kinseyi entry replaced by Belonocnema_kinseyi_GCF
        if genome == 'GCA_010883055.1':
            continue

        row = OrderedDict()

        # Metadata
        row['genome_id'] = genome
        row['species'] = species_map.get(genome, genome)
        row['phylo_order'] = phylo_order_map.get(genome, 999)

        # Check if this genome has a synteny block
        genome_block = locus_blocks[locus_blocks['genome'] == genome]

        if not genome_block.empty:
            block = genome_block.iloc[0]
            row['synteny_pct'] = round(block['num_query_matches'] / total_expected * 100, 1)
            row['num_proteins_found'] = block['num_query_matches']
            row['scaffold'] = block['scaffold']
            row['strand'] = block.get('strand', '.')  # Use '.' if strand not available
            row['start'] = block['start']
            row['end'] = block['end']

            # Get flanking protein details if available
            proteins_found = set()
            if not flanking_df.empty:
                # Calculate search window - use expected locus span if available
                # This ensures we count all flanking genes within the expected region,
                # not just those within the selected block's narrow coordinates
                if expected_span_kb and expected_span_kb > 0:
                    # Use expected span centered on block midpoint
                    block_mid = (block['start'] + block['end']) / 2
                    half_span = (expected_span_kb * 1000) / 2
                    search_start = block_mid - half_span
                    search_end = block_mid + half_span
                else:
                    # Fallback to block coordinates
                    search_start = block['start']
                    search_end = block['end']

                # Build filter conditions using search window
                filter_cond = (
                    (flanking_df['genome'] == genome) &
                    (flanking_df['sseqid'] == block['scaffold']) &
                    (flanking_df['coord_end'] >= search_start) &
                    (flanking_df['coord_start'] <= search_end)
                )
                # Add strand filter only if strand column exists in both
                if 'strand' in flanking_df.columns and 'strand' in block.index:
                    filter_cond = filter_cond & (flanking_df['strand'] == block['strand'])

                block_hits = flanking_df[filter_cond]

                # Map protein IDs to their positions (U1, U2, D1, D2, etc)
                # Query IDs are like "XP_033209222.1|LOC117168016"
                for _, hit in block_hits.iterrows():
                    query_id = hit['qseqid']

                    # Extract protein ID (part before |)
                    protein_id = query_id.split('|')[0] if '|' in query_id else query_id

                    # Check if this protein is in our upstream list
                    if protein_id in upstream_proteins:
                        idx = upstream_proteins.index(protein_id)
                        position_num = len(upstream_proteins) - idx  # U14, U13, ..., U1
                        position_id = f"U{position_num}"
                        proteins_found.add(position_id)
                        proteins_found.add(protein_id)  # Also add for SwissProt lookup

                    # Check if this protein is in our downstream list
                    elif protein_id in downstream_proteins:
                        idx = downstream_proteins.index(protein_id)
                        position_num = idx + 1  # D1, D2, D3, ...
                        position_id = f"D{position_num}"
                        proteins_found.add(position_id)
                        proteins_found.add(protein_id)  # Also add for SwissProt lookup
        else:
            # No synteny block for this genome
            row['synteny_pct'] = 0
            row['num_proteins_found'] = 0
            row['scaffold'] = ""
            row['strand'] = ""
            row['start'] = ""
            row['end'] = ""
            proteins_found = set()

        # Upstream proteins (descending order)
        for i in range(len(upstream_proteins), 0, -1):
            protein_id = upstream_proteins[i-1]

            # Use protein name from faa file if available, otherwise fall back to LOC
            if protein_id in protein_names:
                display_name = protein_names[protein_id]
            else:
                display_name = upstream_names[i-1]

            col_name = f"U{i}_{display_name}"

            if f"U{i}" in proteins_found or protein_id in proteins_found:
                # Look up SwissProt annotation
                annotation = swissprot_map.get((genome, locus_id, protein_id), "no match")
                row[col_name] = annotation
            else:
                row[col_name] = ""

        # Target gene column
        genome_targets = targets_index.get(genome, [])

        if genome_targets:
            # Use extraction metadata keyed by (genome, placed locus = locus_id)
            meta = seq_metadata.get((genome, locus_id), {})
            parts = []
            if meta:
                length = meta.get('length_aa', 0)
                if length > 0:
                    # Just show length - I/F/P status removed as unreliable
                    parts.append(str(length))

            # Only add generic "hit" if no metadata found (not extracted)
            if not parts:
                # No extraction metadata - show generic hit(s)
                parts.extend(["hit"] * len(genome_targets))

            row['TARGET'] = f"{gene_family} [{'; '.join(parts)}]"
        else:
            if not genome_block.empty:
                # Synteny block exists but no target found
                row['TARGET'] = f"{gene_family} [empty]"
            else:
                # No synteny block - report as 0% conservation
                row['TARGET'] = f"{gene_family} [0%]"

        # Downstream proteins (ascending order)
        for i, protein_id in enumerate(downstream_proteins, 1):
            # Use protein name from faa file if available, otherwise fall back to LOC
            if protein_id in protein_names:
                display_name = protein_names[protein_id]
            else:
                display_name = downstream_names[i-1]

            col_name = f"D{i}_{display_name}"

            if f"D{i}" in proteins_found or protein_id in proteins_found:
                # Look up SwissProt annotation
                annotation = swissprot_map.get((genome, locus_id, protein_id), "no match")
                row[col_name] = annotation
            else:
                row[col_name] = ""

        # Recalculate synteny_pct based on actual flanking columns filled
        # Count U/D columns with data (including "no match")
        flanking_cols_with_data = 0
        total_flanking_cols = len(upstream_proteins) + len(downstream_proteins)

        for i in range(len(upstream_proteins), 0, -1):
            protein_id = upstream_proteins[i-1]
            if protein_id in protein_names:
                display_name = protein_names[protein_id]
            else:
                display_name = upstream_names[i-1]
            col_name = f"U{i}_{display_name}"

            val = row.get(col_name, "")
            if val and val != "" and str(val) != "nan":
                flanking_cols_with_data += 1

        for i, protein_id in enumerate(downstream_proteins, 1):
            if protein_id in protein_names:
                display_name = protein_names[protein_id]
            else:
                display_name = downstream_names[i-1]
            col_name = f"D{i}_{display_name}"

            val = row.get(col_name, "")
            if val and val != "" and str(val) != "nan":
                flanking_cols_with_data += 1

        # Update synteny percentage based on flanking columns
        if total_flanking_cols > 0:
            row['synteny_pct'] = round((flanking_cols_with_data / total_flanking_cols) * 100, 1)
            row['num_proteins_found'] = flanking_cols_with_data
        else:
            row['synteny_pct'] = 0
            row['num_proteins_found'] = 0

        matrix_rows.append(row)

    # Create DataFrame and sort
    if matrix_rows:
        matrix_df = pd.DataFrame(matrix_rows)
        matrix_df = matrix_df.sort_values('phylo_order', ascending=True)
        return matrix_df

    return pd.DataFrame()  # Empty if no data

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate locus-specific matrices with flanking proteins and SwissProt annotations",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('--locus-defs', required=True, type=Path,
                        help='Path to locus_definitions.tsv')
    parser.add_argument('--synteny-dir', required=True, type=Path,
                        help='Path to 02_synteny_blocks directory')
    parser.add_argument('--blocks', type=Path,
                        help='Path to synteny_blocks_filtered.tsv (optional)')
    parser.add_argument('--targets', type=Path,
                        help='Path to all_targets_classified.tsv (optional)')
    parser.add_argument('--graded-syntenic', type=Path,
                        help='Optional path to syntenic_targets_graded.tsv (Phase 6 grading output)')
    parser.add_argument('--keep-grades', type=str, default='intact,degraded_fragment',
                        help='Comma-separated list of grades to keep (default: intact,degraded_fragment)')
    parser.add_argument('--swissprot', type=Path,
                        help='Path to genome_specific_swissprot_annotations.tsv (optional)')
    parser.add_argument('--reference-proteins', required=False, type=Path,
                        help='Path to reference protein.faa file (for column headers). Optional - will derive from SwissProt if not provided.')
    parser.add_argument('--extracted-seqs', required=False, type=Path,
                        help='Path to extracted sequences directory (for target length/status metadata)')
    parser.add_argument('--species-map', required=True, type=Path,
                        help='Path to gca_to_species.tsv')
    parser.add_argument('--output-dir', required=True, type=Path,
                        help='Output directory for matrices')

    return parser.parse_args()

def main():
    """Main execution function."""
    args = parse_args()

    print("=" * 80)
    print("STEP 08a: GENERATE LOCUS-SPECIFIC MATRICES")
    print("=" * 80)
    print(f"\nInput files:")
    print(f"  Locus definitions: {args.locus_defs}")
    print(f"  Synteny directory: {args.synteny_dir}")
    print(f"  Blocks: {args.blocks if args.blocks else 'None'}")
    print(f"  Targets: {args.targets if args.targets else 'None'}")
    print(f"  SwissProt: {args.swissprot if args.swissprot else 'None'}")
    print(f"  Reference proteins: {args.reference_proteins}")
    print(f"  Species map: {args.species_map}")
    print(f"  Output directory: {args.output_dir}")

    # Create output directory
    args.output_dir.mkdir(exist_ok=True, parents=True)

    # Load all required data
    print("\n[1] Loading data...")

    # Locus definitions
    loci_df = pd.read_csv(args.locus_defs, sep='\t')
    print(f"  Loaded {len(loci_df)} locus definitions")

    # Filtered synteny blocks
    if args.blocks and args.blocks.exists():
        blocks_df = pd.read_csv(args.blocks, sep='\t')
        print(f"  Loaded {len(blocks_df)} synteny blocks")
    else:
        blocks_df = pd.DataFrame()
        print("  No synteny blocks found")

    # Classified targets
    if args.targets and args.targets.exists():
        targets_df = pd.read_csv(args.targets, sep='\t')
        # Filter to syntenic only for main column
        syntenic_targets = targets_df[targets_df['placement'] == 'synteny']
        # Optional quality grading filter
        graded_path = args.graded_syntenic if args.graded_syntenic else (args.targets.parent / 'syntenic_targets_graded.tsv')
        if graded_path and graded_path.exists():
            try:
                graded_df = pd.read_csv(graded_path, sep='\t')
                keep = {g.strip() for g in str(args.keep_grades).split(',') if g.strip()}
                graded_keep = graded_df[graded_df['grade'].isin(keep)].copy()
                # Composite key for robust match
                def _key(df):
                    return set(zip(df.get('genome', []), df.get('locus_name', []), df.get('query_id', []), df.get('start', []), df.get('end', [])))
                k_keep = _key(graded_keep)
                syntenic_targets = syntenic_targets[syntenic_targets.apply(lambda r: (r.get('genome',''), r.get('locus_name',''), r.get('query_id',''), r.get('start',None), r.get('end',None)) in k_keep, axis=1)]
                print(f"  Loaded {len(syntenic_targets)} syntenic targets (graded kept: {len(graded_keep)})")
            except Exception as e:
                print(f"  Warning: could not apply graded filter ({e}); using unfiltered syntenic")
                print(f"  Loaded {len(syntenic_targets)} syntenic targets")
        else:
            print(f"  Loaded {len(syntenic_targets)} syntenic targets")
    else:
        syntenic_targets = pd.DataFrame()
        print("  No classified targets found")

    # SwissProt annotations
    if args.swissprot and args.swissprot.exists():
        swissprot_df = pd.read_csv(args.swissprot, sep='\t')
        print(f"  Loaded {len(swissprot_df)} SwissProt annotations")
    else:
        swissprot_df = pd.DataFrame()
        print("  No SwissProt annotations found")

    # Species mapping
    species_map, phylo_order_map = load_species_and_phylo(args.species_map)
    print(f"  Loaded species mapping for {len(species_map)} genomes")

    # Load extraction metadata (length and status from FASTA headers)
    if args.extracted_seqs and args.extracted_seqs.exists():
        seq_metadata = load_extracted_seq_metadata(args.extracted_seqs)
        print(f"  Loaded extraction metadata for {len(seq_metadata)} loci")
    else:
        seq_metadata = {}
        print("  No extraction metadata found")

    # Get protein names for better column headers
    if args.reference_proteins and args.reference_proteins.exists():
        protein_names = get_protein_names_from_faa(args.reference_proteins)
        print(f"  Loaded protein names for {len(protein_names)} proteins from reference")
    elif not swissprot_df.empty:
        # Derive from SwissProt annotations (use bk_protein_id and bk_protein_annotation)
        protein_names = {}
        for _, row in swissprot_df.iterrows():
            protein_id = row['bk_protein_id']
            if protein_id not in protein_names and 'bk_protein_annotation' in row:
                # Simplify SwissProt description for column header
                desc = row['bk_protein_annotation']
                if 'RecName: Full=' in desc:
                    name = desc.split('RecName: Full=')[1].split(';')[0].strip()
                elif ' ' in desc:
                    name = desc.split('[')[0].strip() if '[' in desc else desc
                else:
                    name = desc
                # Shorten if too long
                if len(name) > 30:
                    name = name[:27] + '...'
                protein_names[protein_id] = name
        print(f"  Derived protein names for {len(protein_names)} proteins from SwissProt")
    else:
        protein_names = {}
        print(f"  No protein names available - using IDs only")

    # Process each locus
    print("\n[2] Generating matrices...")

    for _, locus_info in loci_df.iterrows():
        locus_id = locus_info['locus_id']

        # Load flanking protein details if available
        flanking_file = args.synteny_dir / locus_id / "flanking_blast_all.tsv"
        if flanking_file.exists():
            flanking_df = pd.read_csv(flanking_file, sep='\t')
        else:
            flanking_df = pd.DataFrame()

        # Create matrix
        matrix_df = create_locus_matrix(
            locus_id, locus_info, blocks_df, syntenic_targets,
            swissprot_df, flanking_df, species_map, phylo_order_map, protein_names,
            args.synteny_dir, seq_metadata
        )

        if not matrix_df.empty:
            # Save matrix
            output_file = args.output_dir / f"{locus_id}_genome_swissprot_matrix.tsv"
            matrix_df.to_csv(output_file, sep='\t', index=False)

            print(f"    Saved: {output_file.name}")
            print(f"    Rows: {len(matrix_df)}")
            print(f"    With synteny: {len(matrix_df[matrix_df['synteny_pct'] > 0])}")

            # Show example target entries
            target_examples = matrix_df[matrix_df['TARGET'].str.contains('[', regex=False, na=False)].head(3)
            if len(target_examples) > 0:
                print(f"    Example target entries:")
                for _, ex in target_examples.iterrows():
                    print(f"      {ex['species']}: {ex['TARGET']}")

    # Summary
    print("\n" + "=" * 80)
    print("LOCUS MATRICES COMPLETE")
    print("=" * 80)
    print(f"\nMatrices saved to: {args.output_dir}")
    print("\nNext step: 08b_generate_summary_matrices.py")

if __name__ == "__main__":
    main()
