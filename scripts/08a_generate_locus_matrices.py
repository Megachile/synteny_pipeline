#!/usr/bin/env python3
"""
Step 08a: Generate locus-specific matrices.

Creates one matrix per locus showing:
- All genomes (rows) sorted by phylogenetic order
- Flanking protein presence/absence with SwissProt annotations
- Target gene status (length + functional/pseudogene)
"""

from pathlib import Path
import pandas as pd
from collections import OrderedDict
import config

def get_protein_names_from_faa():
    """Extract protein names from protein.faa for better column headers."""
    protein_names = {}

    from Bio import SeqIO

    for record in SeqIO.parse(config.BK_PROTEINS_FILE, "fasta"):
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
                    phylo_order_map[genome_id] = 999  # Default for missing

    return species_map, phylo_order_map

def parse_flanking_proteins_from_faa(flanking_file):
    """
    Parse flanking protein file to extract protein IDs and descriptions.

    Returns tuple of:
    - protein_ids: list of protein IDs
    - protein_names: list of protein descriptions
    """
    from Bio import SeqIO

    protein_ids = []
    protein_names = []

    if not flanking_file.exists():
        print(f"    Warning: Flanking file not found: {flanking_file}")
        return [], []

    for record in SeqIO.parse(flanking_file, "fasta"):
        # Header format: >XP_033209112.1|LOC117167954 XP_033209112.1 description [species]
        protein_id = record.id.split('|')[0]  # Get XP_ ID
        protein_ids.append(protein_id)

        # Extract description from full header
        description = record.description
        if ' ' in description:
            # Skip the IDs at the start, get the description part
            parts = description.split(' ', 2)  # Split into 3 parts max
            if len(parts) >= 3:
                desc = parts[2]  # The description part
                # Remove species name in brackets
                if '[' in desc:
                    desc = desc.split('[')[0].strip()
                protein_names.append(desc)
            else:
                protein_names.append(protein_id)
        else:
            protein_names.append(protein_id)

    return protein_ids, protein_names

def create_locus_matrix(locus_id, locus_info, blocks_df, targets_df, swissprot_df, flanking_df, species_map, phylo_order_map, protein_names):
    """Create matrix for one locus."""

    print(f"\n  Processing {locus_id}...")

    # Get locus details
    gene_family = locus_info['gene_family']
    flanking_file = Path(locus_info['flanking_file'])

    # Parse flanking proteins from .faa file
    flanking_protein_ids, flanking_protein_names = parse_flanking_proteins_from_faa(flanking_file)

    # For compatibility, treat all as "upstream" (no distinction in this pipeline)
    upstream_proteins = flanking_protein_ids
    downstream_proteins = []
    upstream_names = flanking_protein_names
    downstream_names = []
    total_expected = len(flanking_protein_ids)

    # Get blocks for this locus
    locus_blocks = blocks_df[blocks_df['locus_id'] == locus_id]

    # Index SwissProt annotations
    swissprot_map = {}
    if not swissprot_df.empty:
        for _, row in swissprot_df.iterrows():
            if row['locus'] == locus_id:
                key = (row['genome'], row['locus'], row['bk_protein'])
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

    # Index target genes
    targets_index = {}
    if not targets_df.empty:
        for _, target in targets_df.iterrows():
            if target.get('parent_locus') == locus_id or target.get('description') == locus_id:
                key = target['genome']
                if key not in targets_index:
                    targets_index[key] = []
                targets_index[key].append(target)

    # Build matrix rows
    matrix_rows = []

    # Get all other genomes (from blocks or from master list)
    all_genomes = set()

    # Add genomes with blocks
    all_genomes.update(locus_blocks['genome'].unique())

    # Add all genomes from species map to ensure complete coverage
    all_genomes.update(species_map.keys())

    for genome in all_genomes:
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
            row['strand'] = block['strand']
            row['start'] = block['start']
            row['end'] = block['end']

            # Get flanking protein details if available
            proteins_found = set()
            if not flanking_df.empty:
                block_hits = flanking_df[
                    (flanking_df['genome'] == genome) &
                    (flanking_df['sseqid'] == block['scaffold']) &
                    (flanking_df['coord_end'] >= block['start']) &
                    (flanking_df['coord_start'] <= block['end']) &
                    (flanking_df['strand'] == block['strand'])
                ]

                # Map protein IDs to their positions (U1, U2, D1, D2, etc)
                # Query IDs in BLAST are like "metacluster211_chr9a_U1"
                for _, hit in block_hits.iterrows():
                    query_id = hit['qseqid']

                    # Extract the position from the query ID
                    # Format is like: "metacluster211_chr9a_U1" or "metacluster211_chr9a_D1"
                    if '_U' in query_id:
                        # Extract position number
                        position_num = query_id.split('_U')[-1]
                        position_id = f"U{position_num}"
                        proteins_found.add(position_id)

                        # Also add the actual protein ID for SwissProt lookup
                        try:
                            idx = int(position_num) - 1
                            if 0 <= idx < len(upstream_proteins):
                                proteins_found.add(upstream_proteins[len(upstream_proteins) - idx - 1])
                        except (ValueError, IndexError):
                            pass

                    elif '_D' in query_id:
                        # Extract position number
                        position_num = query_id.split('_D')[-1]
                        position_id = f"D{position_num}"
                        proteins_found.add(position_id)

                        # Also add the actual protein ID for SwissProt lookup
                        try:
                            idx = int(position_num) - 1
                            if 0 <= idx < len(downstream_proteins):
                                proteins_found.add(downstream_proteins[idx])
                        except (ValueError, IndexError):
                            pass
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
            # Format: "gene_family [301F; 305P]"
            parts = []
            for target in genome_targets:
                length = int(target.get('orf_length_aa', 0))
                status = target.get('status', 'unknown')
                status_letter = 'F' if status == 'functional' else 'P'
                parts.append(f"{length}{status_letter}")

            target_str = f"{gene_family} [{'; '.join(parts)}]"
            row['TARGET'] = target_str
        else:
            if not genome_block.empty:
                row['TARGET'] = f"{gene_family} [not found]"
            else:
                row['TARGET'] = ""  # No block, no target expected

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

        matrix_rows.append(row)

    # Create DataFrame and sort
    if matrix_rows:
        matrix_df = pd.DataFrame(matrix_rows)
        matrix_df = matrix_df.sort_values('phylo_order', ascending=True)
        return matrix_df

    return pd.DataFrame()  # Empty if no data

def main():
    """Main execution function."""
    print("=" * 80)
    print("STEP 08a: GENERATE LOCUS-SPECIFIC MATRICES")
    print("=" * 80)

    # Create output directory
    config.LOCUS_MATRICES_DIR.mkdir(exist_ok=True, parents=True)

    # Load all required data
    print("\n[1] Loading data...")

    # Locus definitions
    loci_df = pd.read_csv(config.LOCI_DEFINITIONS_FILE, sep='\t')
    print(f"  Loaded {len(loci_df)} locus definitions")

    # Filtered synteny blocks
    blocks_file = config.STEP03_FILTERED / "synteny_blocks_filtered.tsv"
    if blocks_file.exists():
        blocks_df = pd.read_csv(blocks_file, sep='\t')
        print(f"  Loaded {len(blocks_df)} synteny blocks")
    else:
        blocks_df = pd.DataFrame()
        print("  No synteny blocks found")

    # Classified targets
    targets_file = config.STEP05_CLASSIFIED / "all_targets_classified.tsv"
    if targets_file.exists():
        targets_df = pd.read_csv(targets_file, sep='\t')
        # Filter to syntenic only for main column
        syntenic_targets = targets_df[targets_df['placement'] == 'synteny']
        print(f"  Loaded {len(syntenic_targets)} syntenic targets")
    else:
        syntenic_targets = pd.DataFrame()
        print("  No classified targets found")

    # SwissProt annotations
    swissprot_file = config.STEP07_SWISSPROT / "genome_specific_swissprot_annotations.tsv"
    if swissprot_file.exists():
        swissprot_df = pd.read_csv(swissprot_file, sep='\t')
        print(f"  Loaded {len(swissprot_df)} SwissProt annotations")
    else:
        swissprot_df = pd.DataFrame()
        print("  No SwissProt annotations found")

    # Species mapping
    species_map, phylo_order_map = load_species_and_phylo()
    print(f"  Loaded species mapping for {len(species_map)} genomes")

    # Get protein names for better column headers
    protein_names = get_protein_names_from_faa()
    print(f"  Loaded protein names for {len(protein_names)} proteins")

    # Process each locus
    print("\n[2] Generating matrices...")

    for _, locus_info in loci_df.iterrows():
        locus_id = locus_info['locus_id']

        # Load flanking protein details if available
        flanking_file = config.STEP02_SYNTENY / locus_id / "flanking_blast_all.tsv"
        if flanking_file.exists():
            flanking_df = pd.read_csv(flanking_file, sep='\t')
        else:
            flanking_df = pd.DataFrame()

        # Create matrix
        matrix_df = create_locus_matrix(
            locus_id, locus_info, blocks_df, syntenic_targets,
            swissprot_df, flanking_df, species_map, phylo_order_map, protein_names
        )

        if not matrix_df.empty:
            # Save matrix
            output_file = config.LOCUS_MATRICES_DIR / f"{locus_id}_genome_swissprot_matrix.tsv"
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
    print(f"\nMatrices saved to: {config.LOCUS_MATRICES_DIR}")
    print("\nNext step: 07b_generate_summary_matrices.py")

if __name__ == "__main__":
    main()