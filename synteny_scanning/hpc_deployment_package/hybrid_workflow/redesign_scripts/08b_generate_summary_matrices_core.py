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
        --output-dir <path/to/gene_type_summaries> \
        [--locus-matrices-dir <path/to/phase8a_v2>] \
        [--exclude-genomes GCA_010883055.1]
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

    Returns dict keyed by (genome_id, locus_name):
      {(genome, locus): {'length_aa': int, 'status': str}}

    Rationale: locus names repeat across genomes (e.g., BK_chr6_a), so using
    only the locus name as a key will overwrite entries from different genomes
    and produce incorrect, repeated values in summaries.
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

        # Traverse locus directories
        for locus_dir in genome_dir.iterdir():
            if not locus_dir.is_dir():
                continue

            # Handle UNPLACED directory (contains parent_locus/block subdirectories)
            if locus_dir.name == 'UNPLACED':
                for parent_locus_dir in locus_dir.iterdir():
                    if not parent_locus_dir.is_dir():
                        continue

                    locus_name = parent_locus_dir.name

                    # UNPLACED sequences are nested one level deeper in block directories
                    # Try direct files first (for backwards compatibility)
                    cds_files = list(parent_locus_dir.glob("*_gene1_cds.fasta"))
                    protein_files = list(parent_locus_dir.glob("*_gene1_protein.fasta"))

                    if not cds_files or not protein_files:
                        # Look in block subdirectories
                        for block_dir in parent_locus_dir.iterdir():
                            if not block_dir.is_dir():
                                continue
                            cds_files = list(block_dir.glob("*_gene1_cds.fasta"))
                            protein_files = list(block_dir.glob("*_gene1_protein.fasta"))
                            if cds_files and protein_files:
                                break

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

                metadata[(genome_id, locus_name)] = {
                    'length_aa': length_aa,
                    'status': status
                }

    return metadata

def create_gene_family_matrix(gene_family, loci_list, all_targets_df, synteny_blocks_df, species_map, phylo_order_map, seq_metadata, locus_matrices_dir, output_dir):
    """Create summary matrix for one gene family."""

    print(f"\n  Processing {gene_family}...")

    # Get all targets for this gene family
    # Fallback: if gene_family column is missing or 'unknown', use ALL targets
    # (we know the family from locus_definitions passed to this function)
    if 'gene_family' in all_targets_df.columns and (all_targets_df['gene_family'] == gene_family).any():
        family_targets = all_targets_df[all_targets_df['gene_family'] == gene_family]
    else:
        # Fallback: Phase 5 didn't set gene_family, so use all targets
        print(f"    Warning: gene_family column missing or all 'unknown', using all {len(all_targets_df)} targets for {gene_family}")
        family_targets = all_targets_df

    # Group targets by genome
    targets_by_genome = defaultdict(list)
    for _, target in family_targets.iterrows():
        targets_by_genome[target['genome']].append(target)

    # Helper: deduplicate targets per genome/category/scaffold by merging overlapping/nearby intervals
    # Keep the best (lowest) e-value per cluster.
    def dedup_targets(rows, unplaceable_col_name):
        MERGE_GAP_BP = 200  # merge hits that overlap or are within 200 bp
        dedup = []
        # Build sortable tuples: (category, scaffold, start, end, evalue, row)
        tuples = []
        for t in rows:
            assigned_to = t.get('assigned_to', t.get('description', 'unknown'))
            placement = t.get('placement', '')
            category = unplaceable_col_name if placement == 'unplaceable' else assigned_to
            try:
                start = int(t.get('start', 0))
                end = int(t.get('end', 0))
            except Exception:
                start = 0; end = 0
            if start and end and start > end:
                start, end = end, start
            scaffold = t.get('scaffold', '')
            try:
                evalue = float(t.get('best_evalue', 'inf'))
            except Exception:
                evalue = float('inf')
            tuples.append((category, scaffold, start, end, evalue, t))

        # Sort for linear clustering
        tuples.sort(key=lambda x: (x[0], x[1], x[2], x[3]))

        # Sweep and cluster per (category, scaffold)
        current_cat = None
        current_scaf = None
        cur_start = None
        cur_end = None
        cur_best = None  # (evalue, row)

        def flush():
            nonlocal cur_start, cur_end, cur_best
            if cur_best is not None:
                dedup.append(cur_best[1])
            cur_start = None
            cur_end = None
            cur_best = None

        for cat, scaf, s, e, ev, row in tuples:
            if (cat != current_cat) or (scaf != current_scaf) or (cur_start is None):
                # new group
                flush()
                current_cat = cat
                current_scaf = scaf
                cur_start, cur_end = s, e
                cur_best = (ev, row)
                continue
            # same (category, scaffold): check overlap or proximity
            if s <= (cur_end + MERGE_GAP_BP):
                # extend cluster
                if e > cur_end:
                    cur_end = e
                # keep best (lowest) e-value
                if ev < (cur_best[0] if cur_best else float('inf')):
                    cur_best = (ev, row)
            else:
                # gap: commit cluster and start a new one
                flush()
                cur_start, cur_end = s, e
                cur_best = (ev, row)

        flush()
        return dedup

    # Read synteny percentages from Phase 8a locus matrices
    # Phase 8a already calculated correct synteny_pct based on flanking columns
    synteny_by_genome_locus = {}  # {genome: {locus: synteny_pct}}

    for locus in loci_list:
        # Load the locus matrix created by Phase 8a
        locus_matrix_file = locus_matrices_dir / f"{locus}_genome_swissprot_matrix.tsv"

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

    # Categories: syntenic loci (from loci_list) + a single collapsed unplaceable category
    unplaceable_col = f"{gene_family}_unplaceable"
    locus_categories = loci_list + [unplaceable_col]

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
        # Deduplicate per genome before summarizing
        genome_targets_all = targets_by_genome.get(genome, [])
        genome_targets = dedup_targets(genome_targets_all, unplaceable_col)

        # Initialize all columns
        for category in locus_categories:
            if category.endswith('_unplaceable'):
                row[category] = ""  # Unplaceable columns don't have synteny percentages
            else:
                row[category] = ""  # Empty by default; filled in if synteny block exists

        # Fill in target counts
        category_counts = defaultdict(list)
        for target in genome_targets:
            assigned_to = target.get('assigned_to', target.get('description', 'unknown'))
            locus_name = target.get('locus_name', target.get('parent_locus', 'unknown'))

            # For unplaceable targets, collapse to single category
            placement = target.get('placement', '')
            confidence = target.get('confidence', '')
            if placement == 'unplaceable':
                assigned_to = unplaceable_col

            # Get metadata from extracted sequences
            if placement == 'unplaceable':
                # For unplaceables, construct unique_tag: scaffold_start_end
                scaf = str(target.get('scaffold'))
                start = str(target.get('start'))
                end = str(target.get('end'))
                unique_tag = f"{scaf}_{start}_{end}"
                lookup_key = (target.get('genome'), unique_tag)
                meta = seq_metadata.get(lookup_key, {})
                length = meta.get('length_aa', 0)
                status = meta.get('status', 'unknown')
            else:
                # Look up extracted sequence metadata strictly by (genome, assigned locus)
                # Do not fall back to parent_locus; if missing, record as generic hit.
                lookup_key = (target.get('genome'), assigned_to)
                meta = seq_metadata.get(lookup_key, {})
                length = meta.get('length_aa', 0)
                status = meta.get('status', 'unknown')

            # Skip targets that failed extraction (both syntenic and unplaceable)
            # Only keep as generic "hit" if no extraction attempted (length=0 and no metadata)
            if length == 0 and meta:
                # Extraction was attempted but failed
                continue

            # Status letters: I=intact, P=pseudogene, F=fragment
            if length > 0:
                if status == 'intact':
                    status_letter = 'I'
                elif status == 'pseudogene':
                    status_letter = 'P'
                elif status == 'fragment':
                    status_letter = 'F'
                else:
                    status_letter = '?'
                target_str = f"{length}{status_letter}"
            else:
                # Unplaceable or failed extraction: record as generic hit
                target_str = "hit"

            category_counts[assigned_to].append(target_str)

        # Format counts for each category with synteny percent
        for category, target_list in category_counts.items():
            if category in row:  # Only if it's a known category
                if target_list:
                    # Get synteny percent if available
                    if (category in loci_list) and genome in synteny_by_genome_locus and category in synteny_by_genome_locus[genome]:
                        synteny_pct = synteny_by_genome_locus[genome][category]
                        row[category] = f"{synteny_pct}% [{'; '.join(target_list)}]"
                    else:
                        row[category] = f"[{'; '.join(target_list)}]"

        # Check for synteny blocks without targets
        # For loci where this genome has a synteny block but no extracted target
        if genome in synteny_by_genome_locus:
            for locus, synteny_pct in synteny_by_genome_locus[genome].items():
                if locus in row and row[locus] == "":
                    # Synteny block exists but no target found/extracted
                    if synteny_pct > 0:
                        row[locus] = f"{synteny_pct}% [empty]"
                    else:
                        # 0% synteny = no syntenic block found at this locus
                        row[locus] = "[not found]"

        # For loci not even attempted (genome not in synteny search), mark as not found
        for locus in loci_list:
            if locus in row and row[locus] == "":
                row[locus] = "[not found]"

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
    parser.add_argument('--locus-matrices-dir', type=Path,
                        help='Directory containing per-locus matrices from Phase 8a (defaults to output-dir)')
    parser.add_argument('--exclude-genomes', type=str, default='GCA_010883055.1',
                        help='Comma-separated genome IDs to exclude (e.g., deprecated BK accession)')
    parser.add_argument('--exclude-genomes-file', type=Path,
                        help='File with genome IDs to exclude (one per line, # comments)')

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
    # Exclude deprecated/renamed genomes
    exclude_set: set[str] = set()
    if args.exclude_genomes:
        for g in [x.strip() for x in args.exclude_genomes.split(',') if x.strip()]:
            exclude_set.add(g)
    # Also read from exclusion file if provided
    if args.exclude_genomes_file and args.exclude_genomes_file.exists():
        with open(args.exclude_genomes_file) as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    exclude_set.add(line)
    # Apply exclusions
    for g in exclude_set:
        species_map.pop(g, None)
        phylo_order_map.pop(g, None)
    if exclude_set:
        print(f"  Excluded {len(exclude_set)} genomes")
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
        locus_matrices_dir = args.locus_matrices_dir if args.locus_matrices_dir else args.output_dir
        matrix_df = create_gene_family_matrix(
            gene_family, family_loci, all_targets_df, synteny_blocks_df,
            species_map, phylo_order_map, seq_metadata, locus_matrices_dir, args.output_dir
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

            # Count syntenic vs unplaceable (collapsed)
            syntenic_count = 0
            unplaceable_count = 0

            for col in family_loci:
                syntenic_count += len(matrix_df[matrix_df[col] != ""])

            unp_col = f"{gene_family}_unplaceable"
            if unp_col in matrix_df.columns:
                unplaceable_count = len(matrix_df[matrix_df[unp_col] != ""]) 

            total_count = syntenic_count + unplaceable_count
            if total_count > 0:
                print(f"    Syntenic: {syntenic_count} ({syntenic_count/total_count*100:.1f}%)")
                print(f"    Unplaceable: {unplaceable_count} ({unplaceable_count/total_count*100:.1f}%)")

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
