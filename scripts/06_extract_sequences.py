#!/usr/bin/env python3
"""
Step 06: Extract sequences for filtered targets using Exonerate.

Reads filtered target list (from Phase 5) and extracts gene structures.

Usage:
    python 06_extract_sequences.py \\
        --syntenic <path/to/syntenic_targets.tsv> \\
        --unplaceable <path/to/unplaceable_targets.tsv> \\
        --query-proteins <path/to/combined_targets.faa> \\
        --genome-fasta-dir <path/to/ragtag_output> \\
        --output-dir <path/to/extracted_sequences> \\
        [--unplaceable-evalue 1e-10] \\
        [--verbose]
"""

from pathlib import Path
import pandas as pd
import sys
import argparse
sys.path.insert(0, str(Path(__file__).parent))
import extract_with_exonerate

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Extract sequences for filtered targets using Exonerate",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('--syntenic', required=True, type=Path,
                        help='Path to syntenic_targets.tsv from Phase 5')
    parser.add_argument('--unplaceable', type=Path,
                        help='Path to unplaceable_targets.tsv from Phase 5 (optional)')
    parser.add_argument('--query-proteins', required=True, type=Path,
                        help='Path to combined_targets.faa with query proteins')
    parser.add_argument('--genome-fasta-dir', required=True, type=Path,
                        help='Directory containing genome FASTA files (e.g., ragtag_output/*/ragtag.scaffold.fasta)')
    parser.add_argument('--output-dir', required=True, type=Path,
                        help='Output directory for extracted sequences')
    parser.add_argument('--unplaceable-evalue', type=float, default=1e-10,
                        help='E-value threshold for unplaceable targets (default: 1e-10)')
    parser.add_argument('--verbose', action='store_true',
                        help='Print verbose extraction progress')

    return parser.parse_args()

def main():
    """Extract sequences for filtered targets."""
    args = parse_args()

    print("=" * 80, flush=True)
    print("STEP 06: EXTRACT SEQUENCES WITH EXONERATE", flush=True)
    print("=" * 80, flush=True)
    print(f"\nInput files:")
    print(f"  Syntenic targets: {args.syntenic}")
    print(f"  Unplaceable targets: {args.unplaceable if args.unplaceable else 'None'}")
    print(f"  Query proteins: {args.query_proteins}")
    print(f"  Genome FASTA dir: {args.genome_fasta_dir}")
    print(f"  Output directory: {args.output_dir}")
    print(f"\nParameters:")
    print(f"  Unplaceable e-value threshold: {args.unplaceable_evalue}")

    # Load filtered targets (syntenic + strong unplaceable)
    print("\n[1] Loading filtered targets...", flush=True)

    if not args.syntenic.exists():
        print(f"  ERROR: Syntenic targets not found: {args.syntenic}", flush=True)
        print("  Run Step 05 (classify_targets.py) first.", flush=True)
        return

    # Load syntenic targets (keep all - validated by synteny)
    syntenic_df = pd.read_csv(args.syntenic, sep='\t')
    print(f"  Loaded {len(syntenic_df)} syntenic targets", flush=True)

    # Load and filter strong unplaceable targets
    if args.unplaceable and args.unplaceable.exists():
        unplaceable_df = pd.read_csv(args.unplaceable, sep='\t')
        print(f"  Loaded {len(unplaceable_df)} unplaceable targets", flush=True)

        # Apply strict e-value filter for unplaceable targets
        strong_unplaceable = unplaceable_df[unplaceable_df['best_evalue'] < args.unplaceable_evalue]
        print(f"  Kept {len(strong_unplaceable)} strong unplaceable targets (e-value < {args.unplaceable_evalue})", flush=True)
        print(f"  Filtered out {len(unplaceable_df) - len(strong_unplaceable)} weak unplaceable targets", flush=True)

        # Combine syntenic + strong unplaceable
        targets_df = pd.concat([syntenic_df, strong_unplaceable], ignore_index=True)
    else:
        targets_df = syntenic_df

    print(f"  Total targets for extraction: {len(targets_df)}", flush=True)

    # Group by genome for efficient processing
    print("\n[2] Grouping targets by genome...", flush=True)
    genome_groups = targets_df.groupby('genome')
    print(f"  {len(genome_groups)} genomes to process", flush=True)

    # Output directory
    args.output_dir.mkdir(exist_ok=True, parents=True)

    # Process each genome
    print("\n[3] Extracting sequences with Exonerate...", flush=True)
    total_extracted = 0
    failed_genomes = []

    for genome_name, genome_targets in genome_groups:
        print(f"\n  {genome_name}: {len(genome_targets)} targets", flush=True)

        # Find genome FASTA
        genome_fasta = args.genome_fasta_dir / genome_name / "ragtag.scaffold.fasta"
        if not genome_fasta.exists():
            print(f"    WARNING: Genome FASTA not found, skipping")
            failed_genomes.append(genome_name)
            continue

        # Group targets by locus for batch extraction
        for parent_locus, locus_targets in genome_targets.groupby('parent_locus'):
            # Process each target
            for _, target in locus_targets.iterrows():
                target_name = target['locus_name']
                expected_query = target['query_id']

                # Create pseudo-block structure for Exonerate
                block = {
                    'block_id': target_name,
                    'scaffold': target['scaffold'],
                    'start': target['start'],
                    'end': target['end'],
                    'strand': target['strand']
                }

                # Output directory for this target
                target_output_dir = args.output_dir / genome_name / target_name
                target_output_dir.mkdir(exist_ok=True, parents=True)

                # Extract SPECIFIC query protein to temp file
                temp_query = target_output_dir / f"query_{expected_query}.faa"
                found_query = False

                from Bio import SeqIO
                with open(args.query_proteins, 'r') as infile:
                    for record in SeqIO.parse(infile, 'fasta'):
                        if expected_query in record.id:
                            with open(temp_query, 'w') as outfile:
                                SeqIO.write(record, outfile, 'fasta')
                            found_query = True
                            break

                if not found_query:
                    print(f"    WARNING: Query protein {expected_query} not found in {args.query_proteins}")
                    continue

                # Extract with Exonerate using specific query
                try:
                    genes = extract_with_exonerate.extract_block_genes(
                        block=block,
                        query_protein_file=temp_query,
                        genome_fasta=genome_fasta,
                        output_dir=target_output_dir
                    )

                    if genes:
                        total_extracted += len(genes)
                        if args.verbose:
                            print(f"    {target_name}: extracted {len(genes)} genes")
                    else:
                        if args.verbose:
                            print(f"    {target_name}: no genes extracted")

                except Exception as e:
                    print(f"    ERROR extracting {target_name}: {e}")

    # Summary
    print("\n" + "=" * 80)
    print("SEQUENCE EXTRACTION COMPLETE")
    print("=" * 80)
    print(f"\nTotal genes extracted: {total_extracted}")
    print(f"Genomes processed: {len(genome_groups) - len(failed_genomes)}/{len(genome_groups)}")
    if failed_genomes:
        print(f"\nFailed genomes (no FASTA): {len(failed_genomes)}")
        for g in failed_genomes[:10]:
            print(f"  - {g}")

    print(f"\nOutput directory: {args.output_dir}")

if __name__ == "__main__":
    main()
