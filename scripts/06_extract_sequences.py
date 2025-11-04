#!/usr/bin/env python3
"""
Step 06: Extract sequences for filtered targets using Exonerate.

Reads filtered target list (from Phase 5) and extracts gene structures.
"""

from pathlib import Path
import pandas as pd
import sys
sys.path.insert(0, str(Path(__file__).parent))
import extract_with_exonerate
import config

def main():
    """Extract sequences for filtered targets."""
    print("=" * 80, flush=True)
    print("STEP 06: EXTRACT SEQUENCES WITH EXONERATE", flush=True)
    print("=" * 80, flush=True)

    # Load filtered targets (syntenic + strong unplaceable)
    print("\n[1] Loading filtered targets...", flush=True)

    # Read from Phase 5 outputs
    syntenic_file = config.STEP05_CLASSIFIED / "syntenic_targets.tsv"
    unplaceable_file = config.STEP05_CLASSIFIED / "unplaceable_targets.tsv"

    if not syntenic_file.exists():
        print(f"  ERROR: Syntenic targets not found: {syntenic_file}", flush=True)
        print("  Run Step 05 (classify_targets.py) first.", flush=True)
        return

    # Load syntenic targets (keep all - validated by synteny)
    syntenic_df = pd.read_csv(syntenic_file, sep='\t')
    print(f"  Loaded {len(syntenic_df)} syntenic targets", flush=True)

    # Load and filter strong unplaceable targets
    if unplaceable_file.exists():
        unplaceable_df = pd.read_csv(unplaceable_file, sep='\t')
        print(f"  Loaded {len(unplaceable_df)} unplaceable targets", flush=True)

        # Apply strict e-value filter for unplaceable targets
        # Syntenic targets are already validated by synteny, but unplaceable
        # targets need very strong BLAST evidence
        UNPLACEABLE_EVALUE_THRESHOLD = 1e-10

        strong_unplaceable = unplaceable_df[unplaceable_df['best_evalue'] < UNPLACEABLE_EVALUE_THRESHOLD]
        print(f"  Kept {len(strong_unplaceable)} strong unplaceable targets (e-value < {UNPLACEABLE_EVALUE_THRESHOLD})", flush=True)
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
    sequences_dir = config.STEP04_TARGETS / "extracted_sequences"
    sequences_dir.mkdir(exist_ok=True, parents=True)

    # Process each genome
    print("\n[3] Extracting sequences with Exonerate...", flush=True)
    total_extracted = 0
    failed_genomes = []

    for genome_name, genome_targets in genome_groups:
        print(f"\n  {genome_name}: {len(genome_targets)} targets", flush=True)

        # Find genome FASTA
        genome_fasta = config.RAGTAG_FASTA_DIR / genome_name / "ragtag.scaffold.fasta"
        if not genome_fasta.exists():
            print(f"    WARNING: Genome FASTA not found, skipping")
            failed_genomes.append(genome_name)
            continue

        # Group targets by locus for batch extraction
        for parent_locus, locus_targets in genome_targets.groupby('parent_locus'):
            # Get query protein for this locus
            locus_id = locus_targets.iloc[0]['parent_locus']
            expected_query = locus_targets.iloc[0]['expected_query']

            # Find query protein file
            query_file = config.STEP04_TARGETS / "combined_targets.faa"
            if not query_file.exists():
                print(f"    ERROR: Query file not found: {query_file}")
                continue

            # Process each target
            for _, target in locus_targets.iterrows():
                target_name = target['locus_name']
                expected_query = target['expected_query']

                # Create pseudo-block structure for Exonerate
                block = {
                    'block_id': target_name,
                    'scaffold': target['scaffold'],
                    'start': target['start'],
                    'end': target['end'],
                    'strand': target['strand']
                }

                # Output directory for this target
                target_output_dir = sequences_dir / genome_name / target_name
                target_output_dir.mkdir(exist_ok=True, parents=True)

                # Extract SPECIFIC query protein to temp file
                temp_query = target_output_dir / f"query_{expected_query}.faa"
                found_query = False

                from Bio import SeqIO
                with open(query_file, 'r') as infile:
                    for record in SeqIO.parse(infile, 'fasta'):
                        if expected_query in record.id:
                            with open(temp_query, 'w') as outfile:
                                SeqIO.write(record, outfile, 'fasta')
                            found_query = True
                            break

                if not found_query:
                    print(f"    WARNING: Query protein {expected_query} not found in {query_file}")
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
                        if config.VERBOSE:
                            print(f"    {target_name}: extracted {len(genes)} genes")
                    else:
                        if config.VERBOSE:
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

    print(f"\nOutput directory: {sequences_dir}")

if __name__ == "__main__":
    main()
