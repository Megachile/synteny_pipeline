#!/usr/bin/env python3
"""
Step 05: Classify target loci as syntenic or unplaceable.

Checks if target loci overlap with synteny blocks to determine placement.

Usage:
    python 05_classify_targets.py \\
        --targets <path/to/all_target_loci.tsv> \\
        --blocks <path/to/synteny_blocks_filtered.tsv> \\
        --output-dir <path/to/classified_targets> \\
        [--min-length 30]
"""

from pathlib import Path
import pandas as pd
from collections import defaultdict
from Bio import SeqIO
import argparse

def check_proximity(target, block, max_distance_kb=500):
    """Check if target gene is near synteny block (within max_distance_kb).

    Synteny blocks are defined by flanking proteins, so the target gene
    should be NEAR the block, not necessarily overlapping it.
    """

    # Must be same genome and scaffold
    if target['genome'] != block['genome'] or target['scaffold'] != block['scaffold']:
        return False

    # Check if target is within max_distance_kb of the block
    # Distance = minimum distance between target and block boundaries
    max_distance_bp = max_distance_kb * 1000

    # If they overlap, distance = 0
    if not (target['end'] < block['start'] or target['start'] > block['end']):
        return True

    # Calculate distance to nearest edge
    if target['end'] < block['start']:
        distance = block['start'] - target['end']
    else:
        distance = target['start'] - block['end']

    return distance <= max_distance_bp

def assess_target_quality(hsp_seq, ref_length, pident, min_length):
    """
    Assess target quality from HSP sequence.

    Args:
        hsp_seq: HSP protein sequence (with gaps)
        ref_length: Reference protein length
        pident: Percent identity
        min_length: Minimum ORF length to consider

    Returns:
        (status, quality_flag, has_stops) tuple
    """
    # Remove gaps
    clean_seq = hsp_seq.replace('-', '') if hsp_seq else ""

    if not clean_seq:
        return "uncertain", "no_sequence", False

    hsp_length = len(clean_seq)
    coverage = (hsp_length / ref_length * 100) if ref_length > 0 else 0

    # Flag 1: Internal stop codons
    has_stops = '*' in clean_seq

    # Flag 2: Severely truncated (< 50% of reference)
    severely_truncated = coverage < 50

    # Flag 3: Very low identity
    low_identity = pident < 30

    # Flag 4: Very short absolute length
    very_short = hsp_length < min_length

    # Decision logic
    if has_stops:
        return "pseudogene", "internal_stop_codon", True
    elif severely_truncated and very_short:
        return "pseudogene", "severe_truncation", False
    elif low_identity and very_short:
        return "pseudogene", "low_quality", False
    elif severely_truncated or low_identity:
        return "uncertain", "degraded_sequence", False
    else:
        return "functional", "passes_basic_checks", False

def load_target_sequences(genome_id, parent_locus, output_dir):
    """Load target sequences from FASTA files."""
    seqs_dir = output_dir / "target_sequences"

    if not seqs_dir.exists():
        return {}

    # Find FASTA file for this genome and locus
    fasta_file = seqs_dir / f"{genome_id}_{parent_locus}_target_proteins.fasta"

    if not fasta_file.exists():
        return {}

    # Load sequences
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = str(record.seq)

    return sequences

def classify_targets(targets_df, blocks_df, targets_dir, bk_proteins_file, min_length):
    """Classify each target as syntenic or unplaceable with quality assessment."""

    # Load reference protein sequences to get lengths
    from Bio import SeqIO
    ref_proteins = {}
    if bk_proteins_file.exists():
        for record in SeqIO.parse(bk_proteins_file, "fasta"):
            ref_proteins[record.id] = len(record.seq)
    
    # Index synteny blocks by (genome, scaffold) for fast lookup
    synteny_index = defaultdict(list)
    for _, block in blocks_df.iterrows():
        key = (block['genome'], block['scaffold'])
        synteny_index[key].append(block)

    classifications = []

    for _, target in targets_df.iterrows():
        # Check for overlap with synteny blocks
        key = (target['genome'], target['scaffold'])
        candidate_blocks = synteny_index.get(key, [])

        assigned_to = None
        for block in candidate_blocks:
            # Check proximity - accept ANY target hit near synteny block
            # (don't filter by locus_id - BLAST e-values can be fickle)
            if check_proximity(target, block):
                assigned_to = block['locus_id']
                break

        if assigned_to:
            classification = {
                **target.to_dict(),
                'placement': 'synteny',
                'assigned_to': assigned_to,
                'description': assigned_to  # For compatibility with matrix script
            }
        else:
            classification = {
                **target.to_dict(),
                'placement': 'unplaceable',
                'assigned_to': f"{target['gene_family']}_unplaceable",
                'description': f"{target['gene_family']}_unplaceable"
            }

        # Load target sequences for quality assessment
        genome_sequences = load_target_sequences(
            target['genome'],
            target['parent_locus'],
            targets_dir
        )

        # Get reference protein length
        query_id = target.get('query_id', target.get('gene_id', ''))
        ref_length = ref_proteins.get(query_id, 300)  # Default 300aa if not found

        # Find matching HSP sequence
        hsp_seq = None
        pident = target.get('pident', 0)
        for seq_id, seq in genome_sequences.items():
            if query_id in seq_id or target['gene_family'] in seq_id:
                hsp_seq = seq
                break

        # Assess quality
        if hsp_seq:
            status, quality_flag, has_stops = assess_target_quality(hsp_seq, ref_length, pident, min_length)
            hsp_length = len(hsp_seq.replace('-', ''))
            coverage_pct = (hsp_length / ref_length * 100) if ref_length > 0 else 0
        else:
            # Fallback to coordinate-based estimate
            orf_length_nt = target['end'] - target['start'] + 1
            hsp_length = orf_length_nt // 3
            status = 'uncertain'
            quality_flag = 'no_hsp_sequence'
            has_stops = False
            coverage_pct = 0
        
        classification['orf_length_aa'] = hsp_length
        classification['status'] = status
        classification['quality_flag'] = quality_flag
        classification['has_stop_codons'] = has_stops
        classification['ref_length'] = ref_length
        classification['coverage_pct'] = round(coverage_pct, 1)
        classification['hsp_sequence'] = hsp_seq if hsp_seq else ""
        classification['genome_id'] = target['genome']  # Alias for matrix script

        classifications.append(classification)

    return pd.DataFrame(classifications)

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Classify target loci as syntenic or unplaceable",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('--targets', required=True, type=Path,
                        help='Path to all_target_loci.tsv from Phase 4')
    parser.add_argument('--blocks', required=True, type=Path,
                        help='Path to synteny_blocks_filtered.tsv from Phase 3')
    parser.add_argument('--output-dir', required=True, type=Path,
                        help='Output directory for classified targets')
    parser.add_argument('--min-length', type=int, default=30,
                        help='Minimum ORF length in amino acids (default: 30)')

    return parser.parse_args()

def main():
    """Main execution function."""
    args = parse_args()

    print("=" * 80, flush=True)
    print("STEP 05: CLASSIFY TARGET LOCI", flush=True)
    print("=" * 80, flush=True)
    print(f"\nInput files:")
    print(f"  Target loci: {args.targets}")
    print(f"  Synteny blocks: {args.blocks}")
    print(f"  Output directory: {args.output_dir}")
    print(f"\nParameters:")
    print(f"  Minimum ORF length: {args.min_length} aa")

    # Create output directory
    args.output_dir.mkdir(exist_ok=True, parents=True)

    # Load synteny blocks
    print("\n[1] Loading synteny blocks...", flush=True)

    if not args.blocks.exists():
        print(f"  ERROR: Filtered blocks not found: {args.blocks}", flush=True)
        print("  Run Step 03 first.", flush=True)
        return

    blocks_df = pd.read_csv(args.blocks, sep='\t')
    print(f"  Loaded {len(blocks_df)} synteny blocks", flush=True)

    # Load target loci
    print("\n[2] Loading target loci...", flush=True)

    if not args.targets.exists():
        print(f"  ERROR: Target loci not found: {args.targets}", flush=True)
        print("  Run Step 04 first.", flush=True)
        return

    targets_df = pd.read_csv(args.targets, sep='\t')
    print(f"  Loaded {len(targets_df)} target loci", flush=True)

    # Infer paths for reference proteins (in same directory as targets file)
    targets_dir = args.targets.parent
    bk_proteins_file = targets_dir.parent / "01_extracted_proteins" / "BK_proteins.faa"

    # Classify targets
    print("\n[3] Classifying targets...", flush=True)
    classified_df = classify_targets(targets_df, blocks_df, targets_dir, bk_proteins_file, args.min_length)

    syntenic = classified_df[classified_df['placement'] == 'synteny']
    unplaceable = classified_df[classified_df['placement'] == 'unplaceable']

    print(f"  Syntenic: {len(syntenic)} ({len(syntenic)/len(classified_df)*100:.1f}%)")
    print(f"  Unplaceable: {len(unplaceable)} ({len(unplaceable)/len(classified_df)*100:.1f}%)")

    # Save classifications
    print("\n[4] Saving classifications...", flush=True)

    # Save all classifications
    all_file = args.output_dir / "all_targets_classified.tsv"
    classified_df.to_csv(all_file, sep='\t', index=False)
    print(f"  Saved all: {all_file.name}", flush=True)

    # Save syntenic only (for matrix generation)
    syntenic_file = args.output_dir / "syntenic_targets.tsv"
    syntenic.to_csv(syntenic_file, sep='\t', index=False)
    print(f"  Saved syntenic: {syntenic_file.name}", flush=True)

    # Save unplaceable
    unplaceable_file = args.output_dir / "unplaceable_targets.tsv"
    unplaceable.to_csv(unplaceable_file, sep='\t', index=False)
    print(f"  Saved unplaceable: {unplaceable_file.name}", flush=True)

    # Summary by locus
    print("\n[5] Summary by locus:")
    for parent_locus in classified_df['parent_locus'].unique():
        locus_targets = classified_df[classified_df['parent_locus'] == parent_locus]
        locus_syntenic = locus_targets[locus_targets['placement'] == 'synteny']
        locus_unplaceable = locus_targets[locus_targets['placement'] == 'unplaceable']

        print(f"\n  {parent_locus}:")
        print(f"    Total targets: {len(locus_targets)}")
        print(f"    Syntenic: {len(locus_syntenic)} ({len(locus_syntenic)/len(locus_targets)*100:.1f}%)")
        print(f"    Unplaceable: {len(locus_unplaceable)} ({len(locus_unplaceable)/len(locus_targets)*100:.1f}%)")

        # Functional status
        functional = locus_targets[locus_targets['status'] == 'functional']
        pseudo = locus_targets[locus_targets['status'] == 'pseudogene']
        print(f"    Functional: {len(functional)} ({len(functional)/len(locus_targets)*100:.1f}%)")
        print(f"    Pseudogenes: {len(pseudo)} ({len(pseudo)/len(locus_targets)*100:.1f}%)")

    # Overall summary
    print("\n" + "=" * 80, flush=True)
    print("CLASSIFICATION COMPLETE", flush=True)
    print("=" * 80, flush=True)
    print(f"\nTotal targets classified: {len(classified_df)}", flush=True)
    print(f"Syntenic: {len(syntenic)} ({len(syntenic)/len(classified_df)*100:.1f}%)", flush=True)
    print(f"Unplaceable: {len(unplaceable)} ({len(unplaceable)/len(classified_df)*100:.1f}%)", flush=True)

    functional = classified_df[classified_df['status'] == 'functional']
    print(f"Functional: {len(functional)} ({len(functional)/len(classified_df)*100:.1f}%)", flush=True)

    print(f"\nOutputs saved to: {args.output_dir}", flush=True)
    print("\nNext step: 06_extract_sequences.py", flush=True)

if __name__ == "__main__":
    main()