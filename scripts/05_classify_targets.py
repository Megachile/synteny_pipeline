#!/usr/bin/env python3
"""
Step 05: Classify target loci as syntenic or unplaceable.

Checks if target loci overlap with synteny blocks to determine placement.
"""

from pathlib import Path
import pandas as pd
from collections import defaultdict
from Bio import SeqIO
import config

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

def assess_target_quality(hsp_seq, ref_length, pident):
    """
    Assess target quality from HSP sequence.
    
    Args:
        hsp_seq: HSP protein sequence (with gaps)
        ref_length: Reference protein length
        pident: Percent identity
    
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
    very_short = hsp_length < 50
    
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

def classify_targets(targets_df, blocks_df):
    """Classify each target as syntenic or unplaceable with quality assessment."""

    # Load reference protein sequences to get lengths
    from Bio import SeqIO
    ref_proteins = {}
    ref_file = config.BK_PROTEINS_FILE
    if ref_file.exists():
        for record in SeqIO.parse(ref_file, "fasta"):
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
            config.STEP04_TARGETS
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
            status, quality_flag, has_stops = assess_target_quality(hsp_seq, ref_length, pident)
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

def main():
    """Main execution function."""
    print("=" * 80)
    print("STEP 05: CLASSIFY TARGET LOCI")
    print("=" * 80)

    # Create output directory
    config.STEP05_CLASSIFIED.mkdir(exist_ok=True, parents=True)

    # Load synteny blocks
    print("\n[1] Loading synteny blocks...")
    blocks_file = config.STEP03_FILTERED / "synteny_blocks_filtered.tsv"

    if not blocks_file.exists():
        print(f"  ERROR: Filtered blocks not found: {blocks_file}")
        print("  Run Step 03 first.")
        return

    blocks_df = pd.read_csv(blocks_file, sep='\t')
    print(f"  Loaded {len(blocks_df)} synteny blocks")

    # Load target loci
    print("\n[2] Loading target loci...")
    targets_file = config.STEP04_TARGETS / "all_target_loci.tsv"

    if not targets_file.exists():
        print(f"  ERROR: Target loci not found: {targets_file}")
        print("  Run Step 04 first.")
        return

    targets_df = pd.read_csv(targets_file, sep='\t')
    print(f"  Loaded {len(targets_df)} target loci")

    # Classify targets
    print("\n[3] Classifying targets...")
    classified_df = classify_targets(targets_df, blocks_df)

    syntenic = classified_df[classified_df['placement'] == 'synteny']
    unplaceable = classified_df[classified_df['placement'] == 'unplaceable']

    print(f"  Syntenic: {len(syntenic)} ({len(syntenic)/len(classified_df)*100:.1f}%)")
    print(f"  Unplaceable: {len(unplaceable)} ({len(unplaceable)/len(classified_df)*100:.1f}%)")

    # Save classifications
    print("\n[4] Saving classifications...")

    # Save all classifications
    all_file = config.STEP05_CLASSIFIED / "all_targets_classified.tsv"
    classified_df.to_csv(all_file, sep='\t', index=False)
    print(f"  Saved all: {all_file.name}")

    # Save syntenic only (for matrix generation)
    syntenic_file = config.STEP05_CLASSIFIED / "syntenic_targets.tsv"
    syntenic.to_csv(syntenic_file, sep='\t', index=False)
    print(f"  Saved syntenic: {syntenic_file.name}")

    # Save unplaceable
    unplaceable_file = config.STEP05_CLASSIFIED / "unplaceable_targets.tsv"
    unplaceable.to_csv(unplaceable_file, sep='\t', index=False)
    print(f"  Saved unplaceable: {unplaceable_file.name}")

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
    print("\n" + "=" * 80)
    print("CLASSIFICATION COMPLETE")
    print("=" * 80)
    print(f"\nTotal targets classified: {len(classified_df)}")
    print(f"Syntenic: {len(syntenic)} ({len(syntenic)/len(classified_df)*100:.1f}%)")
    print(f"Unplaceable: {len(unplaceable)} ({len(unplaceable)/len(classified_df)*100:.1f}%)")

    functional = classified_df[classified_df['status'] == 'functional']
    print(f"Functional: {len(functional)} ({len(functional)/len(classified_df)*100:.1f}%)")

    print(f"\nOutputs saved to: {config.STEP05_CLASSIFIED}")
    print("\nNext step: 06_swissprot_annotation.py")

if __name__ == "__main__":
    main()