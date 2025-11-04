#!/usr/bin/env python3
"""
Step 05: Classify target loci as syntenic or unplaceable.

Checks if target loci overlap with synteny blocks to determine placement.
Filters out tiny fragments (< min-span-kb) to remove likely spurious hits.

Usage:
    python 05_classify_targets.py \\
        --targets <path/to/all_target_loci.tsv> \\
        --blocks <path/to/synteny_blocks_filtered.tsv> \\
        --output-dir <path/to/classified_targets> \\
        [--min-span-kb 0.3] \\
        [--unplaceable-evalue 1e-10] \\
        [--unplaceable-bitscore 100]
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

def classify_targets(targets_df, blocks_df, unplaceable_evalue, unplaceable_bitscore, min_span_kb):
    """Classify each target as syntenic or unplaceable based on proximity to synteny blocks.

    Quality assessment happens AFTER Exonerate extraction (Phase 6).

    Args:
        targets_df: DataFrame of target loci from Phase 4
        blocks_df: DataFrame of synteny blocks from Phase 3
        unplaceable_evalue: E-value threshold for unplaceable targets
        unplaceable_bitscore: Bitscore threshold for unplaceable targets
        min_span_kb: Minimum genomic span in kb to keep (filters tiny fragments)
    """

    # Index synteny blocks by (genome, scaffold) for fast lookup
    synteny_index = defaultdict(list)
    for _, block in blocks_df.iterrows():
        key = (block['genome'], block['scaffold'])
        synteny_index[key].append(block)

    classifications = []
    filtered_tiny = 0

    for _, target in targets_df.iterrows():
        # Filter 1: Remove tiny fragments (likely spurious hits)
        span_kb = target.get('span_kb', 0)
        if span_kb < min_span_kb:
            filtered_tiny += 1
            continue

        # Check for overlap with synteny blocks
        key = (target['genome'], target['scaffold'])
        candidate_blocks = synteny_index.get(key, [])

        assigned_to = None
        for block in candidate_blocks:
            # Check proximity - accept ANY target hit near synteny block
            if check_proximity(target, block):
                assigned_to = block['locus_id']
                break

        if assigned_to:
            # Syntenic hit - keep if passes span filter
            classification = {
                **target.to_dict(),
                'placement': 'synteny',
                'assigned_to': assigned_to,
                'parent_locus': assigned_to,  # For downstream compatibility
                'description': assigned_to
            }
            classifications.append(classification)
        else:
            # Unplaceable hit - apply strict filters
            evalue = target.get('best_evalue', 1.0)
            bitscore = target.get('bitscore', 0) if 'bitscore' in target else target.get('best_bitscore', 0)

            # Only keep strong unplaceable hits
            if evalue <= unplaceable_evalue and bitscore >= unplaceable_bitscore:
                classification = {
                    **target.to_dict(),
                    'placement': 'unplaceable',
                    'assigned_to': f"{target['gene_family']}_unplaceable",
                    'parent_locus': f"{target['gene_family']}_unplaceable",  # For downstream compatibility
                    'description': f"{target['gene_family']}_unplaceable"
                }
                classifications.append(classification)
            # Weak hits are silently filtered (not included in output)

    print(f"  Filtered {filtered_tiny} tiny fragments (< {min_span_kb} kb)")

    return pd.DataFrame(classifications)

def deduplicate_overlapping_targets(classified_df):
    """Deduplicate overlapping target loci.

    When BLAST finds multiple overlapping hits at the same genomic location,
    merge them into a single target, keeping the one with the best e-value.

    This is critical for accurate gene counts - overlapping hits often represent
    the same gene detected multiple times with slightly different boundaries.
    """

    if len(classified_df) == 0:
        return classified_df

    print("\n  Deduplicating overlapping targets...", flush=True)
    initial_count = len(classified_df)

    deduplicated = []
    duplicates_removed = 0

    # Group by (genome, parent_locus) for deduplication
    for (genome, parent_locus), group in classified_df.groupby(['genome', 'parent_locus']):
        if len(group) == 1:
            # No duplicates possible
            deduplicated.append(group.iloc[0].to_dict())
            continue

        # Sort by start position
        group = group.sort_values('start')

        # Check for overlaps
        kept_targets = []
        current_target = group.iloc[0].to_dict()

        for _, next_target in group.iloc[1:].iterrows():
            # Check if next_target overlaps with current_target
            overlap = not (next_target['end'] < current_target['start'] or
                          next_target['start'] > current_target['end'])

            if overlap:
                # Overlapping - merge by keeping the one with best e-value
                duplicates_removed += 1

                if next_target['best_evalue'] < current_target['best_evalue']:
                    # Next target is better, replace current
                    print(f"    Merged {genome}/{parent_locus}: {current_target['locus_name']} "
                          f"→ {next_target['locus_name']} (better e-value)", flush=True)
                    current_target = next_target.to_dict()
                else:
                    # Current target is better, keep it
                    print(f"    Merged {genome}/{parent_locus}: {next_target['locus_name']} "
                          f"→ {current_target['locus_name']} (worse e-value)", flush=True)
            else:
                # No overlap - save current and move to next
                kept_targets.append(current_target)
                current_target = next_target.to_dict()

        # Don't forget the last target
        kept_targets.append(current_target)

        deduplicated.extend(kept_targets)

    deduplicated_df = pd.DataFrame(deduplicated)

    print(f"  Before deduplication: {initial_count} targets", flush=True)
    print(f"  After deduplication: {len(deduplicated_df)} targets", flush=True)
    print(f"  Removed: {duplicates_removed} overlapping duplicates", flush=True)

    return deduplicated_df

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
    parser.add_argument('--min-span-kb', type=float, default=0.3,
                        help='Minimum genomic span in kb to keep (filters tiny fragments, default: 0.3)')
    parser.add_argument('--unplaceable-evalue', type=float, default=1e-10,
                        help='E-value threshold for unplaceable targets (default: 1e-10)')
    parser.add_argument('--unplaceable-bitscore', type=float, default=100,
                        help='Bitscore threshold for unplaceable targets (default: 100)')

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
    print(f"  Minimum span: {args.min_span_kb} kb")
    print(f"  Unplaceable e-value threshold: {args.unplaceable_evalue}")
    print(f"  Unplaceable bitscore threshold: {args.unplaceable_bitscore}")

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
    print(f"  Minimum span filter: {args.min_span_kb} kb")
    print(f"  Unplaceable filters: e-value <= {args.unplaceable_evalue}, bitscore >= {args.unplaceable_bitscore}")

    classified_df = classify_targets(
        targets_df,
        blocks_df,
        args.unplaceable_evalue,
        args.unplaceable_bitscore,
        args.min_span_kb
    )

    if len(classified_df) == 0:
        print("  No targets passed classification filters!", flush=True)
        return

    # Deduplicate overlapping targets
    print("\n[3b] Deduplicating overlapping targets...", flush=True)
    classified_df = deduplicate_overlapping_targets(classified_df)

    if len(classified_df) == 0:
        print("  No targets remained after deduplication!", flush=True)
        return

    syntenic = classified_df[classified_df['placement'] == 'synteny']
    unplaceable = classified_df[classified_df['placement'] == 'unplaceable']

    print(f"\n  Input targets: {len(targets_df)}")
    print(f"  After filtering: {len(classified_df)} ({len(classified_df)/len(targets_df)*100:.1f}%)")
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
    print("\n[5] Summary by assigned locus:")
    for assigned_locus in sorted(classified_df['assigned_to'].unique()):
        locus_targets = classified_df[classified_df['assigned_to'] == assigned_locus]
        locus_syntenic = locus_targets[locus_targets['placement'] == 'synteny']
        locus_unplaceable = locus_targets[locus_targets['placement'] == 'unplaceable']

        print(f"\n  {assigned_locus}:")
        print(f"    Total targets: {len(locus_targets)}")
        print(f"    Syntenic: {len(locus_syntenic)}")
        print(f"    Unplaceable: {len(locus_unplaceable)}")

    # Overall summary
    print("\n" + "=" * 80, flush=True)
    print("CLASSIFICATION COMPLETE", flush=True)
    print("=" * 80, flush=True)
    print(f"\nTotal targets classified: {len(classified_df)}", flush=True)
    print(f"Syntenic: {len(syntenic)} ({len(syntenic)/len(classified_df)*100:.1f}%)", flush=True)
    print(f"Unplaceable: {len(unplaceable)} ({len(unplaceable)/len(classified_df)*100:.1f}%)", flush=True)

    print(f"\nOutputs saved to: {args.output_dir}", flush=True)
    print("\nNext step: 06_extract_sequences.py (Exonerate will extract full gene structures)", flush=True)

if __name__ == "__main__":
    main()