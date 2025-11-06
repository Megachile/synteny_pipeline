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
import os

def check_proximity(target, block, max_distance_kb=100):
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

def check_overlap(target, block):
    """Require overlap between target and block on same genome/scaffold."""
    if target['genome'] != block['genome'] or target['scaffold'] != block['scaffold']:
        return False
    return not (target['end'] < block['start'] or target['start'] > block['end'])

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

def classify_targets(targets_df, blocks_df, unplaceable_evalue, unplaceable_bitscore,
                    allow_nearby=False, proximity_kb=100):
    """Classify each target as syntenic or unplaceable based on proximity to synteny blocks.

    Quality assessment happens AFTER Exonerate extraction (Phase 6).
    """

    # Index synteny blocks by (genome, scaffold) for fast lookup
    synteny_index = defaultdict(list)
    for _, block in blocks_df.iterrows():
        key = (block['genome'], block['scaffold'])
        synteny_index[key].append(block)

    classifications = []
    dropped = []

    for _, target in targets_df.iterrows():
        # Check for overlap with synteny blocks
        key = (target['genome'], target['scaffold'])
        candidate_blocks = synteny_index.get(key, [])

        assigned_to = None
        best_overlap = -1
        best_distance = None
        best_block = None

        # Choose the best block by overlap first, then by nearest distance
        for block in candidate_blocks:
            # Compute overlap in bp
            start, end = int(target['start']), int(target['end'])
            bstart, bend = int(block['start']), int(block['end'])
            overlap = max(0, min(end, bend) - max(start, bstart) + 1)

            if overlap > 0:
                if overlap > best_overlap:
                    best_overlap = overlap
                    best_block = block
                    best_distance = 0
                elif overlap == best_overlap and best_overlap > 0 and best_block is not None:
                    # Tie-breaker: choose block whose center is closest to target center
                    t_center = (start + end) // 2
                    b_center_curr = (bstart + bend) // 2
                    b_center_best = (int(best_block['start']) + int(best_block['end'])) // 2
                    if abs(t_center - b_center_curr) < abs(t_center - b_center_best):
                        best_block = block
                        best_distance = 0
                continue

            # If no overlap, consider proximity when enabled
            if allow_nearby:
                # Distance to nearest edge
                if end < bstart:
                    distance = bstart - end
                else:
                    distance = start - bend
                if distance <= proximity_kb * 1000:
                    if best_distance is None or distance < best_distance:
                        best_distance = distance
                        best_block = block

        if best_block is not None and (best_overlap > 0 or (allow_nearby and best_distance is not None)):
            assigned_to = best_block['locus_id']

        if assigned_to:
            # Syntenic hit - keep regardless of strength
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
            else:
                dropped.append({
                    'locus_name': target.get('locus_name', ''),
                    'scaffold': target.get('scaffold', ''),
                    'strand': target.get('strand', ''),
                    'frame': target.get('frame', 0),
                    'start': target.get('start', 0),
                    'end': target.get('end', 0),
                    'genome': target.get('genome', ''),
                    'query_id': target.get('query_id', ''),
                    'best_evalue': evalue,
                    'bitscore': bitscore,
                    'drop_reason': 'outside_block_and_below_thresholds'
                })
    return pd.DataFrame(classifications), pd.DataFrame(dropped)

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

def cap_targets_by_expected_count(classified_df, locus_defs_path: Path):
    """No-op capping: return input unchanged to remain agnostic on tandem counts."""
    return classified_df

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
    parser.add_argument('--unplaceable-evalue', type=float, default=1e-5,
                        help='E-value threshold for unplaceable targets (default: 1e-5)')
    parser.add_argument('--unplaceable-bitscore', type=float, default=50,
                        help='Bitscore threshold for unplaceable targets (default: 50)')
    parser.add_argument('--allow-nearby', action='store_true',
                        help='Allow targets near (not overlapping) a block to be syntenic')
    parser.add_argument('--proximity-kb', type=int, default=100,
                        help='Max distance in kb for nearby classification (default: 100)')

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
    print(f"  Unplaceable filters: e-value <= {args.unplaceable_evalue}, bitscore >= {args.unplaceable_bitscore}")
    print(f"  Synteny placement: {'overlap-only' if not args.allow_nearby else f'overlap-or-within {args.proximity_kb}kb'}")

    classified_df, dropped_df = classify_targets(
        targets_df,
        blocks_df,
        args.unplaceable_evalue,
        args.unplaceable_bitscore,
        allow_nearby=args.allow_nearby,
        proximity_kb=args.proximity_kb
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

    # Cap per-locus targets using expected counts from locus_definitions.tsv
    locus_defs_candidate = None
    # Prefer sibling of targets file (family root)/locus_definitions.tsv
    try:
        locus_defs_candidate = (args.targets.parent.parent / 'locus_definitions.tsv')
    except Exception:
        locus_defs_candidate = None

    if locus_defs_candidate and locus_defs_candidate.exists():
        print(f"\n[3c] Capping targets by expected locus counts using: {locus_defs_candidate}", flush=True)
        before = len(classified_df)
        classified_df = cap_targets_by_expected_count(classified_df, locus_defs_candidate)
        print(f"  After capping: {len(classified_df)} targets (was {before})", flush=True)

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
    # Save dropped targets, if any
    if 'dropped_df' in locals() and len(dropped_df) > 0:
        dropped_file = args.output_dir / "dropped_targets.tsv"
        dropped_df.to_csv(dropped_file, sep='\t', index=False)
        print(f"  Saved dropped: {dropped_file.name}", flush=True)

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
    if 'dropped_df' in locals():
        print(f"Dropped: {len(dropped_df)}", flush=True)

    print(f"\nOutputs saved to: {args.output_dir}", flush=True)
    print("\nNext step: 06_extract_sequences.py (Exonerate will extract full gene structures)", flush=True)

if __name__ == "__main__":
    main()
