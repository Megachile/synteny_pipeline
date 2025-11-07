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
import math

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
                    allow_nearby=False, proximity_kb=100,
                    pad_kb=0, dynamic_nearby=False, nearby_frac_span=0.3, nearby_upper_kb=200):
    """Classify each target as syntenic or unplaceable based on proximity to synteny blocks.

    Quality assessment happens AFTER Exonerate extraction (Phase 6).
    """

    # Optionally pad blocks for overlap/proximity checks and index by (genome, scaffold)
    synteny_index = defaultdict(list)
    pad_bp = pad_kb * 1000
    for _, block in blocks_df.iterrows():
        b = block.to_dict()
        b['start'] = int(b['start']) - pad_bp
        b['end'] = int(b['end']) + pad_bp
        # Compute span_kb if missing
        b['span_kb'] = b.get('span_kb', (int(b['end']) - int(b['start'])) / 1000)
        key = (b['genome'], b['scaffold'])
        synteny_index[key].append(b)

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
                # Dynamic nearby threshold per block (capped)
                if dynamic_nearby:
                    dyn_kb = max(30, int(block.get('span_kb', 0) * nearby_frac_span))
                    dyn_kb = min(dyn_kb, nearby_upper_kb)
                    thresh_bp = dyn_kb * 1000
                else:
                    thresh_bp = proximity_kb * 1000
                if distance <= thresh_bp:
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

            # Check if bitscore column exists (newer Phase 4 outputs include it)
            has_bitscore = 'bitscore' in target or 'best_bitscore' in target
            bitscore = target.get('best_bitscore', target.get('bitscore', 0)) if has_bitscore else None

            # Only keep strong unplaceable hits
            # If bitscore is available, require both evalue and bitscore thresholds
            # If bitscore is missing (older outputs), use evalue only
            passes_filters = evalue <= unplaceable_evalue
            if has_bitscore:
                passes_filters = passes_filters and bitscore >= unplaceable_bitscore

            if passes_filters:
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
                    'bitscore': bitscore if has_bitscore else 'N/A',
                    'drop_reason': 'outside_block_and_below_thresholds'
                })
    return pd.DataFrame(classifications), pd.DataFrame(dropped)

def compute_distance_debug(targets_df, blocks_df, bins_kb, out_path: Path = None):
    """Compute min distance of each target to any block on same genome/scaffold (0 if overlap).

    Returns a dict with histogram counts and optionally writes a TSV of distances.
    """
    # Build index
    index = defaultdict(list)
    for _, b in blocks_df.iterrows():
        key = (b['genome'], b['scaffold'])
        index[key].append((int(b['start']), int(b['end'])))

    distances = []
    for _, t in targets_df.iterrows():
        key = (t['genome'], t['scaffold'])
        bks = index.get(key, [])
        if not bks:
            continue
        s, e = int(t['start']), int(t['end'])
        md = 10**12
        for bs, be in bks:
            if not (e < bs or s > be):
                md = 0
                break
            d = bs - e if e < bs else s - be
            if d < md:
                md = d
        distances.append({
            'genome': t['genome'],
            'scaffold': t['scaffold'],
            'locus_name': t.get('locus_name', ''),
            'query_id': t.get('query_id', ''),
            'start': s,
            'end': e,
            'min_distance_bp': md
        })

    if not distances:
        return {'total': 0, 'bins': []}

    # Histogram
    bins_bp = [int(k*1000) for k in bins_kb]
    counts = [0]*len(bins_bp)
    total = len(distances)
    for rec in distances:
        d = rec['min_distance_bp']
        for i, b in enumerate(bins_bp):
            if d <= b:
                counts[i] += 1
                break

    # Write TSV if requested
    if out_path:
        df = pd.DataFrame(distances)
        df.to_csv(out_path, sep='\t', index=False)

    return {'total': total, 'bins': list(zip(bins_kb, counts))}

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
                        help='Max distance in kb for nearby classification when not using dynamic (default: 100)')
    parser.add_argument('--pad-kb', type=int, default=0,
                        help='Pad block start/end by this many kb for overlap/proximity checks (default: 0)')
    parser.add_argument('--dynamic-nearby', action='store_true',
                        help='Use per-block dynamic nearby threshold = min(max(frac*span,30kb), nearby-upper-kb)')
    parser.add_argument('--nearby-frac-span', type=float, default=0.3,
                        help='Fraction of block span used for dynamic nearby threshold (default: 0.3)')
    parser.add_argument('--nearby-upper-kb', type=int, default=200,
                        help='Upper cap for dynamic nearby threshold in kb (default: 200)')
    parser.add_argument('--debug-distances', action='store_true',
                        help='Print a histogram of min distances from targets to nearest block and write a TSV file')
    parser.add_argument('--distance-bins-kb', type=str, default='30,50,75,100,150,200,300,500,1000',
                        help='Comma-separated distance bins in kb for debug histogram')

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

    # Optional debug: distances from targets to nearest block (raw, without padding)
    if args.debug_distances:
        try:
            bins_kb = [int(x) for x in args.distance_bins_kb.split(',') if x.strip()]
        except Exception:
            bins_kb = [30,50,75,100,150,200,300,500,1000]
        debug_path = args.output_dir / 'debug_target_block_distances.tsv'
        stats = compute_distance_debug(targets_df, blocks_df, bins_kb, out_path=debug_path)
        print(f"\n[2b] Distance debug (raw, no padding): total targets with candidate blocks: {stats['total']}")
        if stats['total'] > 0:
            for kb, c in stats['bins']:
                pct = (c / stats['total'] * 100) if stats['total'] else 0.0
                print(f"  <= {kb}kb: {c} ({pct:.1f}%)")
            print(f"  Wrote per-target distances to: {debug_path}")

    # Infer paths for reference proteins (in same directory as targets file)
    targets_dir = args.targets.parent
    bk_proteins_file = targets_dir.parent / "01_extracted_proteins" / "BK_proteins.faa"

    # Classify targets
    print("\n[3] Classifying targets...", flush=True)
    print(f"  Unplaceable filters: e-value <= {args.unplaceable_evalue}, bitscore >= {args.unplaceable_bitscore}")
    if args.allow_nearby:
        if args.dynamic_nearby:
            print(f"  Synteny placement: overlap-or-dynamic-nearby (frac={args.nearby_frac_span}, cap={args.nearby_upper_kb}kb), pad={args.pad_kb}kb")
        else:
            print(f"  Synteny placement: overlap-or-within {args.proximity_kb}kb, pad={args.pad_kb}kb")
    else:
        print(f"  Synteny placement: overlap-only, pad={args.pad_kb}kb")

    classified_df, dropped_df = classify_targets(
        targets_df,
        blocks_df,
        args.unplaceable_evalue,
        args.unplaceable_bitscore,
        allow_nearby=args.allow_nearby,
        proximity_kb=args.proximity_kb,
        pad_kb=args.pad_kb,
        dynamic_nearby=args.dynamic_nearby,
        nearby_frac_span=args.nearby_frac_span,
        nearby_upper_kb=args.nearby_upper_kb
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
