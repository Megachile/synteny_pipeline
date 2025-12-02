#!/usr/bin/env python3
"""
Phase 5 (Helixer v3): Classify targets by flanking gene synteny.

Flow:
1. Take Helixer flanking proteins from Phase 4b
2. DIAMOND blastp against full BK proteome
3. Check if BK hits match Phase 1 flanking genes
4. Calculate synteny score

Inputs:
- Phase 4 Helixer targets (all_target_loci.tsv)
- Phase 4b flanking proteins (flanking_genes.tsv, flanking_proteins.faa)
- Phase 1 locus definitions (which flanking genes belong to which locus)
- BK full proteome for DIAMOND search

Outputs:
- syntenic_targets.tsv: Targets with flanking gene matches
- unplaceable_targets.tsv: Targets without synteny support (ENHANCED with detailed info)
- flanking_matches.tsv: Detailed flanking gene comparison
- locus_conflicts.tsv: Targets matching multiple loci

Enhanced unplaceable output includes:
- best_locus_match: Best matching locus even if below threshold
- best_locus_score: Score for best locus
- n_flanking_with_bk_hits: How many flanking genes hit ANY BK protein
- unplaceable_reason: Why classification failed
- all_locus_scores: All locus scores for debugging
"""

from __future__ import annotations

import argparse
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Dict, Set

import pandas as pd


# Default BK proteome location
DEFAULT_BK_PROTEOME = Path("/carc/scratch/projects/emartins/2016456/adam/databases/bkinseyi_proteome.dmnd")


def clean_fasta_for_diamond(input_faa: Path, output_faa: Path) -> int:
    """Clean FASTA file for DIAMOND by replacing invalid characters."""
    valid_aa = set('ACDEFGHIKLMNPQRSTVWYX*')
    n_written = 0

    with open(input_faa) as f_in, open(output_faa, 'w') as f_out:
        current_header = None
        current_seq = []

        for line in f_in:
            line = line.strip()
            if line.startswith('>'):
                if current_header and current_seq:
                    seq = ''.join(current_seq).replace('.', 'X')
                    if len(seq) > 0:
                        f_out.write(f"{current_header}\n{seq}\n")
                        n_written += 1
                current_header = line
                current_seq = []
            else:
                current_seq.append(line)

        if current_header and current_seq:
            seq = ''.join(current_seq).replace('.', 'X')
            if len(seq) > 0:
                f_out.write(f"{current_header}\n{seq}\n")
                n_written += 1

    return n_written


def run_diamond_blastp(
    query_faa: Path,
    db_dmnd: Path,
    output_tsv: Path,
    evalue: float = 1e-3,
    threads: int = 4,
) -> bool:
    """Run DIAMOND blastp against pre-built database."""
    cmd = [
        'diamond', 'blastp',
        '--query', str(query_faa),
        '--db', str(db_dmnd).replace('.dmnd', ''),
        '--out', str(output_tsv),
        '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length', 'evalue', 'bitscore',
        '--evalue', str(evalue),
        '--max-target-seqs', '1',
        '--threads', str(threads),
        '--quiet'
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.returncode == 0


def load_phase1_flanking_genes(phase1_dir: Path) -> Dict[str, Set[str]]:
    """
    Load flanking gene IDs from Phase 1 locus files.

    Returns: {locus_id: set of flanking gene LOC IDs}
    """
    locus_flanking = {}

    locus_defs = phase1_dir / 'locus_definitions.tsv'
    if not locus_defs.exists():
        return {}

    df = pd.read_csv(locus_defs, sep='\t')

    for _, row in df.iterrows():
        locus_id = row['locus_id']
        flanking_file = phase1_dir / f"{locus_id}_flanking.faa"

        if not flanking_file.exists():
            continue

        # Extract XP IDs from flanking protein headers
        flanking_xp_ids = set()
        with open(flanking_file) as f:
            for line in f:
                if line.startswith('>'):
                    # Format: >XP_033208252.1|LOC117167433 U1 ...
                    header = line[1:].split()[0]
                    # Extract XP ID (part before |)
                    xp_id = header.split('|')[0] if '|' in header else header
                    flanking_xp_ids.add(xp_id)

        locus_flanking[locus_id] = flanking_xp_ids

    return locus_flanking


def main():
    parser = argparse.ArgumentParser(
        description="Phase 5 (Helixer v3): Classify targets by flanking gene synteny"
    )
    parser.add_argument("--family", required=True)
    parser.add_argument("--phase4-dir", type=Path, required=True)
    parser.add_argument("--phase4b-dir", type=Path, required=True)
    parser.add_argument("--phase1-dir", type=Path, required=True)
    parser.add_argument("--bk-proteome", type=Path, default=DEFAULT_BK_PROTEOME,
                        help="BK proteome DIAMOND database")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--min-synteny", type=float, default=0.15,
                        help="Minimum synteny score (default: 0.15 = 3/20 flanking matches)")
    parser.add_argument("--evalue", type=float, default=1e-3)
    parser.add_argument("--threads", type=int, default=4)

    args = parser.parse_args()

    print("=" * 80)
    print("PHASE 5 (HELIXER v3): FLANKING GENE SYNTENY")
    print("=" * 80)
    print()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Load targets
    targets_file = args.phase4_dir / "all_target_loci.tsv"
    targets_df = pd.read_csv(targets_file, sep='\t')
    print(f"Loaded {len(targets_df)} targets")

    # Load flanking genes from Phase 4b
    flanking_file = args.phase4b_dir / "flanking_genes.tsv"
    flanking_faa = args.phase4b_dir / "flanking_proteins.faa"

    if not flanking_file.exists():
        print("[ERROR] Phase 4b flanking_genes.tsv not found")
        return 1

    flanking_df = pd.read_csv(flanking_file, sep='\t')
    print(f"Loaded {len(flanking_df)} flanking gene entries")

    # Load Phase 1 flanking gene sets
    locus_flanking = load_phase1_flanking_genes(args.phase1_dir)
    print(f"Loaded {len(locus_flanking)} reference loci:")
    for locus_id, locs in locus_flanking.items():
        print(f"  {locus_id}: {len(locs)} flanking genes")

    # Clean flanking proteins for DIAMOND (remove invalid characters)
    print(f"\nCleaning flanking proteins for DIAMOND...")
    cleaned_faa = args.output_dir / "flanking_proteins_cleaned.faa"
    n_cleaned = clean_fasta_for_diamond(flanking_faa, cleaned_faa)
    print(f"  Cleaned {n_cleaned} sequences")

    # Run DIAMOND against BK proteome
    print(f"Running DIAMOND against BK proteome...")
    diamond_out = args.output_dir / "diamond_vs_bk.tsv"

    success = run_diamond_blastp(
        cleaned_faa, args.bk_proteome, diamond_out,
        evalue=args.evalue, threads=args.threads
    )

    if not success:
        print("[ERROR] DIAMOND failed")
        return 1

    # Parse DIAMOND results
    try:
        hits_df = pd.read_csv(diamond_out, sep='\t', header=None,
                              names=['qseqid', 'sseqid', 'pident', 'length', 'evalue', 'bitscore'])
        print(f"DIAMOND hits: {len(hits_df)}")
    except pd.errors.EmptyDataError:
        hits_df = pd.DataFrame()
        print("DIAMOND hits: 0")

    # Map Helixer flanking → BK XP ID
    # BK proteome headers: XP_033208252.1
    helixer_to_bk = {}
    if not hits_df.empty:
        for _, hit in hits_df.iterrows():
            helixer_id = hit['qseqid']
            bk_xp = hit['sseqid']  # Just XP ID
            helixer_to_bk[helixer_id] = {
                'bk_xp': bk_xp,
                'pident': hit['pident'],
                'evalue': hit['evalue'],
            }

    print(f"Helixer flanking with BK hits: {len(helixer_to_bk)}")

    # Map flanking genes to targets
    flanking_to_target = flanking_df.set_index('flanking_protein_id')['target_gene_id'].to_dict()

    # Count flanking genes per target
    target_flanking_count = flanking_df.groupby('target_gene_id').size().to_dict()

    # Track which flanking genes (per target) have ANY BK hit (not just locus-specific)
    target_flanking_with_bk_hits = defaultdict(set)
    for helixer_id in helixer_to_bk.keys():
        target_id = flanking_to_target.get(helixer_id)
        if target_id:
            target_flanking_with_bk_hits[target_id].add(helixer_id)

    # Calculate synteny for each target
    target_synteny = defaultdict(lambda: defaultdict(list))

    for helixer_id, bk_info in helixer_to_bk.items():
        target_id = flanking_to_target.get(helixer_id)
        if not target_id:
            continue

        bk_xp = bk_info['bk_xp']

        # Check which locus this BK gene belongs to
        for locus_id, locus_xp_ids in locus_flanking.items():
            if bk_xp in locus_xp_ids:
                target_synteny[target_id][locus_id].append({
                    'helixer_flanking': helixer_id,
                    'bk_xp': bk_xp,
                    'pident': bk_info['pident'],
                    'evalue': bk_info['evalue'],
                })

    # Classify targets
    syntenic_rows = []
    unplaceable_rows = []
    match_details = []
    conflict_details = []

    for _, row in targets_df.iterrows():
        target_id = row['target_gene_id']
        n_flanking = target_flanking_count.get(target_id, 0)
        n_flanking_with_bk = len(target_flanking_with_bk_hits.get(target_id, set()))

        # Calculate scores for ALL loci
        locus_scores = {}
        for locus_id, matches in target_synteny[target_id].items():
            n_matches = len(matches)
            score = n_matches / n_flanking if n_flanking > 0 else 0
            locus_scores[locus_id] = {
                'score': score,
                'n_matches': n_matches,
                'matches': matches,
            }

            # Record match details
            for m in matches:
                match_details.append({
                    'target_gene_id': target_id,
                    'locus_id': locus_id,
                    'helixer_flanking': m['helixer_flanking'],
                    'bk_xp': m['bk_xp'],
                    'pident': m['pident'],
                    'evalue': m['evalue'],
                })

        # Find loci meeting threshold
        qualifying_loci = {
            locus_id: info for locus_id, info in locus_scores.items()
            if info['score'] >= args.min_synteny
        }

        # Check for conflicts (multiple loci meet threshold)
        has_conflict = len(qualifying_loci) > 1

        # Find best locus (even if below threshold)
        best_locus = None
        best_score = 0
        best_matches = 0

        for locus_id, info in locus_scores.items():
            if info['score'] > best_score:
                best_score = info['score']
                best_locus = locus_id
                best_matches = info['n_matches']

        row_dict = row.to_dict()
        row_dict['n_flanking_genes'] = n_flanking
        row_dict['n_flanking_with_bk_hits'] = n_flanking_with_bk
        row_dict['n_flanking_matches'] = best_matches
        row_dict['synteny_score'] = round(best_score, 3)
        row_dict['matched_locus'] = best_locus or ''

        # Add conflict info
        row_dict['has_conflict'] = has_conflict
        if has_conflict:
            conflict_loci = sorted(qualifying_loci.keys())
            conflict_scores = [f"{l}:{qualifying_loci[l]['score']:.2f}" for l in conflict_loci]
            row_dict['conflict_loci'] = ';'.join(conflict_loci)
            row_dict['conflict_scores'] = ';'.join(conflict_scores)

            # Record conflict details
            conflict_details.append({
                'target_gene_id': target_id,
                'genome': row.get('genome', ''),
                'n_qualifying_loci': len(qualifying_loci),
                'qualifying_loci': ';'.join(conflict_loci),
                'scores': ';'.join(conflict_scores),
                'assigned_to': best_locus,
            })
        else:
            row_dict['conflict_loci'] = ''
            row_dict['conflict_scores'] = ''

        # Build all_locus_scores string for debugging
        if locus_scores:
            all_scores_str = ';'.join([f"{l}:{info['score']:.3f}" for l, info in sorted(locus_scores.items())])
        else:
            all_scores_str = ''
        row_dict['all_locus_scores'] = all_scores_str

        if best_score >= args.min_synteny:
            row_dict['classification'] = 'syntenic'
            row_dict['locus_id'] = best_locus
            row_dict['unplaceable_reason'] = ''
            syntenic_rows.append(row_dict)
        else:
            row_dict['classification'] = 'unplaceable'
            # Determine reason for being unplaceable
            if n_flanking == 0:
                reason = 'no_flanking_genes'
            elif n_flanking_with_bk == 0:
                reason = 'no_bk_hits'
            elif best_score == 0:
                reason = 'no_locus_matches'
            else:
                reason = f'below_threshold_{best_score:.3f}'
            row_dict['unplaceable_reason'] = reason
            row_dict['best_locus_match'] = best_locus or ''
            row_dict['best_locus_score'] = round(best_score, 3)
            unplaceable_rows.append(row_dict)

    # Write outputs
    syntenic_df = pd.DataFrame(syntenic_rows)
    unplaceable_df = pd.DataFrame(unplaceable_rows)

    # Add columns needed by Phase 8 scripts
    # For syntenic: placement='synteny', assigned_to=locus_id, locus_name=locus_id
    if not syntenic_df.empty:
        syntenic_df['placement'] = 'synteny'
        syntenic_df['assigned_to'] = syntenic_df['locus_id']
        syntenic_df['locus_name'] = syntenic_df['locus_id']

    # For unplaceable: placement='unplaceable', assigned_to='', locus_name=parent_locus
    if not unplaceable_df.empty:
        unplaceable_df['placement'] = 'unplaceable'
        unplaceable_df['assigned_to'] = ''
        unplaceable_df['locus_name'] = unplaceable_df['parent_locus']

    syntenic_df.to_csv(args.output_dir / "syntenic_targets.tsv", sep='\t', index=False)
    unplaceable_df.to_csv(args.output_dir / "unplaceable_targets.tsv", sep='\t', index=False)

    # Create all_targets_classified.tsv (combined) for Phase 8
    all_targets = pd.concat([syntenic_df, unplaceable_df], ignore_index=True)
    all_targets.to_csv(args.output_dir / "all_targets_classified.tsv", sep='\t', index=False)

    if match_details:
        pd.DataFrame(match_details).to_csv(
            args.output_dir / "flanking_matches.tsv", sep='\t', index=False)

    # Write conflict report
    if conflict_details:
        pd.DataFrame(conflict_details).to_csv(
            args.output_dir / "locus_conflicts.tsv", sep='\t', index=False)

    print(f"\n[OUTPUT] Syntenic: {len(syntenic_df)}")
    print(f"[OUTPUT] Unplaceable: {len(unplaceable_df)}")
    print(f"[OUTPUT] Total classified: {len(all_targets)}")

    # Conflict summary
    n_conflicts = len(conflict_details)
    if n_conflicts > 0:
        print(f"\n[WARNING] {n_conflicts} targets have LOCUS CONFLICTS (match multiple loci):")
        for c in conflict_details[:10]:  # Show first 10
            print(f"  {c['target_gene_id']}: {c['qualifying_loci']} -> assigned to {c['assigned_to']}")
        if n_conflicts > 10:
            print(f"  ... and {n_conflicts - 10} more (see locus_conflicts.tsv)")

    if not syntenic_df.empty:
        print(f"\nSyntenic by locus:")
        for locus in syntenic_df['locus_id'].unique():
            n = len(syntenic_df[syntenic_df['locus_id'] == locus])
            n_conflict = len(syntenic_df[(syntenic_df['locus_id'] == locus) & (syntenic_df['has_conflict'] == True)])
            if n_conflict > 0:
                print(f"  {locus}: {n} ({n_conflict} with conflicts)")
            else:
                print(f"  {locus}: {n}")

    # Unplaceable reason summary
    if not unplaceable_df.empty:
        print(f"\nUnplaceable by reason:")
        reason_counts = unplaceable_df['unplaceable_reason'].value_counts()
        for reason, count in reason_counts.items():
            print(f"  {reason}: {count}")

        # Show targets with partial matches (below threshold but has some signal)
        partial_matches = unplaceable_df[unplaceable_df['unplaceable_reason'].str.startswith('below_threshold')]
        if not partial_matches.empty:
            print(f"\nPartial matches (below threshold but has synteny signal):")
            for _, row in partial_matches.head(10).iterrows():
                print(f"  {row['target_gene_id']}: best={row.get('best_locus_match', 'N/A')} "
                      f"score={row.get('best_locus_score', 0):.3f} "
                      f"({row['n_flanking_matches']}/{row['n_flanking_genes']} flanking)")
            if len(partial_matches) > 10:
                print(f"  ... and {len(partial_matches) - 10} more")

    return 0


if __name__ == "__main__":
    exit(main())
