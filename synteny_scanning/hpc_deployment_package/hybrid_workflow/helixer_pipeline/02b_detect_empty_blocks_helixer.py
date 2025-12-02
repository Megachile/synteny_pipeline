#!/usr/bin/env python3
"""
Phase 2b (Helixer): Detect synteny blocks via flanking-first approach.

Searches for BK flanking genes across all Helixer proteomes, then clusters
hits to identify conserved synteny blocks. This is complementary to the
target-first approach in Phase 4/5.

Flow:
1. DIAMOND blastp: BK flanking proteins → all Helixer proteomes
2. Cluster flanking gene hits by scaffold/position
3. Score clusters by flanking gene coverage (vs BK reference)
4. Check if cluster overlaps existing Phase 4 target hits
5. Label blocks as "concordant" (has target) or "empty" (no target)

Inputs:
- Phase 1 locus flanking FASTAs (BK flanking proteins)
- Helixer proteomes for all genomes
- Phase 4 target hits (to identify concordant blocks)

Outputs:
- phase2b_synteny_blocks.tsv: All synteny blocks with concordant/empty status
- phase2b_flanking_details.tsv: Individual flanking gene hits (Helixer IDs)
- locus_coverage_summary.tsv: Per-locus coverage across genomes
"""

from __future__ import annotations

import argparse
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Set, Tuple

import pandas as pd
from Bio import SeqIO


def run_diamond_blastp(
    query_faa: Path,
    db_dmnd: Path,
    output_tsv: Path,
    evalue: float = 1e-5,
    max_targets: int = 5,
    threads: int = 4,
) -> bool:
    """Run DIAMOND blastp."""
    cmd = [
        'diamond', 'blastp',
        '--query', str(query_faa),
        '--db', str(db_dmnd).replace('.dmnd', ''),
        '--out', str(output_tsv),
        '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length', 'evalue', 'bitscore',
        '--evalue', str(evalue),
        '--max-target-seqs', str(max_targets),
        '--threads', str(threads),
        '--quiet'
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.returncode == 0


def build_diamond_db(proteome_faa: Path, db_path: Path) -> bool:
    """Build DIAMOND database from proteome FASTA."""
    db_path.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        'diamond', 'makedb',
        '--in', str(proteome_faa),
        '--db', str(db_path),
        '--quiet'
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.returncode == 0


def parse_helixer_gff3_genes(gff3_path: Path) -> Dict[str, Dict]:
    """Parse Helixer GFF3 to get gene coordinates."""
    genes = {}

    with open(gff3_path) as f:
        for line in f:
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            scaffold, source, feature, start, end, score, strand, phase, attributes = parts

            if feature != 'gene':
                continue

            # Parse attributes
            attr_dict = {}
            for attr in attributes.split(';'):
                if '=' in attr:
                    key, value = attr.split('=', 1)
                    attr_dict[key] = value

            gene_id = attr_dict.get('ID', '')
            if not gene_id:
                continue

            genes[gene_id] = {
                'scaffold': scaffold,
                'strand': strand,
                'start': int(start),
                'end': int(end),
            }

    return genes


def protein_to_gene_id(protein_id: str) -> str:
    """Convert protein ID to gene ID (remove .1 suffix)."""
    if '.' in protein_id:
        parts = protein_id.rsplit('.', 1)
        if parts[1].isdigit():
            return parts[0]
    return protein_id


def load_phase1_flanking_info(phase1_dir: Path) -> Dict[str, Dict]:
    """
    Load Phase 1 locus information.

    Returns: {locus_id: {flanking_file, n_flanking, xp_ids}}
    """
    locus_info = {}

    locus_defs = phase1_dir / 'locus_definitions.tsv'
    if not locus_defs.exists():
        return {}

    df = pd.read_csv(locus_defs, sep='\t')

    for _, row in df.iterrows():
        locus_id = row['locus_id']
        flanking_file = phase1_dir / f"{locus_id}_flanking.faa"

        if not flanking_file.exists():
            continue

        # Count and extract XP IDs from flanking file
        xp_ids = set()
        n_flanking = 0
        for record in SeqIO.parse(flanking_file, 'fasta'):
            n_flanking += 1
            # Header: >XP_033208252.1|LOC117167433 U1 ...
            xp_id = record.id.split('|')[0] if '|' in record.id else record.id
            xp_ids.add(xp_id)

        locus_info[locus_id] = {
            'flanking_file': flanking_file,
            'n_flanking': n_flanking,
            'xp_ids': xp_ids,
            'chromosome': row.get('expected_chromosome', ''),
            'target_gene': row.get('target_gene', ''),
        }

    return locus_info


def cluster_flanking_hits(
    hits: List[Dict],
    gff3_data: Dict[str, Dict],
    max_span_kb: float = 500.0,
) -> List[Dict]:
    """
    Cluster flanking gene hits by scaffold/position.

    Returns list of clusters with score and position info.
    """
    if not hits:
        return []

    # Group by scaffold
    by_scaffold = defaultdict(list)
    for hit in hits:
        gene_id = protein_to_gene_id(hit['sseqid'])
        if gene_id not in gff3_data:
            continue

        gene_info = gff3_data[gene_id]
        by_scaffold[gene_info['scaffold']].append({
            **hit,
            'gene_id': gene_id,
            'start': gene_info['start'],
            'end': gene_info['end'],
        })

    max_span_bp = max_span_kb * 1000
    clusters = []

    for scaffold, scaffold_hits in by_scaffold.items():
        if not scaffold_hits:
            continue

        # Sort by position
        scaffold_hits.sort(key=lambda x: x['start'])

        # Greedy clustering
        current_cluster = [scaffold_hits[0]]
        cluster_start = scaffold_hits[0]['start']
        cluster_end = scaffold_hits[0]['end']

        for hit in scaffold_hits[1:]:
            # Check if this hit extends the cluster too much
            new_end = max(cluster_end, hit['end'])
            new_start = min(cluster_start, hit['start'])

            if (new_end - new_start) <= max_span_bp:
                current_cluster.append(hit)
                cluster_start = new_start
                cluster_end = new_end
            else:
                # Save current cluster and start new one
                if len(current_cluster) >= 2:  # Minimum 2 flanking genes
                    clusters.append({
                        'scaffold': scaffold,
                        'start': cluster_start,
                        'end': cluster_end,
                        'hits': current_cluster,
                    })
                current_cluster = [hit]
                cluster_start = hit['start']
                cluster_end = hit['end']

        # Don't forget the last cluster
        if len(current_cluster) >= 2:
            clusters.append({
                'scaffold': scaffold,
                'start': cluster_start,
                'end': cluster_end,
                'hits': current_cluster,
            })

    return clusters


def check_target_overlap(
    cluster: Dict,
    target_df: pd.DataFrame,
    genome: str,
) -> Tuple[bool, str]:
    """
    Check if a cluster contains any Phase 4 target hit.

    For concordance, we require the target to be INSIDE the cluster
    (no margin) - this ensures both methods found the same region.

    Returns: (has_overlap, overlapping_target_id)
    """
    if target_df.empty:
        return False, ''

    cluster_start = cluster['start']
    cluster_end = cluster['end']

    genome_targets = target_df[
        (target_df['genome'] == genome) &
        (target_df['scaffold'] == cluster['scaffold'])
    ]

    for _, target in genome_targets.iterrows():
        target_start = target['start']
        target_end = target['end']

        # Check if target is INSIDE cluster (strict containment)
        if target_start >= cluster_start and target_end <= cluster_end:
            return True, target['target_gene_id']

    return False, ''


def main():
    parser = argparse.ArgumentParser(
        description="Phase 2b (Helixer): Detect empty synteny blocks"
    )
    parser.add_argument("--family", required=True, help="Gene family name")
    parser.add_argument(
        "--phase1-dir", type=Path, required=True,
        help="Phase 1 directory with locus definitions and flanking FASTAs"
    )
    parser.add_argument(
        "--helixer-dir", type=Path, required=True,
        help="Directory containing Helixer outputs"
    )
    parser.add_argument(
        "--phase4-targets", type=Path, default=None,
        help="Phase 4 targets file (to exclude regions with hits)"
    )
    parser.add_argument(
        "--genome-list", type=Path, default=None,
        help="TSV with genome_id column (optional, auto-detects if not provided)"
    )
    parser.add_argument(
        "--output-dir", type=Path, required=True,
        help="Output directory"
    )
    parser.add_argument(
        "--diamond-db-dir", type=Path, default=None,
        help="Directory for DIAMOND databases"
    )
    parser.add_argument(
        "--evalue", type=float, default=1e-5,
        help="E-value threshold for DIAMOND"
    )
    parser.add_argument(
        "--min-flanking", type=int, default=3,
        help="Minimum flanking genes for a valid block (default: 3)"
    )
    parser.add_argument(
        "--min-score", type=float, default=0.15,
        help="Minimum synteny score (default: 0.15)"
    )
    parser.add_argument(
        "--max-span-kb", type=float, default=500.0,
        help="Maximum span for clustering flanking hits (default: 500kb)"
    )
    parser.add_argument(
        "--threads", type=int, default=4,
        help="DIAMOND threads"
    )

    args = parser.parse_args()

    print("=" * 80)
    print("PHASE 2B (HELIXER): DETECT EMPTY SYNTENY BLOCKS")
    print("=" * 80)
    print()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    diamond_db_dir = args.diamond_db_dir or (args.output_dir / "diamond_dbs")
    diamond_db_dir.mkdir(parents=True, exist_ok=True)
    diamond_results_dir = args.output_dir / "diamond_results"
    diamond_results_dir.mkdir(parents=True, exist_ok=True)

    # Load Phase 1 locus info
    locus_info = load_phase1_flanking_info(args.phase1_dir)
    print(f"Loaded {len(locus_info)} reference loci:")
    for locus_id, info in locus_info.items():
        print(f"  {locus_id}: {info['n_flanking']} flanking genes")

    if not locus_info:
        print("[ERROR] No locus information found")
        return 1

    # Load Phase 4 targets (to exclude)
    target_df = pd.DataFrame()
    if args.phase4_targets and args.phase4_targets.exists():
        target_df = pd.read_csv(args.phase4_targets, sep='\t')
        print(f"Loaded {len(target_df)} Phase 4 targets (will exclude overlaps)")

    # Get genome list
    if args.genome_list and args.genome_list.exists():
        genome_df = pd.read_csv(args.genome_list, sep='\t')
        genomes = genome_df['genome_id'].tolist() if 'genome_id' in genome_df.columns else genome_df.iloc[:, 0].tolist()
    else:
        # Auto-detect from Helixer directory
        genomes = [d.name for d in args.helixer_dir.iterdir()
                   if d.is_dir() and (d / f"{d.name}_helixer_proteins.faa").exists()]

    print(f"Processing {len(genomes)} genomes")

    # Process each locus
    all_empty_blocks = []
    all_hit_details = []
    locus_coverage = []

    for locus_id, info in locus_info.items():
        print(f"\n[{locus_id}] Processing...")

        flanking_faa = info['flanking_file']
        n_reference_flanking = info['n_flanking']

        # Process each genome
        for genome_id in genomes:
            genome_dir = args.helixer_dir / genome_id
            proteome_faa = genome_dir / f"{genome_id}_helixer_proteins.faa"
            gff3_file = genome_dir / f"{genome_id}_helixer.gff3"

            if not proteome_faa.exists() or not gff3_file.exists():
                continue

            # Build/check DIAMOND database
            db_path = diamond_db_dir / f"{genome_id}_helixer"
            if not (db_path.parent / f"{db_path.name}.dmnd").exists():
                if not build_diamond_db(proteome_faa, db_path):
                    continue

            # Run DIAMOND
            diamond_out = diamond_results_dir / f"{locus_id}_{genome_id}.tsv"
            if not diamond_out.exists():
                run_diamond_blastp(
                    flanking_faa, db_path, diamond_out,
                    evalue=args.evalue, threads=args.threads
                )

            # Parse results
            if not diamond_out.exists() or diamond_out.stat().st_size == 0:
                locus_coverage.append({
                    'locus_id': locus_id,
                    'genome': genome_id,
                    'n_flanking_hits': 0,
                    'n_unique_bk_flanking': 0,
                    'coverage': 0,
                    'has_target': genome_id in target_df['genome'].values if not target_df.empty else False,
                    'empty_block_count': 0,
                })
                continue

            hits_df = pd.read_csv(
                diamond_out, sep='\t', header=None,
                names=['qseqid', 'sseqid', 'pident', 'length', 'evalue', 'bitscore']
            )

            # Parse GFF3 for coordinates
            gff3_data = parse_helixer_gff3_genes(gff3_file)

            # Convert to list of dicts
            hits = hits_df.to_dict('records')

            # Get unique BK flanking genes that hit
            unique_bk_flanking = set(h['qseqid'].split('|')[0] if '|' in h['qseqid'] else h['qseqid'] for h in hits)

            # Cluster hits
            clusters = cluster_flanking_hits(hits, gff3_data, args.max_span_kb)

            # Score and filter clusters
            empty_block_count = 0
            for cluster in clusters:
                # Count unique BK flanking genes in cluster
                cluster_bk_ids = set()
                for h in cluster['hits']:
                    bk_id = h['qseqid'].split('|')[0] if '|' in h['qseqid'] else h['qseqid']
                    cluster_bk_ids.add(bk_id)

                n_matches = len(cluster_bk_ids)
                score = n_matches / n_reference_flanking if n_reference_flanking > 0 else 0

                if n_matches < args.min_flanking or score < args.min_score:
                    continue

                # Check for target overlap (but don't filter - just flag it)
                has_overlap, overlapping_target = check_target_overlap(
                    cluster, target_df, genome_id
                )

                # Record ALL synteny blocks (empty OR with target overlap)
                empty_block_count += 1
                all_empty_blocks.append({
                    'locus_id': locus_id,
                    'genome': genome_id,
                    'scaffold': cluster['scaffold'],
                    'start': cluster['start'],
                    'end': cluster['end'],
                    'span_kb': round((cluster['end'] - cluster['start']) / 1000, 2),
                    'n_flanking_matches': n_matches,
                    'n_reference_flanking': n_reference_flanking,
                    'synteny_score': round(score, 3),
                    'flanking_genes_matched': ';'.join(sorted(cluster_bk_ids)),
                    'gene_family': args.family,
                    'has_target_overlap': has_overlap,
                    'overlapping_target_id': overlapping_target if has_overlap else '',
                })

                # Record hit details
                for h in cluster['hits']:
                    all_hit_details.append({
                        'locus_id': locus_id,
                        'genome': genome_id,
                        'scaffold': cluster['scaffold'],
                        'bk_flanking_id': h['qseqid'],
                        'helixer_gene_id': h['gene_id'],
                        'pident': h['pident'],
                        'evalue': h['evalue'],
                    })

            # Record coverage
            locus_coverage.append({
                'locus_id': locus_id,
                'genome': genome_id,
                'n_flanking_hits': len(hits),
                'n_unique_bk_flanking': len(unique_bk_flanking),
                'coverage': round(len(unique_bk_flanking) / n_reference_flanking, 3) if n_reference_flanking > 0 else 0,
                'has_target': genome_id in target_df[target_df['parent_locus'].str.contains(locus_id, na=False)]['genome'].values if not target_df.empty else False,
                'empty_block_count': empty_block_count,
            })

    # Write outputs
    if all_empty_blocks:
        pd.DataFrame(all_empty_blocks).to_csv(
            args.output_dir / "phase2b_synteny_blocks.tsv", sep='\t', index=False
        )

    if all_hit_details:
        pd.DataFrame(all_hit_details).to_csv(
            args.output_dir / "phase2b_flanking_details.tsv", sep='\t', index=False
        )

    if locus_coverage:
        pd.DataFrame(locus_coverage).to_csv(
            args.output_dir / "locus_coverage_summary.tsv", sep='\t', index=False
        )

    # Summary
    print(f"\n[RESULTS]")
    print(f"  Total synteny blocks found: {len(all_empty_blocks)}")

    if all_empty_blocks:
        blocks_df = pd.DataFrame(all_empty_blocks)
        n_empty = len(blocks_df[~blocks_df['has_target_overlap']])
        n_concordant = len(blocks_df[blocks_df['has_target_overlap']])
        print(f"    Empty (no target): {n_empty}")
        print(f"    Concordant (has target): {n_concordant}")

        print(f"\nBlocks by locus:")
        for locus_id in blocks_df['locus_id'].unique():
            locus_df = blocks_df[blocks_df['locus_id'] == locus_id]
            n_total = len(locus_df)
            n_with_target = len(locus_df[locus_df['has_target_overlap']])
            print(f"  {locus_id}: {n_total} blocks ({n_with_target} concordant)")

        print(f"\nBlocks by genome (top 10):")
        genome_counts = blocks_df['genome'].value_counts().head(10)
        for genome, count in genome_counts.items():
            print(f"  {genome}: {count} loci")

    print(f"\n[OUTPUT] {args.output_dir}/phase2b_synteny_blocks.tsv")
    print(f"[OUTPUT] {args.output_dir}/locus_coverage_summary.tsv")

    return 0


if __name__ == "__main__":
    exit(main())
