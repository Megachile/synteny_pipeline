#!/usr/bin/env python3
"""
Phase 5b (Helixer): Validate novel locus candidates.

For unplaceable targets that don't match any known locus flanking pattern,
extract their flanking genes and search for this pattern across all genomes
(including BK/LB references) to validate as a true novel locus.

Uses SAME parameters as Phase 2b for consistency:
- DIAMOND: evalue=1e-5, max_target_seqs=5
- Clustering: max_span=500kb
- Scoring: n_matches / n_reference_flanking

Inputs:
- Phase 5 unplaceable_targets.tsv (to identify novel candidates)
- Helixer proteomes and GFF3 files
- BK/LB reference proteomes (from Phase 1)

Outputs:
- novel_candidates.tsv: Summary of novel candidates
- novel_locus_validation.tsv: Cross-genome validation results
- {candidate}_flanking.faa: Flanking proteins for each candidate
"""

from __future__ import annotations

import argparse
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd


# Match Phase 2b parameters exactly
DIAMOND_EVALUE = 1e-5
DIAMOND_MAX_TARGETS = 5
MAX_SPAN_KB = 500.0
MIN_FLANKING = 3


def run_diamond_blastp(
    query_faa: Path,
    db_path: Path,
    output_tsv: Path,
    threads: int = 4,
) -> bool:
    """Run DIAMOND blastp with Phase 2b parameters."""
    cmd = [
        'diamond', 'blastp',
        '--query', str(query_faa),
        '--db', str(db_path).replace('.dmnd', ''),
        '--out', str(output_tsv),
        '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length', 'evalue', 'bitscore',
        '--evalue', str(DIAMOND_EVALUE),
        '--max-target-seqs', str(DIAMOND_MAX_TARGETS),
        '--threads', str(threads),
        '--quiet'
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.returncode == 0


def build_diamond_db(proteome_faa: Path, db_path: Path) -> bool:
    """Build DIAMOND database."""
    db_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        'diamond', 'makedb',
        '--in', str(proteome_faa),
        '--db', str(db_path),
        '--quiet'
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.returncode == 0


def parse_gff3_genes(gff3_path: Path) -> Dict[str, Dict]:
    """Parse GFF3 to get gene coordinates."""
    genes = {}
    with open(gff3_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9 or parts[2] != 'gene':
                continue

            scaffold, _, _, start, end, _, strand, _, attributes = parts

            attr_dict = {}
            for attr in attributes.split(';'):
                if '=' in attr:
                    key, value = attr.split('=', 1)
                    attr_dict[key] = value

            gene_id = attr_dict.get('ID', '')
            if gene_id:
                genes[gene_id] = {
                    'scaffold': scaffold,
                    'start': int(start),
                    'end': int(end),
                    'strand': strand,
                }
    return genes


def extract_flanking_genes(
    target_scaffold: str,
    target_start: int,
    target_end: int,
    gff3_genes: Dict[str, Dict],
    proteome_faa: Path,
    n_upstream: int = 25,
    n_downstream: int = 25,
) -> Tuple[List[Dict], List[Dict]]:
    """
    Extract flanking genes around a target position.

    Returns (upstream_genes, downstream_genes) lists with gene info.
    """
    # Get genes on same scaffold, sorted by position
    scaffold_genes = [
        {'gene_id': gid, **info}
        for gid, info in gff3_genes.items()
        if info['scaffold'] == target_scaffold
    ]
    scaffold_genes.sort(key=lambda x: x['start'])

    # Find target position in gene order
    target_midpoint = (target_start + target_end) / 2

    # Find index of first gene after target
    insert_idx = 0
    for i, gene in enumerate(scaffold_genes):
        if gene['start'] > target_midpoint:
            insert_idx = i
            break
        insert_idx = i + 1

    # Extract upstream (genes before target)
    upstream = []
    for i in range(max(0, insert_idx - n_upstream), insert_idx):
        gene = scaffold_genes[i]
        # Skip if this IS the target gene
        if gene['start'] <= target_end and gene['end'] >= target_start:
            continue
        upstream.append(gene)

    # Extract downstream (genes after target)
    downstream = []
    for i in range(insert_idx, min(len(scaffold_genes), insert_idx + n_downstream)):
        gene = scaffold_genes[i]
        if gene['start'] <= target_end and gene['end'] >= target_start:
            continue
        downstream.append(gene)

    return upstream, downstream


def create_flanking_fasta(
    upstream_genes: List[Dict],
    downstream_genes: List[Dict],
    proteome_faa: Path,
    output_faa: Path,
    candidate_id: str,
) -> int:
    """
    Create FASTA of flanking proteins with U/D labels.

    Format matches Phase 1: >XP_ID|LOC_ID U{N} ...
    """
    # Load proteome sequences
    sequences = {}
    current_id = None
    current_seq = []

    with open(proteome_faa) as f:
        for line in f:
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]
                # Also store by gene ID (remove .1 suffix)
                if '.' in current_id:
                    gene_id = current_id.rsplit('.', 1)[0]
                    sequences[gene_id] = None  # Will be set below
                current_seq = []
            else:
                current_seq.append(line.strip())
        if current_id:
            sequences[current_id] = ''.join(current_seq)

    # Write flanking FASTA
    n_written = 0
    with open(output_faa, 'w') as f:
        # Upstream genes (U1 is closest to target, UN is furthest)
        for i, gene in enumerate(reversed(upstream_genes)):
            pos_label = f"U{len(upstream_genes) - i}"
            gene_id = gene['gene_id']
            protein_id = f"{gene_id}.1"

            seq = sequences.get(protein_id) or sequences.get(gene_id)
            if seq:
                f.write(f">{protein_id}|{gene_id} {pos_label} {candidate_id}\n")
                f.write(f"{seq}\n")
                n_written += 1

        # Downstream genes (D1 is closest to target)
        for i, gene in enumerate(downstream_genes):
            pos_label = f"D{i + 1}"
            gene_id = gene['gene_id']
            protein_id = f"{gene_id}.1"

            seq = sequences.get(protein_id) or sequences.get(gene_id)
            if seq:
                f.write(f">{protein_id}|{gene_id} {pos_label} {candidate_id}\n")
                f.write(f"{seq}\n")
                n_written += 1

    return n_written


def cluster_hits(
    hits: List[Dict],
    gff3_genes: Dict[str, Dict],
    max_span_kb: float = MAX_SPAN_KB,
) -> List[Dict]:
    """
    Cluster flanking gene hits by scaffold/position.
    Same algorithm as Phase 2b.
    """
    if not hits:
        return []

    # Map protein IDs to gene positions
    def get_gene_pos(protein_id):
        gene_id = protein_id.rsplit('.', 1)[0] if '.' in protein_id else protein_id
        return gff3_genes.get(gene_id, {})

    # Group by scaffold
    by_scaffold = defaultdict(list)
    for hit in hits:
        gene_info = get_gene_pos(hit['sseqid'])
        if gene_info:
            by_scaffold[gene_info['scaffold']].append({
                **hit,
                'start': gene_info['start'],
                'end': gene_info['end'],
            })

    max_span_bp = max_span_kb * 1000
    clusters = []

    for scaffold, scaffold_hits in by_scaffold.items():
        scaffold_hits.sort(key=lambda x: x['start'])

        # Greedy clustering
        current_cluster = [scaffold_hits[0]]
        cluster_start = scaffold_hits[0]['start']
        cluster_end = scaffold_hits[0]['end']

        for hit in scaffold_hits[1:]:
            new_end = max(cluster_end, hit['end'])
            new_start = min(cluster_start, hit['start'])

            if (new_end - new_start) <= max_span_bp:
                current_cluster.append(hit)
                cluster_start = new_start
                cluster_end = new_end
            else:
                if len(current_cluster) >= MIN_FLANKING:
                    clusters.append({
                        'scaffold': scaffold,
                        'start': cluster_start,
                        'end': cluster_end,
                        'hits': current_cluster,
                    })
                current_cluster = [hit]
                cluster_start = hit['start']
                cluster_end = hit['end']

        if len(current_cluster) >= MIN_FLANKING:
            clusters.append({
                'scaffold': scaffold,
                'start': cluster_start,
                'end': cluster_end,
                'hits': current_cluster,
            })

    return clusters


def identify_novel_candidates(
    unplaceable_file: Path,
    min_bitscore: float = 400,
    min_flanking: int = 5,
) -> pd.DataFrame:
    """
    Identify novel locus candidates from unplaceable targets.

    Criteria:
    - High quality hit (bitscore >= threshold)
    - Has flanking genes (not tiny scaffold)
    - no_locus_matches reason (flanking doesn't match any BK locus)
    """
    df = pd.read_csv(unplaceable_file, sep='\t')

    candidates = df[
        (df['best_bitscore'] >= min_bitscore) &
        (df['n_flanking_genes'] >= min_flanking) &
        (df['unplaceable_reason'] == 'no_locus_matches')
    ].copy()

    return candidates


def main():
    parser = argparse.ArgumentParser(
        description="Phase 5b (Helixer): Validate novel locus candidates"
    )
    parser.add_argument("--family", required=True, help="Gene family name")
    parser.add_argument(
        "--phase5-dir", type=Path, required=True,
        help="Phase 5 output directory (with unplaceable_targets.tsv)"
    )
    parser.add_argument(
        "--helixer-dir", type=Path, required=True,
        help="Directory containing Helixer outputs"
    )
    parser.add_argument(
        "--phase1-dir", type=Path, required=True,
        help="Phase 1 directory (for BK/LB reference proteomes)"
    )
    parser.add_argument(
        "--output-dir", type=Path, required=True,
        help="Output directory"
    )
    parser.add_argument(
        "--min-bitscore", type=float, default=400,
        help="Minimum bitscore for novel candidates (default: 400)"
    )
    parser.add_argument(
        "--n-flanking", type=int, default=25,
        help="Number of flanking genes to extract each side (default: 25)"
    )
    parser.add_argument(
        "--threads", type=int, default=4,
        help="DIAMOND threads"
    )

    args = parser.parse_args()

    print("=" * 70)
    print("PHASE 5B (HELIXER): VALIDATE NOVEL LOCUS CANDIDATES")
    print("=" * 70)
    print()
    print(f"Parameters (matching Phase 2b):")
    print(f"  DIAMOND e-value: {DIAMOND_EVALUE}")
    print(f"  Max cluster span: {MAX_SPAN_KB} kb")
    print(f"  Min flanking genes: {MIN_FLANKING}")
    print()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Identify novel candidates
    unplaceable_file = args.phase5_dir / "unplaceable_targets.tsv"
    if not unplaceable_file.exists():
        print(f"[ERROR] Unplaceable targets file not found: {unplaceable_file}")
        return 1

    candidates = identify_novel_candidates(
        unplaceable_file,
        min_bitscore=args.min_bitscore,
    )

    print(f"[1] Identified {len(candidates)} novel locus candidates")

    if len(candidates) == 0:
        print("No novel candidates to validate.")
        return 0

    # Save candidate summary
    candidates.to_csv(args.output_dir / "novel_candidates.tsv", sep='\t', index=False)

    # Collect all proteomes to search
    proteomes = {}

    # Helixer proteomes
    for genome_dir in args.helixer_dir.iterdir():
        if not genome_dir.is_dir():
            continue
        proteome = genome_dir / f"{genome_dir.name}_helixer_proteins.faa"
        gff3 = genome_dir / f"{genome_dir.name}_helixer.gff3"
        if proteome.exists() and gff3.exists():
            proteomes[genome_dir.name] = {
                'proteome': proteome,
                'gff3': gff3,
                'type': 'helixer'
            }

    # BK/LB reference proteomes - check standard locations
    data_dir = args.helixer_dir.parent  # ../data

    # BK proteome (Belonocnema kinseyi)
    bk_proteome = data_dir / "reference" / "protein.faa"
    if bk_proteome.exists():
        proteomes['BK'] = {
            'proteome': bk_proteome,
            'gff3': data_dir / "reference" / "genomic.gff",  # NCBI GFF
            'type': 'reference'
        }

    # LB proteome (Belonocnema latifrons)
    lb_proteome = data_dir / "proteomes" / "Belonocnema_latifrons.faa"
    if lb_proteome.exists():
        proteomes['LB'] = {
            'proteome': lb_proteome,
            'gff3': None,  # No GFF3 for LB
            'type': 'reference'
        }

    print(f"[2] Will search {len(proteomes)} proteomes")
    print(f"    Helixer: {sum(1 for p in proteomes.values() if p['type'] == 'helixer')}")
    print(f"    Reference: {sum(1 for p in proteomes.values() if p['type'] == 'reference')}")

    # Process each candidate
    validation_results = []

    for idx, (_, candidate) in enumerate(candidates.iterrows()):
        candidate_id = f"NOVEL_{candidate['genome']}_{candidate['scaffold']}_{candidate['start']}"
        print(f"\n[3.{idx+1}] Processing: {candidate_id}")

        # Load source genome GFF3
        source_genome = candidate['genome']
        source_gff3 = args.helixer_dir / source_genome / f"{source_genome}_helixer.gff3"
        source_proteome = args.helixer_dir / source_genome / f"{source_genome}_helixer_proteins.faa"

        if not source_gff3.exists():
            print(f"  [SKIP] GFF3 not found for {source_genome}")
            continue

        gff3_genes = parse_gff3_genes(source_gff3)

        # Extract flanking genes
        upstream, downstream = extract_flanking_genes(
            candidate['scaffold'],
            candidate['start'],
            candidate['end'],
            gff3_genes,
            source_proteome,
            n_upstream=args.n_flanking,
            n_downstream=args.n_flanking,
        )

        print(f"  Extracted {len(upstream)} upstream, {len(downstream)} downstream genes")

        if len(upstream) + len(downstream) < MIN_FLANKING:
            print(f"  [SKIP] Not enough flanking genes")
            continue

        # Create flanking FASTA
        flanking_faa = args.output_dir / f"{candidate_id}_flanking.faa"
        n_written = create_flanking_fasta(
            upstream, downstream,
            source_proteome, flanking_faa,
            candidate_id,
        )
        print(f"  Created flanking FASTA with {n_written} proteins")

        n_reference_flanking = n_written

        # Search each proteome
        for genome_id, genome_info in proteomes.items():
            if genome_id == source_genome:
                continue  # Skip source genome

            proteome_faa = genome_info['proteome']

            # Build/find DIAMOND database
            db_dir = args.output_dir / "diamond_dbs"
            db_dir.mkdir(exist_ok=True)
            db_path = db_dir / f"{genome_id}"

            if not (db_path.parent / f"{db_path.name}.dmnd").exists():
                build_diamond_db(proteome_faa, db_path)

            # Run DIAMOND
            diamond_out = args.output_dir / "diamond_results" / f"{candidate_id}_{genome_id}.tsv"
            diamond_out.parent.mkdir(exist_ok=True)

            if not diamond_out.exists():
                run_diamond_blastp(flanking_faa, db_path, diamond_out, args.threads)

            # Parse results
            if not diamond_out.exists() or diamond_out.stat().st_size == 0:
                validation_results.append({
                    'candidate_id': candidate_id,
                    'source_genome': source_genome,
                    'target_genome': genome_id,
                    'genome_type': genome_info['type'],
                    'n_hits': 0,
                    'n_unique_flanking': 0,
                    'synteny_score': 0,
                    'best_cluster_scaffold': '',
                    'best_cluster_score': 0,
                })
                continue

            hits_df = pd.read_csv(
                diamond_out, sep='\t', header=None,
                names=['qseqid', 'sseqid', 'pident', 'length', 'evalue', 'bitscore']
            )

            # Count unique flanking genes that hit
            unique_flanking = set(h.split('|')[0] for h in hits_df['qseqid'])

            # Cluster hits (if we have GFF3)
            best_score = 0
            best_scaffold = ''

            if genome_info['gff3']:
                target_gff3 = parse_gff3_genes(genome_info['gff3'])
                clusters = cluster_hits(hits_df.to_dict('records'), target_gff3)

                for cluster in clusters:
                    cluster_flanking = set(h['qseqid'].split('|')[0] for h in cluster['hits'])
                    score = len(cluster_flanking) / n_reference_flanking if n_reference_flanking > 0 else 0
                    if score > best_score:
                        best_score = score
                        best_scaffold = cluster['scaffold']
            else:
                # No GFF3 - just use hit count as approximation
                best_score = len(unique_flanking) / n_reference_flanking if n_reference_flanking > 0 else 0

            validation_results.append({
                'candidate_id': candidate_id,
                'source_genome': source_genome,
                'target_genome': genome_id,
                'genome_type': genome_info['type'],
                'n_hits': len(hits_df),
                'n_unique_flanking': len(unique_flanking),
                'synteny_score': round(best_score, 3),
                'best_cluster_scaffold': best_scaffold,
                'best_cluster_score': round(best_score, 3),
            })

        # Summary for this candidate
        candidate_results = [r for r in validation_results if r['candidate_id'] == candidate_id]
        n_with_pattern = sum(1 for r in candidate_results if r['synteny_score'] >= 0.15)
        ref_matches = [r for r in candidate_results if r['genome_type'] == 'reference' and r['synteny_score'] >= 0.15]

        print(f"  Found pattern in {n_with_pattern} genomes")
        if ref_matches:
            print(f"  ** Found in reference genome(s): {[r['target_genome'] for r in ref_matches]}")

    # Save validation results
    if validation_results:
        results_df = pd.DataFrame(validation_results)
        results_df.to_csv(args.output_dir / "novel_locus_validation.tsv", sep='\t', index=False)

        # Summary
        print("\n" + "=" * 70)
        print("VALIDATION SUMMARY")
        print("=" * 70)

        for candidate_id in results_df['candidate_id'].unique():
            cdf = results_df[results_df['candidate_id'] == candidate_id]
            n_positive = len(cdf[cdf['synteny_score'] >= 0.15])
            ref_positive = cdf[(cdf['genome_type'] == 'reference') & (cdf['synteny_score'] >= 0.15)]

            print(f"\n{candidate_id}:")
            print(f"  Genomes with pattern (score >= 0.15): {n_positive}/{len(cdf)}")

            if len(ref_positive) > 0:
                print(f"  ** FOUND IN REFERENCE: {list(ref_positive['target_genome'])}")
                print(f"     -> This locus exists in BK/LB! Consider adding to locus definitions.")
            elif n_positive >= 3:
                print(f"  ** VALIDATED: Pattern found in multiple genomes")
                print(f"     -> Novel locus not present in reference.")
            elif n_positive > 0:
                print(f"  ? WEAK: Pattern found in only {n_positive} genome(s)")
            else:
                print(f"  X NOT VALIDATED: Pattern not found in other genomes")
                print(f"     -> Likely species-specific or artifact")

    print(f"\n[OUTPUT] {args.output_dir}")

    return 0


if __name__ == "__main__":
    exit(main())
