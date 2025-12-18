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
    """
    Parse GFF3 to get gene coordinates.

    Handles both:
    - Helixer GFF3: gene IDs match protein IDs directly
    - NCBI GFF: gene IDs are gene-LOC..., protein IDs are XP_... in CDS features

    Returns dict keyed by both gene_id (for Helixer) and protein_id (for NCBI).
    """
    genes = {}
    gene_by_name = {}  # gene= attribute -> gene coordinates

    # First pass: collect gene features
    with open(gff3_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            feature_type = parts[2]
            if feature_type != 'gene':
                continue

            scaffold, _, _, start, end, _, strand, _, attributes = parts

            attr_dict = {}
            for attr in attributes.split(';'):
                if '=' in attr:
                    key, value = attr.split('=', 1)
                    attr_dict[key] = value

            gene_id = attr_dict.get('ID', '')
            gene_name = attr_dict.get('gene', '')  # LOC... number for NCBI

            coords = {
                'scaffold': scaffold,
                'start': int(start),
                'end': int(end),
                'strand': strand,
            }

            if gene_id:
                genes[gene_id] = coords
            if gene_name:
                gene_by_name[gene_name] = coords

    # Second pass: map protein IDs from CDS features to gene coordinates (for NCBI GFFs)
    with open(gff3_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            feature_type = parts[2]
            if feature_type != 'CDS':
                continue

            attributes = parts[8]
            attr_dict = {}
            for attr in attributes.split(';'):
                if '=' in attr:
                    key, value = attr.split('=', 1)
                    attr_dict[key] = value

            protein_id = attr_dict.get('protein_id', '')
            gene_name = attr_dict.get('gene', '')  # LOC... number

            # Map protein_id to gene coordinates
            if protein_id and gene_name and gene_name in gene_by_name:
                # Also add without version suffix for matching
                genes[protein_id] = gene_by_name[gene_name]
                if '.' in protein_id:
                    genes[protein_id.rsplit('.', 1)[0]] = gene_by_name[gene_name]

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
    parser.add_argument(
        "--min-total-genomes", type=int, default=3,
        help="Minimum genomes with synteny block to keep locus (default: 3). "
             "Filters out spurious single-genome hits."
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

    # LB proteome (Leptopilina boulardi) - NCBI reference
    lb_proteome = data_dir / "proteomes" / "GCF_019393585.1_Lboulardi_NCBI.faa"
    lb_gff = data_dir / "reference" / "Leptopilina_boulardi" / "ncbi_dataset" / "data" / "GCF_019393585.1" / "genomic.gff"
    if lb_proteome.exists():
        proteomes['LB'] = {
            'proteome': lb_proteome,
            'gff3': lb_gff if lb_gff.exists() else None,
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

            # Check for DIAMOND database - first alongside proteome, then in output dir
            db_in_genome_dir = proteome_faa.parent / f"{genome_id}_helixer.dmnd"
            db_dir = args.output_dir / "diamond_dbs"
            db_dir.mkdir(exist_ok=True)
            db_in_output_dir = db_dir / f"{genome_id}"

            if db_in_genome_dir.exists():
                db_path = proteome_faa.parent / f"{genome_id}_helixer"
            elif (db_in_output_dir.parent / f"{db_in_output_dir.name}.dmnd").exists():
                db_path = db_in_output_dir
            else:
                db_path = db_in_output_dir
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

    # Cluster novel loci based on mutual flanking gene synteny
    if validation_results:
        # Load species mapping for proper naming
        # gca_to_species.tsv is in hybrid_workflow/data/
        # output_dir is outputs/{family}/phase5b_xxx, need to go up 3 levels
        data_dir = args.output_dir.parent.parent.parent / "data"
        species_mapping = load_species_mapping(data_dir)
        print(f"  Loaded species mapping from {data_dir}: {len(species_mapping)} entries")

        # Load chromosome mapping for CM→chr conversion
        # chromosome_mappings.json is in hpc_deployment_package/
        chr_mapping = load_chromosome_mapping(args.output_dir)
        print(f"  Loaded chromosome mapping: {len(chr_mapping)} entries")

        cluster_novel_loci(
            results_df,
            candidates,
            args.output_dir,
            min_synteny=0.15,
            min_total_genomes=args.min_total_genomes,
            species_mapping=species_mapping,
            chr_mapping=chr_mapping,
        )

    return 0


def load_species_mapping(data_dir: Path) -> dict:
    """Load genome ID to species name mapping."""
    mapping_file = data_dir / "gca_to_species.tsv"
    if not mapping_file.exists():
        return {}

    mapping = {}
    df = pd.read_csv(mapping_file, sep='\t')
    for _, row in df.iterrows():
        mapping[row['accession']] = row['species']
    return mapping


def load_chromosome_mapping(base_dir: Path) -> dict:
    """
    Load chromosome accession to chr# mapping.

    Maps both RefSeq (NC_*) and GenBank (CM*) accessions to chromosome numbers.
    Merges mappings from all genomes (BK, T_remus, etc.) into a single flat dict.
    File location: hpc_deployment_package/chromosome_mappings.json
    """
    import json

    # Try multiple locations going up from base_dir
    candidates = [
        base_dir / "chromosome_mappings.json",
        base_dir.parent / "chromosome_mappings.json",
        base_dir.parent.parent / "chromosome_mappings.json",
        base_dir.parent.parent.parent / "chromosome_mappings.json",
        base_dir.parent.parent.parent.parent / "chromosome_mappings.json",
    ]

    for mapping_file in candidates:
        if mapping_file.exists():
            with open(mapping_file) as f:
                data = json.load(f)

            # Merge all genome mappings into single flat dict
            # Skip keys starting with '_' (comments, metadata)
            merged = {}
            for genome_key, genome_mapping in data.items():
                if isinstance(genome_mapping, dict):
                    for acc, chr_name in genome_mapping.items():
                        if not acc.startswith('_'):
                            merged[acc] = chr_name

            return merged

    return {}


def cluster_novel_loci(
    validation_df: pd.DataFrame,
    candidates_df: pd.DataFrame,
    output_dir: Path,
    min_synteny: float = 0.15,
    min_total_genomes: int = 3,
    species_mapping: dict = None,
    chr_mapping: dict = None,
) -> pd.DataFrame:
    """
    Cluster novel locus candidates based on mutual flanking gene synteny.

    Two candidates are the same locus if:
    - Candidate A has high synteny in genome B
    - Candidate B has high synteny in genome A

    Filtering:
    - Requires synteny block present in >= min_total_genomes
    - This filters out spurious single-genome hits

    This mirrors Phase 2 synteny block detection.

    Also tracks empty blocks (genomes with flanking pattern but no target).
    """
    print("\n" + "=" * 70)
    print("CLUSTERING NOVEL LOCI")
    print("=" * 70)

    # Build lookup: candidate_id -> source_genome
    candidate_genomes = {}
    candidate_bitscores = {}
    for _, row in candidates_df.iterrows():
        cid = f"NOVEL_{row['genome']}_{row['scaffold']}_{row['start']}"
        candidate_genomes[cid] = row['genome']
        candidate_bitscores[cid] = row['best_bitscore']

    # Build synteny matrix: which candidates have synteny to which genomes
    # synteny_to[candidate_id] = {genome: score}
    synteny_to = {}
    for cid in validation_df['candidate_id'].unique():
        cdf = validation_df[validation_df['candidate_id'] == cid]
        synteny_to[cid] = {}
        for _, row in cdf.iterrows():
            if row['synteny_score'] >= min_synteny:
                synteny_to[cid][row['target_genome']] = row['synteny_score']

    # Build graph of mutual synteny connections
    # Two candidates are connected if A has synteny in B's genome AND B has synteny in A's genome
    from collections import defaultdict

    connections = defaultdict(set)
    candidate_ids = list(synteny_to.keys())

    for i, cid_a in enumerate(candidate_ids):
        genome_a = candidate_genomes.get(cid_a)
        if not genome_a:
            continue

        for cid_b in candidate_ids[i+1:]:
            genome_b = candidate_genomes.get(cid_b)
            if not genome_b:
                continue

            # Check mutual synteny
            a_has_b = genome_b in synteny_to.get(cid_a, {})
            b_has_a = genome_a in synteny_to.get(cid_b, {})

            if a_has_b and b_has_a:
                connections[cid_a].add(cid_b)
                connections[cid_b].add(cid_a)

    # Find connected components (each component = one locus)
    visited = set()
    clusters = []

    def dfs(node, cluster):
        if node in visited:
            return
        visited.add(node)
        cluster.add(node)
        for neighbor in connections[node]:
            dfs(neighbor, cluster)

    for cid in candidate_ids:
        if cid not in visited:
            cluster = set()
            dfs(cid, cluster)
            clusters.append(cluster)

    print(f"\nFound {len(clusters)} unique novel loci from {len(candidate_ids)} candidates")

    # For each cluster, determine:
    # - Representative (highest bitscore)
    # - Member genomes (have target)
    # - Empty genomes (have synteny but no target in candidates)
    # - Locus name (from representative's species + scaffold)

    locus_results = []

    # Track base name usage for suffix assignment (a, b, c...)
    base_name_count = {}

    for i, cluster in enumerate(clusters):
        # Get representative (highest bitscore)
        cluster_list = list(cluster)
        rep = max(cluster_list, key=lambda x: candidate_bitscores.get(x, 0))

        # Member genomes (have targets)
        member_genomes = set()
        for cid in cluster:
            genome = candidate_genomes.get(cid)
            if genome:
                member_genomes.add(genome)

        # Genomes with synteny to this locus (from any member)
        all_syntenic_genomes = set()
        for cid in cluster:
            for genome in synteny_to.get(cid, {}):
                all_syntenic_genomes.add(genome)

        # Empty genomes = syntenic but not members
        empty_genomes = all_syntenic_genomes - member_genomes

        # Get representative info for naming
        rep_genome = candidate_genomes.get(rep, 'unknown')

        # Get species name from mapping or use genome ID
        if species_mapping and rep_genome in species_mapping:
            species_name = species_mapping[rep_genome].replace(' ', '_')
        elif '_' in rep_genome and not rep_genome.startswith('GCA'):
            species_name = rep_genome  # Already a species name like Alloxysta_arcuata
        else:
            species_name = rep_genome.replace('GCA_', '').replace('.', '_')

        # Extract scaffold/chromosome from candidate ID (NOVEL_genome_scaffold_position)
        # Parse from representative candidate ID
        rep_parts = rep.split('_')
        rep_scaffold = 'unknown'
        for j, part in enumerate(rep_parts):
            if part.startswith('CM') or part.startswith('NC') or part.startswith('NODE'):
                # Include version suffix if present (e.g., CM021342.1)
                if j + 1 < len(rep_parts) and rep_parts[j + 1] in ['1', 'RagTag']:
                    if rep_parts[j + 1] == 'RagTag':
                        rep_scaffold = f"{part}.1"  # Add .1 for CM codes
                    else:
                        rep_scaffold = f"{part}.{rep_parts[j + 1]}"
                else:
                    rep_scaffold = part
                break

        # Convert scaffold accession to chromosome number using mapping
        # e.g., CM021342.1 → chr5
        chr_name = rep_scaffold
        if chr_mapping:
            # Try with and without version suffix
            chr_name = chr_mapping.get(rep_scaffold,
                       chr_mapping.get(rep_scaffold.split('.')[0] + '.1',
                       rep_scaffold.replace('_RagTag', '').replace('.1', '')))

        # Create base locus name (species_chr)
        base_locus_name = f"{species_name}_{chr_name}"

        # Assign suffix (a, b, c...) for multiple loci on same chromosome
        suffix_idx = base_name_count.get(base_locus_name, 0)
        suffix = chr(ord('a') + suffix_idx)  # 0=a, 1=b, 2=c, etc.
        base_name_count[base_locus_name] = suffix_idx + 1

        locus_name = f"{base_locus_name}_{suffix}"

        locus_results.append({
            'locus_id': f"NOVEL_{i+1}",
            'locus_name': locus_name,
            'representative': rep,
            'n_members': len(member_genomes),
            'n_empty': len(empty_genomes),
            'n_total_genomes': len(member_genomes) + len(empty_genomes),
            'member_genomes': ';'.join(sorted(member_genomes)),
            'empty_genomes': ';'.join(sorted(empty_genomes)),
            'rep_bitscore': candidate_bitscores.get(rep, 0),
        })

        print(f"\n  NOVEL_{i+1}: {locus_name}")
        print(f"    Members: {len(member_genomes)} genomes with targets")
        print(f"    Empty: {len(empty_genomes)} genomes with synteny only")
        print(f"    Representative: {rep} (bitscore {candidate_bitscores.get(rep, 0):.0f})")

    # Filter loci by minimum total genomes (block conservation)
    # A genuine novel locus should have synteny block present in multiple genomes
    # even if the target is only in one (novel insertion vs spurious hit)
    loci_df = pd.DataFrame(locus_results)

    if not loci_df.empty:
        n_before = len(loci_df)
        loci_df = loci_df[loci_df['n_total_genomes'] >= min_total_genomes].copy()
        n_filtered = n_before - len(loci_df)

        if n_filtered > 0:
            print(f"\n[FILTER] Removed {n_filtered} loci with synteny block in < {min_total_genomes} genomes")
            print(f"         (These are likely spurious single-genome hits)")

        # Renumber locus IDs after filtering
        loci_df = loci_df.reset_index(drop=True)
        loci_df['locus_id'] = [f"NOVEL_{i+1}" for i in range(len(loci_df))]

    # Save clustered loci
    loci_df.to_csv(output_dir / "novel_loci_clustered.tsv", sep='\t', index=False)

    print(f"\n[OUTPUT] {output_dir / 'novel_loci_clustered.tsv'}")
    print(f"         {len(loci_df)} novel loci passed filters (block in >= {min_total_genomes} genomes)")

    return loci_df


if __name__ == "__main__":
    exit(main())
