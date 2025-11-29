#!/usr/bin/env python3
"""
Phase 5 (Helixer): Compare flanking genes to BK/LB reference and assign synteny.

For each target gene's flanking genes, run DIAMOND against the BK/LB flanking
gene proteins to find orthologs. Score synteny based on flanking gene matches.

Outputs:
- syntenic_targets.tsv      # Targets with synteny scores
- synteny_details.tsv       # Per-flanking-gene match details
"""

from __future__ import annotations

import argparse
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Set

import pandas as pd


def run_diamond_blastp(
    query_faa: Path,
    db_path: Path,
    output_tsv: Path,
    *,
    evalue: float = 1e-5,
    threads: int = 4,
) -> bool:
    """Run DIAMOND blastp."""
    cmd = [
        "diamond", "blastp",
        "--query", str(query_faa),
        "--db", str(db_path),
        "--out", str(output_tsv),
        "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "evalue", "bitscore",
        "--max-target-seqs", "1",  # Best hit only
        "--evalue", str(evalue),
        "--threads", str(threads),
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.returncode == 0


def build_diamond_db(faa_path: Path, db_path: Path) -> bool:
    """Build DIAMOND database."""
    cmd = [
        "diamond", "makedb",
        "--in", str(faa_path),
        "--db", str(db_path),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.returncode == 0


def load_bk_lb_flanking_genes(phase1_dir: Path) -> Dict[str, Set[str]]:
    """
    Load BK/LB flanking gene protein IDs for each locus.
    Returns dict: locus_id -> set of protein IDs
    """
    locus_flanking = defaultdict(set)

    # Find all flanking FASTA files
    for faa_file in phase1_dir.glob("*_flanking.faa"):
        locus_id = faa_file.stem.replace("_flanking", "")

        with open(faa_file) as f:
            for line in f:
                if line.startswith('>'):
                    protein_id = line[1:].split()[0]
                    locus_flanking[locus_id].add(protein_id)

    return locus_flanking


def calculate_synteny_score(
    target_flanking_matches: List[str],
    reference_flanking: Set[str],
    n_flanking: int,
) -> float:
    """
    Calculate synteny score as fraction of flanking genes that match reference.

    Score = (matched flanking genes) / (total flanking positions checked)
    """
    if n_flanking == 0:
        return 0.0

    matches = sum(1 for m in target_flanking_matches if m in reference_flanking)
    return matches / (n_flanking * 2)  # Both upstream and downstream


def main():
    parser = argparse.ArgumentParser(
        description="Phase 5 (Helixer): Classify targets by synteny"
    )
    parser.add_argument(
        "--family", required=True,
        help="Gene family name"
    )
    parser.add_argument(
        "--phase4-targets", type=Path, required=True,
        help="Phase 4 Helixer all_target_loci.tsv"
    )
    parser.add_argument(
        "--phase4b-flanking", type=Path, required=True,
        help="Phase 4b flanking_genes.tsv"
    )
    parser.add_argument(
        "--phase4b-proteins", type=Path, required=True,
        help="Phase 4b flanking_proteins.faa"
    )
    parser.add_argument(
        "--phase1-dir", type=Path, required=True,
        help="Phase 1 directory with BK/LB flanking FASTAs"
    )
    parser.add_argument(
        "--output-dir", type=Path, required=True,
        help="Output directory"
    )
    parser.add_argument(
        "--n-flanking", type=int, default=10,
        help="Number of flanking genes used (default: 10)"
    )
    parser.add_argument(
        "--min-synteny", type=float, default=0.1,
        help="Minimum synteny score to be classified as syntenic (default: 0.1)"
    )
    parser.add_argument(
        "--threads", type=int, default=4,
        help="DIAMOND threads (default: 4)"
    )

    args = parser.parse_args()

    print("=" * 80)
    print("PHASE 5 (HELIXER): CLASSIFY TARGETS BY SYNTENY")
    print("=" * 80)
    print()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Load Phase 4 targets
    targets_df = pd.read_csv(args.phase4_targets, sep='\t')
    print(f"Loaded {len(targets_df)} targets from Phase 4")

    # Load Phase 4b flanking genes
    flanking_df = pd.read_csv(args.phase4b_flanking, sep='\t')
    print(f"Loaded {len(flanking_df)} flanking records from Phase 4b")

    # Load BK/LB reference flanking genes
    locus_flanking = load_bk_lb_flanking_genes(args.phase1_dir)
    print(f"Loaded flanking genes for {len(locus_flanking)} loci")

    # Combine all reference flanking proteins into one FASTA for DIAMOND DB
    ref_proteins_faa = args.output_dir / "reference_flanking.faa"
    ref_protein_to_locus = {}

    with open(ref_proteins_faa, 'w') as f:
        for locus_id, proteins in locus_flanking.items():
            for faa_file in args.phase1_dir.glob(f"{locus_id}_flanking.faa"):
                with open(faa_file) as src:
                    content = src.read()
                    f.write(content)
                    # Track which locus each protein belongs to
                    for line in content.split('\n'):
                        if line.startswith('>'):
                            protein_id = line[1:].split()[0]
                            ref_protein_to_locus[protein_id] = locus_id

    print(f"Built reference FASTA with {len(ref_protein_to_locus)} proteins")

    # Build DIAMOND database
    ref_db = args.output_dir / "reference_flanking"
    print("Building DIAMOND database...")
    if not build_diamond_db(ref_proteins_faa, ref_db):
        print("[ERROR] Failed to build DIAMOND database")
        return 1

    # Run DIAMOND search
    diamond_out = args.output_dir / "flanking_vs_reference.tsv"
    print("Running DIAMOND search...")
    if not run_diamond_blastp(
        args.phase4b_proteins, ref_db, diamond_out,
        threads=args.threads
    ):
        print("[ERROR] DIAMOND search failed")
        return 1

    # Parse DIAMOND results
    diamond_columns = ['qseqid', 'sseqid', 'pident', 'length', 'evalue', 'bitscore']
    if diamond_out.stat().st_size > 0:
        diamond_df = pd.read_csv(diamond_out, sep='\t', names=diamond_columns, header=None)
        print(f"Found {len(diamond_df)} DIAMOND hits")
    else:
        diamond_df = pd.DataFrame(columns=diamond_columns)
        print("No DIAMOND hits found")

    # Map flanking proteins to their best reference match
    flanking_to_ref = {}
    for _, hit in diamond_df.iterrows():
        qid = hit['qseqid']
        sid = hit['sseqid']
        if qid not in flanking_to_ref:  # Keep best hit only
            flanking_to_ref[qid] = {
                'ref_protein': sid,
                'ref_locus': ref_protein_to_locus.get(sid, 'unknown'),
                'pident': hit['pident'],
                'evalue': hit['evalue'],
            }

    print(f"Mapped {len(flanking_to_ref)} flanking proteins to references")

    # Calculate synteny score per target
    results = []
    details = []

    for _, target in targets_df.iterrows():
        target_gene_id = target['target_gene_id']
        parent_locus = target['parent_locus']

        # Get this target's flanking genes
        target_flanking = flanking_df[flanking_df['target_gene_id'] == target_gene_id]

        # Get reference flanking set for parent locus
        ref_flanking = locus_flanking.get(parent_locus, set())

        # Check matches
        matched_count = 0
        total_count = len(target_flanking)

        for _, flank in target_flanking.iterrows():
            flank_protein = flank['flanking_protein_id']
            match_info = flanking_to_ref.get(flank_protein, {})

            matched_ref = match_info.get('ref_protein', '')
            matched_locus = match_info.get('ref_locus', '')

            # Record detail
            detail = {
                'target_gene_id': target_gene_id,
                'flanking_protein_id': flank_protein,
                'position': flank['position'],
                'matched_ref_protein': matched_ref,
                'matched_ref_locus': matched_locus,
                'is_same_locus': matched_locus == parent_locus,
                'pident': match_info.get('pident', 0),
                'evalue': match_info.get('evalue', 1),
            }
            details.append(detail)

            if matched_locus == parent_locus:
                matched_count += 1

        # Calculate synteny score
        synteny_score = matched_count / total_count if total_count > 0 else 0

        # Classify
        classification = 'syntenic' if synteny_score >= args.min_synteny else 'unplaceable'

        result = target.to_dict()
        result['synteny_score'] = synteny_score
        result['flanking_matches'] = matched_count
        result['flanking_total'] = total_count
        result['classification'] = classification
        results.append(result)

    # Write outputs
    results_df = pd.DataFrame(results)

    # Syntenic targets
    syntenic = results_df[results_df['classification'] == 'syntenic']
    syntenic_file = args.output_dir / "syntenic_targets.tsv"
    syntenic.to_csv(syntenic_file, sep='\t', index=False)
    print(f"\n[OUTPUT] {len(syntenic)} syntenic targets → {syntenic_file}")

    # Unplaceable targets
    unplaceable = results_df[results_df['classification'] == 'unplaceable']
    unplaceable_file = args.output_dir / "unplaceable_targets.tsv"
    unplaceable.to_csv(unplaceable_file, sep='\t', index=False)
    print(f"[OUTPUT] {len(unplaceable)} unplaceable targets → {unplaceable_file}")

    # Details
    details_df = pd.DataFrame(details)
    details_file = args.output_dir / "synteny_details.tsv"
    details_df.to_csv(details_file, sep='\t', index=False)
    print(f"[OUTPUT] {len(details_df)} detail records → {details_file}")

    print("\nPhase 5 (Helixer) complete.")
    return 0


if __name__ == "__main__":
    exit(main())
