#!/usr/bin/env python3
"""
Phase 1.5 (redesign): BK-LB locus unification

Merges LB loci that are clear orthologs of BK loci based on flanking protein
DIAMOND hits. This creates a BK-anchored locus space where:
- BK loci are canonical
- LB loci with matching BK flanking context are treated as aliases (dropped)
- Only true LB-only loci (no BK counterpart) are retained

Outputs:
- Rewrites locus_definitions.tsv with BK-anchored locus space
- Creates bk_lb_locus_mapping.tsv for diagnostics
- Prints summary statistics
"""

from __future__ import annotations

from pathlib import Path
import argparse
import subprocess
import pandas as pd
from Bio import SeqIO
from collections import defaultdict
import tempfile
import sys


def parse_args():
    p = argparse.ArgumentParser(
        description="Phase 1.5: BK-LB locus unification via flanking DIAMOND"
    )
    p.add_argument(
        "--locus-defs",
        type=Path,
        required=True,
        help="Phase1_v2 locus_definitions.tsv"
    )
    p.add_argument(
        "--lb-proteome-db",
        type=Path,
        required=True,
        help="DIAMOND database for LB proteome (e.g., data/proteomes/Belonocnema_latifrons.dmnd)"
    )
    p.add_argument(
        "--diamond-threads",
        type=int,
        default=8,
        help="Number of threads for DIAMOND"
    )
    p.add_argument(
        "--min-anchors",
        type=int,
        default=3,
        help="Minimum anchor hits to consider LB locus as BK alias (default 3)"
    )
    p.add_argument(
        "--strong-evalue",
        type=float,
        default=1e-20,
        help="E-value threshold for strong 2-anchor matches (default 1e-20)"
    )
    return p.parse_args()


def extract_protein_id(header: str) -> str:
    """Extract protein ID from FASTA header."""
    # Handle formats like "XP_033208249.1|LOC117167432" or "XP_033208249.1"
    token = header.split()[0]
    if "|" in token:
        return token.split("|")[0]
    return token


def build_lb_protein_to_locus_map(lb_loci: pd.DataFrame) -> dict[str, str]:
    """Map LB flanking protein IDs to their locus."""
    lb_protein_to_locus = {}

    for _, row in lb_loci.iterrows():
        locus_id = row['locus_id']
        flanking_file = Path(row['flanking_file'])

        if not flanking_file.exists():
            print(f"  WARNING: Flanking file not found for {locus_id}: {flanking_file}")
            continue

        for rec in SeqIO.parse(str(flanking_file), 'fasta'):
            prot_id = extract_protein_id(rec.id)
            lb_protein_to_locus[prot_id] = locus_id

    return lb_protein_to_locus


def build_bk_flanking_fasta(bk_loci: pd.DataFrame, output_fasta: Path) -> dict[str, str]:
    """
    Build combined BK flanking FASTA with headers encoding locus and anchor position.

    Header format: >BK_chr2_a|U1|XP_033208249.1

    Returns: Map from query_id -> (locus_id, anchor_id)
    """
    query_to_locus = {}
    records = []

    for _, row in bk_loci.iterrows():
        locus_id = row['locus_id']
        flanking_file = Path(row['flanking_file'])

        if not flanking_file.exists():
            print(f"  WARNING: Flanking file not found for {locus_id}: {flanking_file}")
            continue

        # Track position in flanking file (U1, U2, ..., I1, ..., D1, D2, ...)
        anchor_idx = 0
        for rec in SeqIO.parse(str(flanking_file), 'fasta'):
            prot_id = extract_protein_id(rec.id)
            anchor_label = f"A{anchor_idx}"  # Generic anchor label
            query_id = f"{locus_id}|{anchor_label}|{prot_id}"

            query_to_locus[query_id] = (locus_id, prot_id)

            # Create new record with modified header
            rec.id = query_id
            rec.description = ""
            records.append(rec)
            anchor_idx += 1

    SeqIO.write(records, str(output_fasta), 'fasta')
    return query_to_locus


def run_diamond(
    query_fasta: Path,
    db: Path,
    output_tsv: Path,
    threads: int,
    evalue: float = 1e-5
) -> bool:
    """Run DIAMOND blastp."""
    cmd = [
        'diamond', 'blastp',
        '--query', str(query_fasta),
        '--db', str(db),
        '--out', str(output_tsv),
        '--outfmt', '6', 'qseqid', 'sseqid', 'bitscore', 'evalue', 'pident',
        '--max-target-seqs', '50',
        '--evalue', str(evalue),
        '--threads', str(threads),
        '--quiet'
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.returncode == 0


def parse_diamond_output(
    diamond_tsv: Path,
    query_to_locus: dict[str, tuple[str, str]],
    lb_protein_to_locus: dict[str, str]
) -> dict[str, dict[str, list]]:
    """
    Parse DIAMOND output and map BK loci to candidate LB loci.

    Returns: bk_to_lb_hits[bk_locus_id][lb_locus_id] = [hit_dicts]
    """
    bk_to_lb_hits = defaultdict(lambda: defaultdict(list))

    with open(diamond_tsv) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 4:
                continue

            query_id = parts[0]
            subject_id = parts[1]
            bitscore = float(parts[2])
            evalue = float(parts[3])

            # Extract BK locus and anchor from query
            if query_id not in query_to_locus:
                continue
            bk_locus_id, bk_anchor_id = query_to_locus[query_id]

            # Extract LB protein and find its locus
            lb_protein = extract_protein_id(subject_id)
            lb_locus_id = lb_protein_to_locus.get(lb_protein)

            if not lb_locus_id:
                continue

            bk_to_lb_hits[bk_locus_id][lb_locus_id].append({
                'bk_locus': bk_locus_id,
                'lb_locus': lb_locus_id,
                'bk_anchor': bk_anchor_id,
                'lb_protein': lb_protein,
                'bitscore': bitscore,
                'evalue': evalue
            })

    return bk_to_lb_hits


def identify_aliases(
    bk_to_lb_hits: dict,
    min_anchors: int,
    strong_evalue: float
) -> dict[str, str]:
    """
    Identify LB loci that are aliases of BK loci.

    Returns: aliases[lb_locus_id] = bk_locus_id
    """
    aliases = {}
    alias_details = []

    for bk_locus, lb_hits in bk_to_lb_hits.items():
        for lb_locus, hits in lb_hits.items():
            # Count unique BK anchors
            unique_anchors = {h['bk_anchor'] for h in hits}
            n_anchors = len(unique_anchors)

            # Get best e-value
            best_evalue = min(h['evalue'] for h in hits)
            total_bitscore = sum(h['bitscore'] for h in hits)

            # Criteria: ≥3 anchors OR ≥2 anchors with strong e-value
            is_alias = (n_anchors >= min_anchors) or (n_anchors >= 2 and best_evalue <= strong_evalue)

            if is_alias:
                # If LB locus already assigned to different BK locus, pick better one
                if lb_locus in aliases:
                    prev_bk = aliases[lb_locus]
                    prev_hits = bk_to_lb_hits[prev_bk][lb_locus]
                    prev_anchors = len({h['bk_anchor'] for h in prev_hits})
                    prev_bitscore = sum(h['bitscore'] for h in prev_hits)

                    # Keep the one with more anchors, or better bitscore if tied
                    if n_anchors > prev_anchors or (n_anchors == prev_anchors and total_bitscore > prev_bitscore):
                        aliases[lb_locus] = bk_locus
                else:
                    aliases[lb_locus] = bk_locus

                alias_details.append({
                    'bk_locus': bk_locus,
                    'lb_locus': lb_locus,
                    'n_anchors': n_anchors,
                    'best_evalue': best_evalue,
                    'total_bitscore': total_bitscore,
                    'anchor_list': ','.join(sorted(unique_anchors))
                })

    return aliases, alias_details


def main():
    args = parse_args()

    print("=" * 80)
    print("PHASE 1.5: BK-LB LOCUS UNIFICATION")
    print("=" * 80)
    print(f"Locus definitions: {args.locus_defs}")
    print(f"LB proteome DB: {args.lb_proteome_db}")
    print()

    # Load locus definitions
    df = pd.read_csv(args.locus_defs, sep='\t')

    # Split BK and LB
    bk = df[df['genome'] == 'BK'].copy()
    lb = df[df['genome'] == 'LB'].copy()

    print(f"[1] Loaded {len(df)} loci:")
    print(f"    BK: {len(bk)} loci")
    print(f"    LB: {len(lb)} loci")
    print()

    if len(lb) == 0:
        print("No LB loci found - nothing to merge. Exiting.")
        return

    if len(bk) == 0:
        print("No BK loci found - LB loci will remain as-is. Exiting.")
        return

    # Build LB protein → locus map
    print("[2] Building LB protein → locus map...")
    lb_protein_to_locus = build_lb_protein_to_locus_map(lb)
    print(f"    Mapped {len(lb_protein_to_locus)} LB proteins to loci")
    print()

    # Build combined BK flanking FASTA
    print("[3] Building combined BK flanking FASTA...")
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        bk_fasta = tmpdir / "bk_flanking_combined.faa"
        diamond_out = tmpdir / "bk_vs_lb_flanking.tsv"

        query_to_locus = build_bk_flanking_fasta(bk, bk_fasta)
        print(f"    Created {bk_fasta} with {len(query_to_locus)} BK anchors")
        print()

        # Run DIAMOND
        print("[4] Running DIAMOND: BK flanking → LB proteome...")
        if not run_diamond(bk_fasta, args.lb_proteome_db, diamond_out, args.diamond_threads):
            print("ERROR: DIAMOND failed")
            sys.exit(1)
        print(f"    DIAMOND complete: {diamond_out}")
        print()

        # Parse DIAMOND output
        print("[5] Mapping DIAMOND hits to LB loci...")
        bk_to_lb_hits = parse_diamond_output(diamond_out, query_to_locus, lb_protein_to_locus)

        total_hits = sum(len(lb_dict) for lb_dict in bk_to_lb_hits.values())
        print(f"    Found {total_hits} BK→LB locus connections")
        print()

        # Identify aliases
        print("[6] Identifying LB loci that are BK aliases...")
        print(f"    Criteria: ≥{args.min_anchors} anchors OR ≥2 anchors with e-value ≤{args.strong_evalue}")
        aliases, alias_details = identify_aliases(
            bk_to_lb_hits,
            args.min_anchors,
            args.strong_evalue
        )
        print(f"    Found {len(aliases)} LB loci that are BK aliases")
        print()

        # Save alias details for diagnostics
        if alias_details:
            alias_tsv = args.locus_defs.parent / "bk_lb_locus_mapping.tsv"
            pd.DataFrame(alias_details).to_csv(alias_tsv, sep='\t', index=False)
            print(f"    Saved mapping details: {alias_tsv}")
            print()

    # Rewrite locus_definitions.tsv
    print("[7] Rewriting locus_definitions.tsv (BK-anchored)...")

    # Add lb_aliases column to BK loci
    if 'lb_aliases' not in df.columns:
        df['lb_aliases'] = ''

    bk_alias_map = defaultdict(list)
    for lb_locus, bk_locus in aliases.items():
        bk_alias_map[bk_locus].append(lb_locus)

    for bk_locus, alias_list in bk_alias_map.items():
        mask = df['locus_id'] == bk_locus
        if mask.any():
            df.loc[mask, 'lb_aliases'] = ','.join(sorted(alias_list))

    # Drop LB alias rows
    alias_lb_ids = set(aliases.keys())
    df_merged = df[~df['locus_id'].isin(alias_lb_ids)].copy()

    # Save
    df_merged.to_csv(args.locus_defs, sep='\t', index=False)
    print(f"    Saved merged locus definitions: {args.locus_defs}")
    print()

    # Summary
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)

    bk_final = df_merged[df_merged['genome'] == 'BK']
    lb_final = df_merged[df_merged['genome'] == 'LB']

    print(f"BK loci (canonical):        {len(bk_final)}")
    print(f"  With LB aliases:          {len([1 for _, row in bk_final.iterrows() if row.get('lb_aliases')])}")
    print(f"LB-only loci (no BK match): {len(lb_final)}")
    print(f"LB loci merged as aliases:  {len(alias_lb_ids)}")
    print()
    print(f"Total loci in merged space: {len(df_merged)} (was {len(df)})")
    print()
    print("BK-anchored locus space created successfully!")
    print("=" * 80)


if __name__ == "__main__":
    main()
