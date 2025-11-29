#!/usr/bin/env python3
"""
Phase 7 (Helixer): Annotate flanking genes with SwissProt.

Takes the Helixer flanking proteins from Phase 4b/5 and runs DIAMOND against
SwissProt to get functional annotations. These annotations are displayed in
the Phase 8 matrices to visually validate synteny.

Inputs:
- Phase 4b flanking proteins (flanking_proteins.faa)
- Phase 5 flanking matches (flanking_matches.tsv) - tells us which flanking
  genes matched which BK proteins
- SwissProt DIAMOND database

Outputs:
- flanking_swissprot_annotations.tsv: SwissProt annotations for all flanking genes
- flanking_annotations_by_target.tsv: Annotations organized by target for Phase 8

The key output columns used by Phase 8:
- helixer_flanking: The Helixer gene ID
- bk_xp: The BK protein it matched (from Phase 5)
- swissprot_name: Short functional name (e.g., "Cytochrome P450 3A31")
- swissprot_pident: Percent identity to SwissProt hit
"""

from __future__ import annotations

import argparse
import re
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, Optional

import pandas as pd


DEFAULT_SWISSPROT_DB = Path("/carc/scratch/projects/emartins/2016456/adam/databases/uniprot_sprot.dmnd")


def run_diamond_swissprot(
    query_fasta: Path,
    output_tsv: Path,
    swissprot_db: Path,
    evalue: float = 1e-5,
    threads: int = 8,
) -> bool:
    """Run DIAMOND blastp against SwissProt database."""
    cmd = [
        'diamond', 'blastp',
        '--query', str(query_fasta),
        '--db', str(swissprot_db).replace('.dmnd', ''),
        '--out', str(output_tsv),
        '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length', 'evalue', 'bitscore', 'stitle',
        '--evalue', str(evalue),
        '--max-target-seqs', '1',
        '--threads', str(threads),
        '--quiet'
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  [ERROR] DIAMOND failed: {result.stderr}")
        return False
    return True


def parse_swissprot_title(stitle: str) -> Dict[str, str]:
    """
    Parse SwissProt hit title to extract gene name and description.

    Example titles:
    - "sp|P29981|CP4C1_BLADI Cytochrome P450 4C1 OS=Blaberus discoidalis..."
    - "RecName: Full=Cytochrome P450 3A31; AltName: Full=CYPIIIA31"
    """
    result = {
        'swissprot_id': '',
        'swissprot_gene': '',
        'swissprot_name': '',
        'swissprot_organism': '',
    }

    if not stitle:
        return result

    # Try to parse sp|ACCESSION|GENE_ORGANISM format
    if stitle.startswith('sp|') or stitle.startswith('tr|'):
        parts = stitle.split('|')
        if len(parts) >= 3:
            result['swissprot_id'] = parts[1]
            gene_rest = parts[2]
            if '_' in gene_rest:
                gene_parts = gene_rest.split(' ', 1)
                if '_' in gene_parts[0]:
                    result['swissprot_gene'] = gene_parts[0].split('_')[0]
                if len(gene_parts) > 1:
                    desc = gene_parts[1]
                    # Extract name before OS=
                    if ' OS=' in desc:
                        result['swissprot_name'] = desc.split(' OS=')[0].strip()
                        org_rest = desc.split(' OS=')[1]
                        if ' OX=' in org_rest:
                            result['swissprot_organism'] = org_rest.split(' OX=')[0].strip()
                        else:
                            result['swissprot_organism'] = org_rest.strip()
                    else:
                        result['swissprot_name'] = desc.strip()

    # Try RecName format
    elif 'RecName:' in stitle:
        match = re.search(r'RecName: Full=([^;]+)', stitle)
        if match:
            result['swissprot_name'] = match.group(1).strip()

    # Fallback: just use the title
    else:
        result['swissprot_name'] = stitle[:50].strip()

    return result


def shorten_name(name: str, max_len: int = 40) -> str:
    """Shorten SwissProt name for display in matrix."""
    if not name:
        return ''

    # Remove common prefixes
    name = re.sub(r'^(Probable|Putative|Uncharacterized)\s+', '', name, flags=re.IGNORECASE)

    # Truncate if too long
    if len(name) > max_len:
        name = name[:max_len-3] + '...'

    return name


def main():
    parser = argparse.ArgumentParser(
        description="Phase 7 (Helixer): Annotate flanking genes with SwissProt"
    )
    parser.add_argument("--family", required=True, help="Gene family name")
    parser.add_argument("--phase4b-dir", type=Path, required=True,
                        help="Phase 4b directory with flanking_proteins.faa")
    parser.add_argument("--phase5-dir", type=Path, required=True,
                        help="Phase 5 directory with flanking_matches.tsv")
    parser.add_argument("--swissprot-db", type=Path, default=DEFAULT_SWISSPROT_DB,
                        help="SwissProt DIAMOND database")
    parser.add_argument("--output-dir", type=Path, required=True,
                        help="Output directory")
    parser.add_argument("--evalue", type=float, default=1e-5,
                        help="E-value threshold (default: 1e-5)")
    parser.add_argument("--threads", type=int, default=8,
                        help="DIAMOND threads (default: 8)")

    args = parser.parse_args()

    print("=" * 80)
    print("PHASE 7 (HELIXER): ANNOTATE FLANKING GENES WITH SWISSPROT")
    print("=" * 80)
    print()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Check inputs
    flanking_faa = args.phase4b_dir / "flanking_proteins.faa"
    if not flanking_faa.exists():
        print(f"[ERROR] Flanking proteins not found: {flanking_faa}")
        return 1

    flanking_matches_file = args.phase5_dir / "flanking_matches.tsv"
    if not flanking_matches_file.exists():
        print(f"[WARNING] flanking_matches.tsv not found, will annotate all flanking genes")
        flanking_matches_df = pd.DataFrame()
    else:
        flanking_matches_df = pd.read_csv(flanking_matches_file, sep='\t')
        print(f"Loaded {len(flanking_matches_df)} flanking matches from Phase 5")

    if not args.swissprot_db.exists():
        print(f"[ERROR] SwissProt database not found: {args.swissprot_db}")
        return 1

    # Run DIAMOND against SwissProt
    print(f"\nRunning DIAMOND against SwissProt...")
    diamond_out = args.output_dir / "diamond_swissprot.tsv"

    success = run_diamond_swissprot(
        flanking_faa, diamond_out, args.swissprot_db,
        evalue=args.evalue, threads=args.threads
    )

    if not success:
        print("[ERROR] DIAMOND failed")
        return 1

    # Parse results
    try:
        hits_df = pd.read_csv(
            diamond_out, sep='\t', header=None,
            names=['qseqid', 'sseqid', 'pident', 'length', 'evalue', 'bitscore', 'stitle']
        )
        print(f"DIAMOND hits: {len(hits_df)}")
    except pd.errors.EmptyDataError:
        hits_df = pd.DataFrame()
        print("DIAMOND hits: 0")

    # Build annotation lookup
    annotations = {}
    for _, hit in hits_df.iterrows():
        flanking_id = hit['qseqid']
        parsed = parse_swissprot_title(hit['stitle'])

        annotations[flanking_id] = {
            'helixer_flanking': flanking_id,
            'swissprot_id': parsed['swissprot_id'],
            'swissprot_gene': parsed['swissprot_gene'],
            'swissprot_name': parsed['swissprot_name'],
            'swissprot_name_short': shorten_name(parsed['swissprot_name']),
            'swissprot_organism': parsed['swissprot_organism'],
            'swissprot_pident': round(hit['pident'], 1),
            'swissprot_evalue': hit['evalue'],
            'swissprot_bitscore': round(hit['bitscore'], 1),
        }

    print(f"Annotated {len(annotations)} flanking genes")

    # Save all annotations
    if annotations:
        annot_df = pd.DataFrame(list(annotations.values()))
        annot_df.to_csv(
            args.output_dir / "flanking_swissprot_annotations.tsv",
            sep='\t', index=False
        )

    # Merge with flanking matches for Phase 8
    if not flanking_matches_df.empty:
        # Add SwissProt annotations to flanking matches
        enhanced_matches = []
        for _, row in flanking_matches_df.iterrows():
            match_dict = row.to_dict()
            flanking_id = row['helixer_flanking']

            if flanking_id in annotations:
                annot = annotations[flanking_id]
                match_dict['swissprot_name'] = annot['swissprot_name_short']
                match_dict['swissprot_pident'] = annot['swissprot_pident']
            else:
                match_dict['swissprot_name'] = ''
                match_dict['swissprot_pident'] = 0

            enhanced_matches.append(match_dict)

        enhanced_df = pd.DataFrame(enhanced_matches)
        enhanced_df.to_csv(
            args.output_dir / "flanking_matches_annotated.tsv",
            sep='\t', index=False
        )
        print(f"\n[OUTPUT] flanking_matches_annotated.tsv: {len(enhanced_df)} rows")

        # Show sample annotations
        n_with_annot = len([m for m in enhanced_matches if m['swissprot_name']])
        print(f"  {n_with_annot} flanking matches have SwissProt annotations")

        if n_with_annot > 0:
            print("\nSample annotations:")
            for m in enhanced_matches[:5]:
                if m['swissprot_name']:
                    print(f"  {m['helixer_flanking'][:40]} -> {m['bk_xp']}")
                    print(f"    SwissProt: {m['swissprot_name']} ({m['swissprot_pident']}%)")

    print(f"\n[OUTPUT] flanking_swissprot_annotations.tsv: {len(annotations)} annotations")
    print("\nPhase 7 complete.")
    return 0


if __name__ == "__main__":
    exit(main())
