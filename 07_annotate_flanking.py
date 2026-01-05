#!/usr/bin/env python3
"""
Phase 7 (Helixer): Annotate flanking genes using NR annotations from FASTA headers.

SIMPLIFIED VERSION: Uses annotations already present in Helixer proteome FASTA headers
(from DIAMOND NR search) instead of running a separate SwissProt search.

The Helixer proteome headers contain NR annotations like:
  >Genome_Scaffold_Gene.1 Cytochrome P450 3A31
  >Genome_Scaffold_Gene.1 hypothetical protein
  >Genome_Scaffold_Gene.1 Uncharacterized protein LOC117175202

This script:
1. Loads flanking proteins from Phase 4b/5 and Phase 2b
2. Parses NR annotations directly from FASTA headers
3. Adds chromosome placement info (on BK-mapped chromosome vs unplaced)
4. Outputs annotated tables for Phase 8

No DIAMOND search needed - annotations are already there!

Inputs:
- Phase 4b flanking proteins (flanking_proteins.faa)
- Phase 5 flanking matches (flanking_matches.tsv)
- Phase 2b flanking details (phase2b_flanking_details.tsv)
- Helixer proteomes (for NR annotations in headers)

Outputs:
- flanking_annotations.tsv: All flanking genes with NR annotations + chromosome info
- flanking_matches_annotated.tsv: Phase 5 matches with annotations
- phase2b_flanking_annotated.tsv: Phase 2b flanking with annotations
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, Set, Tuple

import pandas as pd
from Bio import SeqIO

from constants import (
    FILE_PHASE2B_FLANKING_DETAILS,
    FILE_EMPTY_BLOCK_DETAILS_LEGACY,
)


def parse_nr_annotation(description: str) -> Tuple[str, str]:
    """
    Parse NR annotation from FASTA description line.

    Input: "Genome_Scaffold_Gene.1 Cytochrome P450 3A31"
    Returns: (protein_id, annotation)

    Also handles:
    - "hypothetical protein" -> "hypothetical"
    - "Uncharacterized protein LOC..." -> "uncharacterized LOC..."
    - Long descriptions -> truncated to 60 chars
    """
    parts = description.split(None, 1)  # Split on first whitespace
    protein_id = parts[0]

    if len(parts) < 2:
        annotation = ""
    else:
        annotation = parts[1].strip()

    # Clean up common patterns
    annotation = re.sub(r'^hypothetical protein$', 'hypothetical', annotation, flags=re.I)
    annotation = re.sub(r'^Uncharacterized protein ', 'unchar ', annotation, flags=re.I)
    annotation = re.sub(r'^Putative ', '', annotation, flags=re.I)
    annotation = re.sub(r'^Probable ', '', annotation, flags=re.I)

    # Truncate if too long
    if len(annotation) > 60:
        annotation = annotation[:57] + '...'

    return protein_id, annotation


def is_on_chromosome(scaffold: str) -> bool:
    """
    Check if a scaffold is a BK-mapped chromosome.

    BK chromosomes are named CM02* (e.g., CM021338.1)
    Unplaced scaffolds have other names (e.g., JAJGSE010000001.1)
    """
    # CM02 scaffolds are BK-mapped chromosomes
    return scaffold.startswith('CM02') if scaffold else False


def extract_scaffold_from_gene_id(gene_id: str) -> str:
    """
    Extract scaffold name from Helixer gene ID.

    Format: Genome_Scaffold_GeneNum.1
    Example: Belonocnema_kinseyi_GCF_CM021338.1_000001.1 -> CM021338.1
             Aulacidea_tavakolii_JAJGSE010000001.1_RagTag_000001.1 -> JAJGSE010000001.1_RagTag

    This is tricky because genome names can have underscores.
    We look for scaffold patterns: CM*, JAJGSE*, scaffold*, etc.
    """
    # Remove .1 suffix if present
    if gene_id.endswith('.1'):
        gene_id = gene_id[:-2]

    parts = gene_id.split('_')

    # Find the scaffold part - it's usually the part that looks like CM02*, JAJGSE*, etc.
    for i, part in enumerate(parts):
        if part.startswith('CM02') or part.startswith('JAJGSE') or part.startswith('scaffold'):
            # Include RagTag suffix if present
            scaffold_parts = [part]
            if i + 1 < len(parts) and parts[i + 1] == 'RagTag':
                scaffold_parts.append('RagTag')
            return '_'.join(scaffold_parts)

    # Fallback: can't determine scaffold
    return ''


def load_proteome_annotations(helixer_dir: Path) -> Dict[str, Dict[str, str]]:
    """
    Load NR annotations from all Helixer proteome FASTA files.

    Returns dict: protein_id -> {annotation, scaffold, on_chromosome}
    """
    annotations = {}

    for genome_dir in helixer_dir.iterdir():
        if not genome_dir.is_dir():
            continue

        # Find proteome file
        proteome_file = None
        for pattern in [f"{genome_dir.name}_helixer_proteins.faa",
                        f"{genome_dir.name}_proteins.faa"]:
            candidate = genome_dir / pattern
            if candidate.exists():
                proteome_file = candidate
                break

        if not proteome_file:
            continue

        # Parse all protein headers
        for record in SeqIO.parse(proteome_file, 'fasta'):
            protein_id, annotation = parse_nr_annotation(record.description)
            scaffold = extract_scaffold_from_gene_id(protein_id)

            annotations[protein_id] = {
                'nr_annotation': annotation,
                'scaffold': scaffold,
                'on_chromosome': is_on_chromosome(scaffold),
            }

            # Also store without .1 suffix for matching Phase 2b IDs
            if protein_id.endswith('.1'):
                base_id = protein_id[:-2]
                annotations[base_id] = annotations[protein_id]

    return annotations


def main():
    parser = argparse.ArgumentParser(
        description="Phase 7 (Helixer): Annotate flanking genes using NR annotations from headers"
    )
    parser.add_argument("--family", required=True, help="Gene family name")
    parser.add_argument("--phase4b-dir", type=Path, required=True,
                        help="Phase 4b directory with flanking_proteins.faa")
    parser.add_argument("--phase5-dir", type=Path, required=True,
                        help="Phase 5 directory with flanking_matches.tsv")
    parser.add_argument("--phase2b-dir", type=Path, default=None,
                        help="Phase 2b directory with phase2b_flanking_details.tsv (optional)")
    parser.add_argument("--phase5b-dir", type=Path, default=None,
                        help="Phase 5b directory with novel loci flanking FAAs (optional)")
    parser.add_argument("--helixer-dir", type=Path, required=True,
                        help="Helixer proteomes directory (contains NR annotations in headers)")
    parser.add_argument("--output-dir", type=Path, required=True,
                        help="Output directory")

    args = parser.parse_args()

    print("=" * 80)
    print("PHASE 7 (HELIXER): ANNOTATE FLANKING GENES (NR FROM HEADERS)")
    print("=" * 80)
    print()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # =========================================================================
    # 1. Load all NR annotations from Helixer proteomes
    # =========================================================================
    print("Loading NR annotations from Helixer proteome headers...")
    all_annotations = load_proteome_annotations(args.helixer_dir)
    print(f"  Loaded annotations for {len(all_annotations)} proteins")

    # Count chromosome placement
    on_chrom = sum(1 for a in all_annotations.values() if a['on_chromosome'])
    print(f"  On BK-mapped chromosomes: {on_chrom}")
    print(f"  Unplaced: {len(all_annotations) - on_chrom}")

    # =========================================================================
    # 2. Load Phase 4b/5 flanking proteins
    # =========================================================================
    phase5_flanking_ids = set()

    flanking_faa = args.phase4b_dir / "flanking_proteins.faa"
    if flanking_faa.exists():
        for record in SeqIO.parse(flanking_faa, 'fasta'):
            phase5_flanking_ids.add(record.id)
        print(f"\nPhase 4b/5 flanking proteins: {len(phase5_flanking_ids)}")
    else:
        print(f"[WARNING] Phase 4b flanking proteins not found: {flanking_faa}")

    # Load flanking matches
    flanking_matches_df = pd.DataFrame()
    flanking_matches_file = args.phase5_dir / "flanking_matches.tsv"
    if flanking_matches_file.exists():
        flanking_matches_df = pd.read_csv(flanking_matches_file, sep='\t')
        print(f"Phase 5 flanking matches: {len(flanking_matches_df)}")

    # =========================================================================
    # 3. Load Phase 2b flanking proteins
    # =========================================================================
    phase2b_flanking_ids = set()
    phase2b_details_df = pd.DataFrame()

    if args.phase2b_dir:
        phase2b_file = args.phase2b_dir / FILE_PHASE2B_FLANKING_DETAILS
        if not phase2b_file.exists():
            phase2b_file = args.phase2b_dir / FILE_EMPTY_BLOCK_DETAILS_LEGACY

        if phase2b_file.exists():
            phase2b_details_df = pd.read_csv(phase2b_file, sep='\t')
            phase2b_flanking_ids = set(phase2b_details_df['helixer_gene_id'].dropna().unique())
            print(f"Phase 2b flanking proteins: {len(phase2b_flanking_ids)}")

    # =========================================================================
    # 4. Load Phase 5b novel loci flanking
    # =========================================================================
    phase5b_flanking_ids = set()
    novel_loci_flanking_info = []

    if args.phase5b_dir and args.phase5b_dir.exists():
        novel_faa_files = list(args.phase5b_dir.glob("NOVEL_*_flanking.faa"))
        print(f"\nPhase 5b novel loci flanking files: {len(novel_faa_files)}")

        for faa_file in novel_faa_files:
            candidate_id = faa_file.stem.replace('_flanking', '')

            for record in SeqIO.parse(faa_file, 'fasta'):
                protein_id = record.id.split('|')[0]
                phase5b_flanking_ids.add(protein_id)
                novel_loci_flanking_info.append({
                    'helixer_gene_id': protein_id,
                    'candidate_id': candidate_id,
                    'full_header': record.description,
                })

        print(f"  Loaded {len(phase5b_flanking_ids)} unique flanking proteins from novel loci")

    # =========================================================================
    # 5. Build unified annotation table
    # =========================================================================
    all_flanking_ids = phase5_flanking_ids | phase2b_flanking_ids | phase5b_flanking_ids
    print(f"\nTotal unique flanking proteins: {len(all_flanking_ids)}")

    annotation_records = []
    matched = 0

    for fid in all_flanking_ids:
        # Determine source
        sources = []
        if fid in phase5_flanking_ids:
            sources.append('phase5')
        if fid in phase2b_flanking_ids:
            sources.append('phase2b')
        if fid in phase5b_flanking_ids:
            sources.append('phase5b')

        # Get annotation
        if fid in all_annotations:
            annot = all_annotations[fid]
            matched += 1
        else:
            annot = {'nr_annotation': '', 'scaffold': '', 'on_chromosome': False}

        # Use swissprot_name for backward compatibility with Phase 8
        # (even though it's NR annotation, not SwissProt)
        annotation_records.append({
            'helixer_flanking': fid,
            'source': '+'.join(sources),
            'swissprot_name': annot['nr_annotation'],  # Phase 8 expects this column name
            'swissprot_pident': 100.0 if annot['nr_annotation'] else 0,  # Placeholder for Phase 8
            'nr_annotation': annot['nr_annotation'],
            'scaffold': annot['scaffold'],
            'on_chromosome': annot['on_chromosome'],
        })

    print(f"  Matched annotations: {matched}/{len(all_flanking_ids)}")

    # =========================================================================
    # 6. Save outputs
    # =========================================================================

    # All annotations
    if annotation_records:
        annot_df = pd.DataFrame(annotation_records)
        annot_df.to_csv(
            args.output_dir / "flanking_annotations.tsv",
            sep='\t', index=False
        )
        print(f"\n[OUTPUT] flanking_annotations.tsv: {len(annot_df)} rows")

        # Stats
        n_annotated = annot_df['nr_annotation'].apply(lambda x: bool(x and x != 'hypothetical')).sum()
        n_on_chrom = annot_df['on_chromosome'].sum()
        print(f"  With functional annotation: {n_annotated}")
        print(f"  On BK-mapped chromosome: {n_on_chrom}")

    # Phase 5 matches with annotations
    if not flanking_matches_df.empty:
        annot_lookup = {r['helixer_flanking']: r for r in annotation_records}

        enhanced_matches = []
        for _, row in flanking_matches_df.iterrows():
            match_dict = row.to_dict()
            flanking_id = row['helixer_flanking']

            if flanking_id in annot_lookup:
                annot = annot_lookup[flanking_id]
                match_dict['swissprot_name'] = annot['swissprot_name']
                match_dict['swissprot_pident'] = annot['swissprot_pident']
                match_dict['nr_annotation'] = annot['nr_annotation']
                match_dict['on_chromosome'] = annot['on_chromosome']
            else:
                match_dict['swissprot_name'] = ''
                match_dict['swissprot_pident'] = 0
                match_dict['nr_annotation'] = ''
                match_dict['on_chromosome'] = False

            enhanced_matches.append(match_dict)

        enhanced_df = pd.DataFrame(enhanced_matches)
        enhanced_df.to_csv(
            args.output_dir / "flanking_matches_annotated.tsv",
            sep='\t', index=False
        )
        print(f"[OUTPUT] flanking_matches_annotated.tsv: {len(enhanced_df)} rows")

    # Phase 2b flanking with annotations
    if not phase2b_details_df.empty:
        annot_lookup = {r['helixer_flanking']: r for r in annotation_records}

        phase2b_annotated = []
        for _, row in phase2b_details_df.iterrows():
            record = row.to_dict()
            helixer_id = row['helixer_gene_id']

            if helixer_id in annot_lookup:
                annot = annot_lookup[helixer_id]
                record['swissprot_name'] = annot['swissprot_name']
                record['swissprot_pident'] = annot['swissprot_pident']
                record['nr_annotation'] = annot['nr_annotation']
                record['on_chromosome'] = annot['on_chromosome']
            else:
                record['swissprot_name'] = ''
                record['swissprot_pident'] = 0
                record['nr_annotation'] = ''
                record['on_chromosome'] = False

            phase2b_annotated.append(record)

        phase2b_annot_df = pd.DataFrame(phase2b_annotated)
        phase2b_annot_df.to_csv(
            args.output_dir / "phase2b_flanking_annotated.tsv",
            sep='\t', index=False
        )
        print(f"[OUTPUT] phase2b_flanking_annotated.tsv: {len(phase2b_annot_df)} rows")

    # Phase 5b novel loci flanking with annotations
    if novel_loci_flanking_info:
        annot_lookup = {r['helixer_flanking']: r for r in annotation_records}

        phase5b_annotated = []
        for record in novel_loci_flanking_info:
            helixer_id = record['helixer_gene_id']
            out_record = record.copy()

            if helixer_id in annot_lookup:
                annot = annot_lookup[helixer_id]
                out_record['swissprot_name'] = annot['swissprot_name']
                out_record['swissprot_pident'] = annot['swissprot_pident']
                out_record['nr_annotation'] = annot['nr_annotation']
                out_record['on_chromosome'] = annot['on_chromosome']
            else:
                out_record['swissprot_name'] = ''
                out_record['swissprot_pident'] = 0
                out_record['nr_annotation'] = ''
                out_record['on_chromosome'] = False

            phase5b_annotated.append(out_record)

        phase5b_annot_df = pd.DataFrame(phase5b_annotated)
        phase5b_annot_df.to_csv(
            args.output_dir / "novel_loci_flanking_annotated.tsv",
            sep='\t', index=False
        )
        print(f"[OUTPUT] novel_loci_flanking_annotated.tsv: {len(phase5b_annot_df)} rows")

    print("\nPhase 7 complete.")
    return 0


if __name__ == "__main__":
    exit(main())
