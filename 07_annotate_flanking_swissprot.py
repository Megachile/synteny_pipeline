#!/usr/bin/env python3
"""
Phase 7 (Helixer): Annotate flanking genes with SwissProt.

UNIONS flanking proteins from both detection methods:
- Phase 4b/5 (target-first): Flanking genes extracted around target hits
- Phase 2b (flanking-first): Flanking genes found by searching BK flanking proteins

Then runs DIAMOND against SwissProt to get functional annotations for ALL
Helixer flanking proteins. Phase 8 displays these annotations in matrices.

Inputs:
- Phase 4b flanking proteins (flanking_proteins.faa)
- Phase 5 flanking matches (flanking_matches.tsv)
- Phase 2b flanking details (phase2b_flanking_details.tsv)
- Helixer proteomes (to extract Phase 2b proteins)
- SwissProt DIAMOND database

Outputs:
- flanking_swissprot_annotations.tsv: SwissProt annotations for all flanking genes
- flanking_matches_annotated.tsv: Phase 5 matches with SwissProt annotations
- phase2b_flanking_annotated.tsv: Phase 2b flanking with SwissProt annotations
- combined_flanking_proteins.faa: Union FASTA for reference

The key output columns used by Phase 8:
- helixer_flanking: The Helixer gene ID
- bk_xp: The BK protein it matched
- swissprot_name: Short functional name (e.g., "Cytochrome P450 3A31")
- swissprot_pident: Percent identity to SwissProt hit
- source: "phase5" or "phase2b"
"""

from __future__ import annotations

import argparse
import re
import subprocess
from pathlib import Path
from typing import Dict, Set

import pandas as pd
from Bio import SeqIO

from constants import (
    FILE_PHASE2B_FLANKING_DETAILS,
    FILE_EMPTY_BLOCK_DETAILS_LEGACY,
)


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


def load_fasta_ids(fasta_path: Path) -> Set[str]:
    """Load protein IDs from a FASTA file."""
    ids = set()
    if fasta_path.exists():
        for record in SeqIO.parse(fasta_path, 'fasta'):
            ids.add(record.id)
    return ids


def extract_proteins_from_helixer(
    protein_ids: Set[str],
    helixer_dir: Path,
    output_fasta: Path,
) -> Dict[str, str]:
    """
    Extract specific proteins from Helixer proteomes.

    Handles ID format differences:
    - Phase 2b stores: Genome_Scaffold_Number (without .1)
    - FASTA has: Genome_Scaffold_Number.1 (with .1)

    Returns dict of protein_id -> sequence (using original input IDs as keys).
    """
    # Group protein IDs by genome
    genome_proteins = {}
    for pid in protein_ids:
        # Phase 2b IDs are like: Aulacidea_tavakolii_CM021343.1_RagTag_000668
        # Try to match genome by directory name
        for genome_dir in helixer_dir.iterdir():
            if genome_dir.is_dir() and pid.startswith(genome_dir.name + '_'):
                if genome_dir.name not in genome_proteins:
                    genome_proteins[genome_dir.name] = set()
                genome_proteins[genome_dir.name].add(pid)
                break

    sequences = {}

    for genome, pids in genome_proteins.items():
        proteome_file = helixer_dir / genome / f"{genome}_helixer_proteins.faa"
        if not proteome_file.exists():
            continue

        # Build lookup: base_id -> full_id with .1 suffix
        # Phase 2b ID: Genome_Scaffold_Number
        # FASTA ID: Genome_Scaffold_Number.1
        pids_with_suffix = {f"{pid}.1" for pid in pids}
        pid_map = {f"{pid}.1": pid for pid in pids}  # Map FASTA ID back to original

        for record in SeqIO.parse(proteome_file, 'fasta'):
            if record.id in pids_with_suffix:
                # Clean sequence for DIAMOND
                seq = str(record.seq).replace('.', 'X').replace('*', '')
                # Use original Phase 2b ID as key (without .1)
                original_id = pid_map[record.id]
                sequences[original_id] = seq

    # Write output FASTA
    with open(output_fasta, 'w') as f:
        for pid, seq in sequences.items():
            f.write(f">{pid}\n{seq}\n")

    return sequences


def main():
    parser = argparse.ArgumentParser(
        description="Phase 7 (Helixer): Annotate flanking genes with SwissProt"
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
    parser.add_argument("--helixer-dir", type=Path, default=None,
                        help="Helixer proteomes directory (required if --phase2b-dir provided)")
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

    # =========================================================================
    # 1. Load Phase 4b/5 flanking proteins (target-first method)
    # =========================================================================
    phase5_flanking_ids = set()
    phase5_sequences = {}

    flanking_faa = args.phase4b_dir / "flanking_proteins.faa"
    if flanking_faa.exists():
        for record in SeqIO.parse(flanking_faa, 'fasta'):
            phase5_flanking_ids.add(record.id)
            seq = str(record.seq).replace('.', 'X').replace('*', '')
            phase5_sequences[record.id] = seq
        print(f"Phase 4b/5 flanking proteins: {len(phase5_flanking_ids)}")
    else:
        print(f"[WARNING] Phase 4b flanking proteins not found: {flanking_faa}")

    # Load flanking matches for BK mapping
    flanking_matches_df = pd.DataFrame()
    flanking_matches_file = args.phase5_dir / "flanking_matches.tsv"
    if flanking_matches_file.exists():
        flanking_matches_df = pd.read_csv(flanking_matches_file, sep='\t')
        print(f"Phase 5 flanking matches: {len(flanking_matches_df)}")

    # =========================================================================
    # 2. Load Phase 2b flanking proteins (flanking-first method)
    # =========================================================================
    phase2b_flanking_ids = set()
    phase2b_details_df = pd.DataFrame()

    if args.phase2b_dir:
        # Try new filename first, then legacy filename for backward compatibility
        phase2b_file = args.phase2b_dir / FILE_PHASE2B_FLANKING_DETAILS
        if not phase2b_file.exists():
            phase2b_file = args.phase2b_dir / FILE_EMPTY_BLOCK_DETAILS_LEGACY

        if phase2b_file.exists():
            phase2b_details_df = pd.read_csv(phase2b_file, sep='\t')
            phase2b_flanking_ids = set(phase2b_details_df['helixer_gene_id'].dropna().unique())
            print(f"Phase 2b flanking proteins: {len(phase2b_flanking_ids)}")
        else:
            print(f"[WARNING] Phase 2b flanking details not found")

    # =========================================================================
    # 2b. Load Phase 5b novel loci flanking proteins
    # =========================================================================
    phase5b_flanking_ids = set()
    phase5b_sequences = {}
    novel_loci_flanking_info = []  # Track which candidate each flanking belongs to

    if args.phase5b_dir and args.phase5b_dir.exists():
        # Look for NOVEL_*_flanking.faa files
        novel_faa_files = list(args.phase5b_dir.glob("NOVEL_*_flanking.faa"))
        print(f"\nPhase 5b novel loci flanking files: {len(novel_faa_files)}")

        for faa_file in novel_faa_files:
            # Extract candidate_id from filename: NOVEL_{...}_flanking.faa
            candidate_id = faa_file.stem.replace('_flanking', '')

            for record in SeqIO.parse(faa_file, 'fasta'):
                # ID format: Genome_Scaffold_Gene.1|Genome_Scaffold_Gene U{pos} NOVEL_{...}
                protein_id = record.id.split('|')[0]  # Get first part before |
                phase5b_flanking_ids.add(protein_id)
                seq = str(record.seq).replace('.', 'X').replace('*', '')
                phase5b_sequences[protein_id] = seq

                # Track info for output
                novel_loci_flanking_info.append({
                    'helixer_gene_id': protein_id,
                    'candidate_id': candidate_id,
                    'full_header': record.description,
                })

        print(f"  Loaded {len(phase5b_flanking_ids)} unique flanking proteins from novel loci")

    # =========================================================================
    # 3. Union flanking proteins and identify what needs extraction
    # =========================================================================
    all_flanking_ids = phase5_flanking_ids | phase2b_flanking_ids | phase5b_flanking_ids
    print(f"\nTotal unique flanking proteins (union): {len(all_flanking_ids)}")

    # Find Phase 2b proteins not already in Phase 5 FASTA
    need_extraction = phase2b_flanking_ids - phase5_flanking_ids
    print(f"  From Phase 5 only: {len(phase5_flanking_ids - phase2b_flanking_ids)}")
    print(f"  From Phase 2b only: {len(need_extraction)}")
    print(f"  In both (overlap): {len(phase5_flanking_ids & phase2b_flanking_ids)}")

    # =========================================================================
    # 4. Extract Phase 2b proteins from Helixer proteomes
    # =========================================================================
    phase2b_sequences = {}
    if need_extraction and args.helixer_dir:
        print(f"\nExtracting {len(need_extraction)} Phase 2b proteins from Helixer proteomes...")
        phase2b_fasta = args.output_dir / "phase2b_flanking_proteins.faa"
        phase2b_sequences = extract_proteins_from_helixer(
            need_extraction, args.helixer_dir, phase2b_fasta
        )
        print(f"  Extracted: {len(phase2b_sequences)}")
        missing = need_extraction - set(phase2b_sequences.keys())
        if missing:
            print(f"  [WARNING] Could not extract: {len(missing)}")
    elif need_extraction:
        print(f"[WARNING] Need --helixer-dir to extract Phase 2b proteins")

    # =========================================================================
    # 5. Combine all flanking proteins into single FASTA
    # =========================================================================
    all_sequences = {**phase5_sequences, **phase2b_sequences, **phase5b_sequences}
    combined_fasta = args.output_dir / "combined_flanking_proteins.faa"

    with open(combined_fasta, 'w') as f:
        for pid, seq in all_sequences.items():
            f.write(f">{pid}\n{seq}\n")
    print(f"\nCombined FASTA: {len(all_sequences)} proteins")

    # =========================================================================
    # 6. Run DIAMOND against SwissProt
    # =========================================================================
    if not args.swissprot_db.exists():
        print(f"[ERROR] SwissProt database not found: {args.swissprot_db}")
        return 1

    print(f"\nRunning DIAMOND against SwissProt...")
    diamond_out = args.output_dir / "diamond_swissprot.tsv"

    success = run_diamond_swissprot(
        combined_fasta, diamond_out, args.swissprot_db,
        evalue=args.evalue, threads=args.threads
    )

    if not success:
        print("[ERROR] DIAMOND failed")
        return 1

    # =========================================================================
    # 7. Parse SwissProt results
    # =========================================================================
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

        # Track source
        sources = []
        if flanking_id in phase5_flanking_ids:
            sources.append('phase5')
        if flanking_id in phase2b_flanking_ids:
            sources.append('phase2b')
        if flanking_id in phase5b_flanking_ids:
            sources.append('phase5b')
        source = '+'.join(sources) if sources else 'unknown'

        annotations[flanking_id] = {
            'helixer_flanking': flanking_id,
            'source': source,
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

    # =========================================================================
    # 8. Save outputs
    # =========================================================================

    # All annotations
    if annotations:
        annot_df = pd.DataFrame(list(annotations.values()))
        annot_df.to_csv(
            args.output_dir / "flanking_swissprot_annotations.tsv",
            sep='\t', index=False
        )
        print(f"\n[OUTPUT] flanking_swissprot_annotations.tsv: {len(annot_df)} annotations")

        # Count by source
        source_counts = annot_df['source'].value_counts()
        for src, cnt in source_counts.items():
            print(f"  {src}: {cnt}")

    # Phase 5 matches with annotations
    if not flanking_matches_df.empty:
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
        print(f"[OUTPUT] flanking_matches_annotated.tsv: {len(enhanced_df)} rows")

        n_with_annot = len([m for m in enhanced_matches if m['swissprot_name']])
        print(f"  {n_with_annot} have SwissProt annotations")

    # Phase 2b flanking with annotations
    if not phase2b_details_df.empty:
        phase2b_annotated = []
        for _, row in phase2b_details_df.iterrows():
            record = row.to_dict()
            helixer_id = row['helixer_gene_id']

            if helixer_id in annotations:
                annot = annotations[helixer_id]
                record['swissprot_name'] = annot['swissprot_name_short']
                record['swissprot_pident'] = annot['swissprot_pident']
            else:
                record['swissprot_name'] = ''
                record['swissprot_pident'] = 0

            phase2b_annotated.append(record)

        phase2b_annot_df = pd.DataFrame(phase2b_annotated)
        phase2b_annot_df.to_csv(
            args.output_dir / "phase2b_flanking_annotated.tsv",
            sep='\t', index=False
        )
        print(f"[OUTPUT] phase2b_flanking_annotated.tsv: {len(phase2b_annot_df)} rows")

        n_with_annot = len([r for r in phase2b_annotated if r['swissprot_name']])
        print(f"  {n_with_annot} have SwissProt annotations")

    # Phase 5b novel loci flanking with annotations
    if novel_loci_flanking_info:
        phase5b_annotated = []
        for record in novel_loci_flanking_info:
            helixer_id = record['helixer_gene_id']
            out_record = record.copy()

            if helixer_id in annotations:
                annot = annotations[helixer_id]
                out_record['swissprot_name'] = annot['swissprot_name_short']
                out_record['swissprot_pident'] = annot['swissprot_pident']
                out_record['swissprot_id'] = annot['swissprot_id']
            else:
                out_record['swissprot_name'] = ''
                out_record['swissprot_pident'] = 0
                out_record['swissprot_id'] = ''

            phase5b_annotated.append(out_record)

        phase5b_annot_df = pd.DataFrame(phase5b_annotated)
        phase5b_annot_df.to_csv(
            args.output_dir / "novel_loci_flanking_annotated.tsv",
            sep='\t', index=False
        )
        print(f"[OUTPUT] novel_loci_flanking_annotated.tsv: {len(phase5b_annot_df)} rows")

        n_with_annot = len([r for r in phase5b_annotated if r['swissprot_name']])
        print(f"  {n_with_annot} have SwissProt annotations")

    print("\nPhase 7 complete.")
    return 0


if __name__ == "__main__":
    exit(main())
