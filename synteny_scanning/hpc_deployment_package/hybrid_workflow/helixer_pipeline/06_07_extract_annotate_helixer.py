#!/usr/bin/env python3
"""
Phase 6-7 (Helixer): Extract target sequences and annotate with SwissProt.

Since Helixer already predicts proteins, we can directly extract target sequences
from the Helixer proteomes and run SwissProt annotation in a single step.

Inputs:
- Phase 4 Helixer targets (all_target_loci.tsv)
- Phase 5 Helixer classification (syntenic_targets.tsv, unplaceable_targets.tsv)
- Helixer proteomes ({genome}_helixer_proteins.faa)

Outputs:
- target_proteins.faa: Combined protein sequences for all targets
- swissprot_annotations.tsv: SwissProt top hits for each target
- annotated_targets.tsv: Targets with SwissProt annotations merged
"""

from __future__ import annotations

import argparse
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, List, Optional, Set

import pandas as pd
from Bio import SeqIO


def find_proteome_file(genome: str, helixer_dir: Path) -> Optional[Path]:
    """Find the Helixer proteome file for a genome."""
    genome_dir = helixer_dir / genome
    if not genome_dir.exists():
        return None

    proteome = genome_dir / f"{genome}_helixer_proteins.faa"
    if proteome.exists():
        return proteome

    return None


def extract_target_sequences(
    targets_df: pd.DataFrame,
    helixer_dir: Path,
    output_fasta: Path,
) -> Dict[str, str]:
    """
    Extract protein sequences for all targets from Helixer proteomes.

    Returns dict mapping target_gene_id to sequence.
    """
    # Group targets by genome
    genome_targets = targets_df.groupby('genome')['target_gene_id'].apply(set).to_dict()

    sequences = {}
    missing = []

    for genome, target_ids in genome_targets.items():
        proteome_file = find_proteome_file(genome, helixer_dir)
        if proteome_file is None:
            print(f"  [WARNING] No proteome found for {genome}")
            missing.extend(target_ids)
            continue

        # Load proteome and extract matching sequences
        found_in_genome = 0
        for record in SeqIO.parse(proteome_file, 'fasta'):
            # Try exact match first
            if record.id in target_ids:
                sequences[record.id] = str(record.seq)
                found_in_genome += 1
            # Also try matching without version suffix (e.g., .1, .2)
            elif record.id.rsplit('.', 1)[0] in target_ids:
                base_id = record.id.rsplit('.', 1)[0]
                sequences[base_id] = str(record.seq)
                found_in_genome += 1

        found_ids = set(sequences.keys())
        genome_missing = target_ids - found_ids
        if genome_missing:
            # Try with .1 suffix
            for target_id in list(genome_missing):
                if f"{target_id}.1" in sequences:
                    sequences[target_id] = sequences[f"{target_id}.1"]
                    genome_missing.discard(target_id)

        if genome_missing:
            print(f"  [WARNING] {genome}: {len(genome_missing)} targets not found in proteome")
            missing.extend(genome_missing)
        else:
            print(f"  {genome}: extracted {found_in_genome} target sequences")

    # Write combined FASTA (clean sequences for DIAMOND compatibility)
    with open(output_fasta, 'w') as f:
        for target_id, seq in sequences.items():
            # Replace invalid characters: . -> X, * -> (remove)
            clean_seq = seq.replace('.', 'X').replace('*', '')
            f.write(f">{target_id}\n{clean_seq}\n")

    print(f"\nExtracted {len(sequences)} sequences, {len(missing)} missing")
    return sequences


def run_diamond_swissprot(
    query_fasta: Path,
    output_xml: Path,
    swissprot_db: Path,
    evalue: float = 1e-5,
    threads: int = 16,
) -> bool:
    """Run DIAMOND blastp against SwissProt database."""
    cmd = [
        'diamond', 'blastp',
        '--query', str(query_fasta),
        '--db', str(swissprot_db),
        '--outfmt', '5',  # XML format
        '--evalue', str(evalue),
        '--max-target-seqs', '1',
        '--threads', str(threads),
        '--quiet'
    ]

    print(f"Running DIAMOND against SwissProt...")
    with open(output_xml, 'w') as f:
        result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, text=True)

    if result.returncode != 0:
        print(f"  [ERROR] DIAMOND failed: {result.stderr}")
        return False

    return True


def parse_swissprot_xml(xml_file: Path) -> Dict[str, Dict]:
    """Parse SwissProt BLAST XML results."""
    annotations = {}

    try:
        tree = ET.parse(xml_file)
        root = tree.getroot()

        for iteration in root.findall('.//Iteration'):
            query_def = iteration.find('Iteration_query-def').text
            query_id = query_def.split()[0]

            # Check for no hits
            message = iteration.find('Iteration_message')
            if message is not None and message.text == 'No hits found':
                annotations[query_id] = {
                    'target_gene_id': query_id,
                    'swissprot_id': '',
                    'swissprot_gene': '',
                    'swissprot_description': '',
                    'swissprot_organism': '',
                    'percent_identity': 0.0,
                    'alignment_length': 0,
                    'evalue': 1.0,
                    'bitscore': 0.0,
                }
                continue

            # Get first hit
            hit = iteration.find('.//Hit')
            if hit is None:
                annotations[query_id] = {
                    'target_gene_id': query_id,
                    'swissprot_id': '',
                    'swissprot_gene': '',
                    'swissprot_description': '',
                    'swissprot_organism': '',
                    'percent_identity': 0.0,
                    'alignment_length': 0,
                    'evalue': 1.0,
                    'bitscore': 0.0,
                }
                continue

            hit_id = hit.find('Hit_id').text
            hit_def = hit.find('Hit_def').text if hit.find('Hit_def') is not None else ''

            # Parse SwissProt format: sp|ACCESSION|GENE_ORGANISM Description
            swissprot_id = ''
            swissprot_gene = ''
            swissprot_organism = ''
            description = hit_def

            if '|' in hit_id:
                parts = hit_id.split('|')
                if len(parts) >= 3:
                    swissprot_id = parts[1]
                    gene_org = parts[2]
                    if '_' in gene_org:
                        swissprot_gene = gene_org.split('_')[0]
                        swissprot_organism = gene_org.split('_', 1)[1] if '_' in gene_org else ''

            # Parse description for OS= organism
            if 'OS=' in hit_def:
                desc_parts = hit_def.split(' OS=')
                description = desc_parts[0]
                if len(desc_parts) > 1:
                    org_parts = desc_parts[1].split(' ')
                    swissprot_organism = ' '.join(org_parts[:2])  # Usually "Genus species"

            # Get HSP statistics
            hsp = hit.find('.//Hsp')
            if hsp is not None:
                identity = float(hsp.find('Hsp_identity').text)
                align_len = int(hsp.find('Hsp_align-len').text)
                evalue = float(hsp.find('Hsp_evalue').text)
                bitscore = float(hsp.find('Hsp_bit-score').text)
                percent_identity = (identity / align_len) * 100 if align_len > 0 else 0
            else:
                percent_identity = 0.0
                align_len = 0
                evalue = 1.0
                bitscore = 0.0

            annotations[query_id] = {
                'target_gene_id': query_id,
                'swissprot_id': swissprot_id,
                'swissprot_gene': swissprot_gene,
                'swissprot_description': description,
                'swissprot_organism': swissprot_organism,
                'percent_identity': round(percent_identity, 1),
                'alignment_length': align_len,
                'evalue': evalue,
                'bitscore': round(bitscore, 1),
            }

    except ET.ParseError as e:
        print(f"  [ERROR] Failed to parse XML: {e}")

    return annotations


def main():
    parser = argparse.ArgumentParser(
        description="Phase 6-7 (Helixer): Extract and annotate target proteins"
    )
    parser.add_argument(
        "--family", required=True,
        help="Gene family name"
    )
    parser.add_argument(
        "--phase4-dir", type=Path, required=True,
        help="Phase 4 Helixer output directory"
    )
    parser.add_argument(
        "--phase5-dir", type=Path, required=True,
        help="Phase 5 Helixer output directory"
    )
    parser.add_argument(
        "--helixer-dir", type=Path, required=True,
        help="Helixer proteomes directory (ragtag_output)"
    )
    parser.add_argument(
        "--swissprot-db", type=Path, required=True,
        help="SwissProt DIAMOND database (.dmnd)"
    )
    parser.add_argument(
        "--output-dir", type=Path, required=True,
        help="Output directory"
    )
    parser.add_argument(
        "--evalue", type=float, default=1e-5,
        help="E-value threshold for SwissProt search (default: 1e-5)"
    )
    parser.add_argument(
        "--threads", type=int, default=16,
        help="Number of threads for DIAMOND (default: 16)"
    )

    args = parser.parse_args()

    print("=" * 80)
    print("PHASE 6-7 (HELIXER): EXTRACT AND ANNOTATE TARGET PROTEINS")
    print("=" * 80)
    print()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Load Phase 4 targets
    phase4_targets = args.phase4_dir / "all_target_loci.tsv"
    if not phase4_targets.exists():
        print(f"[ERROR] Phase 4 targets not found: {phase4_targets}")
        return 1

    targets_df = pd.read_csv(phase4_targets, sep='\t')
    print(f"Loaded {len(targets_df)} targets from Phase 4")

    # Load Phase 5 classifications if available
    syntenic_file = args.phase5_dir / "syntenic_targets.tsv"
    unplaceable_file = args.phase5_dir / "unplaceable_targets.tsv"

    classification_map = {}
    locus_map = {}

    if syntenic_file.exists() and syntenic_file.stat().st_size > 0:
        try:
            syntenic_df = pd.read_csv(syntenic_file, sep='\t')
            for _, row in syntenic_df.iterrows():
                classification_map[row['target_gene_id']] = 'syntenic'
                if 'locus_id' in row:
                    locus_map[row['target_gene_id']] = row['locus_id']
            print(f"  {len(syntenic_df)} syntenic targets")
        except pd.errors.EmptyDataError:
            print("  0 syntenic targets (empty file)")

    if unplaceable_file.exists() and unplaceable_file.stat().st_size > 0:
        try:
            unplaceable_df = pd.read_csv(unplaceable_file, sep='\t')
            for _, row in unplaceable_df.iterrows():
                if row['target_gene_id'] not in classification_map:
                    classification_map[row['target_gene_id']] = 'unplaceable'
            print(f"  {len(unplaceable_df)} unplaceable targets")
        except pd.errors.EmptyDataError:
            print("  0 unplaceable targets (empty file)")

    # Extract target sequences
    print("\nExtracting target protein sequences...")
    target_fasta = args.output_dir / "target_proteins.faa"
    sequences = extract_target_sequences(targets_df, args.helixer_dir, target_fasta)

    if not sequences:
        print("[ERROR] No sequences extracted")
        return 1

    # Run SwissProt annotation
    swissprot_xml = args.output_dir / "swissprot_results.xml"
    success = run_diamond_swissprot(
        target_fasta, swissprot_xml, args.swissprot_db,
        evalue=args.evalue, threads=args.threads
    )

    if not success:
        print("[ERROR] SwissProt search failed")
        return 1

    # Parse annotations
    print("\nParsing SwissProt annotations...")
    annotations = parse_swissprot_xml(swissprot_xml)
    print(f"  Parsed {len(annotations)} annotations")

    # Count hits
    n_with_hits = sum(1 for a in annotations.values() if a['swissprot_id'])
    print(f"  {n_with_hits} targets have SwissProt matches")

    # Save annotations
    annotations_df = pd.DataFrame(list(annotations.values()))
    annotations_file = args.output_dir / "swissprot_annotations.tsv"
    annotations_df.to_csv(annotations_file, sep='\t', index=False)
    print(f"\n[OUTPUT] SwissProt annotations: {annotations_file}")

    # Merge with target data
    merged_df = targets_df.copy()

    # Add classification
    merged_df['classification'] = merged_df['target_gene_id'].map(classification_map).fillna('unclassified')
    merged_df['assigned_locus'] = merged_df['target_gene_id'].map(locus_map).fillna('')

    # Merge SwissProt annotations
    if not annotations_df.empty:
        annotations_subset = annotations_df[[
            'target_gene_id', 'swissprot_id', 'swissprot_gene',
            'swissprot_description', 'percent_identity', 'evalue'
        ]]
        merged_df = merged_df.merge(annotations_subset, on='target_gene_id', how='left')

    # Save merged output
    merged_file = args.output_dir / "annotated_targets.tsv"
    merged_df.to_csv(merged_file, sep='\t', index=False)
    print(f"[OUTPUT] Annotated targets: {merged_file}")

    # Summary by genome
    print("\nSummary by genome:")
    for genome in targets_df['genome'].unique():
        genome_df = merged_df[merged_df['genome'] == genome]
        n_total = len(genome_df)
        n_syntenic = len(genome_df[genome_df['classification'] == 'syntenic'])
        n_annotated = len(genome_df[genome_df['swissprot_id'].notna() & (genome_df['swissprot_id'] != '')])
        print(f"  {genome}: {n_total} targets, {n_syntenic} syntenic, {n_annotated} with SwissProt")

    print("\nPhase 6-7 (Helixer) complete.")
    return 0


if __name__ == "__main__":
    exit(main())
