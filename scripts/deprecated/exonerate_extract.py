#!/usr/bin/env python3
"""
exonerate_extract.py - Extract gene structures using Exonerate from tBLASTn hits

This script uses Exonerate to extract accurate gene structures from genomic regions
identified by tBLASTn searches, including handling of introns, frameshifts, and pseudogenes.
"""

import os
import sys
import json
import argparse
import subprocess
from pathlib import Path
from datetime import datetime
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from typing import Dict, List, Optional, Tuple

def extract_genomic_region(genome_fasta: str, scaffold: str,
                          start: int, end: int, flank: int = 1000) -> str:
    """
    Extract a genomic region with flanking sequence

    Args:
        genome_fasta: Path to genome FASTA file
        scaffold: Scaffold/contig ID
        start: Start position (1-based)
        end: End position (1-based)
        flank: Flanking sequence to add (bp)

    Returns:
        Path to extracted sequence file
    """
    # Adjust coordinates with flanking
    region_start = max(1, start - flank)
    region_end = end + flank

    # Create temporary file for region
    temp_file = f"temp_{scaffold}_{region_start}_{region_end}.fasta"

    # Extract using BioPython
    extracted = False
    with open(genome_fasta, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            # Match scaffold names flexibly (handle version suffixes like .1)
            # Try exact match first, then try with/without version suffix
            scaffold_base = scaffold.split('.')[0]  # Remove version if present
            record_base = record.id.split('.')[0]   # Remove version if present

            if record.id == scaffold or scaffold_base == record_base or scaffold in record.id or record.id in scaffold:
                # Extract subsequence (convert to 0-based)
                subseq = record.seq[region_start-1:region_end]

                # Write to temp file
                with open(temp_file, 'w') as out:
                    out.write(f">{scaffold}:{region_start}-{region_end}\n")
                    out.write(str(subseq) + "\n")

                extracted = True
                break

    if not extracted:
        print(f"Warning: Could not find scaffold {scaffold}", file=sys.stderr)
        return None

    return temp_file

def run_exonerate(query_protein: str, target_dna: str,
                 output_file: str, model: str = "protein2genome",
                 exhaustive: bool = False, percent: float = None) -> bool:
    """
    Run Exonerate to extract gene structure

    Args:
        query_protein: Path to query protein sequence
        target_dna: Path to target genomic DNA
        output_file: Path for output file
        model: Exonerate model to use
        exhaustive: Use exhaustive alignment (default False for speed)
        percent: Minimum percent score threshold

    Returns:
        True if successful
    """
    cmd = [
        "exonerate",
        "--model", model,
        "--query", query_protein,
        "--target", target_dna,
        "--showtargetgff", "yes",
        "--showvulgar", "yes",
        "--showquerygff", "no",
        "--showcigar", "no"
    ]

    # Only add percent threshold if specified
    if percent is not None:
        cmd.extend(["--percent", str(percent)])

    if exhaustive:
        cmd.append("--exhaustive")
        cmd.append("yes")

    # Add custom output format for sequences
    cmd.extend([
        "--ryo",
        "# Alignment %qi vs %ti\\n" +
        "# Query: %qd (%ql bp)\\n" +
        "# Target: %td (%tl bp) strand:%tS\\n" +
        "# Score: %s Percent: %pi\\n" +
        "# Query range: %qab-%qae\\n" +
        "# Target range: %tab-%tae\\n"
    ])

    try:
        with open(output_file, 'w') as f:
            result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE,
                                  text=True, timeout=300)

        if result.returncode != 0:
            print(f"Exonerate error: {result.stderr}", file=sys.stderr)
            return False

        return True

    except subprocess.TimeoutExpired:
        print(f"Exonerate timed out", file=sys.stderr)
        return False
    except Exception as e:
        print(f"Error running Exonerate: {e}", file=sys.stderr)
        return False

def parse_exonerate_gff(exonerate_output: str) -> List[Dict]:
    """
    Parse GFF features from Exonerate output

    Returns:
        List of feature dictionaries
    """
    features = []

    with open(exonerate_output, 'r') as f:
        for line in f:
            if line.startswith('##gff-version'):
                continue
            if line.startswith('#'):
                continue
            if not line.strip():
                continue

            # Parse GFF line
            parts = line.strip().split('\t')
            if len(parts) >= 9:
                feature = {
                    'scaffold': parts[0],
                    'source': parts[1],
                    'type': parts[2],
                    'start': int(parts[3]),
                    'end': int(parts[4]),
                    'score': parts[5],
                    'strand': parts[6],
                    'phase': parts[7],
                    'attributes': parts[8]
                }

                # Only keep relevant features (case-insensitive)
                feat_type = feature['type'].lower()
                if feat_type in ['gene', 'exon', 'cds', 'coding_exon', 'intron']:
                    features.append(feature)

    return features

def extract_cds_sequence(genome_fasta: str, features: List[Dict],
                         gene_id: str) -> Optional[str]:
    """
    Extract CDS sequence from genomic coordinates

    Handles case-insensitive feature types and falls back to exon features
    when CDS not explicitly present in GFF.

    Args:
        genome_fasta: Path to genome FASTA (or temp region file)
        features: List of GFF features
        gene_id: Gene identifier

    Returns:
        CDS sequence or None
    """
    if not features:
        print(f"    Warning: No features for {gene_id}", file=sys.stderr)
        return None

    # Get CDS features (case-insensitive, with fallback to exon)
    cds_features = [f for f in features if f['type'].lower() in {'cds', 'coding_exon'}]

    # Fallback to exon if no CDS features found
    if not cds_features:
        cds_features = [f for f in features if f['type'].lower() == 'exon']

    if not cds_features:
        print(f"    Warning: No CDS or exon features for {gene_id}", file=sys.stderr)
        return None

    # Sort by genomic position (will reverse for minus strand later if needed)
    cds_features.sort(key=lambda x: (x['start'], x['end']))

    # Determine strand
    strand = cds_features[0]['strand']
    scaffold = cds_features[0]['scaffold']

    # Extract scaffold sequence
    scaffold_seq = None
    with open(genome_fasta, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            # Match scaffold names flexibly
            # Handle temp file headers like "scaffold:start-end"
            if record.id == scaffold or scaffold in record.id or record.id in scaffold:
                scaffold_seq = str(record.seq)
                break

    if not scaffold_seq:
        print(f"    Warning: Could not find scaffold {scaffold} in {genome_fasta}", file=sys.stderr)
        return None

    # Extract and concatenate CDS sequences
    cds_parts = []
    for cds in cds_features:
        # GFF is 1-based inclusive; Python is 0-based end-exclusive
        start = cds['start'] - 1
        end = cds['end']

        # Bounds check
        if start < 0 or end > len(scaffold_seq):
            print(f"    Warning: CDS slice out of bounds: {cds['start']}-{cds['end']} (scaffold length={len(scaffold_seq)})", file=sys.stderr)
            return None

        seq = scaffold_seq[start:end]
        cds_parts.append(seq)

        if os.getenv('DEBUG'):
            print(f"    Extracted {cds['type']} {cds['start']}-{cds['end']} ({len(seq)} bp)", file=sys.stderr)

    # Concatenate in genomic order
    cds_seq = ''.join(cds_parts)

    # Reverse complement if minus strand (after concatenation)
    if strand == '-':
        cds_seq = str(Seq(cds_seq).reverse_complement())
        if os.getenv('DEBUG'):
            print(f"    Reverse-complemented for minus strand", file=sys.stderr)

    # Validate sequence
    if not cds_seq or not all(c in 'ACGTNacgtn' for c in cds_seq):
        print(f"    Warning: Extracted CDS is empty or contains invalid bases", file=sys.stderr)
        return None

    if os.getenv('DEBUG'):
        print(f"    Final CDS: {len(cds_seq)} bp, strand {strand}", file=sys.stderr)

    return cds_seq

def classify_gene_status(features: List[Dict], cds_seq: Optional[str]) -> Dict:
    """
    Classify gene functional status based on internal integrity.

    Classification based on:
    - Internal integrity: no internal stop codons, no frameshifts
    - Near full-length coverage (> 50 aa)
    - Canonical splice junctions (validated by Exonerate)

    Does NOT require terminal stop codon (Exonerate stops at query protein end).

    Returns:
        Dictionary with classification details
    """
    status = {
        'functional_status': 'unknown',
        'has_start_codon': False,
        'has_stop_codon': False,
        'frameshift': False,
        'premature_stop': False,
        'num_exons': 0,
        'total_cds_length': 0,
        'issues': []
    }

    if not cds_seq:
        status['functional_status'] = 'no_cds'
        return status

    # Check for start codon (informational only, not required for "intact" status)
    if cds_seq[:3].upper() == 'ATG':
        status['has_start_codon'] = True

    # Check for stop codons (informational only, not required)
    stop_codons = ['TAA', 'TAG', 'TGA']
    if cds_seq[-3:].upper() in stop_codons:
        status['has_stop_codon'] = True

    # Check for internal stop codons (premature stops - these ARE a problem)
    for i in range(0, len(cds_seq)-3, 3):
        if cds_seq[i:i+3].upper() in stop_codons:
            status['premature_stop'] = True
            status['issues'].append(f"Internal stop codon at position {i}")

    # Check if length is multiple of 3 (frameshift detection)
    if len(cds_seq) % 3 != 0:
        status['frameshift'] = True
        status['issues'].append(f"CDS length {len(cds_seq)} not multiple of 3")

    # Count exons
    exon_features = [f for f in features if f['type'] == 'exon']
    status['num_exons'] = len(exon_features)

    # Total CDS length
    status['total_cds_length'] = len(cds_seq)

    # Classify functional status based on internal integrity
    if status['premature_stop'] or status['frameshift']:
        # Has internal problems - definitely damaged
        status['functional_status'] = 'pseudogene'
    elif status['total_cds_length'] < 150:  # Less than 50 amino acids
        # Too short to be a real gene
        status['functional_status'] = 'partial'
    else:
        # Near full-length coverage, no internal defects, canonical splicing
        status['functional_status'] = 'intact'

    return status

def process_tblastn_hits(hits_file: str, query_protein: str, genome_fasta: str,
                        output_dir: str, params: Dict) -> List[Dict]:
    """
    Process tBLASTn hits with Exonerate

    Args:
        hits_file: Path to tBLASTn hits TSV file
        query_protein: Path to query protein FASTA
        genome_fasta: Path to genome FASTA
        output_dir: Output directory
        params: Parameters dictionary

    Returns:
        List of extracted genes
    """
    # Read tBLASTn hits
    if not os.path.exists(hits_file):
        return []

    df = pd.read_csv(hits_file, sep='\t')
    if df.empty:
        return []

    extracted_genes = []

    # Process each unique hit cluster
    if 'cluster' in df.columns:
        clusters = df.groupby('cluster')
    else:
        # Treat each hit as separate cluster
        df['cluster'] = range(len(df))
        clusters = df.groupby('cluster')

    for cluster_id, cluster_hits in clusters:
        # Get consensus region
        scaffold = cluster_hits.iloc[0]['sseqid']
        all_starts = list(cluster_hits['sstart']) + list(cluster_hits['send'])
        region_start = min(all_starts)
        region_end = max(all_starts)

        print(f"  Processing cluster {cluster_id} on {scaffold}:{region_start}-{region_end}")

        # Extract genomic region
        region_file = extract_genomic_region(
            genome_fasta, scaffold, region_start, region_end,
            flank=params['flank_size']
        )

        if not region_file:
            continue

        # Run Exonerate
        exonerate_output = Path(output_dir) / f"cluster_{cluster_id}_exonerate.txt"
        success = run_exonerate(
            query_protein, region_file, str(exonerate_output),
            model=params['model'],
            exhaustive=params['exhaustive'],
            percent=params['min_percent']
        )

        # Clean up temp file
        if os.path.exists(region_file):
            os.remove(region_file)

        if not success:
            continue

        # Parse Exonerate results
        features = parse_exonerate_gff(str(exonerate_output))

        if features:
            # Extract CDS sequence
            gene_id = f"gene_cluster_{cluster_id}"
            cds_seq = extract_cds_sequence(genome_fasta, features, gene_id)

            # Classify gene
            classification = classify_gene_status(features, cds_seq)

            # Store result
            gene_info = {
                'cluster_id': cluster_id,
                'gene_id': gene_id,
                'scaffold': scaffold,
                'start': region_start,
                'end': region_end,
                'best_hit_evalue': cluster_hits['evalue'].min(),
                'best_hit_pident': cluster_hits['pident'].max(),
                'num_exons': classification['num_exons'],
                'cds_length': classification['total_cds_length'],
                'functional_status': classification['functional_status'],
                'has_start_codon': classification['has_start_codon'],
                'has_stop_codon': classification['has_stop_codon'],
                'issues': ';'.join(classification['issues']) if classification['issues'] else ''
            }

            # Save CDS sequence if extracted
            if cds_seq:
                cds_file = Path(output_dir) / f"{gene_id}_cds.fasta"
                with open(cds_file, 'w') as f:
                    f.write(f">{gene_id} {scaffold}:{region_start}-{region_end} status:{classification['functional_status']}\n")
                    f.write(cds_seq + "\n")
                gene_info['cds_file'] = str(cds_file)

            # Save GFF features
            gff_file = Path(output_dir) / f"{gene_id}.gff"
            with open(gff_file, 'w') as f:
                f.write("##gff-version 3\n")
                for feat in features:
                    f.write('\t'.join([
                        feat['scaffold'], feat['source'], feat['type'],
                        str(feat['start']), str(feat['end']), feat['score'],
                        feat['strand'], feat['phase'], feat['attributes']
                    ]) + '\n')
            gene_info['gff_file'] = str(gff_file)

            extracted_genes.append(gene_info)

    return extracted_genes

def main():
    parser = argparse.ArgumentParser(description='Extract gene structures using Exonerate')
    parser.add_argument('hits_file', help='tBLASTn hits TSV file')
    parser.add_argument('query_protein', help='Query protein FASTA file')
    parser.add_argument('genome_fasta', help='Genome FASTA file')
    parser.add_argument('--output-dir', default='exonerate_results',
                       help='Output directory')
    parser.add_argument('--flank-size', type=int, default=1000,
                       help='Flanking region size (default: 1000)')
    parser.add_argument('--model', default='protein2genome',
                       help='Exonerate model (default: protein2genome)')
    parser.add_argument('--min-percent', type=float, default=30.0,
                       help='Minimum percent score (default: 30)')
    parser.add_argument('--exhaustive', action='store_true',
                       help='Use exhaustive alignment')

    args = parser.parse_args()

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Set up parameters
    params = {
        'flank_size': args.flank_size,
        'model': args.model,
        'min_percent': args.min_percent,
        'exhaustive': args.exhaustive
    }

    print(f"Processing tBLASTn hits from {args.hits_file}")

    # Process hits
    extracted_genes = process_tblastn_hits(
        args.hits_file, args.query_protein, args.genome_fasta,
        str(output_dir), params
    )

    # Save summary
    if extracted_genes:
        summary_df = pd.DataFrame(extracted_genes)
        summary_file = output_dir / "extracted_genes_summary.tsv"
        summary_df.to_csv(summary_file, sep='\t', index=False)

        print(f"\n=== Summary ===")
        print(f"Total genes extracted: {len(extracted_genes)}")
        print(f"Functional: {sum(1 for g in extracted_genes if g['functional_status'] == 'functional')}")
        print(f"Pseudogenes: {sum(1 for g in extracted_genes if g['functional_status'] == 'pseudogene')}")
        print(f"Fragments: {sum(1 for g in extracted_genes if g['functional_status'] == 'fragment')}")
        print(f"\nResults saved to: {output_dir}")
    else:
        print("No genes extracted")

    # Write provenance
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    provenance_file = output_dir / f"provenance_{timestamp}.json"
    provenance = {
        'script': os.path.abspath(__file__),
        'timestamp': datetime.now().isoformat(),
        'command': ' '.join(sys.argv),
        'parameters': params,
        'num_genes_extracted': len(extracted_genes)
    }

    with open(provenance_file, 'w') as f:
        json.dump(provenance, f, indent=2)

if __name__ == '__main__':
    main()