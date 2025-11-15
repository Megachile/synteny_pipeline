#!/usr/bin/env python3
"""
Step 06: Extract sequences for filtered targets using Exonerate.

Consolidated script - all Exonerate functions included directly.

Reads filtered target list (from Phase 5) and extracts gene structures.

Usage:
    python 06_extract_sequences.py \\
        --syntenic <path/to/syntenic_targets.tsv> \\
        --unplaceable <path/to/unplaceable_targets.tsv> \\
        --query-proteins <path/to/combined_targets.faa> \\
        --genome-fasta-dir <path/to/ragtag_output> \\
        --output-dir <path/to/extracted_sequences> \\
        [--unplaceable-evalue 1e-10] \\
        [--verbose]
"""

from pathlib import Path
import pandas as pd
import sys
import argparse
import shutil
import subprocess
import tempfile
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from typing import Dict, List, Optional
from dataclasses import dataclass

# =============================================================================
# EXONERATE FUNCTIONS (consolidated from exonerate_extract.py)
# =============================================================================

def extract_genomic_region(genome_fasta: str, scaffold: str,
                          start: int, end: int, flank: int = 1000) -> Optional[str]:
    """
    Extract a genomic region with flanking sequence.

    Args:
        genome_fasta: Path to genome FASTA file
        scaffold: Scaffold/contig ID
        start: Start position (1-based)
        end: End position (1-based)
        flank: Flanking sequence to add (bp)

    Returns:
        Path to extracted sequence file or None if failed
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
            scaffold_base = scaffold.split('.')[0]
            record_base = record.id.split('.')[0]

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
                 output_file: str, model: str = "protein2genome") -> bool:
    """
    Run Exonerate to extract gene structure.

    Args:
        query_protein: Path to query protein sequence
        target_dna: Path to target genomic DNA
        output_file: Path for output file
        model: Exonerate model to use

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

    # Add custom output format
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
    """Parse GFF features from Exonerate output."""
    features = []

    with open(exonerate_output, 'r') as f:
        for line in f:
            if line.startswith('##gff-version'):
                continue
            if line.startswith('#'):
                continue
            if not line.strip():
                continue

            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            feature = {
                'seqname': parts[0],
                'source': parts[1],
                'type': parts[2],
                'start': int(parts[3]),
                'end': int(parts[4]),
                'score': parts[5],
                'strand': parts[6],
                'frame': parts[7],
                'attributes': parts[8]
            }

            features.append(feature)

    return features

def extract_cds_sequence(genome_file: str, features: List[Dict]) -> Optional[str]:
    """Extract CDS sequence from genomic coordinates."""
    if not features:
        return None

    # Get CDS features
    cds_features = sorted(
        [f for f in features if f['type'].lower() in ['cds', 'coding_exon', 'exon']],
        key=lambda x: x['start']
    )

    if not cds_features:
        return None

    # Load genomic sequence
    seqname = cds_features[0]['seqname']
    genomic_seq = None

    with open(genome_file, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            if seqname in record.id or record.id in seqname:
                genomic_seq = str(record.seq)
                break

    if not genomic_seq:
        return None

    # Extract and concatenate CDS regions
    cds_seq = ""
    for feature in cds_features:
        start = feature['start'] - 1  # Convert to 0-based
        end = feature['end']
        exon_seq = genomic_seq[start:end]
        cds_seq += exon_seq

    # Reverse complement if on minus strand
    if cds_features[0]['strand'] == '-':
        cds_seq = str(Seq(cds_seq).reverse_complement())

    return cds_seq

def classify_gene_status(features: List[Dict], cds_seq: Optional[str], query_length_aa: int) -> Dict:
    """
    Classify gene functional status with query-length-aware classification.

    Classification based on:
    - Internal integrity: no internal stop codons, no frameshifts
    - Coverage relative to query protein (>= 90% of query length = intact)
    - Canonical splice junctions (validated by Exonerate)

    Args:
        features: Exonerate GFF features
        cds_seq: Extracted CDS nucleotide sequence
        query_length_aa: Length of query protein in amino acids

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
        'extracted_length_aa': 0,
        'query_length_aa': query_length_aa,
        'coverage_pct': 0.0,
        'issues': []
    }

    if not cds_seq:
        status['functional_status'] = 'no_cds'
        return status

    # Calculate extracted protein length
    extracted_length_aa = len(cds_seq) // 3
    status['extracted_length_aa'] = extracted_length_aa
    status['coverage_pct'] = (extracted_length_aa / query_length_aa * 100) if query_length_aa > 0 else 0

    # Check for start codon (informational only)
    if cds_seq[:3].upper() == 'ATG':
        status['has_start_codon'] = True

    # Check for stop codons (informational only)
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

    # Classify functional status with query-length awareness
    if status['premature_stop'] or status['frameshift']:
        # Has internal problems - definitely damaged
        status['functional_status'] = 'pseudogene'
    elif extracted_length_aa >= (query_length_aa * 0.9):
        # Near-full-length (>= 90% of query) with no internal defects
        # This is essentially intact even if terminal stop is missing
        status['functional_status'] = 'intact'
    else:
        # Too short - likely a fragment
        status['functional_status'] = 'fragment'

    return status

# =============================================================================
# EXTRACTION WITH ADAPTIVE WINDOWING (consolidated from extract_with_exonerate.py)
# =============================================================================

def extract_block_genes(block, query_protein_file, genome_fasta, output_dir):
    """
    Extract gene structures from a synteny block using Exonerate.

    Uses adaptive windowing: starts with small flanking region and expands
    progressively until a complete gene is found or max window reached.

    Args:
        block: Dict with keys: scaffold, start, end, strand, block_id
        query_protein_file: Path to query protein file
        genome_fasta: Path to genome FASTA file
        output_dir: Output directory for this block

    Returns:
        List of extracted gene dictionaries with sequences
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)

    # Adaptive windowing strategy
    flank_sizes = [0, 1000, 5000, 10000, 15000, 20000, 50000, 100000, 200000]

    # Get query protein length for classification
    query_length_aa = 0
    try:
        for record in SeqIO.parse(query_protein_file, 'fasta'):
            query_length_aa = len(record.seq)
            break
    except:
        query_length_aa = 169  # Default fallback

    best_features = []
    best_region_file = None
    best_exonerate_output = None
    best_flank = 0
    best_coverage = 0

    for flank in flank_sizes:
        # Extract genomic region with current flanking size
        region_file = extract_genomic_region(
            genome_fasta=str(genome_fasta),
            scaffold=block['scaffold'],
            start=block['start'],
            end=block['end'],
            flank=flank
        )

        if not region_file:
            continue

        # Run Exonerate
        exonerate_output = output_dir / f"{block['block_id']}_exonerate_flank{flank}.txt"
        success = run_exonerate(
            query_protein=str(query_protein_file),
            target_dna=region_file,
            output_file=str(exonerate_output),
            model="protein2genome"
        )

        if not success:
            Path(region_file).unlink(missing_ok=True)
            continue

        # Parse Exonerate results
        features = parse_exonerate_gff(str(exonerate_output))

        if not features:
            Path(region_file).unlink(missing_ok=True)
            continue

        # Calculate coverage
        cds_features = [f for f in features if f['type'].lower() in ['cds', 'coding_exon', 'exon']]
        total_cds_length = sum(f['end'] - f['start'] + 1 for f in cds_features) if cds_features else 0
        coverage = total_cds_length / (query_length_aa * 3) if query_length_aa > 0 else 0

        # Check for completeness
        is_complete = coverage >= 0.90  # 90% coverage threshold

        # Keep track of best result
        if coverage > best_coverage:
            if best_region_file and best_region_file != region_file:
                Path(best_region_file).unlink(missing_ok=True)

            best_features = features
            best_region_file = region_file
            best_exonerate_output = exonerate_output
            best_flank = flank
            best_coverage = coverage

        # If we found a complete gene, stop searching
        if is_complete:
            break

        # Clean up this region file if not the best
        if region_file != best_region_file:
            Path(region_file).unlink(missing_ok=True)

    # If we didn't find anything, return empty
    if not best_features:
        return []

    # Group features by gene (gene_id from GFF attributes)
    genes_by_id = defaultdict(lambda: {'features': [], 'cds_features': []})

    for feature in best_features:
        if feature['type'] == 'gene':
            # Extract gene_id from attributes
            attrs = feature['attributes']
            gene_id = None
            for attr in attrs.split(';'):
                attr = attr.strip()
                if attr.startswith('gene_id'):
                    gene_id = attr.split()[1]
                    break

            if gene_id:
                genes_by_id[gene_id]['gene_info'] = feature
        elif feature['type'] in ['cds', 'exon', 'coding_exon']:
            # Add to first gene (Exonerate typically outputs one gene per run)
            gene_id = '1'
            genes_by_id[gene_id]['cds_features'].append(feature)
            genes_by_id[gene_id]['features'].append(feature)

    # If no genes found by gene_id, create a single gene from all CDS features
    if not genes_by_id:
        cds_features = [f for f in best_features if f['type'] in ['cds', 'exon', 'coding_exon']]
        if cds_features:
            genes_by_id['1'] = {
                'cds_features': cds_features,
                'features': cds_features,
                'strand': cds_features[0]['strand']
            }

    # Extract CDS and classify each gene
    extracted_genes = []

    for gene_id, gene_data in genes_by_id.items():
        cds_features = gene_data['cds_features']

        if not cds_features:
            continue

        # Extract CDS sequence
        cds_seq = extract_cds_sequence(
            genome_file=best_region_file,
            features=cds_features
        )

        # Classify functional status with query-length awareness
        classification = classify_gene_status(cds_features, cds_seq, query_length_aa)

        # Extract full genomic sequence (gene with introns)
        genomic_seq = None
        if cds_features:
            gene_start = min(f['start'] for f in cds_features)
            gene_end = max(f['end'] for f in cds_features)

            with open(best_region_file, 'r') as f:
                for record in SeqIO.parse(f, 'fasta'):
                    genomic_seq = str(record.seq)[gene_start-1:gene_end]
                    if cds_features[0]['strand'] == '-':
                        genomic_seq = str(Seq(genomic_seq).reverse_complement())
                    break

        # Translate CDS to protein
        protein_seq = None
        if cds_seq and len(cds_seq) % 3 == 0:
            try:
                protein_seq = str(Seq(cds_seq).translate(to_stop=False))
            except:
                protein_seq = None

        # Save sequences
        if cds_seq:
            # Save CDS sequence
            cds_file = output_dir / f"{block['block_id']}_gene{gene_id}_cds.fasta"
            with open(cds_file, 'w') as f:
                header = f">gene{gene_id} {block['scaffold']}:{block['start']}-{block['end']} strand:{gene_data.get('strand', '+')} status:{classification['functional_status']} exons:{len(cds_features)}"
                f.write(header + "\n")
                f.write(cds_seq + "\n")

        if genomic_seq:
            # Save genomic sequence
            genomic_file = output_dir / f"{block['block_id']}_gene{gene_id}_genomic.fasta"
            with open(genomic_file, 'w') as f:
                header = f">gene{gene_id} {block['scaffold']}:{block['start']}-{block['end']} strand:{gene_data.get('strand', '+')} type:genomic_with_introns length:{len(genomic_seq)}bp"
                f.write(header + "\n")
                f.write(genomic_seq + "\n")

        if protein_seq:
            # Save protein sequence
            protein_file = output_dir / f"{block['block_id']}_gene{gene_id}_protein.fasta"
            with open(protein_file, 'w') as f:
                header = f">gene{gene_id} {block['scaffold']}:{block['start']}-{block['end']} strand:{gene_data.get('strand', '+')} length:{len(protein_seq)}aa"
                f.write(header + "\n")
                f.write(protein_seq + "\n")

        extracted_genes.append({
            'gene_id': gene_id,
            'scaffold': block['scaffold'],
            'start': block['start'],
            'end': block['end'],
            'strand': gene_data.get('strand', '+'),
            'functional_status': classification['functional_status'],
            'cds_seq': cds_seq,
            'protein_seq': protein_seq,
            'genomic_seq': genomic_seq,
            'num_exons': len(cds_features),
            'best_flank': best_flank
        })

    # Clean up temp files
    if best_region_file:
        Path(best_region_file).unlink(missing_ok=True)

    return extracted_genes

# =============================================================================
# DEDUPLICATION
# =============================================================================

def deduplicate_extracted_sequences(targets_df, output_dir):
    """
    Deduplicate extracted sequences within each (genome, parent_locus) group.

    When multiple BLAST seeds lead to identical extracted proteins, keep only one.
    """
    duplicates_removed = 0
    # Log file to record dedup events for downstream QC
    dedup_log = output_dir / 'dedup_events.tsv'
    if not dedup_log.exists():
        try:
            with open(dedup_log, 'w') as w:
                w.write('\t'.join([
                    'genome', 'parent_locus', 'kept_locus', 'removed_locus',
                    'kept_evalue', 'removed_evalue'
                ]) + '\n')
        except Exception:
            pass

    for (genome, parent_locus), group in targets_df.groupby(['genome', 'parent_locus']):
        if len(group) == 1:
            continue

        sequences = {}
        target_info = {}

        for _, target in group.iterrows():
            target_name = target.get('locus_name', target.get('locus_id', target.get('parent_locus', 'unknown')))
            extraction_dir = output_dir / genome / target_name

            protein_files = list(extraction_dir.glob("*_gene1_protein.fasta"))
            if not protein_files:
                continue

            protein_file = protein_files[0]

            try:
                for record in SeqIO.parse(protein_file, 'fasta'):
                    seq = str(record.seq)
                    sequences[target_name] = seq
                    target_info[target_name] = {
                        'evalue': target['best_evalue'],
                        'extraction_dir': extraction_dir
                    }
                    break
            except Exception as e:
                print(f"    WARNING: Could not read {protein_file}: {e}")
                continue

        if len(sequences) <= 1:
            continue

        target_names = list(sequences.keys())
        removed_targets = set()

        for i, target1 in enumerate(target_names):
            if target1 in removed_targets:
                continue

            seq1 = sequences[target1]

            for target2 in target_names[i+1:]:
                if target2 in removed_targets:
                    continue

                seq2 = sequences[target2]

                if seq1 == seq2:
                    # Keep the one with better e-value
                    if target_info[target1]['evalue'] <= target_info[target2]['evalue']:
                        to_remove = target2
                        to_keep = target1
                    else:
                        to_remove = target1
                        to_keep = target2

                    print(f"    Duplicate in {genome}/{parent_locus}:")
                    print(f"      Keeping:  {to_keep} (e-value: {target_info[to_keep]['evalue']:.2e})")
                    print(f"      Removing: {to_remove} (e-value: {target_info[to_remove]['evalue']:.2e})")

                    # Append to log for manual review
                    try:
                        with open(dedup_log, 'a') as w:
                            w.write('\t'.join([
                                str(genome),
                                str(parent_locus),
                                str(to_keep),
                                str(to_remove),
                                f"{target_info[to_keep]['evalue']:.3e}" if target_info[to_keep]['evalue'] != float('inf') else 'inf',
                                f"{target_info[to_remove]['evalue']:.3e}" if target_info[to_remove]['evalue'] != float('inf') else 'inf',
                            ]) + '\n')
                    except Exception:
                        pass

                    extraction_dir = target_info[to_remove]['extraction_dir']
                    if extraction_dir.exists():
                        shutil.rmtree(extraction_dir)
                        duplicates_removed += 1

                    removed_targets.add(to_remove)

    return duplicates_removed

# =============================================================================
# PRE-EXTRACTION DEDUPLICATION (avoid double Exonerate on identical sites)
# =============================================================================

def _interval_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    a1, a2 = (a_start, a_end) if a_start <= a_end else (a_end, a_start)
    b1, b2 = (b_start, b_end) if b_start <= b_end else (b_end, b_start)
    return max(0, min(a2, b2) - max(a1, b1) + 1)

def deduplicate_targets_by_overlap(targets_df):
    """
    Collapse overlapping targets per (genome, assigned_to, scaffold) keeping best e-value.

    This is conservative: only overlapping intervals are merged (no distance-based merging),
    so tandem duplicates separated by gaps remain intact.
    """
    if targets_df is None or len(targets_df) == 0:
        return targets_df, []

    df = targets_df.copy()
    if 'assigned_to' not in df.columns:
        # If classification did not provide assigned_to, do nothing
        return df, []

    kept_idx = []
    removed_rows = []

    # Group by placed locus
    for (genome, assigned, scaffold), grp in df.groupby(['genome', 'assigned_to', 'scaffold']):
        # Build sortable tuples
        items = []
        for ix, r in grp.iterrows():
            try:
                s = int(r['start']); e = int(r['end'])
                ev = float(r.get('best_evalue', 'inf'))
            except Exception:
                s = int(r.get('start', 0) or 0); e = int(r.get('end', 0) or 0)
                ev = float('inf')
            items.append((ix, min(s,e), max(s,e), ev))
        # Sort by start
        items.sort(key=lambda x: (x[1], x[2]))

        # Sweep clusters by overlap only
        clusters = []
        cur = []
        cur_s = cur_e = None
        for ix, s, e, ev in items:
            if not cur:
                cur = [(ix, s, e, ev)]; cur_s = s; cur_e = e
            else:
                if s <= cur_e:  # overlap
                    cur.append((ix, s, e, ev)); cur_e = max(cur_e, e)
                else:
                    clusters.append(cur); cur = [(ix, s, e, ev)]; cur_s = s; cur_e = e
        if cur:
            clusters.append(cur)

        # Keep best by lowest e-value within each overlapping cluster
        for cl in clusters:
            best = min(cl, key=lambda t: t[3])
            kept_idx.append(best[0])
            for t in cl:
                if t[0] != best[0]:
                    removed_rows.append(df.loc[t[0]].to_dict())

    dedup_df = df.loc[sorted(set(kept_idx))].copy()
    return dedup_df, removed_rows

# =============================================================================
# MAIN WORKFLOW
# =============================================================================

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Extract sequences for filtered targets using Exonerate",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('--syntenic', required=True, type=Path,
                        help='Path to syntenic_targets.tsv from Phase 5')
    parser.add_argument('--unplaceable', type=Path,
                        help='Path to unplaceable_targets.tsv from Phase 5 (optional)')
    parser.add_argument('--query-proteins', required=True, type=Path,
                        help='Path to combined_targets.faa with query proteins')
    parser.add_argument('--genome-fasta-dir', required=True, type=Path,
                        help='Directory containing genome FASTA files')
    parser.add_argument('--output-dir', required=True, type=Path,
                        help='Output directory for extracted sequences')
    parser.add_argument('--unplaceable-evalue', type=float, default=1e-10,
                        help='E-value threshold for unplaceable targets (default: 1e-10)')
    parser.add_argument('--verbose', action='store_true',
                        help='Print verbose extraction progress')

    return parser.parse_args()

def main():
    """Extract sequences for filtered targets."""
    args = parse_args()

    print("=" * 80, flush=True)
    print("STEP 06: EXTRACT SEQUENCES WITH EXONERATE (CONSOLIDATED)", flush=True)
    print("=" * 80, flush=True)
    print(f"\nInput files:")
    print(f"  Syntenic targets: {args.syntenic}")
    print(f"  Unplaceable targets: {args.unplaceable if args.unplaceable else 'None'}")
    print(f"  Query proteins: {args.query_proteins}")
    print(f"  Genome FASTA dir: {args.genome_fasta_dir}")
    print(f"  Output directory: {args.output_dir}")
    print(f"\nParameters:")
    print(f"  Unplaceable e-value threshold: {args.unplaceable_evalue}")
    print(f"  Classification: Query-length-aware (>= 90% = intact)")

    # Load filtered targets
    print("\n[1] Loading filtered targets...", flush=True)

    if not args.syntenic.exists():
        print(f"  ERROR: Syntenic targets not found: {args.syntenic}", flush=True)
        return

    syntenic_df = pd.read_csv(args.syntenic, sep='\t')
    print(f"  Loaded {len(syntenic_df)} syntenic targets", flush=True)

    if args.unplaceable and args.unplaceable.exists():
        unplaceable_df = pd.read_csv(args.unplaceable, sep='\t')
        print(f"  Loaded {len(unplaceable_df)} unplaceable targets", flush=True)

        strong_unplaceable = unplaceable_df[unplaceable_df['best_evalue'] < args.unplaceable_evalue]
        print(f"  Kept {len(strong_unplaceable)} strong unplaceable targets", flush=True)

        targets_df = pd.concat([syntenic_df, strong_unplaceable], ignore_index=True)
    else:
        targets_df = syntenic_df

    print(f"  Total targets for extraction: {len(targets_df)}", flush=True)

    # Pre-extraction dedup on syntenic targets only (do not touch unplaceables)
    try:
        syn_only = targets_df[targets_df.get('placement', '') == 'synteny'].copy()
        rest = targets_df[targets_df.get('placement', '') != 'synteny'].copy()
        syn_dedup, removed = deduplicate_targets_by_overlap(syn_only)
        if removed:
            # Log removed overlaps
            logp = args.output_dir / 'dedup_pre_phase6.tsv'
            with open(logp, 'w') as w:
                w.write('\t'.join(syn_only.columns) + '\n')
                for r in removed:
                    w.write('\t'.join(str(r.get(c, '')) for c in syn_only.columns) + '\n')
            print(f"  Pre-extraction dedup removed {len(removed)} overlapping synteny targets (logged to {logp})", flush=True)
        targets_df = pd.concat([syn_dedup, rest], ignore_index=True)
        print(f"  Targets after overlap-dedup: {len(targets_df)}", flush=True)
    except Exception as e:
        print(f"  WARNING: pre-extraction dedup failed: {e}", flush=True)

    # Group by genome
    print("\n[2] Grouping targets by genome...", flush=True)
    genome_groups = targets_df.groupby('genome')
    print(f"  {len(genome_groups)} genomes to process", flush=True)

    args.output_dir.mkdir(exist_ok=True, parents=True)

    # Process each genome
    print("\n[3] Extracting sequences with Exonerate...", flush=True)
    total_extracted = 0
    failed_genomes = []

    for genome_name, genome_targets in genome_groups:
        print(f"\n  {genome_name}: {len(genome_targets)} targets", flush=True)

        genome_fasta = args.genome_fasta_dir / genome_name / "ragtag.scaffold.fasta"
        if not genome_fasta.exists():
            print(f"    WARNING: Genome FASTA not found, skipping")
            failed_genomes.append(genome_name)
            continue

        for parent_locus, locus_targets in genome_targets.groupby('parent_locus'):
            for _, target in locus_targets.iterrows():
                expected_query = target.get('query_id', None)

                # Prefer assigned block coordinates for syntenic targets;
                # fall back to the target hit coordinates otherwise.
                assigned_scaf = target.get('assigned_block_scaffold') if 'assigned_block_scaffold' in target else None
                assigned_start = target.get('assigned_block_start') if 'assigned_block_start' in target else None
                assigned_end = target.get('assigned_block_end') if 'assigned_block_end' in target else None
                assigned_id = target.get('assigned_block_id') if 'assigned_block_id' in target else None
                assigned_to = target.get('assigned_to', None)

                use_assigned = (
                    str(target.get('placement', '')) == 'synteny' and
                    assigned_scaf not in (None, '', 'nan') and
                    assigned_start not in (None, '', 'nan') and
                    assigned_end not in (None, '', 'nan')
                )

                if use_assigned:
                    # Syntenic targets: use assigned_to as directory name (NOT parent_locus)
                    target_name = assigned_to if assigned_to not in (None, '', 'nan') else target.get('locus_name', target.get('locus_id', target.get('parent_locus', 'unknown')))
                    # Always extract using the target interval, but keep block_id from the assigned block
                    block = {
                        'block_id': assigned_id if assigned_id not in (None, '', 'nan') else target_name,
                        'scaffold': target['scaffold'],
                        'start': target['start'],
                        'end': target['end'],
                        'strand': target.get('strand', '+')
                    }
                    # Canonical syntenic outputs go under genome/assigned_to
                    target_output_dir = args.output_dir / genome_name / target_name
                else:
                    # Unplaceable targets: no assigned locus, store in flat UNPLACED directory
                    target_name = target.get('locus_name', target.get('locus_id', target.get('parent_locus', 'unknown')))
                    block = {
                        'block_id': target_name,
                        'scaffold': target['scaffold'],
                        'start': target['start'],
                        'end': target['end'],
                        'strand': target['strand']
                    }
                    # Unplaceables go in genome/UNPLACED/unique_tag (no parent_locus subdivision)
                    up_scaf = str(target.get('scaffold'))
                    up_start = str(target.get('start'))
                    up_end = str(target.get('end'))
                    unique_tag = f"{up_scaf}_{up_start}_{up_end}"
                    target_output_dir = args.output_dir / genome_name / 'UNPLACED' / unique_tag
                target_output_dir.mkdir(exist_ok=True, parents=True)

                # Extract specific query protein
                temp_query = target_output_dir / f"query_{expected_query}.faa"
                found_query = False

                # Strategy:
                # 1) If expected_query provided, try to pull from combined_targets.faa
                # 2) Fallback: use first locus-specific target from Phase 1 directory
                def _write_first_record(src_fa: Path) -> bool:
                    try:
                        for rec in SeqIO.parse(str(src_fa), 'fasta'):
                            with open(temp_query, 'w') as outfile:
                                SeqIO.write(rec, outfile, 'fasta')
                            return True
                    except Exception:
                        return False
                    return False

                if expected_query:
                    with open(args.query_proteins, 'r') as infile:
                        for record in SeqIO.parse(infile, 'fasta'):
                            if expected_query in record.id:
                                with open(temp_query, 'w') as outfile:
                                    SeqIO.write(record, outfile, 'fasta')
                                found_query = True
                                break

                if not found_query:
                    # Infer family base dir from combined_targets path: outputs/<family>/phase4*/combined_targets.faa
                    base_dir = args.query_proteins.parent.parent if args.query_proteins.parent.name.startswith('phase4') else args.query_proteins.parent
                    # Try phase1_v2 then phase1
                    locus = target.get('parent_locus', target.get('locus_id', ''))
                    if locus:
                        cand1 = base_dir / 'phase1_v2' / locus / f"{locus}_targets.faa"
                        cand2 = base_dir / 'phase1' / locus / f"{locus}_targets.faa"
                        for cand in (cand1, cand2):
                            if cand.exists() and _write_first_record(cand):
                                found_query = True
                                break

                if not found_query:
                    msg_q = expected_query if expected_query else '(none)'
                    print(f"    WARNING: Could not resolve query protein (query={msg_q}); skipping target {target_name}")
                    continue

                # Extract with Exonerate
                try:
                    genes = extract_block_genes(
                        block=block,
                        query_protein_file=temp_query,
                        genome_fasta=genome_fasta,
                        output_dir=target_output_dir
                    )

                    if genes:
                        total_extracted += len(genes)
                        if args.verbose:
                            status_summary = genes[0]['functional_status']
                            print(f"    {target_name}: {len(genes)} genes ({status_summary})")
                    else:
                        if args.verbose:
                            print(f"    {target_name}: no genes extracted")

                except Exception as e:
                    print(f"    ERROR extracting {target_name}: {e}")

    # Deduplicate
    print("\n[4] Deduplicating extracted sequences...", flush=True)
    duplicates_removed = deduplicate_extracted_sequences(targets_df, args.output_dir)
    print(f"  Removed {duplicates_removed} duplicate extractions", flush=True)

    # Aggregate all protein sequences into single file
    print("\n[5] Aggregating extracted proteins...", flush=True)
    all_proteins_file = args.output_dir / "all_extracted_genes.faa"
    protein_count = 0

    with open(all_proteins_file, 'w') as outfile:
        # Find all *_protein.fasta files in genome subdirectories
        for genome_dir in args.output_dir.iterdir():
            if genome_dir.is_dir():
                for target_dir in genome_dir.iterdir():
                    if target_dir.is_dir():
                        for protein_file in target_dir.glob("*_protein.fasta"):
                            # Read and write each protein sequence
                            for record in SeqIO.parse(protein_file, 'fasta'):
                                SeqIO.write(record, outfile, 'fasta')
                                protein_count += 1

    print(f"  Wrote {protein_count} protein sequences to {all_proteins_file.name}", flush=True)

    # Summary
    print("\n" + "=" * 80)
    print("SEQUENCE EXTRACTION COMPLETE")
    print("=" * 80)
    print(f"\nTotal genes extracted: {total_extracted}")
    print(f"Duplicates removed: {duplicates_removed}")
    print(f"Unique genes: {total_extracted - duplicates_removed}")
    print(f"Aggregated proteins: {protein_count}")
    print(f"Genomes processed: {len(genome_groups) - len(failed_genomes)}/{len(genome_groups)}")
    print(f"\nOutput directory: {args.output_dir}")
    print(f"Aggregated file: {all_proteins_file}")

if __name__ == "__main__":
    main()
