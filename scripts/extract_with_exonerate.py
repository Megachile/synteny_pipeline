#!/usr/bin/env python3
"""
Wrapper to extract gene structures from synteny blocks using Exonerate.

Replaces BLAST HSP extraction with proper gene structure extraction.
"""

from pathlib import Path
import tempfile
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import exonerate_extract


def extract_block_genes(block, query_protein_file, genome_fasta, output_dir):
    """
    Extract gene structures from a synteny block using Exonerate.

    Uses adaptive windowing: starts with small flanking region and expands
    progressively until a complete gene is found or max window reached.

    Args:
        block: Dict with keys: scaffold, start, end, strand, hits
        query_protein_file: Path to query flanking proteins
        genome_fasta: Path to genome FASTA file
        output_dir: Output directory for this block

    Returns:
        List of extracted gene dictionaries with sequences
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)

    # Adaptive windowing strategy - try progressively larger flanking regions
    # Start small (fast for single-exon genes), expand until complete or max reached
    flank_sizes = [0, 1000, 5000, 10000, 15000, 20000, 50000, 100000, 200000]  # 0 to 200kb (400kb total window)

    # Get query protein length for coverage calculation
    query_length = 0
    try:
        for record in SeqIO.parse(query_protein_file, 'fasta'):
            query_length = len(record.seq)
            break
    except:
        query_length = 174  # Default to BK ferritin length if can't read

    best_features = []
    best_region_file = None
    best_exonerate_output = None
    best_flank = 0
    best_coverage = 0

    for flank in flank_sizes:
        # Extract genomic region with current flanking size
        region_file = exonerate_extract.extract_genomic_region(
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
        success = exonerate_extract.run_exonerate(
            query_protein=str(query_protein_file),
            target_dna=region_file,
            output_file=str(exonerate_output),
            model="protein2genome"
        )

        if not success:
            Path(region_file).unlink(missing_ok=True)
            continue

        # Parse Exonerate results
        features = exonerate_extract.parse_exonerate_gff(str(exonerate_output))

        if not features:
            Path(region_file).unlink(missing_ok=True)
            continue

        # Calculate coverage from CDS features
        cds_features = [f for f in features if f['type'].lower() in ['cds', 'coding_exon', 'exon']]
        total_cds_length = sum(f['end'] - f['start'] + 1 for f in cds_features) if cds_features else 0
        coverage = total_cds_length / (query_length * 3) if query_length > 0 else 0  # Divide by 3 for nucleotides

        # Check for completeness
        is_complete = coverage >= 0.90  # 90% coverage threshold

        # Keep track of best result so far
        if coverage > best_coverage:
            # Clean up previous best
            if best_region_file and best_region_file != region_file:
                Path(best_region_file).unlink(missing_ok=True)

            best_features = features
            best_region_file = region_file
            best_exonerate_output = exonerate_output
            best_flank = flank
            best_coverage = coverage
        else:
            # This result was worse, clean it up
            Path(region_file).unlink(missing_ok=True)

        # If complete, stop expanding
        if is_complete:
            print(f"      Complete gene at {flank}bp flanking (coverage: {coverage:.1%})", flush=True)
            break

    # Use best result found
    if not best_features or best_coverage == 0:
        if best_region_file:
            Path(best_region_file).unlink(missing_ok=True)
        return []

    # Report result
    if best_coverage < 0.90:
        print(f"      Partial gene at {best_flank}bp flanking (coverage: {best_coverage:.1%})", flush=True)

    region_file = best_region_file
    features = best_features

    if not features:
        # Clean up temp file if no features
        Path(region_file).unlink(missing_ok=True)
        return []

    # Group features by gene
    # First, find all "gene" features with gene_id in attributes
    gene_features_list = [f for f in features if f['type'].lower() == 'gene']

    # For each gene, collect all associated CDS/exon features
    genes_by_id = {}

    if gene_features_list:
        # We have explicit gene features - group by gene_id
        for gene_feat in gene_features_list:
            # Parse gene ID from attributes (format: "gene_id <id> ; ...")
            if 'gene_id' in gene_feat['attributes']:
                gene_id = gene_feat['attributes'].split('gene_id')[1].split(';')[0].strip()

                # Collect all CDS/exon features within this gene's span
                gene_start = gene_feat['start']
                gene_end = gene_feat['end']
                gene_strand = gene_feat['strand']

                # Initialize gene entry with the gene feature
                genes_by_id[gene_id] = {
                    'gene_feature': gene_feat,
                    'cds_features': [],
                    'exon_features': [],
                    'start': gene_start,
                    'end': gene_end,
                    'strand': gene_strand
                }

                # Collect CDS and exon features within this gene
                for feat in features:
                    if feat['type'].lower() in ['cds', 'coding_exon']:
                        # Check if this CDS is within the gene span
                        if feat['start'] >= gene_start and feat['end'] <= gene_end:
                            genes_by_id[gene_id]['cds_features'].append(feat)
                    elif feat['type'].lower() == 'exon':
                        if feat['start'] >= gene_start and feat['end'] <= gene_end:
                            genes_by_id[gene_id]['exon_features'].append(feat)
    else:
        # No explicit gene features - treat all features as one gene
        # This handles cases where Exonerate only outputs CDS/exon without gene wrapper
        gene_id = '1'
        cds_features = [f for f in features if f['type'].lower() in ['cds', 'coding_exon', 'exon']]

        if cds_features:
            genes_by_id[gene_id] = {
                'gene_feature': None,
                'cds_features': cds_features,
                'exon_features': [],
                'start': min(f['start'] for f in cds_features),
                'end': max(f['end'] for f in cds_features),
                'strand': cds_features[0]['strand']
            }

    # Extract CDS and classify each gene
    extracted_genes = []

    for gene_id, gene_data in genes_by_id.items():
        # Use CDS features for extraction
        cds_features = gene_data['cds_features']
        if not cds_features:
            cds_features = gene_data['exon_features']  # Fallback to exons

        # Extract CDS sequence from temp region file (GFF coords are relative to it)
        cds_seq = exonerate_extract.extract_cds_sequence(
            genome_fasta=region_file,  # Use temp region file, not full genome
            features=cds_features,
            gene_id=gene_id
        )

        # Classify functional status
        classification = exonerate_extract.classify_gene_status(cds_features, cds_seq)

        # Extract full genomic sequence (gene with introns)
        genomic_seq = None
        try:
            with open(region_file, 'r') as f:
                for record in SeqIO.parse(f, 'fasta'):
                    # Extract gene span from temp region
                    gene_start = gene_data['start'] - 1  # GFF 1-based to Python 0-based
                    gene_end = gene_data['end']
                    genomic_seq = str(record.seq[gene_start:gene_end])

                    # Reverse complement if minus strand
                    if gene_data['strand'] == '-':
                        genomic_seq = str(Seq(genomic_seq).reverse_complement())
                    break
        except Exception as e:
            print(f"    Warning: Could not extract genomic sequence for {gene_id}: {e}", flush=True)

        # Save sequences
        if cds_seq:
            # Save CDS sequence
            cds_file = output_dir / f"{block['block_id']}_gene{gene_id}_cds.fasta"
            with open(cds_file, 'w') as f:
                header = f">gene{gene_id} {block['scaffold']}:{block['start']}-{block['end']} strand:{gene_data['strand']} status:{classification['functional_status']} exons:{len(cds_features)}"
                f.write(header + "\n")
                f.write(cds_seq + "\n")

            # Save genomic sequence (with introns)
            if genomic_seq:
                genomic_file = output_dir / f"{block['block_id']}_gene{gene_id}_genomic.fasta"
                with open(genomic_file, 'w') as f:
                    header = f">gene{gene_id} {block['scaffold']}:{block['start']}-{block['end']} strand:{gene_data['strand']} type:genomic_with_introns length:{len(genomic_seq)}bp"
                    f.write(header + "\n")
                    f.write(genomic_seq + "\n")

            # Translate to protein
            try:
                protein_seq = str(Seq(cds_seq).translate(to_stop=True))

                # Save protein sequence
                protein_file = output_dir / f"{block['block_id']}_gene{gene_id}_protein.fasta"
                with open(protein_file, 'w') as f:
                    header = f">gene{gene_id} {block['scaffold']}:{block['start']}-{block['end']} strand:{gene_data['strand']} length:{len(protein_seq)}aa"
                    f.write(header + "\n")
                    f.write(protein_seq + "\n")

                extracted_genes.append({
                    'gene_id': f"gene{gene_id}",
                    'scaffold': block['scaffold'],
                    'start': block['start'],
                    'end': block['end'],
                    'strand': gene_data['strand'],
                    'gene_span': gene_data['end'] - gene_data['start'] + 1,
                    'cds_length': len(cds_seq),
                    'protein_length': len(protein_seq),
                    'functional_status': classification['functional_status'],
                    'has_start_codon': classification['has_start_codon'],
                    'has_stop_codon': classification['has_stop_codon'],
                    'num_exons': len(cds_features),
                    'cds_seq': cds_seq,
                    'genomic_seq': genomic_seq,
                    'protein_seq': protein_seq
                })
            except Exception as e:
                print(f"    Warning: Could not translate gene{gene_id}: {e}", flush=True)

    # Clean up temp genomic region file after all genes extracted
    Path(region_file).unlink(missing_ok=True)

    return extracted_genes
