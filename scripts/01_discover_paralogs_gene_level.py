#!/usr/bin/env python3
"""
Phase 1: Gene-level paralog discovery with isoform deduplication.

Key improvements:
1. Maps proteins to LOC IDs to handle isoforms properly
2. Extracts flanking at GENE level (not protein level)
3. Detects tandem duplications among input LOCs
4. Detects tandem duplications among paralogs
5. Names loci systematically by chromosome

Usage:
    python 01_discover_paralogs_gene_level.py <LOC_ID(s)> <output_dir>

    LOC_ID(s): Single LOC or comma-separated LOCs
               Single: LOC117167432
               Multiple: LOC117167432,LOC117167433,LOC117167434

When multiple LOCs are provided, the script detects if any are tandem duplications
(adjacent on the same chromosome) and groups them into clusters. One representative
from each cluster is used for paralog discovery.
"""

import sys
import re
import subprocess
from pathlib import Path
from collections import defaultdict
import pandas as pd
from Bio import SeqIO
import json
from datetime import datetime

# Configuration
LANDMARKS = {
    'BK': {
        'proteome': 'data/proteomes/GCF_010883055.1_Bkinseyi_NCBI.faa',
        'gff': 'data/genomes/GCF_010883055.1_Bkinseyi.gff',
        'db': 'data/proteomes/GCF_010883055.1_Bkinseyi_NCBI.dmnd'
    },
    'LB': {
        'proteome': 'data/proteomes/GCF_019393585.1_Lboulardi_NCBI.faa',
        'gff': 'data/genomes/GCF_019393585.1_Lboulardi.gff',
        'db': 'data/proteomes/GCF_019393585.1_Lboulardi_NCBI.dmnd'
    },
    'TR': {
        'proteome': 'data/proteomes/GCA_020615435.1.faa',
        'gff': 'data/genomes/GCA_020615435.1.gff3',  # BRAKER3 GFF3
        'db': 'data/proteomes/GCA_020615435.1.dmnd'
    },
    'DR': {
        'proteome': 'data/proteomes/GCA_030998225.1.faa',
        'gff': 'data/genomes/GCA_030998225.1.gff3',  # BRAKER3 GFF3
        'db': 'data/proteomes/GCA_030998225.1.dmnd'
    }
}

FLANKING_DISTANCE_KB = 1000  # Genomic distance (kb) to extract flanking genes
EVALUE_THRESHOLD = 1e-3
MIN_PIDENT = 40.0  # Minimum percent identity to avoid distant paralogs
MIN_SYNTENY_MATCHES = 9  # Out of 80 flanking genes


class GenomicMapper:
    """Handles protein-to-gene mapping for different annotation types."""

    def __init__(self, genome_code, gff_file=None, proteome_file=None):
        self.genome_code = genome_code
        self.protein_to_gene = {}
        self.gene_to_proteins = defaultdict(set)
        self.gene_positions = {}
        self.chr_genes = defaultdict(list)
        self.gene_order = []  # For BRAKER3 - ordered list of genes

        if gff_file and Path(gff_file).exists():
            self._parse_gff(gff_file)
        else:
            # For BRAKER3 genomes, use protein ID pattern and order
            self._use_braker_pattern(proteome_file)

    def _parse_gff(self, gff_file):
        """Parse GFF to extract gene/protein relationships (handles NCBI and BRAKER3)."""
        print(f"  Parsing GFF for {self.genome_code}...")

        with open(gff_file) as f:
            for line in f:
                if line.startswith('#'):
                    continue

                parts = line.split('\t')
                if len(parts) < 9:
                    continue

                feature_type = parts[2]
                chrom = parts[0]
                start = int(parts[3])
                attributes = parts[8]

                # Extract gene positions (both NCBI and BRAKER3)
                if feature_type == 'gene':
                    # NCBI format: Name=LOC117167432
                    ncbi_match = re.search(r'Name=(LOC\d+)', attributes)
                    # BRAKER3 format: ID=g10440;
                    braker_match = re.search(r'ID=(g\d+);', attributes)

                    if ncbi_match:
                        gene_id = ncbi_match.group(1)
                        self.chr_genes[chrom].append((start, gene_id))
                        self.gene_positions[gene_id] = (chrom, start)
                    elif braker_match:
                        gene_id = braker_match.group(1)
                        self.chr_genes[chrom].append((start, gene_id))
                        self.gene_positions[gene_id] = (chrom, start)

                # Extract protein-gene mapping (NCBI only - BRAKER3 inferred from proteome)
                if feature_type == 'CDS' and 'protein_id=' in attributes:
                    protein_match = re.search(r'protein_id=([^;\s]+)', attributes)
                    gene_match = re.search(r'gene=(LOC\d+)', attributes)

                    if protein_match and gene_match:
                        protein_id = protein_match.group(1)
                        gene_id = gene_match.group(1)
                        self.protein_to_gene[protein_id] = gene_id
                        self.gene_to_proteins[gene_id].add(protein_id)

                # For BRAKER3 mRNA features, map transcript to gene
                if feature_type == 'mRNA' and 'Parent=' in attributes:
                    # ID=g1.t1;Parent=g1;
                    transcript_match = re.search(r'ID=(g\d+\.t\d+);', attributes)
                    parent_match = re.search(r'Parent=(g\d+);', attributes)

                    if transcript_match and parent_match:
                        transcript_id = transcript_match.group(1)
                        gene_id = parent_match.group(1)
                        self.protein_to_gene[transcript_id] = gene_id
                        self.gene_to_proteins[gene_id].add(transcript_id)

        # Sort genes by position
        for chrom in self.chr_genes:
            self.chr_genes[chrom].sort()

        print(f"    Mapped {len(self.protein_to_gene)} proteins to {len(self.gene_to_proteins)} genes")

    def _use_braker_pattern(self, proteome_file):
        """For BRAKER3 annotations, extract gene from protein ID pattern."""
        print(f"  Using BRAKER3 pattern matching for {self.genome_code}")

        if proteome_file and Path(proteome_file).exists():
            seen_genes = set()

            # Parse proteome to get gene order
            for record in SeqIO.parse(proteome_file, 'fasta'):
                protein_id = record.id

                # Extract gene from protein ID (g123.t1 -> g123)
                match = re.match(r'(g\d+)\.t\d+', protein_id)
                if match:
                    gene_id = match.group(1)
                    self.protein_to_gene[protein_id] = gene_id
                    self.gene_to_proteins[gene_id].add(protein_id)

                    # Track gene order (first occurrence)
                    if gene_id not in seen_genes:
                        seen_genes.add(gene_id)
                        self.gene_order.append(gene_id)
                        # Use position in list as pseudo-coordinate
                        self.gene_positions[gene_id] = ('scaffold', len(self.gene_order))

            print(f"    Mapped {len(self.protein_to_gene)} proteins to {len(self.gene_order)} genes (BRAKER3)")

    def get_gene_for_protein(self, protein_id):
        """Get gene ID for a protein, handling BRAKER3 format."""
        if protein_id in self.protein_to_gene:
            return self.protein_to_gene[protein_id]

        # Try BRAKER3 pattern (g123.t1 -> g123)
        match = re.match(r'(g\d+)\.t\d+', protein_id)
        if match:
            return match.group(1)

        return None

    def get_flanking_genes(self, target_gene, flanking_distance_kb=1000):
        """Get flanking genes within genomic distance, excluding target itself.

        Args:
            target_gene: Gene ID to find flanking genes around
            flanking_distance_kb: Maximum distance in kb from target gene

        Returns:
            List of gene IDs within genomic distance of target
        """

        if target_gene not in self.gene_positions:
            return []

        chrom, target_pos = self.gene_positions[target_gene]
        flanking_genes = []

        # Check if NCBI or BRAKER3
        if self.gene_order:  # BRAKER3
            # For BRAKER3, get all genes on same scaffold within distance
            for gene in self.gene_order:
                if gene == target_gene:
                    continue
                if gene in self.gene_positions:
                    g_chrom, g_pos = self.gene_positions[gene]
                    if g_chrom == chrom:
                        distance_kb = abs(g_pos - target_pos) / 1000
                        if distance_kb <= flanking_distance_kb:
                            flanking_genes.append(gene)
        else:  # NCBI
            # For NCBI, chr_genes has [(pos, gene_id), ...]
            if chrom in self.chr_genes:
                for pos, gene_id in self.chr_genes[chrom]:
                    if gene_id == target_gene:
                        continue
                    distance_kb = abs(pos - target_pos) / 1000
                    if distance_kb <= flanking_distance_kb:
                        flanking_genes.append(gene_id)

        return flanking_genes

    def get_representative_protein(self, gene_id):
        """Get one representative protein for a gene."""
        if gene_id in self.gene_to_proteins:
            # Return first alphabetically for consistency
            return sorted(self.gene_to_proteins[gene_id])[0]

        # For BRAKER3, might just be g123.t1
        if re.match(r'g\d+$', gene_id):
            return f"{gene_id}.t1"

        return None


def find_paralogs_in_genome(query_protein_seq, genome_code, mapper):
    """Find all paralogs of query in a genome using DIAMOND."""

    genome_info = LANDMARKS[genome_code]
    # Create unique temp files to avoid race conditions in parallel array jobs
    import os
    pid = os.getpid()
    temp_query = Path(f"temp_query_{genome_code}_{pid}.faa")
    temp_output = Path(f"temp_hits_{genome_code}_{pid}.tsv")

    # Write query
    with open(temp_query, 'w') as f:
        SeqIO.write([query_protein_seq], f, 'fasta')

    # Run DIAMOND
    print(f"\n  Running DIAMOND search in {genome_code}...")
    subprocess.run([
        'diamond', 'blastp',
        '--query', str(temp_query),
        '--db', genome_info['db'],
        '--out', str(temp_output),
        '--outfmt', '6',
        '--evalue', str(EVALUE_THRESHOLD),
        '--max-target-seqs', '100',
        '--threads', '2'
    ], check=True, capture_output=True)

    # Parse hits
    hits = pd.read_csv(temp_output, sep='\t', header=None,
                      names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch',
                            'gapopen', 'qstart', 'qend', 'sstart', 'send',
                            'evalue', 'bitscore'])

    # Filter by percent identity
    hits = hits[hits['pident'] >= MIN_PIDENT]

    # Map to genes and deduplicate
    hits['target_gene'] = hits['sseqid'].apply(mapper.get_gene_for_protein)
    unique_genes = hits[hits['target_gene'].notna()]['target_gene'].unique()

    print(f"    Found {len(hits)} protein hits (≥{MIN_PIDENT}% pident) -> {len(unique_genes)} unique genes")

    # Clean up
    temp_query.unlink()
    temp_output.unlink()

    return hits, unique_genes


def detect_input_tandem_clusters(input_locs, mapper, max_distance_kb=50):
    """
    Detect tandem clusters among input LOCs.

    Returns list of clusters, where each cluster is a list of tandem LOCs.
    Non-tandem LOCs are returned as single-element clusters.
    """

    # Get positions for all input LOCs
    positions = {}
    for loc in input_locs:
        if loc in mapper.gene_positions:
            chrom, start = mapper.gene_positions[loc]
            positions[loc] = {'chrom': chrom, 'start': start}
        else:
            print(f"  Warning: {loc} not found in BK genome")

    # Sort by chromosome and position
    sorted_locs = sorted(positions.keys(),
                        key=lambda x: (positions[x]['chrom'], positions[x]['start']))

    # Cluster tandem genes
    clusters = []
    current_cluster = []

    for i, loc in enumerate(sorted_locs):
        if not current_cluster:
            current_cluster.append(loc)
        else:
            prev_loc = current_cluster[-1]

            # Check if on same chromosome and within max distance
            # Using start-to-start distance since we don't have gene lengths
            if (positions[loc]['chrom'] == positions[prev_loc]['chrom'] and
                abs(positions[loc]['start'] - positions[prev_loc]['start']) < max_distance_kb * 1000):
                current_cluster.append(loc)
            else:
                # Start new cluster
                clusters.append(current_cluster)
                current_cluster = [loc]

    # Add final cluster
    if current_cluster:
        clusters.append(current_cluster)

    # Add any LOCs that weren't in positions as single-element clusters
    for loc in input_locs:
        if loc not in positions:
            clusters.append([loc])

    return clusters


def detect_tandems(target_genes, flanking_genes):
    """Detect tandem duplications where targets appear in flanking."""

    tandems = []
    target_set = set(target_genes)

    for gene in flanking_genes:
        if gene in target_set:
            tandems.append(gene)

    return tandems


def extract_flanking_proteins(flanking_genes, mapper, proteome_file):
    """Extract protein sequences for flanking genes."""

    flanking_proteins = {}
    proteome = SeqIO.to_dict(SeqIO.parse(proteome_file, 'fasta'))

    for gene in flanking_genes:
        protein_id = mapper.get_representative_protein(gene)
        if protein_id and protein_id in proteome:
            flanking_proteins[gene] = proteome[protein_id]
        elif f"{protein_id}.1" in proteome:  # Try with version
            flanking_proteins[gene] = proteome[f"{protein_id}.1"]

    return flanking_proteins


def name_locus(genome_code, chromosome, locus_number):
    """Create systematic locus name."""

    # Use proper chromosome mappings
    if genome_code == 'BK':
        # BK has 10 chromosomes
        bk_chr_map = {
            'NC_046657.1': 'chr1', 'NC_046658.1': 'chr2', 'NC_046659.1': 'chr3',
            'NC_046660.1': 'chr4', 'NC_046661.1': 'chr5', 'NC_046662.1': 'chr6',
            'NC_046663.1': 'chr7', 'NC_046664.1': 'chr8', 'NC_046665.1': 'chr9',
            'NC_046666.1': 'chr10'
        }
        chr_name = bk_chr_map.get(chromosome, f"unk{chromosome[-4:]}")
    elif chromosome.startswith('NW_'):
        # Unplaced scaffolds (like in LB) - use last 4 digits of accession
        chr_name = f"scf{chromosome.split('.')[-2][-4:]}"
    elif chromosome.startswith('scaffold'):
        # BRAKER3 scaffolds
        chr_name = chromosome.replace('scaffold_', 'scf')[:10]
    else:
        # Fallback
        chr_name = chromosome.replace('_', '')[:10]

    # Convert locus number to letter (a, b, c...)
    locus_letter = chr(ord('a') + locus_number)

    return f"{genome_code}_{chr_name}_{locus_letter}"


def main():
    if len(sys.argv) < 3:
        print("Usage: python 01_discover_paralogs_gene_level.py <LOC_ID(s)> <output_dir> [gene_family]")
        print("  LOC_ID(s): Single LOC or comma-separated LOCs (e.g., 'LOC117167432' or 'LOC117167432,LOC117167433')")
        print("  gene_family: Optional gene family name for batch processing")
        sys.exit(1)

    target_input = sys.argv[1]
    output_dir = Path(sys.argv[2])
    gene_family = sys.argv[3] if len(sys.argv) > 3 else "unknown"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Parse input LOCs (may be comma-separated)
    input_locs = [loc.strip() for loc in target_input.split(',')]

    print("="*80)
    print("GENE-LEVEL PARALOG DISCOVERY WITH SYNTENY")
    print("="*80)
    print(f"\nInput LOCs: {len(input_locs)}")
    for loc in input_locs:
        print(f"  - {loc}")
    print(f"Output: {output_dir}")

    # Step 1: Build mappers for all genomes
    print("\n[1] Building genome mappers...")
    mappers = {}
    for genome_code, info in LANDMARKS.items():
        mappers[genome_code] = GenomicMapper(genome_code, info.get('gff'), info.get('proteome'))

    # Step 1b: Detect tandem clusters among input LOCs
    print("\n[1b] Detecting tandem clusters in input LOCs...")
    bk_mapper = mappers['BK']
    tandem_clusters = detect_input_tandem_clusters(input_locs, bk_mapper)

    print(f"  Found {len(tandem_clusters)} clusters:")
    for i, cluster in enumerate(tandem_clusters):
        if len(cluster) > 1:
            print(f"    Cluster {i+1} (TANDEM): {', '.join(cluster)}")
        else:
            print(f"    Cluster {i+1}: {cluster[0]}")

    # Use first LOC from each cluster as representative
    representative_locs = [cluster[0] for cluster in tandem_clusters]
    print(f"\n  Using {len(representative_locs)} representatives for paralog discovery")

    # Save tandem clustering info
    cluster_info = {
        'input_locs': input_locs,
        'clusters': [[loc for loc in cluster] for cluster in tandem_clusters],
        'representatives': representative_locs,
        'tandem_cluster_count': sum(1 for c in tandem_clusters if len(c) > 1)
    }
    with open(output_dir / "input_tandem_clusters.json", 'w') as f:
        json.dump(cluster_info, f, indent=2)

    # Step 2: Find query protein in BK for first representative
    print("\n[2] Finding query protein...")
    query_protein = None
    query_gene = None

    bk_proteome = SeqIO.to_dict(SeqIO.parse(LANDMARKS['BK']['proteome'], 'fasta'))

    # Use first representative LOC for paralog discovery
    target_id = representative_locs[0]
    if target_id.startswith('LOC'):
        # Find a protein for this LOC
        for protein_id in mappers['BK'].gene_to_proteins.get(target_id, []):
            if protein_id in bk_proteome:
                query_protein = bk_proteome[protein_id]
                query_gene = target_id
                print(f"  Using {protein_id} as representative for {target_id}")
                break
    else:
        # Direct protein ID
        if target_id in bk_proteome:
            query_protein = bk_proteome[target_id]
            query_gene = mappers['BK'].get_gene_for_protein(target_id)
            print(f"  Found {target_id} -> gene {query_gene}")

    if not query_protein:
        print(f"ERROR: Could not find {target_id} in BK proteome")
        sys.exit(1)

    # Step 3: Find paralogs in all landmark genomes
    print("\n[3] Finding paralogs in landmark genomes...")
    all_paralogs = {}

    for genome_code in LANDMARKS:
        hits, unique_genes = find_paralogs_in_genome(query_protein, genome_code, mappers[genome_code])
        all_paralogs[genome_code] = {
            'hits': hits,
            'unique_genes': list(unique_genes),
            'gene_count': len(unique_genes)
        }

    # Step 4: Process each unique locus
    print("\n[4] Processing unique loci...")
    unique_loci = []
    locus_counter = defaultdict(int)

    for genome_code, paralog_data in all_paralogs.items():
        print(f"\n  {genome_code}: {paralog_data['gene_count']} unique genes")

        for target_gene in paralog_data['unique_genes']:
            # Get flanking genes
            flanking_genes = mappers[genome_code].get_flanking_genes(target_gene, FLANKING_DISTANCE_KB)

            if not flanking_genes:
                print(f"    Warning: No flanking genes found for {target_gene}")
                continue

            # Detect tandems
            tandems = detect_tandems(paralog_data['unique_genes'], flanking_genes)

            # Get chromosome
            if target_gene in mappers[genome_code].gene_positions:
                chrom = mappers[genome_code].gene_positions[target_gene][0]
            else:
                chrom = 'unknown'

            # Create locus entry
            locus_name = name_locus(genome_code, chrom, locus_counter[genome_code])
            locus_counter[genome_code] += 1

            # Extract flanking proteins
            flanking_proteins = extract_flanking_proteins(
                flanking_genes, mappers[genome_code], LANDMARKS[genome_code]['proteome']
            )

            # Save flanking proteins
            flanking_file = output_dir / f"{locus_name}_flanking.faa"
            with open(flanking_file, 'w') as f:
                for gene, seq in flanking_proteins.items():
                    seq.id = f"{seq.id}|{gene}"
                    SeqIO.write([seq], f, 'fasta')

            locus = {
                'locus_id': locus_name,
                'genome': genome_code,
                'chromosome': chrom,
                'target_gene': target_gene,
                'flanking_count': len(flanking_genes),
                'tandem_count': len(tandems),
                'is_tandem': len(tandems) > 0,
                'flanking_file': str(flanking_file),
                'gene_family': gene_family
            }

            unique_loci.append(locus)

            print(f"    {locus_name}:")
            print(f"      Target: {target_gene}")
            print(f"      Flanking: {len(flanking_genes)} genes -> {len(flanking_proteins)} proteins")
            print(f"      Tandems: {len(tandems)}")

    # Step 5: Save results
    print("\n[5] Saving results...")

    # Save locus summary
    locus_df = pd.DataFrame(unique_loci)
    locus_df.to_csv(output_dir / "unique_loci.tsv", sep='\t', index=False)

    # Save detailed paralog info
    with open(output_dir / "paralog_details.json", 'w') as f:
        # Convert DataFrames to dicts for JSON serialization
        for genome in all_paralogs:
            all_paralogs[genome]['hits'] = all_paralogs[genome]['hits'].to_dict('records')
        json.dump(all_paralogs, f, indent=2)

    print(f"\nResults saved to {output_dir}/")
    print(f"  - unique_loci.tsv: {len(unique_loci)} loci across {len(LANDMARKS)} genomes")
    print(f"  - paralog_details.json: Detailed BLAST results")
    print(f"  - *_flanking.faa: Flanking proteins for each locus")

    # Summary
    print("\n" + "="*80)
    print("SUMMARY")
    print("-"*80)
    for genome in LANDMARKS:
        genome_loci = [l for l in unique_loci if l['genome'] == genome]
        tandem_loci = [l for l in genome_loci if l['is_tandem']]
        print(f"{genome}: {len(genome_loci)} loci ({len(tandem_loci)} tandem arrays)")

    return unique_loci


if __name__ == "__main__":
    main()