#!/usr/bin/env python3
"""
Phase 1: Gene-level paralog discovery with isoform deduplication.

Key improvements:
1. Maps proteins to LOC IDs to handle isoforms properly
2. Extracts flanking at GENE level (not protein level) in ORDERED upstream/downstream format
3. Detects tandem duplications among input LOCs
4. Detects tandem duplications among paralogs
5. Names loci systematically by chromosome

Usage:
    python 01_discover_paralogs_gene_level.py \\
        --loc-ids LOC117167432,LOC117167433 \\
        --output-dir outputs/ferritin_phase1 \\
        --gene-family ferritin_MC102 \\
        --flanking-distance 1000 \\
        --evalue 1e-3 \\
        --min-pident 40

When multiple LOCs are provided, the script detects if any are tandem duplications
(adjacent on the same chromosome) and groups them into clusters. One representative
from each cluster is used for paralog discovery.

**NEW**: Flanking genes are now saved in genomic order with U/D labels:
  - U1 = closest upstream gene
  - U2 = second closest upstream gene
  - D1 = closest downstream gene
  - D2 = second closest downstream gene
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
import argparse
import xml.etree.ElementTree as ET

# Configuration
# Import landmark paths from data_paths module (handles per-species directory structure)
from data_paths import get_landmark_dict
LANDMARKS = get_landmark_dict()

# Genome-specific flanking distance caps to normalize gene counts
# BK: 10.4 genes/Mb, LB: 42.5 genes/Mb (4.1x difference)
# Target: ~35 flanking genes per locus for fair synteny comparison
FLANKING_DISTANCE_KB = {
    'BK': 1500,  # 1.5 Mb → ~31 genes expected
    'LB': 400    # 400 kb → ~34 genes expected
}
MAX_FLANKING_GENES = 50  # Hard cap to prevent outliers
EVALUE_THRESHOLD = 1e-3
MIN_PIDENT = 40.0  # Minimum percent identity to avoid distant paralogs


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
                    end = int(parts[4])  # Get end position
                    # NCBI format: Name=LOC117167432
                    ncbi_match = re.search(r'Name=(LOC\d+)', attributes)
                    # BRAKER3 format: ID=g10440;
                    braker_match = re.search(r'ID=(g\d+);', attributes)

                    if ncbi_match:
                        gene_id = ncbi_match.group(1)
                        self.chr_genes[chrom].append((start, gene_id))
                        self.gene_positions[gene_id] = (chrom, start, end)
                    elif braker_match:
                        gene_id = braker_match.group(1)
                        self.chr_genes[chrom].append((start, gene_id))
                        self.gene_positions[gene_id] = (chrom, start, end)

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
                        # Use position in list as pseudo-coordinate (start, end)
                        pos = len(self.gene_order)
                        self.gene_positions[gene_id] = ('scaffold', pos, pos)

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

    def get_flanking_genes(self, target_gene, flanking_distance_kb=1000, max_genes=None):
        """Get flanking genes split into upstream and downstream, sorted by position.

        Args:
            target_gene: Gene ID to find flanking genes around
            flanking_distance_kb: Maximum distance in kb from target gene
            max_genes: Maximum total flanking genes to return (splits evenly between upstream/downstream)

        Returns:
            Tuple of (upstream_genes, downstream_genes)
            - upstream_genes: List sorted by position (closest to target first)
            - downstream_genes: List sorted by position (closest to target first)
        """

        if target_gene not in self.gene_positions:
            return [], []

        target_pos_tuple = self.gene_positions[target_gene]
        chrom = target_pos_tuple[0]
        target_pos = target_pos_tuple[1]
        candidates = []  # (distance, position, gene_id) - track distance for sorting

        # Check if NCBI or BRAKER3
        if self.gene_order:  # BRAKER3
            # For BRAKER3, get all genes on same scaffold within distance
            for gene in self.gene_order:
                if gene == target_gene:
                    continue
                if gene in self.gene_positions:
                    g_pos_tuple = self.gene_positions[gene]
                    g_chrom = g_pos_tuple[0]
                    g_pos = g_pos_tuple[1]
                    if g_chrom == chrom:
                        distance_kb = abs(g_pos - target_pos) / 1000
                        if distance_kb <= flanking_distance_kb:
                            candidates.append((distance_kb, g_pos, gene))
        else:  # NCBI
            # For NCBI, chr_genes has [(pos, gene_id), ...]
            if chrom in self.chr_genes:
                for pos, gene_id in self.chr_genes[chrom]:
                    if gene_id == target_gene:
                        continue
                    distance_kb = abs(pos - target_pos) / 1000
                    if distance_kb <= flanking_distance_kb:
                        candidates.append((distance_kb, pos, gene_id))

        # Sort by position and split into upstream/downstream
        upstream = [(dist, pos, gene) for dist, pos, gene in candidates if pos < target_pos]
        downstream = [(dist, pos, gene) for dist, pos, gene in candidates if pos >= target_pos]

        # Sort upstream by distance (closest first) then extract genes
        upstream.sort(key=lambda x: x[0])
        upstream_genes = [gene for dist, pos, gene in upstream]

        # Sort downstream by distance (closest first) then extract genes
        downstream.sort(key=lambda x: x[0])
        downstream_genes = [gene for dist, pos, gene in downstream]

        # Apply max_genes cap if specified
        if max_genes is not None:
            max_per_side = max_genes // 2
            upstream_genes = upstream_genes[:max_per_side]
            downstream_genes = downstream_genes[:max_per_side]

        return upstream_genes, downstream_genes

    def get_flanking_genes_for_cluster(self, cluster_genes, flanking_distance_kb=1000, max_genes=None, exclude_genes=None,
                                       include_internal_nonfamily=True, max_internal_genes=None):
        """Get flanking genes around a tandem cluster treated as a unit.

        Defines cluster boundaries by the leftmost and rightmost gene positions.
        - Upstream genes: same chromosome with positions < left boundary
        - Downstream genes: same chromosome with positions > right boundary
        - Internal genes (optional): same chromosome with positions between left and right boundaries,
          excluding any cluster members and any genes listed in exclude_genes (e.g., same-family genes)

        Cluster member genes AND all genes in exclude_genes are excluded from the returned lists.

        Args:
            cluster_genes: Genes defining the cluster boundaries
            flanking_distance_kb: Max distance for flanking genes (applies to upstream/downstream only)
            max_genes: Max number of flanking genes total (apportioned evenly to upstream/downstream)
            exclude_genes: Additional genes to exclude (e.g., all family members)
            include_internal_nonfamily: When True, include non-family genes that lie inside the cluster span
            max_internal_genes: Optional cap on number of internal genes to return (closest-to-center first)

        Returns (upstream_genes, downstream_genes, internal_genes) where each list is ordered by proximity
        to the relevant boundary (up/down) or to the cluster center (internal).
        """
        # Filter to genes with known positions and collect chrom/positions
        positions = {}
        for g in cluster_genes:
            if g in self.gene_positions:
                pos_tuple = self.gene_positions[g]
                chrom = pos_tuple[0]
                pos = pos_tuple[1]
                positions[g] = (chrom, pos)

        if not positions:
            return [], []

        # Ensure all genes are on same chromosome; if not, fallback to per-gene flanking
        chroms = {c for c, _ in positions.values()}
        if len(chroms) != 1:
            # Fallback: choose the majority chromosome
            chrom = max(chroms, key=lambda c: sum(1 for (cc, _) in positions.values() if cc == c))
        else:
            chrom = next(iter(chroms))

        left = min(pos for (c, pos) in positions.values() if c == chrom)
        right = max(pos for (c, pos) in positions.values() if c == chrom)

        # Collect candidate flanking genes
        candidates = []  # (distance_to_boundary_kb, position, gene_id, side)
        internal_candidates = []  # (distance_to_center_kb, position, gene_id)

        # Build exclusion set from both cluster members and additional excludes
        exclude_set = set(cluster_genes)
        if exclude_genes:
            exclude_set.update(exclude_genes)

        def maybe_add(gene_id, pos):
            if gene_id in exclude_set:
                return
            if pos < left:
                dist_kb = (left - pos) / 1000
                if dist_kb <= flanking_distance_kb:
                    candidates.append((dist_kb, pos, gene_id, 'U'))
            elif pos > right:
                dist_kb = (pos - right) / 1000
                if dist_kb <= flanking_distance_kb:
                    candidates.append((dist_kb, pos, gene_id, 'D'))
            elif include_internal_nonfamily and left <= pos <= right:
                # Internal non-family gene between cluster members
                # Order internal genes by proximity to cluster center (closest first)
                center = (left + right) // 2
                dist_center_kb = abs(pos - center) / 1000
                internal_candidates.append((dist_center_kb, pos, gene_id))

        if self.gene_order:  # BRAKER3
            for gene in self.gene_order:
                if gene in self.gene_positions:
                    g_pos_tuple = self.gene_positions[gene]
                    g_chrom = g_pos_tuple[0]
                    g_pos = g_pos_tuple[1]
                    if g_chrom == chrom:
                        maybe_add(gene, g_pos)
        else:  # NCBI
            if chrom in self.chr_genes:
                for pos, gene_id in self.chr_genes[chrom]:
                    maybe_add(gene_id, pos)

        # Split and sort each side by increasing distance (closest first)
        upstream = [(d, p, g) for (d, p, g, s) in candidates if s == 'U']
        downstream = [(d, p, g) for (d, p, g, s) in candidates if s == 'D']
        upstream.sort(key=lambda x: x[0])
        downstream.sort(key=lambda x: x[0])

        upstream_genes = [g for d, p, g in upstream]
        downstream_genes = [g for d, p, g in downstream]

        # Internal: sort by distance to center (closest first)
        internal_candidates.sort(key=lambda x: x[0])
        internal_genes = [g for d, p, g in internal_candidates]

        # Optional cap on internal genes
        if max_internal_genes is not None and max_internal_genes >= 0:
            internal_genes = internal_genes[:max_internal_genes]

        # Apply max cap evenly between sides
        if max_genes is not None:
            max_per_side = max_genes // 2
            upstream_genes = upstream_genes[:max_per_side]
            downstream_genes = downstream_genes[:max_per_side]

        return upstream_genes, downstream_genes, internal_genes

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


def run_tblastn(query_file, genome_db, output_xml, evalue="1e-5", max_targets=50):
    """Run tBLASTn search (copied from Phase 3)."""
    cmd = [
        'tblastn',
        '-query', str(query_file),
        '-db', str(genome_db),
        '-outfmt', '5',  # XML format
        '-evalue', str(evalue),
        '-max_target_seqs', str(max_targets),
        '-num_threads', '2'
    ]

    with open(output_xml, 'w') as f:
        result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, text=True)

    return result.returncode == 0


def parse_blast_xml(xml_file):
    """Parse BLAST XML and extract hits (copied from Phase 3)."""
    hits = []

    try:
        tree = ET.parse(xml_file)
        root = tree.getroot()

        for iteration in root.findall('.//Iteration'):
            query_id = iteration.find('Iteration_query-def').text.split()[0]

            for hit in iteration.findall('.//Hit'):
                hit_def = hit.find('Hit_def').text
                hit_accession = hit.find('Hit_accession').text

                for hsp in hit.findall('.//Hsp'):
                    evalue = float(hsp.find('Hsp_evalue').text)
                    bitscore = float(hsp.find('Hsp_bit-score').text)
                    identity = int(hsp.find('Hsp_identity').text)
                    align_length = int(hsp.find('Hsp_align-len').text)

                    sstart = int(hsp.find('Hsp_hit-from').text)
                    send = int(hsp.find('Hsp_hit-to').text)

                    # Determine strand
                    if sstart < send:
                        strand = '+'
                        start = sstart
                        end = send
                    else:
                        strand = '-'
                        start = send
                        end = sstart

                    hits.append({
                        'qseqid': query_id,
                        'sseqid': hit_accession,
                        'scaffold_desc': hit_def,
                        'strand': strand,
                        'coord_start': start,
                        'coord_end': end,
                        'evalue': evalue,
                        'bitscore': bitscore,
                        'pident': identity / align_length * 100,
                        'length': align_length
                    })

    except Exception as e:
        print(f"    Parse error: {e}")

    return hits


def find_paralogs_with_tblastn(query_protein_seq, genome_code, mapper):
    """
    Find paralogs using tblastn-first approach (for BRAKER3 genomes).

    This avoids 600kb coordinate mismatch by:
    1. Running tblastn on genome FIRST → finds actual target location
    2. Building locus around THAT location (±500kb window)
    3. Parsing BRAKER3 GFF for flanking genes within window

    Returns:
        tuple: (hits_list, unique_genes_set) compatible with DIAMOND approach
    """

    genome_info = LANDMARKS[genome_code]
    genome_db = genome_info['genome_db']
    gff_file = genome_info['gff']

    # Create temp query file
    import os
    pid = os.getpid()
    temp_query = Path(f"temp_tblastn_query_{genome_code}_{pid}.faa")
    temp_xml = Path(f"temp_tblastn_{genome_code}_{pid}.xml")

    try:
        # Write query
        with open(temp_query, 'w') as f:
            SeqIO.write([query_protein_seq], f, 'fasta')

        # Run tblastn
        print(f"    Running tblastn on genome database...")
        if not run_tblastn(temp_query, genome_db, temp_xml, evalue="1e-5", max_targets=50):
            print(f"    tblastn failed")
            return [], set()

        # Parse hits
        hits = parse_blast_xml(temp_xml)
        if not hits:
            print(f"    No tblastn hits found")
            return [], set()

        print(f"    Found {len(hits)} tblastn hits")

        # Filter by coverage and identity
        filtered_hits = []
        for hit in hits:
            # Relaxed thresholds for tblastn (fragmented across exons)
            if hit['pident'] >= 35.0:  # 35% identity
                filtered_hits.append(hit)

        if not filtered_hits:
            print(f"    No hits passed filters")
            return [], set()

        print(f"    After filtering (≥35% identity): {len(filtered_hits)} hits")

        # Group by scaffold and take best hit per scaffold
        from collections import defaultdict
        scaffold_hits = defaultdict(list)
        for hit in filtered_hits:
            scaffold_hits[hit['sseqid']].append(hit)

        best_hits = []
        for scaffold, shits in scaffold_hits.items():
            best_hit = max(shits, key=lambda h: h['bitscore'])
            best_hits.append(best_hit)

        print(f"    Best hits on {len(best_hits)} scaffolds")

        # For each target location, find flanking genes in GFF
        unique_genes = set()
        all_hit_records = []

        for hit in best_hits:
            scaffold = hit['sseqid']
            start = hit['coord_start']
            end = hit['coord_end']
            strand = hit['strand']

            # Define ±500kb flanking window
            window_bp = 500 * 1000
            window_start = max(1, start - window_bp)
            window_end = end + window_bp

            print(f"\n    Target at {scaffold}:{start}-{end} ({strand})")
            print(f"      Flanking window: {window_start}-{window_end} (±500kb)")

            # Parse GFF for genes in window
            genes_in_window = []
            with open(gff_file) as f:
                for line in f:
                    if line.startswith('#'):
                        continue

                    parts = line.strip().split('\t')
                    if len(parts) < 9:
                        continue

                    gff_scaffold = parts[0]
                    feature_type = parts[2]
                    gene_start = int(parts[3])
                    gene_end = int(parts[4])
                    gene_strand = parts[6]
                    attributes = parts[8]

                    # Check scaffold match (handle long names with prefix matching)
                    # BLAST may return "CM036341" while GFF has "CM036341.1_Telenomus..."
                    scaffold_match = (
                        gff_scaffold == scaffold or
                        gff_scaffold.startswith(scaffold + '_') or
                        gff_scaffold.startswith(scaffold + '.')
                    )
                    if not scaffold_match:
                        continue

                    # Only process genes within window
                    if gene_end < window_start or gene_start > window_end:
                        continue

                    # Extract gene ID (BRAKER3: ID=g10440;)
                    if feature_type == 'gene':
                        match = re.search(r'ID=(g\d+);', attributes)
                        if match:
                            gene_id = match.group(1)
                            gene_mid = (gene_start + gene_end) // 2
                            target_mid = (start + end) // 2
                            distance = abs(gene_mid - target_mid)

                            genes_in_window.append({
                                'gene_id': gene_id,
                                'start': gene_start,
                                'end': gene_end,
                                'strand': gene_strand,
                                'distance': distance
                            })

            if not genes_in_window:
                print(f"      Warning: No genes found in window, skipping")
                continue

            print(f"      Found {len(genes_in_window)} genes in window")

            # Use first gene as the "target gene" for this locus
            # (maintains compatibility with downstream pipeline)
            genes_in_window.sort(key=lambda g: g['distance'])
            target_gene = genes_in_window[0]['gene_id']

            unique_genes.add(target_gene)

            # Create hit record (mimic DIAMOND format)
            hit_record = {
                'qseqid': str(query_protein_seq.id),
                'sseqid': target_gene,
                'pident': hit['pident'],
                'evalue': hit['evalue'],
                'bitscore': hit['bitscore'],
                'target_gene': target_gene,
                'tblastn_coords': f"{scaffold}:{start}-{end}"
            }
            all_hit_records.append(hit_record)

        print(f"    Final: {len(unique_genes)} unique genes (loci)")

        return all_hit_records, unique_genes

    finally:
        # Cleanup
        temp_query.unlink(missing_ok=True)
        temp_xml.unlink(missing_ok=True)


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
            pos = mapper.gene_positions[loc]
            chrom = pos[0]
            start = pos[1]
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
            # Use start-to-start distance since we don't store gene end positions
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

    # Convert locus number to letter (a, b, c..., aa, bb, cc...)
    if locus_number < 26:
        locus_letter = chr(ord('a') + locus_number)
    else:
        # For numbers >= 26, use aa, bb, cc, etc.
        letter_index = (locus_number - 26) % 26
        repeat_count = ((locus_number - 26) // 26) + 2
        locus_letter = chr(ord('a') + letter_index) * repeat_count

    return f"{genome_code}_{chr_name}_{locus_letter}"


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Gene-level paralog discovery with ordered upstream/downstream flanking genes",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('--loc-ids', required=True, type=str,
                        help='LOC ID(s) - single or comma-separated (e.g., LOC117167432,LOC117167433)')
    parser.add_argument('--output-dir', required=True, type=Path,
                        help='Output directory for locus definitions and flanking proteins')
    # Note: flanking-distance is now genome-specific (hardcoded in script constants)
    # BK: 1500 kb, LB: 400 kb to normalize gene counts
    parser.add_argument('--evalue', type=float, default=EVALUE_THRESHOLD,
                        help=f'E-value threshold for DIAMOND search (default: {EVALUE_THRESHOLD})')
    parser.add_argument('--min-pident', type=float, default=MIN_PIDENT,
                        help=f'Minimum percent identity for paralogs (default: {MIN_PIDENT})')
    parser.add_argument('--gene-family', type=str, default='unknown',
                        help='Gene family name for classification (default: unknown)')

    return parser.parse_args()

def main():
    args = parse_args()

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Parse input LOCs (may be comma-separated)
    input_locs = [loc.strip() for loc in args.loc_ids.split(',')]

    # Update global parameters from arguments
    global EVALUE_THRESHOLD, MIN_PIDENT
    EVALUE_THRESHOLD = args.evalue
    MIN_PIDENT = args.min_pident

    print("="*80)
    print("GENE-LEVEL PARALOG DISCOVERY WITH SYNTENY")
    print("="*80)
    print(f"\nParameters:")
    print(f"  Flanking distance (genome-specific):")
    for genome, dist_kb in FLANKING_DISTANCE_KB.items():
        print(f"    {genome}: {dist_kb} kb")
    print(f"  Max flanking genes: {MAX_FLANKING_GENES}")
    print(f"  E-value threshold: {EVALUE_THRESHOLD}")
    print(f"  Min percent identity: {MIN_PIDENT}%")
    print(f"\nInput LOCs: {len(input_locs)}")
    for loc in input_locs:
        print(f"  - {loc}")
    print(f"Output: {args.output_dir}")

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
    with open(args.output_dir / "input_tandem_clusters.json", 'w') as f:
        json.dump(cluster_info, f, indent=2)

    # Step 2: Find query proteins in BK for ALL representatives
    print("\n[2] Finding query proteins...")

    bk_proteome = SeqIO.to_dict(SeqIO.parse(LANDMARKS['BK']['proteome'], 'fasta'))
    query_proteins = []  # List of (protein, gene_id) tuples

    for target_id in representative_locs:
        found = False
        if target_id.startswith('LOC'):
            # Find a protein for this LOC
            for protein_id in mappers['BK'].gene_to_proteins.get(target_id, []):
                if protein_id in bk_proteome:
                    query_proteins.append((bk_proteome[protein_id], target_id))
                    print(f"  Using {protein_id} as representative for {target_id}")
                    found = True
                    break
        else:
            # Direct protein ID
            if target_id in bk_proteome:
                gene_id = mappers['BK'].get_gene_for_protein(target_id)
                query_proteins.append((bk_proteome[target_id], gene_id))
                print(f"  Found {target_id} -> gene {gene_id}")
                found = True

        if not found:
            print(f"  Warning: Could not find {target_id} in BK proteome, skipping")

    if not query_proteins:
        print(f"ERROR: No query proteins found")
        sys.exit(1)

    print(f"\n  Total query proteins: {len(query_proteins)}")

    # Step 3: Find paralogs in all landmark genomes using ALL query proteins
    print("\n[3] Finding paralogs in landmark genomes (searching with all query proteins)...")
    all_paralogs = {}

    for genome_code in LANDMARKS:
        print(f"\n  {genome_code}:")

        # Accumulate results across all queries
        combined_hits = []
        combined_genes = set()

        for query_protein, query_gene_id in query_proteins:
            print(f"    Query: {query_gene_id}")

            # CRITICAL: Use different approaches based on genome type
            if genome_code in ['TR', 'DR']:
                # BRAKER3 genomes: Use tblastn-first to avoid coordinate mismatch
                hits, unique_genes = find_paralogs_with_tblastn(query_protein, genome_code, mappers[genome_code])
            else:
                # NCBI genomes (BK, LB): Use DIAMOND on proteome (ground truth)
                hits, unique_genes = find_paralogs_in_genome(query_protein, genome_code, mappers[genome_code])

            combined_hits.extend(hits)
            combined_genes.update(unique_genes)

        print(f"    Total after combining all queries: {len(combined_genes)} unique genes")

        all_paralogs[genome_code] = {
            'hits': combined_hits,
            'unique_genes': list(combined_genes),
            'gene_count': len(combined_genes)
        }

    # Step 4: Process loci at cluster level (deduplicate tandems ahead of time)
    print("\n[4] Processing loci at tandem-cluster level…")
    unique_loci = []
    locus_counter = defaultdict(int)

    # Map target genes to input cluster metadata for reporting
    input_tandem_map = {}
    for cluster in tandem_clusters:
        for gene in cluster:
            input_tandem_map[gene] = {
                'is_input_tandem': len(cluster) > 1,
                'cluster_size': len(cluster),
                'cluster_members': cluster
            }

    for genome_code, paralog_data in all_paralogs.items():
        genes = [g for g in paralog_data['unique_genes'] if g in mappers[genome_code].gene_positions]
        print(f"\n  {genome_code}: {len(genes)} target genes with positions")
        if not genes:
            continue

        # Build per-gene flanking map without excluding family genes (for mutual-flanking detection)
        distance_cap = FLANKING_DISTANCE_KB.get(genome_code, 1000)
        per_gene_flanking = {}
        debug_flanking_rows = []
        for g in genes:
            up, dn = mappers[genome_code].get_flanking_genes(
                g, flanking_distance_kb=distance_cap, max_genes=MAX_FLANKING_GENES
            )
            per_gene_flanking[g] = set(up + dn)
            pos_tuple = mappers[genome_code].gene_positions.get(g, ("unknown", 0, 0))
            chrom = pos_tuple[0]
            pos = pos_tuple[1]
            debug_flanking_rows.append({
                'gene': g,
                'chromosome': chrom,
                'position': pos,
                'upstream_count': len(up),
                'downstream_count': len(dn),
                'flanking_total': len(up) + len(dn),
                'upstream_genes': ','.join(up),
                'downstream_genes': ','.join(dn)
            })

        # Determine mutual-flanking clusters (undirected components)
        visited = set()
        mf_clusters = []
        for g in genes:
            if g in visited:
                continue
            comp = set([g])
            stack = [g]
            visited.add(g)
            while stack:
                cur = stack.pop()
                for h in genes:
                    if h == cur or h in visited:
                        continue
                    if (h in per_gene_flanking.get(cur, set())) and (cur in per_gene_flanking.get(h, set())):
                        visited.add(h)
                        comp.add(h)
                        stack.append(h)
            mf_clusters.append(sorted(list(comp), key=lambda x: mappers[genome_code].gene_positions.get(x, ('', 0))[1]))

        print(f"    Mutual-flanking clusters: {len(mf_clusters)} (will emit one locus per cluster)")

        # Debug: compute mutual edges for reporting
        mutual_edges = set()
        for g in genes:
            for h in genes:
                if h <= g:
                    continue
                if (h in per_gene_flanking.get(g, set())) and (g in per_gene_flanking.get(h, set())):
                    mutual_edges.add((g, h))
        print(f"    Mutual links (edges): {len(mutual_edges)}")

        # For Phase 2 hand-off, exclude ALL family genes from flanking
        exclude_family = set(genes)

        cluster_debug_rows = []
        for idx, cluster in enumerate(mf_clusters, start=1):
            chrom = mappers[genome_code].gene_positions[cluster[0]][0]
            rep_gene = cluster[0]

            upstream_genes, downstream_genes, internal_genes = mappers[genome_code].get_flanking_genes_for_cluster(
                cluster,
                flanking_distance_kb=distance_cap,
                max_genes=MAX_FLANKING_GENES,
                exclude_genes=exclude_family,
                include_internal_nonfamily=True,
                max_internal_genes=None
            )

            # Extract proteins for cleaned flanking
            upstream_proteins = extract_flanking_proteins(
                upstream_genes, mappers[genome_code], LANDMARKS[genome_code]['proteome']
            )
            downstream_proteins = extract_flanking_proteins(
                downstream_genes, mappers[genome_code], LANDMARKS[genome_code]['proteome']
            )
            internal_proteins = extract_flanking_proteins(
                internal_genes, mappers[genome_code], LANDMARKS[genome_code]['proteome']
            )

            # Extract target protein(s) for the cluster
            target_proteins = extract_flanking_proteins(
                cluster, mappers[genome_code], LANDMARKS[genome_code]['proteome']
            )

            # Name and write
            locus_name = name_locus(genome_code, chrom, locus_counter[genome_code])
            locus_counter[genome_code] += 1
            flanking_file = args.output_dir / f"{locus_name}_flanking.faa"
            with open(flanking_file, 'w') as f:
                for i, gene in enumerate(upstream_genes, 1):
                    if gene in upstream_proteins:
                        seq = upstream_proteins[gene]
                        seq.id = f"{seq.id}|{gene}"
                        seq.description = f"U{i} {seq.description}"
                        SeqIO.write([seq], f, 'fasta')
                # Write internal (between boundary) non-family genes if present
                for i, gene in enumerate(internal_genes, 1):
                    if gene in internal_proteins:
                        seq = internal_proteins[gene]
                        seq.id = f"{seq.id}|{gene}"
                        seq.description = f"I{i} {seq.description}"
                        SeqIO.write([seq], f, 'fasta')
                for i, gene in enumerate(downstream_genes, 1):
                    if gene in downstream_proteins:
                        seq = downstream_proteins[gene]
                        seq.id = f"{seq.id}|{gene}"
                        seq.description = f"D{i} {seq.description}"
                        SeqIO.write([seq], f, 'fasta')

            # Write target proteins (for Phase 4 BLAST)
            target_dir = args.output_dir / locus_name
            target_dir.mkdir(exist_ok=True)
            target_file = target_dir / f"{locus_name}_targets.faa"
            with open(target_file, 'w') as f:
                for gene in cluster:
                    if gene in target_proteins:
                        seq = target_proteins[gene]
                        seq.id = f"{seq.id}|{gene}"
                        seq.description = f"TARGET {seq.description}"
                        SeqIO.write([seq], f, 'fasta')

            input_meta = None
            for g in cluster:
                if g in input_tandem_map:
                    input_meta = input_tandem_map[g]
                    break
            if input_meta is None:
                input_meta = {'is_input_tandem': False, 'cluster_size': 1, 'cluster_members': [rep_gene]}

            target_pos = mappers[genome_code].gene_positions.get(rep_gene)
            target_start = target_pos[1] if target_pos and len(target_pos) >= 2 else None
            target_end = target_pos[2] if target_pos and len(target_pos) >= 3 else None

            locus = {
                'locus_id': locus_name,
                'gene_family': args.gene_family,
                'genome': genome_code,
                'chromosome': chrom,
                'target_gene': rep_gene,
                'target_start': target_start,
                'target_end': target_end,
                'upstream_count': len(upstream_genes),
                'downstream_count': len(downstream_genes),
                'flanking_count': len(upstream_genes) + len(downstream_genes),
                'tandem_count': max(0, len(cluster) - 1),
                'is_tandem': len(cluster) > 1,
                'input_is_tandem': input_meta['is_input_tandem'],
                'input_cluster_size': input_meta['cluster_size'],
                'input_cluster_members': ','.join(input_meta['cluster_members']),
                'cluster_size': len(cluster),
                'cluster_members': ','.join(cluster),
                'flanking_file': str(flanking_file)
            }

            unique_loci.append(locus)

            print(f"    {locus_name} (cluster size {len(cluster)}): rep={rep_gene}")
            print(f"      Cleaned flanking (excl family): {len(upstream_genes)}U + {len(internal_genes)}I + {len(downstream_genes)}D = {len(upstream_genes)+len(internal_genes)+len(downstream_genes)} genes")

            # Add cluster debug row
            # Compute approximate boundary from member positions
            positions = [mappers[genome_code].gene_positions.get(g, (chrom, 0))[1] for g in cluster]
            left = min(positions) if positions else 0
            right = max(positions) if positions else 0
            cluster_debug_rows.append({
                'cluster_id': f"{genome_code}_cluster_{idx:03d}",
                'chromosome': chrom,
                'member_count': len(cluster),
                'members': ','.join(cluster),
                'boundary_start': left,
                'boundary_end': right,
                'cleaned_upstream_count': len(upstream_genes),
                'cleaned_internal_count': len(internal_genes),
                'cleaned_downstream_count': len(downstream_genes),
                'cleaned_flanking_total': len(upstream_genes) + len(internal_genes) + len(downstream_genes),
                'locus_id': locus_name,
                'flanking_file': str(flanking_file)
            })

        # Write debug files for this genome
        try:
            import pandas as pd  # already imported at top
            dbg_dir = args.output_dir
            # Per-gene flanking map
            if debug_flanking_rows:
                pd.DataFrame(debug_flanking_rows).to_csv(
                    dbg_dir / f"debug_flanking_{genome_code}.tsv", sep='\t', index=False
                )
            # Mutual edges
            if mutual_edges:
                import csv
                with open(dbg_dir / f"debug_mutual_edges_{genome_code}.tsv", 'w', newline='') as f:
                    w = csv.writer(f, delimiter='\t')
                    w.writerow(['gene_a', 'gene_b'])
                    for a, b in sorted(mutual_edges):
                        w.writerow([a, b])
            # Cluster summary
            if cluster_debug_rows:
                pd.DataFrame(cluster_debug_rows).to_csv(
                    dbg_dir / f"debug_clusters_{genome_code}.tsv", sep='\t', index=False
                )
        except Exception as e:
            print(f"    DEBUG write failed for {genome_code}: {e}")

    # Step 5: Save results
    print("\n[5] Saving results...")

    # Save locus summary (used by downstream pipeline scripts)
    locus_df = pd.DataFrame(unique_loci)
    locus_df.to_csv(args.output_dir / "locus_definitions.tsv", sep='\t', index=False)

    # Save detailed paralog info
    with open(args.output_dir / "paralog_details.json", 'w') as f:
        # Convert DataFrames to dicts for JSON serialization
        for genome in all_paralogs:
            hits = all_paralogs[genome]['hits']
            # Handle both DataFrame (DIAMOND) and list (tblastn) formats
            if isinstance(hits, pd.DataFrame):
                all_paralogs[genome]['hits'] = hits.to_dict('records')
            # else: already a list, keep as is
        json.dump(all_paralogs, f, indent=2)

    print(f"\nResults saved to {args.output_dir}/")
    print(f"  - locus_definitions.tsv: {len(unique_loci)} loci across {len(LANDMARKS)} genomes")
    print(f"  - paralog_details.json: Detailed BLAST results")
    print(f"  - *_flanking.faa: Flanking proteins for each locus (ORDERED: U1...Un, D1...Dn)")

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
