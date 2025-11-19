#!/usr/bin/env python3
"""
Phase 4 (redesign): tBLASTn target detection with greedy query-coverage clustering.

Key differences vs legacy:
- Build combined multi-query FASTA from Phase1_v2 per-locus targets.
- Run tBLASTn once per genome.
- Cluster HSPs per (query, scaffold, strand) using greedy UNION query-coverage
  expansion with a paralog threshold (~120%) and a 200 kb gap cap.
- Emit one row per gene-level hit with genomic coordinates; no extra locus
  clustering is performed here (locus reasoning happens in Phase 5).

Outputs (per family):
- all_target_loci.tsv        # per-(query,gene)-level hits
- combined_targets.faa       # combined query protein sequences
- blast_xml/*.xml            # raw BLAST XML per genome
"""

from __future__ import annotations

import argparse
import subprocess
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path
from typing import Dict, List

import pandas as pd
from Bio import SeqIO

try:
    from synteny_utils import normalize_scaffold
except Exception:  # fallback when running redesign in isolation
    def normalize_scaffold(name: str) -> str:
        return str(name).split()[0]


# Greedy clustering parameters (mirroring tblastn_recovery_test v4)
PARALOG_COVERAGE_THRESHOLD = 1.2   # >120% union coverage → distinct paralog
MIN_CLUSTER_COVERAGE = 0.5         # keep clusters with ≥50% coverage
MAX_GAP_KB = 20.0                  # Split if gap > 20kb (separates tandem copies)
QUERY_OVERLAP_THRESHOLD = 0.5      # Query overlap >50% of shorter HSP → tandem copy
MAX_CLUSTER_SPAN_KB = 500.0        # Hard split if cluster span > 500kb (prevents mega-gene blobs)

# Family-level heuristic limit on how large a merged same-query cluster
# is allowed to be. Estimated from BK gene lengths during calibration.
FAMILY_MAX_MERGE_SPAN_KB: float | None = None


def run_tblastn(
    query_file: Path,
    genome_db: Path,
    output_xml: Path,
    *,
    evalue: str,
    max_targets: int,
    threads: int,
) -> bool:
    """Run tBLASTn for a combined multi-query FASTA against one genome DB."""
    cmd = [
        "tblastn",
        "-query", str(query_file),
        "-db", str(genome_db),
        "-outfmt", "5",
        "-evalue", str(evalue),
        "-max_target_seqs", str(max_targets),
        "-num_threads", str(threads),
    ]
    with open(output_xml, "w") as out:
        res = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True)
    return res.returncode == 0


def parse_blast_xml_with_query_coords(xml_file: Path) -> List[Dict]:
    """Parse BLAST XML and extract full HSP information including query coords."""
    hits: List[Dict] = []

    try:
        tree = ET.parse(xml_file)
        root = tree.getroot()

        for iteration in root.findall(".//Iteration"):
            qdef = iteration.find("Iteration_query-def")
            qlen = iteration.find("Iteration_query-len")
            if qdef is None or qlen is None:
                continue

            query_def = qdef.text or ""
            query_id = query_def.split()[0]
            try:
                query_length = int(qlen.text)
            except Exception:
                query_length = 0

            for hit in iteration.findall(".//Hit"):
                hit_def_el = hit.find("Hit_def")
                hit_acc_el = hit.find("Hit_accession")
                if hit_acc_el is None:
                    continue
                hit_def = (hit_def_el.text if hit_def_el is not None else "") or ""
                hit_acc = hit_acc_el.text or ""

                for hsp in hit.findall(".//Hsp"):
                    try:
                        q_start = int(hsp.find("Hsp_query-from").text)
                        q_end = int(hsp.find("Hsp_query-to").text)
                        h_start = int(hsp.find("Hsp_hit-from").text)
                        h_end = int(hsp.find("Hsp_hit-to").text)
                        evalue = float(hsp.find("Hsp_evalue").text)
                        bitscore = float(hsp.find("Hsp_bit-score").text)
                    except Exception:
                        continue

                    if h_start <= h_end:
                        strand = "+"
                        g_start = h_start
                        g_end = h_end
                    else:
                        strand = "-"
                        g_start = h_end
                        g_end = h_start

                    if query_length > 0:
                        q_cov = (q_end - q_start + 1) / float(query_length)
                    else:
                        q_cov = 0.0

                    hits.append(
                        {
                            "query_id": query_id,
                            "query_length": query_length,
                            "query_start": q_start,
                            "query_end": q_end,
                            "query_coverage": q_cov,
                            "scaffold": normalize_scaffold(hit_acc),
                            "scaffold_raw": hit_acc,
                            "scaffold_desc": hit_def,
                            "strand": strand,
                            "genomic_start": g_start,
                            "genomic_end": g_end,
                            "evalue": evalue,
                            "bitscore": bitscore,
                        }
                    )
    except Exception as e:
        print(f"  WARNING: error parsing BLAST XML {xml_file}: {e}")

    return hits


def filter_hsps_by_evalue(hsps: List[Dict], max_evalue: float) -> List[Dict]:
    """Filter HSPs by E-value threshold."""
    return [h for h in hsps if h["evalue"] <= max_evalue]


def calculate_union_query_coverage(hsps: List[Dict]) -> float:
    """Calculate UNION of query regions covered by HSPs for one query."""
    if not hsps:
        return 0.0

    qlen = hsps[0].get("query_length", 0) or 0
    if qlen <= 0:
        return 0.0

    regions = [(h["query_start"], h["query_end"]) for h in hsps]
    regions.sort()

    merged: List[tuple[int, int]] = []
    for start, end in regions:
        if merged and start <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))

    total = sum(end - start + 1 for start, end in merged)
    return total / float(qlen)


def calculate_cumulative_query_coverage(hsps: List[Dict]) -> float:
    """
    Calculate cumulative coverage across potentially multiple queries.

    For tandem duplicates, HSPs from different queries may appear in the same
    cluster. We calculate the sum of each query's individual coverage:

    cumulative_cov = coverage(queryA) + coverage(queryB) + ...

    Example: If queryA covers 40% and queryB covers 60%, cumulative = 1.0
    This represents "1.0 query equivalents" worth of coverage.
    """
    if not hsps:
        return 0.0

    # Group HSPs by query_id and calculate coverage for each query
    query_coverages: Dict[str, float] = {}

    for query_id in set(h["query_id"] for h in hsps):
        query_hsps = [h for h in hsps if h["query_id"] == query_id]
        query_coverages[query_id] = calculate_union_query_coverage(query_hsps)

    # Return sum of all query coverages
    return sum(query_coverages.values())


def summarize_cluster_multi_query(hsps: List[Dict]) -> Dict:
    """
    Summarize a cluster that may contain HSPs from multiple queries (tandem duplicates).

    For multi-query clusters, we pick the query with the best bitscore as the
    representative, but report all contributing queries.
    """
    g_start = min(h["genomic_start"] for h in hsps)
    g_end = max(h["genomic_end"] for h in hsps)
    span_kb = (g_end - g_start) / 1000.0

    # Calculate cumulative coverage across all queries
    cumulative_cov = calculate_cumulative_query_coverage(hsps)

    # Pick the query with the best bitscore as representative
    best_hsp = max(hsps, key=lambda h: h["bitscore"])
    representative_query = best_hsp["query_id"]

    # Get query-specific stats for the representative
    rep_hsps = [h for h in hsps if h["query_id"] == representative_query]
    q_start_min = min(h["query_start"] for h in rep_hsps)
    q_end_max = max(h["query_end"] for h in rep_hsps)
    qlen = rep_hsps[0].get("query_length", 0) or 0
    if qlen > 0:
        q_span_cov = (q_end_max - q_start_min + 1) / float(qlen)
    else:
        q_span_cov = 0.0

    best_evalue = min(h["evalue"] for h in hsps)
    best_bitscore = max(h["bitscore"] for h in hsps)

    # List all contributing queries
    contributing_queries = sorted(set(h["query_id"] for h in hsps))

    return {
        "query_id": representative_query,
        "scaffold": hsps[0]["scaffold"],
        "strand": hsps[0]["strand"],
        "genomic_start": g_start,
        "genomic_end": g_end,
        "genomic_span_kb": span_kb,
        "query_start_min": q_start_min,
        "query_end_max": q_end_max,
        "query_span_coverage": q_span_cov,
        "combined_query_coverage": cumulative_cov,  # Now cumulative across queries!
        "num_hsps": len(hsps),
        "num_queries": len(contributing_queries),
        "contributing_queries": ",".join(contributing_queries),
        "best_evalue": best_evalue,
        "best_bitscore": best_bitscore,
    }


def summarize_cluster(hsps: List[Dict]) -> Dict:
    """Legacy function - redirect to multi-query version."""
    return summarize_cluster_multi_query(hsps)


def cluster_hsps_per_query(
    hsps: List[Dict],
    query_id: str,
    paralog_coverage_threshold: float = PARALOG_COVERAGE_THRESHOLD,
    min_cluster_coverage: float = MIN_CLUSTER_COVERAGE,
    query_overlap_threshold: float = QUERY_OVERLAP_THRESHOLD,
    max_cluster_span_kb: float = MAX_CLUSTER_SPAN_KB,
    max_gap_kb: float = MAX_GAP_KB,
) -> List[Dict]:
    """
    Cluster HSPs for a SINGLE query into gene candidates.

    This handles multi-exon genes for one query by merging HSPs with:
    - Gap < MAX_GAP_KB (exons of same gene)
    - Query overlap < query_overlap_threshold (minor overlaps allowed)
    - Total coverage for THIS query <= paralog_coverage_threshold

    Returns list of gene candidates from this query.
    """
    if not hsps:
        return []

    df = pd.DataFrame(hsps)
    out: List[Dict] = []

    # Group by (scaffold, strand) for this query
    for (scaffold, strand), group in df.groupby(["scaffold", "strand"]):
        sg_sorted = group.sort_values("genomic_start").reset_index(drop=True)
        used: set[int] = set()

        for start_idx in range(len(sg_sorted)):
            if start_idx in used:
                continue

            current = [sg_sorted.iloc[start_idx].to_dict()]
            used.add(start_idx)

            for next_idx in range(start_idx + 1, len(sg_sorted)):
                if next_idx in used:
                    continue

                nxt = sg_sorted.iloc[next_idx].to_dict()
                gap_bp = nxt["genomic_start"] - current[-1]["genomic_end"]
                gap_kb = gap_bp / 1000.0

                # Stop if gap too large (different genes in tandem array)
                if gap_kb > max_gap_kb:
                    break

                # Check if this HSP overlaps SUBSTANTIALLY in QUERY SPACE with any HSP already in cluster
                # Substantial overlap (>threshold) indicates different genomic copy (tandem duplicate)
                # Minor overlap (<threshold) is allowed (adjacent exons, alignment artifacts)
                nxt_q_start = nxt["query_start"]
                nxt_q_end = nxt["query_end"]
                nxt_q_len = nxt_q_end - nxt_q_start + 1
                has_substantial_query_overlap = False

                for existing_hsp in current:
                    ex_q_start = existing_hsp["query_start"]
                    ex_q_end = existing_hsp["query_end"]
                    ex_q_len = ex_q_end - ex_q_start + 1

                    # Check for overlap in query coordinates
                    overlap_start = max(nxt_q_start, ex_q_start)
                    overlap_end = min(nxt_q_end, ex_q_end)

                    if overlap_end >= overlap_start:
                        # Calculate overlap as percentage of shorter HSP
                        overlap_len = overlap_end - overlap_start + 1
                        min_hsp_len = min(nxt_q_len, ex_q_len)
                        overlap_pct = overlap_len / min_hsp_len if min_hsp_len > 0 else 0

                        if overlap_pct > query_overlap_threshold:
                            # Substantial overlap → different gene copies (tandem duplicates)
                            has_substantial_query_overlap = True
                            break
                        # else: minor overlap allowed → continue merging

                if has_substantial_query_overlap:
                    # This HSP substantially overlaps same query region as existing HSP
                    # → Different genomic copy → Start new gene
                    break

                # Safety valve: Check if adding this HSP would create unreasonably large cluster
                # This prevents "mega-gene blobs" from repetitive/unstructured proteins (e.g., Glutenin)
                # Only applies when other filters (gap, query overlap, coverage) fail to split
                cluster_start = min(h["genomic_start"] for h in current)
                cluster_end = max(h["genomic_end"] for h in current)
                potential_span = max(nxt["genomic_end"], cluster_end) - min(nxt["genomic_start"], cluster_start)
                potential_span_kb = potential_span / 1000.0

                if potential_span_kb > MAX_CLUSTER_SPAN_KB:
                    # Cluster would span >500kb - force split to prevent mega-genes
                    # This is a safety valve for pathological cases where normal filters fail
                    break

                # Check coverage for THIS query
                test_cluster = current + [nxt]
                test_cov = calculate_union_query_coverage(test_cluster)

                # Stop if we exceed threshold (different tandem copy)
                if test_cov > paralog_coverage_threshold:
                    break

                current.append(nxt)
                used.add(next_idx)

            # Keep if meets minimum coverage
            cov = calculate_union_query_coverage(current)
            if cov >= min_cluster_coverage:
                g_start = min(h["genomic_start"] for h in current)
                g_end = max(h["genomic_end"] for h in current)
                span_kb = (g_end - g_start) / 1000.0

                # Filter out abnormally long clusters that likely bridge multiple tandem copies
                # If coverage >95% but span >10kb, it's probably bridging tandem genes
                if cov > 0.95 and span_kb > 10.0:
                    # Skip this mega-cluster - it's a bridge across tandem array
                    continue

                out.append({
                    "query_id": query_id,
                    "scaffold": current[0]["scaffold"],
                    "strand": current[0]["strand"],
                    "genomic_start": g_start,
                    "genomic_end": g_end,
                    "query_coverage": cov,
                    "num_hsps": len(current),
                    "hsps": current,
                })

    return out


def deduplicate_genes_across_queries(gene_candidates: List[Dict]) -> List[Dict]:
    """
    Merge gene candidates from different queries that overlap genomically.

    This handles the case where 7-8 different queries all find the same gene copy.
    Any genomic overlap means they're the same gene - merge them.
    """
    if not gene_candidates:
        return []

    # Sort by position
    candidates = sorted(gene_candidates, key=lambda x: (x["scaffold"], x["strand"], x["genomic_start"]))

    merged = []
    current_group = [candidates[0]]

    for candidate in candidates[1:]:
        # Check if on same scaffold/strand as current group
        if (candidate["scaffold"] != current_group[0]["scaffold"] or
            candidate["strand"] != current_group[0]["strand"]):
            # Different location - finalize current group
            merged.append(current_group)
            current_group = [candidate]
            continue

        # Calculate the extent of the current merged group
        group_start = min(g["genomic_start"] for g in current_group)
        group_end = max(g["genomic_end"] for g in current_group)
        group_span = group_end - group_start

        # Calculate overlap with merged group extent
        overlap_start = max(group_start, candidate["genomic_start"])
        overlap_end = min(group_end, candidate["genomic_end"])
        overlap_bp = max(0, overlap_end - overlap_start)

        # Calculate gap from group extent
        gap_bp = candidate["genomic_start"] - group_end

        # Candidate span
        cand_span = candidate["genomic_end"] - candidate["genomic_start"]
        min_span = min(group_span, cand_span)

        # Merge if significant overlap (>50% of smaller span)
        # This means they're the same gene from different queries
        if min_span > 0 and overlap_bp > 0.5 * min_span:
            current_group.append(candidate)
        else:
            # No significant overlap - finalize current group, start new
            merged.append(current_group)
            current_group = [candidate]

    # Don't forget last group
    if current_group:
        merged.append(current_group)

    # Summarize each merged group
    final_genes = []
    for group in merged:
        all_hsps = []
        for gene in group:
            all_hsps.extend(gene["hsps"])

        final_genes.append(summarize_cluster_multi_query(all_hsps))

    return final_genes


def merge_proximal_gene_candidates_same_query(
    gene_candidates: List[Dict],
    merge_gap_kb: float,
) -> List[Dict]:
    """
    Post-processing merge: combine same-query candidates that are genomically proximate.

    This addresses cases where the per-query greedy clustering created multiple nearby
    fragments for the same query within a locus. We merge candidates when:
        - Same query_id
        - Same scaffold and strand
        - Gap between candidates <= merge_gap_kb

    Coverage and coordinates are recomputed from the underlying HSPs.
    """
    if not gene_candidates or merge_gap_kb <= 0:
        return gene_candidates

    max_merge_span_kb = FAMILY_MAX_MERGE_SPAN_KB

    grouped: Dict[tuple[str, str, str], List[Dict]] = defaultdict(list)
    for cand in gene_candidates:
        key = (cand["query_id"], cand["scaffold"], cand["strand"])
        grouped[key].append(cand)

    merged_all: List[Dict] = []

    for (query_id, scaffold, strand), cands in grouped.items():
        cands_sorted = sorted(cands, key=lambda c: c["genomic_start"])
        i = 0
        n = len(cands_sorted)

        while i < n:
            current_group = [cands_sorted[i]]
            current_end = cands_sorted[i]["genomic_end"]
            group_start = cands_sorted[i]["genomic_start"]
            j = i + 1

            while j < n:
                nxt = cands_sorted[j]
                gap_bp = nxt["genomic_start"] - current_end
                gap_kb = gap_bp / 1000.0

                if gap_kb > merge_gap_kb:
                    break

                # Check that merged span would not exceed family-level limit
                potential_start = min(group_start, nxt["genomic_start"])
                potential_end = max(current_end, nxt["genomic_end"])
                potential_span_kb = (potential_end - potential_start) / 1000.0
                if max_merge_span_kb is not None and potential_span_kb > max_merge_span_kb:
                    break

                current_group.append(nxt)
                current_end = max(current_end, nxt["genomic_end"])
                group_start = min(group_start, nxt["genomic_start"])
                j += 1

            # Build merged candidate from all HSPs in current_group
            all_hsps = []
            for cand in current_group:
                all_hsps.extend(cand["hsps"])

            if not all_hsps:
                i = j
                continue

            g_start = min(h["genomic_start"] for h in all_hsps)
            g_end = max(h["genomic_end"] for h in all_hsps)
            cov = calculate_union_query_coverage(all_hsps)

            merged_all.append(
                {
                    "query_id": query_id,
                    "scaffold": scaffold,
                    "strand": strand,
                    "genomic_start": g_start,
                    "genomic_end": g_end,
                    "query_coverage": cov,
                    "num_hsps": len(all_hsps),
                    "hsps": all_hsps,
                }
            )

            i = j

    return merged_all


def cluster_hsps_greedy(
    hsps: List[Dict],
    paralog_coverage_threshold: float = PARALOG_COVERAGE_THRESHOLD,
    min_cluster_coverage: float = MIN_CLUSTER_COVERAGE,
    query_overlap_threshold: float = QUERY_OVERLAP_THRESHOLD,
    max_gap_kb: float = MAX_GAP_KB,
    merge_gap_kb: float | None = None,
) -> List[Dict]:
    """
    Two-step clustering: per-query, then cross-query deduplication.

    Step 1: For each query, cluster its HSPs into gene candidates
            - Handles multi-exon genes within one query
            - Uses gap threshold and per-query coverage
            - Uses query overlap threshold to distinguish tandem copies

    Step 2: Merge gene candidates from different queries that overlap
            - Handles tandem duplicates (7-8 queries finding same gene)
            - Uses genomic overlap to deduplicate

    This cleanly separates two concerns and avoids cumulative coverage confusion.
    """
    if not hsps:
        return []

    # Step 1: Per-query clustering
    all_gene_candidates = []
    for query_id in set(h["query_id"] for h in hsps):
        query_hsps = [h for h in hsps if h["query_id"] == query_id]
        gene_candidates = cluster_hsps_per_query(
            query_hsps,
            query_id,
            paralog_coverage_threshold,
            min_cluster_coverage,
            query_overlap_threshold,
            MAX_CLUSTER_SPAN_KB,
            max_gap_kb,
        )
        all_gene_candidates.extend(gene_candidates)

    # Optional Step 1b: merge genomically proximate candidates from SAME query.
    if merge_gap_kb is not None and merge_gap_kb > 0:
        all_gene_candidates = merge_proximal_gene_candidates_same_query(
            all_gene_candidates,
            merge_gap_kb,
        )

    # Step 2: Cross-query deduplication
    final_genes = deduplicate_genes_across_queries(all_gene_candidates)

    return final_genes


def _normalize_genome_id(name: str) -> str:
    """Normalize genome ID to match Phase 2/3 conventions."""
    if name.startswith("GCA_") or name.startswith("GCF_"):
        parts = name.split("_")
        if len(parts) >= 2:
            return f"{parts[0]}_{parts[1]}"
    return name


def count_phase1_bk_targets(phase1_dir: Path) -> int:
    """Count BK targets from Phase 1 per-locus target files (ground truth)."""
    total = 0
    for locus_dir in phase1_dir.glob("BK_*"):
        if locus_dir.is_dir():
            for faa_file in locus_dir.glob("*_targets.faa"):
                with open(faa_file) as f:
                    total += sum(1 for line in f if line.startswith('>'))
    return total


def get_phase1_bk_genes(phase1_dir: Path) -> tuple[list[dict], int]:
    """
    Extract BK ground truth at GENE level from Phase 1 locus_definitions.tsv.

    Uses gene_chromosomes / gene_starts / gene_ends columns added by
    add_gene_coordinates.py to build one record per gene.

    Returns:
        (genes, total_gene_count) where each gene is:
        {
            'locus_id': str,
            'gene_id': str,
            'chromosome': str,
            'start': int,
            'end': int,
        }
    """
    locus_defs = phase1_dir / "locus_definitions.tsv"
    if not locus_defs.exists():
        return [], 0

    genes: list[dict] = []
    total_genes = 0

    with open(locus_defs) as f:
        header = next(f).strip().split('\t')

        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < len(header):
                continue

            row = dict(zip(header, parts))
            locus_id = row.get('locus_id', '')
            if not locus_id.startswith('BK_'):
                continue

            # Require per-gene coordinates – if missing, skip this locus for gene-level validation.
            if 'gene_starts' not in row or 'gene_ends' not in row:
                continue
            if not row['gene_starts'] or not row['gene_ends']:
                continue

            cluster_members = [
                g.strip()
                for g in row.get('cluster_members', '').split(',')
                if g.strip()
            ]
            gene_chroms = [
                c.strip()
                for c in row.get('gene_chromosomes', '').split(',')
                if c.strip()
            ]

            try:
                gene_starts = [int(x) for x in row['gene_starts'].split(',') if x]
                gene_ends = [int(x) for x in row['gene_ends'].split(',') if x]
            except ValueError:
                # Malformed coordinates – skip this locus.
                continue

            # Determine how many genes we can safely use from this row.
            n_coords = min(len(gene_starts), len(gene_ends))
            if cluster_members:
                n_coords = min(n_coords, len(cluster_members))
            if gene_chroms:
                n_coords = min(n_coords, len(gene_chroms))

            if n_coords == 0:
                continue

            chrom_fallback = row.get('chromosome', '')

            for i in range(n_coords):
                gene_id = (
                    cluster_members[i]
                    if i < len(cluster_members)
                    else f"{locus_id}_gene{i+1}"
                )
                chrom = (
                    gene_chroms[i]
                    if gene_chroms and i < len(gene_chroms)
                    else chrom_fallback
                )

                start = gene_starts[i]
                end = gene_ends[i]

                genes.append(
                    {
                        'locus_id': locus_id,
                        'gene_id': gene_id,
                        'chromosome': chrom,
                        'start': start,
                        'end': end,
                    }
                )
                total_genes += 1

    return genes, total_genes


def get_phase1_bk_loci(phase1_dir: Path) -> tuple[list[dict], int]:
    """
    Extract ground truth BK loci from Phase 1 locus_definitions.tsv.

    Returns:
        (loci, total_gene_count) where loci is list of dicts with:
        {locus_id, chromosome, start, end, expected_genes}
    """
    locus_defs = phase1_dir / "locus_definitions.tsv"
    if not locus_defs.exists():
        return [], 0

    loci = []
    total_genes = 0

    with open(locus_defs) as f:
        header = next(f).strip().split('\t')
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < len(header):
                continue

            row = dict(zip(header, parts))
            if row.get('locus_id', '').startswith('BK_'):
                expected_genes = int(row['cluster_size'])
                loci.append({
                    'locus_id': row['locus_id'],
                    'chromosome': row['chromosome'],
                    'start': int(row['locus_span_start']),
                    'end': int(row['locus_span_end']),
                    'expected_genes': expected_genes
                })
                total_genes += expected_genes

    return loci, total_genes


def validate_detections_vs_groundtruth(detected_clusters: list[dict], ground_truth_coords: list[dict]) -> dict:
    """
    Legacy locus-level validation used by older calibration flows.

    Kept for backwards compatibility; new factorial calibration uses
    validate_detections_vs_groundtruth_per_gene for gene-level scoring.
    """
    if not ground_truth_coords:
        return {
            'true_positives': 0,
            'false_negatives': 0,
            'false_positives': len(detected_clusters),
            'total_detected': len(detected_clusters),
            'total_ground_truth_genes': 0,
            'matched_loci': [],
            'missed_loci': [],
            'extra_genes': detected_clusters
        }

    # Helper to normalize chromosome names (strip .1 suffix)
    def normalize_chrom(chrom: str) -> str:
        return chrom.split('.')[0] if chrom else ''

    # Track which ground truth loci were found
    matched_gt_loci = set()
    matched_detections = []
    unmatched_detections = []

    for det in detected_clusters:
        det_chrom = normalize_chrom(det.get('scaffold', ''))
        det_start = det.get('genomic_start', 0)
        det_end = det.get('genomic_end', 0)

        # Check if this detection overlaps any ground truth locus
        matched = False
        for gt in ground_truth_coords:
            gt_chrom = normalize_chrom(gt['chromosome'])
            if det_chrom == gt_chrom:
                # Check for overlap (allowing some slop)
                if not (det_end < gt['start'] - 50000 or det_start > gt['end'] + 50000):
                    matched_gt_loci.add(gt['locus_id'])
                    matched_detections.append(det)
                    matched = True
                    break

        if not matched:
            unmatched_detections.append(det)

    # Calculate metrics
    total_gt_genes = sum(gt['expected_genes'] for gt in ground_truth_coords)
    all_gt_loci = set(gt['locus_id'] for gt in ground_truth_coords)
    missed_loci = all_gt_loci - matched_gt_loci

    # True positives: detected genes that overlap ground truth
    # False negatives: ground truth genes NOT found (may include partial loci)
    tp_count = len(matched_detections)
    fn_count = total_gt_genes - tp_count
    fp_count = len(unmatched_detections)

    return {
        'true_positives': tp_count,
        'false_negatives': fn_count,
        'false_positives': fp_count,
        'total_detected': len(detected_clusters),
        'total_ground_truth_genes': total_gt_genes,
        'matched_loci': list(matched_gt_loci),
        'missed_loci': list(missed_loci),
        'extra_genes': unmatched_detections
    }


def validate_detections_vs_groundtruth_per_gene(
    detected_clusters: list[dict],
    ground_truth_genes: list[dict],
) -> dict:
    """
    Compare detected gene coordinates to ground truth at GENE level.

    Each ground truth gene is counted at most once:
        - If ≥1 detection overlaps the gene → 1 true positive.
        - Additional detections overlapping the same gene are counted as false positives.

    Detections that do not overlap any ground truth gene are also false positives.

    Returns:
        {
            'true_positives': int,   # Ground truth genes with ≥1 detection
            'false_negatives': int,  # Ground truth genes not detected
            'false_positives': int,  # Detected genes not matching any ground truth
            'total_detected': int,
            'total_ground_truth_genes': int,
            'matched_loci': list,    # Loci containing ≥1 detected gene
            'missed_loci': list,     # Loci where all genes were missed
            'extra_genes': list,     # Detection dicts classified as false positives
            'matched_genes': list,   # Ground truth gene dicts with ≥1 detection
            'missed_genes': list,    # Ground truth gene dicts with 0 detections
        }
    """
    if not ground_truth_genes:
        return {
            'true_positives': 0,
            'false_negatives': 0,
            'false_positives': len(detected_clusters),
            'total_detected': len(detected_clusters),
            'total_ground_truth_genes': 0,
            'matched_loci': [],
            'missed_loci': [],
            'extra_genes': detected_clusters
        }

    # Helper to normalize chromosome names (strip .1 suffix)
    def normalize_chrom(chrom: str) -> str:
        return chrom.split('.')[0] if chrom else ''

    total_detected = len(detected_clusters)
    total_gt_genes = len(ground_truth_genes)

    # Map detections to the single best-overlapping ground truth gene (if any).
    det_to_gene_idx: list[int | None] = [None] * total_detected

    for i, det in enumerate(detected_clusters):
        det_chrom = normalize_chrom(det.get('scaffold', ''))
        det_start = det.get('genomic_start', 0)
        det_end = det.get('genomic_end', 0)

        best_idx: int | None = None
        best_overlap = 0

        for j, gt in enumerate(ground_truth_genes):
            gt_chrom = normalize_chrom(gt['chromosome'])
            if det_chrom != gt_chrom:
                continue

            gt_start = gt['start']
            gt_end = gt['end']

            overlap_start = max(det_start, gt_start)
            overlap_end = min(det_end, gt_end)
            if overlap_end < overlap_start:
                continue

            overlap_len = overlap_end - overlap_start + 1
            if overlap_len > best_overlap:
                best_overlap = overlap_len
                best_idx = j

        det_to_gene_idx[i] = best_idx

    # Count detections per gene.
    gene_hit_counts = [0] * total_gt_genes
    for idx in det_to_gene_idx:
        if idx is not None:
            gene_hit_counts[idx] += 1

    tp_genes = sum(1 for c in gene_hit_counts if c > 0)
    fn_genes = total_gt_genes - tp_genes

    # Any detection that is unassigned, or a second+ hit to the same gene, is a false positive.
    extra_detections: list[dict] = []
    first_hit_taken = [False] * total_gt_genes

    for det, idx in zip(detected_clusters, det_to_gene_idx):
        if idx is None:
            extra_detections.append(det)
        else:
            if not first_hit_taken[idx]:
                first_hit_taken[idx] = True
            else:
                extra_detections.append(det)

    fp_count = len(extra_detections)

    # Derive matched/missed loci and gene lists for reporting.
    matched_genes = [g for g, c in zip(ground_truth_genes, gene_hit_counts) if c > 0]
    missed_genes = [g for g, c in zip(ground_truth_genes, gene_hit_counts) if c == 0]

    matched_loci = sorted({g['locus_id'] for g in matched_genes})
    missed_loci = sorted({g['locus_id'] for g in missed_genes})

    return {
        'true_positives': tp_genes,
        'false_negatives': fn_genes,
        'false_positives': fp_count,
        'total_detected': total_detected,
        'total_ground_truth_genes': total_gt_genes,
        'matched_loci': matched_loci,
        'missed_loci': missed_loci,
        'extra_genes': extra_detections,
        'matched_genes': matched_genes,
        'missed_genes': missed_genes,
    }


def find_bk_blast_xml(blast_xml_dir: Path) -> Path | None:
    """Find BK BLAST XML file in the blast_xml directory."""
    for xml_file in blast_xml_dir.glob("*.xml"):
        if 'Belonocnema_kinseyi' in xml_file.name or 'BK' in xml_file.name:
            return xml_file
    return None


def calibrate_threshold(
    phase1_dir: Path,
    blast_xml_dir: Path,
    default_threshold: float = QUERY_OVERLAP_THRESHOLD,
    max_evalue: float = 1e-5,
) -> tuple[float, float, float]:
    """Auto-calibrate query overlap threshold, min coverage, AND E-value using BK ground truth.

    Tests multiple combinations of parameters on BK genome and selects the one that best
    matches Phase 1 BK gene count. If recovery is poor (<80%), retries with
    lower min_cluster_coverage to handle fragmentary HSPs (unstructured proteins).

    Returns:
        (calibrated_query_overlap_threshold, calibrated_min_coverage, calibrated_evalue)
    """
    print()
    print("=" * 80)
    print("AUTO-CALIBRATING CLUSTERING PARAMETERS")
    print("=" * 80)

    # Count Phase 1 BK targets (ground truth)
    phase1_bk_count = count_phase1_bk_targets(phase1_dir)

    if phase1_bk_count == 0:
        print("  No BK genes in Phase 1 - using default parameters")
        print(f"  → Defaults: query_overlap={default_threshold}, min_coverage={MIN_CLUSTER_COVERAGE}, evalue={max_evalue}")
        return (default_threshold, MIN_CLUSTER_COVERAGE, max_evalue, PARALOG_COVERAGE_THRESHOLD, MAX_GAP_KB)

    print(f"  Phase 1 BK ground truth: {phase1_bk_count} genes")

    # Find BK BLAST XML
    bk_xml = find_bk_blast_xml(blast_xml_dir)
    if not bk_xml:
        print("  No BK BLAST results found - using default parameters")
        print(f"  → Defaults: query_overlap={default_threshold}, min_coverage={MIN_CLUSTER_COVERAGE}, evalue={max_evalue}")
        return (default_threshold, MIN_CLUSTER_COVERAGE, max_evalue, PARALOG_COVERAGE_THRESHOLD, MAX_GAP_KB)

    print(f"  BK BLAST XML: {bk_xml.name}")
    print()

    # Parse BLAST XML once (contains all HSPs with their E-values)
    all_hsps = parse_blast_xml_with_query_coords(bk_xml)

    # STEP 1: Calibrate E-value threshold
    # Test multiple E-value cutoffs to find the one that best matches Phase 1
    evalue_thresholds = [1e-5, 1e-10, 1e-15, 1e-20, 1e-30]
    test_thresholds = [0.1, 0.2, 0.3, 0.5, 0.7]

    print("  STEP 1: Calibrating E-value threshold")
    print("  " + "-" * 76)

    evalue_results = []
    for evalue in evalue_thresholds:
        # Filter HSPs by this E-value
        hsps = filter_hsps_by_evalue(all_hsps, evalue)

        if not hsps:
            print(f"    E-value ≤ {evalue:.0e}: 0 HSPs (too strict)")
            continue

        # Test query overlap thresholds with this E-value
        best_for_evalue = None
        best_diff_for_evalue = float('inf')

        for threshold in test_thresholds:
            clusters = cluster_hsps_greedy(
                hsps,
                query_overlap_threshold=threshold,
                min_cluster_coverage=MIN_CLUSTER_COVERAGE,
            )

            bk_count = len(clusters)
            diff = abs(bk_count - phase1_bk_count)
            recovery_pct = (bk_count / phase1_bk_count) * 100 if phase1_bk_count > 0 else 0

            if diff < best_diff_for_evalue:
                best_diff_for_evalue = diff
                best_for_evalue = (threshold, MIN_CLUSTER_COVERAGE, evalue, bk_count, diff, recovery_pct)

        if best_for_evalue:
            evalue_results.append(best_for_evalue)
            threshold, min_cov, eval_val, count, diff, rec_pct = best_for_evalue
            print(f"    E-value ≤ {eval_val:.0e}: {len(hsps)} HSPs → best: {count} genes ({rec_pct:.0f}% recovery, ±{diff})")

    if not evalue_results:
        print("  ERROR: No valid E-value thresholds found")
        return (default_threshold, MIN_CLUSTER_COVERAGE, max_evalue)

    # Select best E-value (minimize difference from ground truth)
    best = min(evalue_results, key=lambda x: x[4])
    best_threshold, best_min_coverage, best_evalue, best_count, best_diff, best_recovery_pct = best

    print()
    print(f"  → Best E-value: {best_evalue:.0e} ({best_count} genes, {best_recovery_pct:.0f}% recovery)")
    print()

    # STEP 1.5: Calibrate paralog_coverage and gap thresholds
    print(f"  STEP 1.5: Calibrating deduplication & gap parameters (E-value={best_evalue:.0e})")
    print("  " + "-" * 76)

    hsps = filter_hsps_by_evalue(all_hsps, best_evalue)

    # Test paralog_coverage_threshold (affects deduplication)
    paralog_thresholds = [0.1, 0.2, 0.3, 0.5, 0.7]
    # Test max_gap_kb (affects HSP merging)
    gap_thresholds = [10, 20, 50, 100]

    param_results = []
    for paralog_thresh in paralog_thresholds:
        for gap_thresh in gap_thresholds:
            # Test with default query_overlap
            clusters = cluster_hsps_greedy(
                hsps,
                paralog_coverage_threshold=paralog_thresh,
                query_overlap_threshold=QUERY_OVERLAP_THRESHOLD,
                min_cluster_coverage=MIN_CLUSTER_COVERAGE,
                max_gap_kb=gap_thresh,
            )

            bk_count = len(clusters)
            diff = abs(bk_count - phase1_bk_count)
            recovery_pct = (bk_count / phase1_bk_count) * 100 if phase1_bk_count > 0 else 0

            param_results.append((paralog_thresh, gap_thresh, bk_count, diff, recovery_pct))

    # Select best parameter combination
    best_params = min(param_results, key=lambda x: x[3])
    best_paralog_thresh, best_gap_thresh, best_param_count, best_param_diff, best_param_recovery = best_params

    print(f"    Tested {len(param_results)} combinations")
    print(f"    Best: paralog_coverage={best_paralog_thresh:.1f}, gap={best_gap_thresh:.0f}kb")
    print(f"          → {best_param_count} genes ({best_param_recovery:.0f}% recovery, ±{best_param_diff})")
    print()

    # STEP 2: Fine-tune query overlap threshold with all calibrated parameters
    print(f"  STEP 2: Fine-tuning query overlap")
    print(f"         (E-value={best_evalue:.0e}, paralog={best_paralog_thresh:.1f}, gap={best_gap_thresh:.0f}kb)")
    print("  " + "-" * 76)

    results = []

    for threshold in test_thresholds:
        clusters = cluster_hsps_greedy(
            hsps,
            paralog_coverage_threshold=best_paralog_thresh,
            query_overlap_threshold=threshold,
            min_cluster_coverage=MIN_CLUSTER_COVERAGE,
            max_gap_kb=best_gap_thresh,
        )

        bk_count = len(clusters)
        diff = abs(bk_count - phase1_bk_count)
        recovery_pct = (bk_count / phase1_bk_count) * 100 if phase1_bk_count > 0 else 0
        results.append((threshold, MIN_CLUSTER_COVERAGE, best_evalue, best_paralog_thresh, best_gap_thresh, bk_count, diff, recovery_pct))

        status = "✓" if diff == 0 else f"±{diff}"
        print(f"    Threshold {threshold:.1f}: {bk_count} genes ({recovery_pct:.0f}% recovery, {status})")

    # Select best parameters
    best = min(results, key=lambda x: x[6])
    best_threshold, best_min_coverage, best_evalue, best_paralog_thresh, best_gap_thresh, best_count, best_diff, best_recovery_pct = best

    # STAGED CALIBRATION: Progressively relax min_coverage only when needed
    # Stage 1 (above): min_coverage=0.5 (normal, conservative)
    # Stage 2 (below): min_coverage=0.3 if recovery <80% (fragmentary HSPs)
    # Stage 3 (below): min_coverage=0.1 if still <80% (extreme cases like Glutenin)

    # STEP 3: If recovery is poor (<80%), retry with lower min_cluster_coverage
    if best_recovery_pct < 80:
        print()
        print(f"  STEP 3: Low recovery ({best_recovery_pct:.0f}%) - testing min_coverage=0.3")
        print(f"          (handles fragmentary HSPs from unstructured proteins)")
        print("  " + "-" * 76)

        stage2_results = []
        for threshold in test_thresholds:
            clusters = cluster_hsps_greedy(
                hsps,  # Already filtered by best_evalue
                paralog_coverage_threshold=best_paralog_thresh,
                query_overlap_threshold=threshold,
                min_cluster_coverage=0.3,
                max_gap_kb=best_gap_thresh,
            )

            bk_count = len(clusters)
            diff = abs(bk_count - phase1_bk_count)
            recovery_pct = (bk_count / phase1_bk_count) * 100 if phase1_bk_count > 0 else 0
            stage2_results.append((threshold, 0.3, best_evalue, best_paralog_thresh, best_gap_thresh, bk_count, diff, recovery_pct))

            status = "✓" if diff == 0 else f"±{diff}"
            print(f"    Threshold {threshold:.1f}: {bk_count} genes ({recovery_pct:.0f}% recovery, {status})")

        # Select best from stage 2
        best_stage2 = min(stage2_results, key=lambda x: x[6])

        # Use stage 2 if it improves recovery
        if best_stage2[6] < best[6]:
            best = best_stage2
            print()
            print(f"  ✓ min_coverage=0.3 improved recovery!")

        # STEP 4: If STILL <80%, try even lower min_coverage=0.1 (extreme cases)
        stage2_recovery = best[7]
        if stage2_recovery < 80:
            print()
            print(f"  STEP 4: Still low recovery ({stage2_recovery:.0f}%) - testing min_coverage=0.1")
            print(f"          (extreme fallback for highly fragmentary proteins)")
            print("  " + "-" * 76)

            stage3_results = []
            for threshold in test_thresholds:
                clusters = cluster_hsps_greedy(
                    hsps,  # Already filtered by best_evalue
                    paralog_coverage_threshold=best_paralog_thresh,
                    query_overlap_threshold=threshold,
                    min_cluster_coverage=0.1,
                    max_gap_kb=best_gap_thresh,
                )

                bk_count = len(clusters)
                diff = abs(bk_count - phase1_bk_count)
                recovery_pct = (bk_count / phase1_bk_count) * 100 if phase1_bk_count > 0 else 0
                stage3_results.append((threshold, 0.1, best_evalue, best_paralog_thresh, best_gap_thresh, bk_count, diff, recovery_pct))

                status = "✓" if diff == 0 else f"±{diff}"
                print(f"    Threshold {threshold:.1f}: {bk_count} genes ({recovery_pct:.0f}% recovery, {status})")

            # Select best from stage 3
            best_stage3 = min(stage3_results, key=lambda x: x[6])

            # Use stage 3 if it improves recovery (AND doesn't oversplit by >20%)
            # Oversplit check prevents false positives from excessive fragmentation
            if best_stage3[6] < best[6] and best_stage3[7] <= 120:
                best = best_stage3
                print()
                print(f"  ✓ min_coverage=0.1 improved recovery!")
            else:
                if best_stage3[7] > 120:
                    print()
                    print(f"  ⚠️  min_coverage=0.1 oversplits ({best_stage3[7]:.0f}% recovery) - keeping min_coverage={best[1]:.1f}")

    calibrated_threshold = best[0]
    calibrated_min_coverage = best[1]
    calibrated_evalue = best[2]
    calibrated_paralog_thresh = best[3]
    calibrated_gap_thresh = best[4]

    print()
    print("=" * 80)
    print("CALIBRATION COMPLETE")
    print("=" * 80)
    print(f"  query_overlap_threshold: {calibrated_threshold:.1f}")
    print(f"  min_cluster_coverage: {calibrated_min_coverage:.1f}")
    print(f"  evalue_threshold: {calibrated_evalue:.0e}")
    print(f"  paralog_coverage_threshold: {calibrated_paralog_thresh:.1f}")
    print(f"  max_gap_kb: {calibrated_gap_thresh:.0f}")
    print(f"  BK recovery: {best[5]}/{phase1_bk_count} genes ({best[7]:.0f}%)")
    print("=" * 80)
    print()

    return (calibrated_threshold, calibrated_min_coverage, calibrated_evalue, calibrated_paralog_thresh, calibrated_gap_thresh)


def calibrate_threshold_factorial(
    phase1_dir: Path,
    blast_xml_dir: Path,
    default_threshold: float = QUERY_OVERLAP_THRESHOLD,
    max_evalue: float = 1e-5,
) -> tuple[float, float, float, float, float, float]:
    """
    FACTORIAL calibration: Test ALL parameter combinations to find global optimum.

    Returns: (query_overlap, min_coverage, evalue, paralog_coverage, max_gap_kb)
    """
    print()
    print("=" * 80)
    print("FACTORIAL PARAMETER CALIBRATION (GLOBAL OPTIMUM)")
    print("=" * 80)

    # Get Phase 1 BK ground truth from locus_definitions.tsv
    ground_truth_genes, total_gt_genes = get_phase1_bk_genes(phase1_dir)

    if total_gt_genes == 0:
        print("  No BK genes with per-gene coordinates in Phase 1 - using defaults")
        return (
            default_threshold,
            MIN_CLUSTER_COVERAGE,
            max_evalue,
            PARALOG_COVERAGE_THRESHOLD,
            MAX_GAP_KB,
            0.0,
        )

    # Estimate a family-level cap on merged same-query cluster spans from BK gene lengths.
    spans_kb = [
        (g['end'] - g['start']) / 1000.0
        for g in ground_truth_genes
        if g['end'] > g['start']
    ]
    if spans_kb:
        spans_kb_sorted = sorted(spans_kb)
        median_span = spans_kb_sorted[len(spans_kb_sorted) // 2]
        max_span = spans_kb_sorted[-1]
        # Allow merges up to 1.5× the longest known gene or 3× the median, whichever is larger.
        # This permits large-intron genes (like RL) while guarding against fusing entire tandem arrays.
        global FAMILY_MAX_MERGE_SPAN_KB
        FAMILY_MAX_MERGE_SPAN_KB = max(max_span * 1.5, median_span * 3.0)
    else:
        FAMILY_MAX_MERGE_SPAN_KB = None

    num_loci = len({g['locus_id'] for g in ground_truth_genes})
    print(f"  Phase 1 BK ground truth: {total_gt_genes} genes across {num_loci} loci")

    # Find BK BLAST XML
    bk_xml = find_bk_blast_xml(blast_xml_dir)
    if not bk_xml:
        print("  No BK BLAST results found - using defaults")
        return (
            default_threshold,
            MIN_CLUSTER_COVERAGE,
            max_evalue,
            PARALOG_COVERAGE_THRESHOLD,
            MAX_GAP_KB,
            0.0,
        )

    print(f"  BK BLAST XML: {bk_xml.name}")
    print()

    # Parse all HSPs once
    all_hsps = parse_blast_xml_with_query_coords(bk_xml)

    # Define parameter ranges
    evalue_vals = [1e-5, 1e-10, 1e-15, 1e-20, 1e-30]
    paralog_vals = [0.1, 0.3, 0.5, 0.7]
    gap_vals = [10, 20, 50, 100, 150, 200]
    # Same-query merge distance (0 = disabled)
    merge_gap_vals = [0.0, 10.0, 30.0, 60.0]
    query_overlap_vals = [0.1, 0.2, 0.5]
    min_coverage_vals = [0.5, 0.3, 0.1]

    total_combos = (
        len(evalue_vals)
        * len(paralog_vals)
        * len(gap_vals)
        * len(merge_gap_vals)
        * len(query_overlap_vals)
        * len(min_coverage_vals)
    )

    print(f"  Testing {total_combos} parameter combinations:")
    print(f"    E-values: {len(evalue_vals)}, Paralog: {len(paralog_vals)}, Gap: {len(gap_vals)}")
    print(f"    Merge gap: {len(merge_gap_vals)}, Query overlap: {len(query_overlap_vals)}, Min coverage: {len(min_coverage_vals)}")
    print()

    best_params = None
    best_score = float("inf")
    best_validation = None
    best_recall = 0.0
    best_fp_ratio = 0.0

    tested = 0
    for evalue in evalue_vals:
        # Filter HSPs once per E-value
        hsps = filter_hsps_by_evalue(all_hsps, evalue)
        if not hsps:
            continue

        for paralog_thresh in paralog_vals:
            for gap_thresh in gap_vals:
                for merge_gap in merge_gap_vals:
                    for query_thresh in query_overlap_vals:
                        for min_cov in min_coverage_vals:
                            tested += 1

                            # Test this combination
                            clusters = cluster_hsps_greedy(
                                hsps,
                                paralog_coverage_threshold=paralog_thresh,
                                min_cluster_coverage=min_cov,
                                query_overlap_threshold=query_thresh,
                                max_gap_kb=gap_thresh,
                                merge_gap_kb=merge_gap if merge_gap > 0 else None,
                            )

                            # Validate against ground truth gene coordinates
                            validation = validate_detections_vs_groundtruth_per_gene(
                                clusters,
                                ground_truth_genes,
                            )
                            tp = validation["true_positives"]
                            fp = validation["false_positives"]

                            # Recall and FP ratio relative to BK ground truth.
                            if total_gt_genes > 0:
                                recall = tp / total_gt_genes
                                fp_ratio = fp / total_gt_genes
                            else:
                                recall = 0.0
                                fp_ratio = 0.0

                            # FP-aware objective (mirrors prototype_high_fp_calibration.py):
                            # - Strongly penalize recall < 90%
                            # - Within the acceptable window, trade off FP ratio against distance from 100% recall.
                            if total_gt_genes == 0:
                                score = float("inf")
                            elif recall < 0.9:
                                score = 1e6 + (0.9 - recall) * 1e3
                            else:
                                score = fp_ratio * 100.0 + abs(1.0 - recall) * 10.0

                            if score < best_score:
                                best_score = score
                                best_validation = validation
                                best_recall = recall
                                best_fp_ratio = fp_ratio
                                best_params = (
                                    query_thresh,
                                    min_cov,
                                    evalue,
                                    paralog_thresh,
                                    gap_thresh,
                                    merge_gap,
                                )

    if best_params is None:
        print("  ERROR: No valid parameter combinations found")
        return (
            default_threshold,
            MIN_CLUSTER_COVERAGE,
            max_evalue,
            PARALOG_COVERAGE_THRESHOLD,
            MAX_GAP_KB,
            0.0,
        )

    query_thresh, min_cov, evalue, paralog_thresh, gap_thresh, merge_gap = best_params

    print(f"  Tested {tested} combinations")
    print()
    print("=" * 80)
    print("CALIBRATION COMPLETE (GLOBAL OPTIMUM)")
    print("=" * 80)
    print(f"  query_overlap_threshold: {query_thresh:.1f}")
    print(f"  min_cluster_coverage: {min_cov:.1f}")
    print(f"  evalue_threshold: {evalue:.0e}")
    print(f"  paralog_coverage_threshold: {paralog_thresh:.1f}")
    print(f"  max_gap_kb: {gap_thresh:.0f}")
    print(f"  merge_gap_kb: {merge_gap:.0f}")
    print()
    print("  VALIDATION METRICS:")
    print(f"    Ground truth: {total_gt_genes} genes across {num_loci} loci")
    print(f"    True positives: {best_validation['true_positives']} genes found")
    print(f"    False negatives: {best_validation['false_negatives']} genes missed")
    print(f"    False positives: {best_validation['false_positives']} extra genes")
    if total_gt_genes > 0:
        print(f"    FP ratio (FP/GT): {best_fp_ratio:.3f}")
    print(f"    Recall: {best_recall*100:.1f}% of ground truth genes detected")
    if best_validation['missed_loci']:
        print(f"    Missed loci: {', '.join(best_validation['missed_loci'])}")
    print("=" * 80)
    print()

    return best_params


def calibrate_threshold_factorial_per_locus(
    phase1_dir: Path,
    blast_xml_dir: Path,
    query_to_locus: Dict[str, str],
    default_threshold: float = QUERY_OVERLAP_THRESHOLD,
    max_evalue: float = 1e-5,
) -> tuple[tuple[float, float, float, float, float, float], Dict[str, tuple[float, float, float, float, float, float]]]:
    """
    FACTORIAL calibration per BK locus.

    Returns:
        (global_best_params, params_by_locus)

    Each params tuple is:
        (query_overlap, min_coverage, evalue, paralog_coverage, max_gap_kb, merge_gap_kb)
    """
    # Get Phase 1 BK ground truth genes
    ground_truth_genes, total_gt_genes = get_phase1_bk_genes(phase1_dir)

    if total_gt_genes == 0:
        return (
            (
                default_threshold,
                MIN_CLUSTER_COVERAGE,
                max_evalue,
                PARALOG_COVERAGE_THRESHOLD,
                MAX_GAP_KB,
                0.0,
            ),
            {},
        )

    # FAMILY_MAX_MERGE_SPAN_KB is a family-level parameter: set from all BK genes.
    spans_kb = [
        (g["end"] - g["start"]) / 1000.0
        for g in ground_truth_genes
        if g["end"] > g["start"]
    ]
    if spans_kb:
        spans_kb_sorted = sorted(spans_kb)
        median_span = spans_kb_sorted[len(spans_kb_sorted) // 2]
        max_span = spans_kb_sorted[-1]
        global FAMILY_MAX_MERGE_SPAN_KB
        FAMILY_MAX_MERGE_SPAN_KB = max(max_span * 1.5, median_span * 3.0)
    else:
        FAMILY_MAX_MERGE_SPAN_KB = None

    # Find BK BLAST XML
    bk_xml = find_bk_blast_xml(blast_xml_dir)
    if not bk_xml:
        return (
            (
                default_threshold,
                MIN_CLUSTER_COVERAGE,
                max_evalue,
                PARALOG_COVERAGE_THRESHOLD,
                MAX_GAP_KB,
                0.0,
            ),
            {},
        )

    # Parse all HSPs once
    all_hsps = parse_blast_xml_with_query_coords(bk_xml)

    # Map BK genes by locus
    genes_by_locus: Dict[str, List[dict]] = {}
    for g in ground_truth_genes:
        lid = g["locus_id"]
        genes_by_locus.setdefault(lid, []).append(g)

    # Map HSPs to locus via query_to_locus
    hsps_by_locus: Dict[str, List[Dict]] = {}
    for h in all_hsps:
        qid = h.get("query_id", "")
        base_qid = qid.split("|")[0]
        locus = query_to_locus.get(qid) or query_to_locus.get(base_qid)
        if not locus:
            continue
        if not locus.startswith("BK_"):
            continue
        hsps_by_locus.setdefault(locus, []).append(h)

    # Inner factorial calibration for a single (genes, hsps) subset
    def _calibrate_for_subset(
        genes_subset: List[dict],
        hsps_subset: List[Dict],
    ) -> tuple[float, float, float, float, float, float]:
        total_gt = len(genes_subset)
        if total_gt == 0 or not hsps_subset:
            return (
                default_threshold,
                MIN_CLUSTER_COVERAGE,
                max_evalue,
                PARALOG_COVERAGE_THRESHOLD,
                MAX_GAP_KB,
                0.0,
            )

        evalue_vals = [1e-5, 1e-10, 1e-15, 1e-20, 1e-30]
        paralog_vals = [0.1, 0.3, 0.5, 0.7]
        gap_vals = [10, 20, 50, 100, 150, 200]
        merge_gap_vals = [0.0, 10.0, 30.0, 60.0]
        query_overlap_vals = [0.1, 0.2, 0.5]
        min_coverage_vals = [0.5, 0.3, 0.1]

        best_params_local = None
        best_score_local = float("inf")

        for evalue in evalue_vals:
            loc_hsps = filter_hsps_by_evalue(hsps_subset, evalue)
            if not loc_hsps:
                continue

            for paralog_thresh in paralog_vals:
                for gap_thresh in gap_vals:
                    for merge_gap in merge_gap_vals:
                        for query_thresh in query_overlap_vals:
                            for min_cov in min_coverage_vals:
                                clusters = cluster_hsps_greedy(
                                    loc_hsps,
                                    paralog_coverage_threshold=paralog_thresh,
                                    min_cluster_coverage=min_cov,
                                    query_overlap_threshold=query_thresh,
                                    max_gap_kb=gap_thresh,
                                    merge_gap_kb=merge_gap if merge_gap > 0 else None,
                                )

                                validation = validate_detections_vs_groundtruth_per_gene(
                                    clusters,
                                    genes_subset,
                                )
                                tp = validation["true_positives"]
                                fp = validation["false_positives"]

                                if total_gt > 0:
                                    recall = tp / total_gt
                                    fp_ratio = fp / total_gt
                                else:
                                    recall = 0.0
                                    fp_ratio = 0.0

                                if total_gt == 0:
                                    score = float("inf")
                                elif recall < 0.9:
                                    score = 1e6 + (0.9 - recall) * 1e3
                                else:
                                    score = fp_ratio * 100.0 + abs(1.0 - recall) * 10.0

                                if score < best_score_local:
                                    best_score_local = score
                                    best_params_local = (
                                        query_thresh,
                                        min_cov,
                                        evalue,
                                        paralog_thresh,
                                        gap_thresh,
                                        merge_gap,
                                    )

        if best_params_local is None:
            return (
                default_threshold,
                MIN_CLUSTER_COVERAGE,
                max_evalue,
                PARALOG_COVERAGE_THRESHOLD,
                MAX_GAP_KB,
                0.0,
            )
        return best_params_local

    # Global best (fallback/default): use all BK genes and all BK HSPs
    global_best = _calibrate_for_subset(ground_truth_genes, all_hsps)

    # Per-locus parameters
    params_by_locus: Dict[str, tuple[float, float, float, float, float, float]] = {}
    for locus_id, genes_subset in genes_by_locus.items():
        hsps_subset = hsps_by_locus.get(locus_id)
        if not hsps_subset:
            continue
        params_by_locus[locus_id] = _calibrate_for_subset(genes_subset, hsps_subset)

    return global_best, params_by_locus


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Phase 4 (redesign): detect targets with coverage-driven clustering"
    )
    p.add_argument(
        "--locus-defs",
        type=Path,
        required=True,
        help="Phase1_v2 locus_definitions.tsv",
    )
    p.add_argument(
        "--phase1-dir",
        type=Path,
        required=True,
        help="Phase1_v2 directory containing <locus>/*_targets.faa",
    )
    p.add_argument(
        "--genome-db-dir",
        type=Path,
        required=True,
        help="Directory with *.nhr BLAST DBs",
    )
    p.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Phase4_v2 output directory",
    )
    p.add_argument("--evalue", type=str, default="1e-15")
    p.add_argument("--max-targets", type=int, default=10000)
    p.add_argument("--threads", type=int, default=16)
    p.add_argument(
        "--query-overlap-threshold",
        type=float,
        default=None,
        help="Query overlap threshold (0-1). Overlap >threshold → tandem copy. If not specified, auto-calibrates using BK ground truth.",
    )
    p.add_argument(
        "--disable-calibration",
        action="store_true",
        help="Disable auto-calibration and use fixed threshold (default: auto-calibrate)",
    )
    # Retained for CLI compatibility; algorithm now uses coverage-driven rules.
    p.add_argument(
        "--base-min-split-kb",
        type=int,
        default=40,
        help="(unused) retained for backward-compatible CLI",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    blast_dir = args.output_dir / "blast_xml"
    blast_dir.mkdir(exist_ok=True)

    print("=" * 80)
    print("PHASE 4 (REDESIGN): TARGET DETECTION (COVERAGE-DRIVEN)")
    print("=" * 80)
    print(f"Locus defs: {args.locus_defs}")
    print(f"Phase1 dir: {args.phase1_dir}")
    print(f"Genome DBs: {args.genome_db_dir}")
    print(f"Output dir: {args.output_dir}")
    print()

    loci_df = pd.read_csv(args.locus_defs, sep="\t")

    # Optional gene_family annotation (same for all loci in a family run)
    gene_family = None
    if "gene_family" in loci_df.columns and not loci_df.empty:
        try:
            gene_family = str(loci_df["gene_family"].iloc[0])
        except Exception:
            gene_family = None

    # Build combined targets FASTA and mapping query_id -> parent locus_id
    combined = args.output_dir / "combined_targets.faa"
    query_to_locus: Dict[str, str] = {}
    records: List[SeqIO.SeqRecord] = []

    for _, row in loci_df.iterrows():
        locus_id = str(row["locus_id"])
        target_fa = args.phase1_dir / locus_id / f"{locus_id}_targets.faa"
        if not target_fa.exists():
            continue
        for rec in SeqIO.parse(str(target_fa), "fasta"):
            qid = rec.id.split()[0]
            query_to_locus[qid] = locus_id
            records.append(rec)

    if not records:
        print("ERROR: no target proteins found in phase1_v2 locus directories")
        return

    SeqIO.write(records, str(combined), "fasta")
    print(f"[1] Combined {len(records)} target proteins into {combined.name}")
    print(f"    Unique queries: {len({r.id.split()[0] for r in records})}")
    print()

    # Discover genome databases
    genome_dbs: Dict[str, Dict[str, object]] = {}
    for db_file in sorted(args.genome_db_dir.glob("*.nhr")):
        genome_name = db_file.stem
        genome_dbs[genome_name] = {
            "db": args.genome_db_dir / genome_name,
            "genome_id": _normalize_genome_id(genome_name),
        }

    if not genome_dbs:
        print(f"ERROR: no *.nhr BLAST DBs found in {args.genome_db_dir}")
        return

    print(f"[2] Found {len(genome_dbs)} genome databases")
    print()

    # Defaults for clustering parameters (may be overridden by calibration)
    final_threshold = QUERY_OVERLAP_THRESHOLD
    final_min_coverage = MIN_CLUSTER_COVERAGE
    final_evalue = float(args.evalue)
    final_paralog_thresh = PARALOG_COVERAGE_THRESHOLD
    final_gap_thresh = MAX_GAP_KB
    final_merge_gap = 0.0
    per_locus_params: Dict[str, tuple[float, float, float, float, float, float]] = {}

    # Determine query overlap threshold and min_coverage (calibrate or use fixed)
    if args.query_overlap_threshold is not None:
        # User specified threshold explicitly
        final_threshold = args.query_overlap_threshold
        final_min_coverage = MIN_CLUSTER_COVERAGE
        print(f"Using user-specified parameters: query_overlap={final_threshold}, min_coverage={final_min_coverage}")
    elif args.disable_calibration:
        # Calibration disabled, use defaults
        final_threshold = QUERY_OVERLAP_THRESHOLD
        final_min_coverage = MIN_CLUSTER_COVERAGE
        print(f"Auto-calibration disabled, using defaults: query_overlap={final_threshold}, min_coverage={final_min_coverage}")
    else:
        # Auto-calibrate using BK ground truth
        # First, run BLAST for BK genome if needed
        bk_genome = None
        for genome_name in genome_dbs:
            if 'Belonocnema_kinseyi' in genome_name or 'BK' in genome_name:
                bk_genome = genome_name
                break

        if bk_genome:
            bk_db = genome_dbs[bk_genome]["db"]
            bk_xml = blast_dir / f"{bk_genome}.xml"

            if not bk_xml.exists():
                print(f"[2.5] Running tBLASTn for BK genome (for calibration)...", end="")
                ok = run_tblastn(
                    combined,
                    bk_db,  # type: ignore[arg-type]
                    bk_xml,
                    evalue=args.evalue,
                    max_targets=args.max_targets,
                    threads=args.threads,
                )
                print(" done" if ok else " FAILED")

        # Calibrate clustering parameters using FP-aware factorial search.
        # This returns a global default AND per-locus parameters for BK loci.
        global_best, per_locus_params = calibrate_threshold_factorial_per_locus(
            args.phase1_dir,
            blast_dir,
            query_to_locus,
            default_threshold=QUERY_OVERLAP_THRESHOLD,
            max_evalue=float(args.evalue),
        )
        (
            final_threshold,
            final_min_coverage,
            final_evalue,
            final_paralog_thresh,
            final_gap_thresh,
            final_merge_gap,
        ) = global_best

    all_rows: List[Dict] = []

    # Run tBLASTn per genome and cluster HSPs into gene-level hits
    for genome_name, meta in genome_dbs.items():
        genome_db = meta["db"]
        genome_id = str(meta["genome_id"])
        xml_path = blast_dir / f"{genome_name}.xml"

        if xml_path.exists():
            print(f"[3] {genome_name}: using existing BLAST XML")
        else:
            print(f"[3] {genome_name}: tBLASTn combined targets...", end="")
            ok = run_tblastn(
                combined,
                genome_db,  # type: ignore[arg-type]
                xml_path,
                evalue=args.evalue,
                max_targets=args.max_targets,
                threads=args.threads,
            )
            print(" done" if ok else " FAILED")
            if not ok:
                continue

        hsps_all = parse_blast_xml_with_query_coords(xml_path)
        if not hsps_all:
            print(f"    {genome_name}: no HSPs found")
            continue

        # Partition HSPs by parent locus (via query_to_locus) and cluster per locus.
        locus_to_hsps: Dict[str, List[Dict]] = {}
        for h in hsps_all:
            qid = h.get("query_id", "")
            base_qid = qid.split("|")[0]
            parent_locus = query_to_locus.get(qid) or query_to_locus.get(base_qid)
            if not parent_locus:
                continue
            locus_to_hsps.setdefault(parent_locus, []).append(h)

        clusters: List[Dict] = []
        for parent_locus, hsps in locus_to_hsps.items():
            (
                locus_thresh,
                locus_min_cov,
                locus_evalue,
                locus_paralog,
                locus_gap,
                locus_merge_gap,
            ) = per_locus_params.get(
                parent_locus,
                (
                    final_threshold,
                    final_min_coverage,
                    final_evalue,
                    final_paralog_thresh,
                    final_gap_thresh,
                    final_merge_gap,
                ),
            )

            hsps_filt = filter_hsps_by_evalue(hsps, locus_evalue)
            if not hsps_filt:
                continue

            loc_clusters = cluster_hsps_greedy(
                hsps_filt,
                paralog_coverage_threshold=locus_paralog,
                min_cluster_coverage=locus_min_cov,
                query_overlap_threshold=locus_thresh,
                max_gap_kb=locus_gap,
                merge_gap_kb=locus_merge_gap if locus_merge_gap > 0 else None,
            )
            clusters.extend(loc_clusters)

        if not clusters:
            print(f"    {genome_name}: no clusters above coverage threshold")
            continue

        print(
            f"    {genome_name}: {len(hsps_all)} HSPs → {len(clusters)} gene-level hits"
        )

        # Convert clusters to rows, attaching parent_locus via query_to_locus
        for c in clusters:
            qid = c["query_id"]
            parent_locus = query_to_locus.get(qid) or query_to_locus.get(
                qid.split("|")[0], ""
            )
            if not parent_locus:
                continue

            row = {
                "locus_id": parent_locus,
                "parent_locus": parent_locus,
                "genome": genome_id,
                "scaffold": c["scaffold"],
                "strand": c["strand"],
                "start": c["genomic_start"],
                "end": c["genomic_end"],
                "span_kb": c["genomic_span_kb"],
                "num_hsps": c["num_hsps"],
                "best_evalue": c["best_evalue"],
                "best_bitscore": c["best_bitscore"],
                "query_id": qid,
                "query_start_min": c["query_start_min"],
                "query_end_max": c["query_end_max"],
                "query_span_coverage": c["query_span_coverage"],
                "combined_query_coverage": c["combined_query_coverage"],
            }
            if gene_family is not None:
                row["gene_family"] = gene_family
            all_rows.append(row)

    out_tsv = args.output_dir / "all_target_loci.tsv"
    if all_rows:
        out_df = pd.DataFrame(all_rows)
        out_df.to_csv(out_tsv, sep="\t", index=False)
        print()
        print(
            f"[4] Saved {len(out_df)} target hits to {out_tsv} "
            f"({out_df['genome'].nunique()} genomes)"
        )
    else:
        print()
        print("[4] No target hits passed clustering thresholds; nothing written")

    # Validate BK gene recovery against Phase 1
    print()
    print("=" * 80)
    print("BK GENE RECOVERY VALIDATION")
    print("=" * 80)

    # Count Phase 1 BK targets
    phase1_bk_count = 0
    phase1_dir = args.phase1_dir
    for locus_dir in phase1_dir.glob("BK_*"):
        if locus_dir.is_dir():
            for faa_file in locus_dir.glob("*_targets.faa"):
                with open(faa_file) as f:
                    phase1_bk_count += sum(1 for line in f if line.startswith('>'))

    # Count Phase 4 BK targets
    phase4_bk_count = 0
    if all_rows:
        phase4_bk_count = sum(1 for row in all_rows if 'Belonocnema_kinseyi' in row['genome'])

    print(f"Phase 1 BK targets (ground truth): {phase1_bk_count}")
    print(f"Phase 4 BK targets (detected): {phase4_bk_count}")

    if phase4_bk_count < phase1_bk_count:
        missing = phase1_bk_count - phase4_bk_count
        print(f"⚠️  WARNING: Phase 4 failed to recover {missing} BK genes from Phase 1!")
        print(f"   This may indicate clustering parameters need adjustment.")
    elif phase4_bk_count > phase1_bk_count:
        extra = phase4_bk_count - phase1_bk_count
        print(f"ℹ️  Phase 4 found {extra} additional BK targets beyond Phase 1")
        print(f"   This may indicate unannotated genes or over-splitting.")
    else:
        print(f"✓ Perfect BK recovery: {phase4_bk_count}/{phase1_bk_count}")

    print()
    print("=" * 80)
    print("PHASE 4 (REDESIGN) COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    main()
