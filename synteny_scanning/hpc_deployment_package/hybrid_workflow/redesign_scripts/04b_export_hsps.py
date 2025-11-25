#!/usr/bin/env python3
"""
Phase 4b: Export individual HSPs with cluster assignments.

This script re-processes existing Phase 4 BLAST XML outputs to create a detailed
HSP-level output file for use by Phase 6 HSP-guided extraction.

Outputs:
- all_hsps.tsv: Every HSP with its coordinates and assigned target/cluster ID

Usage:
    python 04b_export_hsps.py --phase4-dir outputs/family/phase4_v2 --phase1-dir outputs/family/phase1_v2

Can also be run as part of Phase 4 by importing export_hsps_from_clusters().
"""

from __future__ import annotations

import argparse
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

try:
    from synteny_utils import normalize_scaffold
except Exception:
    def normalize_scaffold(name: str) -> str:
        return str(name).split()[0]


# Import clustering functions from Phase 4
# These are duplicated here for standalone operation
PARALOG_COVERAGE_THRESHOLD = 1.2
MIN_CLUSTER_COVERAGE = 0.5
MAX_GAP_KB = 20.0
QUERY_OVERLAP_THRESHOLD = 0.5
MAX_CLUSTER_SPAN_KB = 500.0


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


def cluster_hsps_and_track(
    hsps: List[Dict],
    genome_name: str,
    query_to_locus: Dict[str, str],
    max_evalue: float = 1e-5,
) -> tuple[List[Dict], List[Dict]]:
    """
    Cluster HSPs into gene candidates and track which HSPs belong to which cluster.

    Returns:
        (clusters, hsp_assignments)

        clusters: List of cluster dicts with genomic coordinates
        hsp_assignments: List of HSP dicts with added 'target_id' field
    """
    # Filter by E-value
    hsps_filtered = [h for h in hsps if h["evalue"] <= max_evalue]
    if not hsps_filtered:
        return [], []

    # Group by (scaffold, strand)
    by_location: Dict[tuple[str, str], List[Dict]] = defaultdict(list)
    for h in hsps_filtered:
        key = (h["scaffold"], h["strand"])
        by_location[key].append(h)

    clusters = []
    hsp_assignments = []
    cluster_counter = 0

    for (scaffold, strand), location_hsps in by_location.items():
        # Sort by genomic position
        sorted_hsps = sorted(location_hsps, key=lambda x: x["genomic_start"])

        # Simple greedy clustering by genomic proximity
        current_cluster: List[Dict] = []
        current_end = 0

        for hsp in sorted_hsps:
            gap_kb = (hsp["genomic_start"] - current_end) / 1000.0 if current_end > 0 else 0

            if not current_cluster or gap_kb <= MAX_GAP_KB:
                # Add to current cluster
                current_cluster.append(hsp)
                current_end = max(current_end, hsp["genomic_end"])
            else:
                # Finalize current cluster and start new one
                if current_cluster:
                    cluster_id = f"{genome_name}_{scaffold}_{min(h['genomic_start'] for h in current_cluster)}_{max(h['genomic_end'] for h in current_cluster)}"

                    # Get parent locus from first HSP's query
                    first_query = current_cluster[0]["query_id"]
                    parent_locus = query_to_locus.get(first_query) or query_to_locus.get(first_query.split("|")[0], "unknown")

                    clusters.append({
                        "target_id": cluster_id,
                        "genome": genome_name,
                        "scaffold": scaffold,
                        "strand": strand,
                        "start": min(h["genomic_start"] for h in current_cluster),
                        "end": max(h["genomic_end"] for h in current_cluster),
                        "num_hsps": len(current_cluster),
                        "parent_locus": parent_locus,
                    })

                    # Assign HSPs to this cluster
                    for h in current_cluster:
                        h_copy = h.copy()
                        h_copy["target_id"] = cluster_id
                        h_copy["genome"] = genome_name
                        hsp_assignments.append(h_copy)

                    cluster_counter += 1

                # Start new cluster
                current_cluster = [hsp]
                current_end = hsp["genomic_end"]

        # Don't forget last cluster
        if current_cluster:
            cluster_id = f"{genome_name}_{scaffold}_{min(h['genomic_start'] for h in current_cluster)}_{max(h['genomic_end'] for h in current_cluster)}"

            first_query = current_cluster[0]["query_id"]
            parent_locus = query_to_locus.get(first_query) or query_to_locus.get(first_query.split("|")[0], "unknown")

            clusters.append({
                "target_id": cluster_id,
                "genome": genome_name,
                "scaffold": scaffold,
                "strand": strand,
                "start": min(h["genomic_start"] for h in current_cluster),
                "end": max(h["genomic_end"] for h in current_cluster),
                "num_hsps": len(current_cluster),
                "parent_locus": parent_locus,
            })

            for h in current_cluster:
                h_copy = h.copy()
                h_copy["target_id"] = cluster_id
                h_copy["genome"] = genome_name
                hsp_assignments.append(h_copy)

    return clusters, hsp_assignments


def export_hsps_from_phase4(
    phase4_dir: Path,
    phase1_dir: Path,
    output_file: Optional[Path] = None,
) -> Path:
    """
    Re-process Phase 4 BLAST XML files and export all HSPs with cluster assignments.

    Args:
        phase4_dir: Phase 4 output directory containing blast_xml/
        phase1_dir: Phase 1 directory for building query_to_locus mapping
        output_file: Output TSV path (default: phase4_dir/all_hsps.tsv)

    Returns:
        Path to output file
    """
    blast_xml_dir = phase4_dir / "blast_xml"
    if not blast_xml_dir.exists():
        raise FileNotFoundError(f"No blast_xml directory found in {phase4_dir}")

    if output_file is None:
        output_file = phase4_dir / "all_hsps.tsv"

    # Build query_to_locus mapping from Phase 1
    query_to_locus: Dict[str, str] = {}
    locus_defs = phase1_dir / "locus_definitions.tsv"

    if locus_defs.exists():
        loci_df = pd.read_csv(locus_defs, sep="\t")
        for _, row in loci_df.iterrows():
            locus_id = str(row["locus_id"])
            target_fa = phase1_dir / locus_id / f"{locus_id}_targets.faa"
            if target_fa.exists():
                with open(target_fa) as f:
                    for line in f:
                        if line.startswith(">"):
                            qid = line[1:].split()[0]
                            query_to_locus[qid] = locus_id

    print(f"Loaded {len(query_to_locus)} query-to-locus mappings from Phase 1")

    # Process each XML file
    all_hsps = []
    all_clusters = []

    xml_files = list(blast_xml_dir.glob("*.xml"))
    print(f"Processing {len(xml_files)} BLAST XML files...")

    for xml_file in xml_files:
        genome_name = xml_file.stem
        print(f"  {genome_name}...", end="")

        hsps = parse_blast_xml_with_query_coords(xml_file)
        if not hsps:
            print(" no HSPs")
            continue

        clusters, hsp_assignments = cluster_hsps_and_track(
            hsps, genome_name, query_to_locus
        )

        all_hsps.extend(hsp_assignments)
        all_clusters.extend(clusters)

        print(f" {len(hsps)} HSPs -> {len(clusters)} clusters")

    # Sort HSPs by genome, scaffold, position
    all_hsps.sort(key=lambda x: (x["genome"], x["scaffold"], x["genomic_start"]))

    # Write output
    output_columns = [
        "genome",
        "scaffold",
        "strand",
        "genomic_start",
        "genomic_end",
        "query_id",
        "query_start",
        "query_end",
        "query_length",
        "evalue",
        "bitscore",
        "target_id",
    ]

    df = pd.DataFrame(all_hsps)
    df = df[output_columns]
    df.to_csv(output_file, sep="\t", index=False)

    print()
    print(f"Wrote {len(all_hsps)} HSPs to {output_file}")
    print(f"  {len(all_clusters)} unique target clusters")
    print(f"  {df['genome'].nunique()} genomes")

    return output_file


def main():
    parser = argparse.ArgumentParser(
        description="Export Phase 4 HSPs with cluster assignments for Phase 6"
    )
    parser.add_argument(
        "--phase4-dir",
        type=Path,
        required=True,
        help="Phase 4 output directory containing blast_xml/",
    )
    parser.add_argument(
        "--phase1-dir",
        type=Path,
        required=True,
        help="Phase 1 directory containing locus_definitions.tsv",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output TSV file (default: phase4_dir/all_hsps.tsv)",
    )

    args = parser.parse_args()

    export_hsps_from_phase4(
        phase4_dir=args.phase4_dir,
        phase1_dir=args.phase1_dir,
        output_file=args.output,
    )


if __name__ == "__main__":
    main()
