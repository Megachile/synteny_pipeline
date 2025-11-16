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

from pathlib import Path
import argparse
import subprocess
import xml.etree.ElementTree as ET
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
                if gap_kb > MAX_GAP_KB:
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


def cluster_hsps_greedy(
    hsps: List[Dict],
    paralog_coverage_threshold: float = PARALOG_COVERAGE_THRESHOLD,
    min_cluster_coverage: float = MIN_CLUSTER_COVERAGE,
    query_overlap_threshold: float = QUERY_OVERLAP_THRESHOLD,
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
            query_overlap_threshold
        )
        all_gene_candidates.extend(gene_candidates)

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
) -> float:
    """Auto-calibrate query overlap threshold using BK ground truth.

    Tests multiple thresholds on BK genome and selects the one that best
    matches Phase 1 BK gene count.

    Returns:
        Calibrated threshold (or default if calibration fails)
    """
    print()
    print("=" * 80)
    print("AUTO-CALIBRATING QUERY OVERLAP THRESHOLD")
    print("=" * 80)

    # Count Phase 1 BK targets (ground truth)
    phase1_bk_count = count_phase1_bk_targets(phase1_dir)

    if phase1_bk_count == 0:
        print("  No BK genes in Phase 1 - using default threshold")
        print(f"  → Default threshold: {default_threshold}")
        return default_threshold

    print(f"  Phase 1 BK ground truth: {phase1_bk_count} genes")

    # Find BK BLAST XML
    bk_xml = find_bk_blast_xml(blast_xml_dir)
    if not bk_xml:
        print("  No BK BLAST results found - using default threshold")
        print(f"  → Default threshold: {default_threshold}")
        return default_threshold

    print(f"  BK BLAST XML: {bk_xml.name}")
    print()

    # Test thresholds
    test_thresholds = [0.1, 0.2, 0.3, 0.5, 0.7]
    results = []

    print("  Testing thresholds on BK genome:")
    for threshold in test_thresholds:
        # Parse BLAST XML
        hsps = parse_blast_xml_with_query_coords(bk_xml)

        # Cluster with this threshold
        clusters = cluster_hsps_greedy(
            hsps,
            query_overlap_threshold=threshold,
        )

        bk_count = len(clusters)
        diff = abs(bk_count - phase1_bk_count)
        results.append((threshold, bk_count, diff))

        status = "✓" if diff == 0 else f"±{diff}"
        print(f"    Threshold {threshold:.1f}: {bk_count} BK genes ({status})")

    # Select best threshold (minimum difference from Phase 1)
    best = min(results, key=lambda x: x[2])
    calibrated_threshold = best[0]

    print()
    print(f"  → Calibrated threshold: {calibrated_threshold:.1f}")
    print(f"    (matches Phase 1 with {best[1]} genes, diff={best[2]})")
    print()

    return calibrated_threshold


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
    p.add_argument("--evalue", type=str, default="1e-5")
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

    # Determine query overlap threshold (calibrate or use fixed)
    if args.query_overlap_threshold is not None:
        # User specified threshold explicitly
        final_threshold = args.query_overlap_threshold
        print(f"Using user-specified threshold: {final_threshold}")
    elif args.disable_calibration:
        # Calibration disabled, use default
        final_threshold = QUERY_OVERLAP_THRESHOLD
        print(f"Auto-calibration disabled, using default threshold: {final_threshold}")
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

        # Calibrate threshold
        final_threshold = calibrate_threshold(
            args.phase1_dir,
            blast_dir,
            default_threshold=QUERY_OVERLAP_THRESHOLD,
        )

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

        hsps = parse_blast_xml_with_query_coords(xml_path)
        if not hsps:
            print(f"    {genome_name}: no HSPs found")
            continue

        clusters = cluster_hsps_greedy(
            hsps,
            paralog_coverage_threshold=PARALOG_COVERAGE_THRESHOLD,
            min_cluster_coverage=MIN_CLUSTER_COVERAGE,
            query_overlap_threshold=final_threshold,
        )

        if not clusters:
            print(f"    {genome_name}: no clusters above coverage threshold")
            continue

        print(
            f"    {genome_name}: {len(hsps)} HSPs → {len(clusters)} gene-level hits"
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
