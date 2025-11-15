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
MAX_GAP_KB = 200.0                 # always split if gap > 200 kb


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


def summarize_cluster(hsps: List[Dict]) -> Dict:
    """Summarize a cluster of HSPs belonging to one query/scaffold/strand."""
    g_start = min(h["genomic_start"] for h in hsps)
    g_end = max(h["genomic_end"] for h in hsps)
    span_kb = (g_end - g_start) / 1000.0

    combined_cov = calculate_union_query_coverage(hsps)

    q_start_min = min(h["query_start"] for h in hsps)
    q_end_max = max(h["query_end"] for h in hsps)
    qlen = hsps[0].get("query_length", 0) or 0
    if qlen > 0:
        q_span_cov = (q_end_max - q_start_min + 1) / float(qlen)
    else:
        q_span_cov = 0.0

    best_evalue = min(h["evalue"] for h in hsps)
    best_bitscore = max(h["bitscore"] for h in hsps)

    return {
        "query_id": hsps[0]["query_id"],
        "scaffold": hsps[0]["scaffold"],
        "strand": hsps[0]["strand"],
        "genomic_start": g_start,
        "genomic_end": g_end,
        "genomic_span_kb": span_kb,
        "query_start_min": q_start_min,
        "query_end_max": q_end_max,
        "query_span_coverage": q_span_cov,
        "combined_query_coverage": combined_cov,
        "num_hsps": len(hsps),
        "best_evalue": best_evalue,
        "best_bitscore": best_bitscore,
    }


def cluster_hsps_greedy(
    hsps: List[Dict],
    paralog_coverage_threshold: float = PARALOG_COVERAGE_THRESHOLD,
    min_cluster_coverage: float = MIN_CLUSTER_COVERAGE,
) -> List[Dict]:
    """
    Cluster HSPs using greedy expansion based on UNION query coverage.

    For each (query_id, scaffold, strand):
      1. Sort HSPs by genomic_start.
      2. Greedily add consecutive HSPs as long as:
         - combined coverage <= paralog_coverage_threshold
         - gap between last and next HSP <= MAX_GAP_KB
      3. When either condition fails, finalize current cluster and start a new one.
      4. Keep clusters with combined coverage >= min_cluster_coverage.

    Note: this is per-query clustering. If multiple queries hit the same gene,
    that shows up as multiple rows with overlapping coordinates and different
    query_ids; cross-query dedup can be handled later (Phase 5).
    """
    if not hsps:
        return []

    df = pd.DataFrame(hsps)
    out: List[Dict] = []

    for _, q_group in df.groupby("query_id"):
        for (scaffold, strand), sg in q_group.groupby(["scaffold", "strand"]):
            sg_sorted = sg.sort_values("genomic_start").reset_index(drop=True)
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

                    test_cluster = current + [nxt]
                    test_cov = calculate_union_query_coverage(test_cluster)

                    if test_cov > paralog_coverage_threshold:
                        break
                    if gap_kb > MAX_GAP_KB:
                        break

                    current.append(nxt)
                    used.add(next_idx)

                cov = calculate_union_query_coverage(current)
                if cov >= min_cluster_coverage:
                    out.append(summarize_cluster(current))

    return out


def _normalize_genome_id(name: str) -> str:
    """Normalize genome ID to match Phase 2/3 conventions."""
    if name.startswith("GCA_") or name.startswith("GCF_"):
        parts = name.split("_")
        if len(parts) >= 2:
            return f"{parts[0]}_{parts[1]}"
    return name


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

    print()
    print("=" * 80)
    print("PHASE 4 (REDESIGN) COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    main()
