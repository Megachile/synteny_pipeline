#!/usr/bin/env python3
"""
Validate Phase 4 HSP clustering against Phase 1 gene coordinates.

Now updated for LOCUS-AGNOSTIC Phase 4 output:
- Phase 4 outputs gene positions without locus assignment
- Validation checks if each ground truth gene is recovered by position overlap
- Reports recovery rate and extra discoveries

This validates that "tandem duplicates are being divided correctly,
rather than oversplit, not deduplicated by overlap, etc."
"""

import os
import sys
import pandas as pd
from pathlib import Path
from collections import defaultdict

def get_families_with_phase4(phase4_suffix="phase4_v2"):
    """Find all families with both phase1 and phase4 outputs."""
    outputs_dir = Path("outputs")
    families = []
    for family_dir in outputs_dir.iterdir():
        if not family_dir.is_dir():
            continue
        phase1_file = family_dir / "phase1_v2" / "locus_definitions.tsv"
        phase4_file = family_dir / phase4_suffix / "all_target_loci.tsv"
        if phase1_file.exists() and phase4_file.exists():
            families.append(family_dir.name)
    return sorted(families)

def load_phase1_bk_genes(family):
    """Load individual BK gene coordinates from Phase 1."""
    coords_file = Path(f"outputs/{family}/phase1_v2/cluster_gene_coordinates.tsv")
    if not coords_file.exists():
        return pd.DataFrame()

    df = pd.read_csv(coords_file, sep='\t')
    # Filter to BK loci only
    bk_genes = df[df['locus_id'].str.startswith('BK_')].copy()
    return bk_genes

def load_phase4_bk_hits(family, phase4_suffix="phase4_v2"):
    """Load Phase 4 hits for BK genome."""
    phase4_file = Path(f"outputs/{family}/{phase4_suffix}/all_target_loci.tsv")
    if not phase4_file.exists():
        return pd.DataFrame()

    df = pd.read_csv(phase4_file, sep='\t')
    # Filter to BK genome (Belonocnema_kinseyi_GCF)
    bk_hits = df[df['genome'] == 'Belonocnema_kinseyi_GCF'].copy()
    return bk_hits

def normalize_chrom(chrom):
    """Remove version suffix from chromosome names (e.g., NC_046658.1 -> NC_046658)."""
    if pd.isna(chrom):
        return chrom
    if '.' in str(chrom):
        return str(chrom).rsplit('.', 1)[0]
    return str(chrom)

def check_overlap(gene_start, gene_end, hit_start, hit_end, min_overlap_bp=100):
    """
    Check if a Phase 4 hit overlaps with a ground truth gene.

    Uses minimum bp overlap rather than fraction because:
    - Genes span large regions (including introns)
    - Phase 4 finds coding exons which may be a small fraction of gene span
    - Any significant overlap indicates the same gene was found
    """
    overlap_start = max(gene_start, hit_start)
    overlap_end = min(gene_end, hit_end)

    if overlap_start >= overlap_end:
        return False, 0

    overlap_len = overlap_end - overlap_start
    gene_len = gene_end - gene_start

    overlap_frac = overlap_len / gene_len if gene_len > 0 else 0
    return overlap_len >= min_overlap_bp, overlap_frac

def validate_gene_recovery(bk_genes, bk_hits, family):
    """
    Validate per-gene recovery.

    Returns list of per-gene results with recovery status.
    """
    results = []

    for _, gene in bk_genes.iterrows():
        locus_id = gene['locus_id']
        gene_id = gene['gene_id']
        chrom = gene['chromosome']
        chrom_normalized = normalize_chrom(chrom)
        gene_start = gene['start']
        gene_end = gene['end']

        # Find overlapping Phase 4 hits
        overlapping_hits = []
        if len(bk_hits) > 0:
            for _, hit in bk_hits.iterrows():
                hit_chrom_normalized = normalize_chrom(hit['scaffold'])
                if hit_chrom_normalized == chrom_normalized:
                    has_overlap, overlap_frac = check_overlap(
                        gene_start, gene_end, hit['start'], hit['end']
                    )
                    if has_overlap:
                        overlapping_hits.append({
                            'start': hit['start'],
                            'end': hit['end'],
                            'overlap_frac': overlap_frac
                        })

        # Determine recovery status
        recovered = len(overlapping_hits) > 0
        num_hits = len(overlapping_hits)

        # Status: RECOVERED, SPLIT (multiple hits for one gene), or MISSED
        if num_hits == 0:
            status = 'MISSED'
        elif num_hits == 1:
            status = 'RECOVERED'
        else:
            status = 'SPLIT'  # Multiple Phase 4 hits cover one gene

        results.append({
            'family': family,
            'locus_id': locus_id,
            'gene_id': gene_id,
            'chromosome': chrom,
            'gene_start': gene_start,
            'gene_end': gene_end,
            'gene_span_bp': gene_end - gene_start,
            'recovered': recovered,
            'num_hits': num_hits,
            'status': status,
            'hit_positions': '; '.join([f"{h['start']}-{h['end']}" for h in overlapping_hits]) if overlapping_hits else 'NONE'
        })

    return results

def count_extra_hits(bk_genes, bk_hits, family):
    """
    Count Phase 4 hits that don't overlap any ground truth gene.
    These represent additional discoveries (possibly pseudogenes or unannotated genes).
    """
    if len(bk_hits) == 0:
        return 0, []

    extra_hits = []
    for _, hit in bk_hits.iterrows():
        hit_chrom = normalize_chrom(hit['scaffold'])
        hit_start = hit['start']
        hit_end = hit['end']

        # Check if this hit overlaps any ground truth gene
        overlaps_gene = False
        for _, gene in bk_genes.iterrows():
            gene_chrom = normalize_chrom(gene['chromosome'])
            if gene_chrom == hit_chrom:
                has_overlap, _ = check_overlap(
                    gene['start'], gene['end'], hit_start, hit_end, min_overlap_bp=100
                )
                if has_overlap:
                    overlaps_gene = True
                    break

        if not overlaps_gene:
            extra_hits.append({
                'family': family,
                'scaffold': hit['scaffold'],
                'start': hit_start,
                'end': hit_end,
                'span_kb': hit['span_kb']
            })

    return len(extra_hits), extra_hits

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Validate Phase 4 clustering')
    parser.add_argument('--phase4-suffix', default='phase4_v2',
                        help='Phase 4 output directory suffix (default: phase4_v2)')
    args = parser.parse_args()

    print("=" * 80)
    print("PHASE 4 CLUSTERING VALIDATION (LOCUS-AGNOSTIC)")
    print(f"Checking Phase 4 suffix: {args.phase4_suffix}")
    print("Validating per-gene recovery from ground truth coordinates")
    print("=" * 80)
    print()

    families = get_families_with_phase4(args.phase4_suffix)
    print(f"Found {len(families)} families with Phase 4 outputs")
    print()

    all_gene_results = []
    all_extra_hits = []
    family_summaries = []

    for family in families:
        bk_genes = load_phase1_bk_genes(family)
        bk_hits = load_phase4_bk_hits(family, args.phase4_suffix)

        if len(bk_genes) == 0:
            continue

        # Validate per-gene recovery
        gene_results = validate_gene_recovery(bk_genes, bk_hits, family)
        all_gene_results.extend(gene_results)

        # Count extra discoveries
        num_extra, extra_hits = count_extra_hits(bk_genes, bk_hits, family)
        all_extra_hits.extend(extra_hits)

        # Family summary
        recovered = sum(1 for r in gene_results if r['recovered'])
        total = len(gene_results)
        recovery_rate = recovered / total if total > 0 else 0

        family_summaries.append({
            'family': family,
            'ground_truth_genes': total,
            'recovered': recovered,
            'missed': total - recovered,
            'recovery_rate': recovery_rate,
            'phase4_hits': len(bk_hits),
            'extra_hits': num_extra,
            'ratio': len(bk_hits) / total if total > 0 else 0
        })

    results_df = pd.DataFrame(all_gene_results)
    summary_df = pd.DataFrame(family_summaries)

    # Overall summary
    print("=" * 80)
    print("OVERALL VALIDATION SUMMARY")
    print("=" * 80)

    total_genes = len(results_df)
    recovered = sum(results_df['recovered'])
    missed = total_genes - recovered

    print(f"Total ground truth BK genes: {total_genes}")
    print(f"Recovered by Phase 4: {recovered} ({100*recovered/total_genes:.1f}%)")
    print(f"Missed: {missed} ({100*missed/total_genes:.1f}%)")
    print()

    # Status breakdown
    status_counts = results_df['status'].value_counts()
    print("Recovery status breakdown:")
    for status in ['RECOVERED', 'SPLIT', 'MISSED']:
        count = status_counts.get(status, 0)
        pct = 100 * count / total_genes
        print(f"  {status}: {count} ({pct:.1f}%)")
    print()

    # Phase 4 hit counts
    total_phase4_hits = summary_df['phase4_hits'].sum()
    total_extra = len(all_extra_hits)
    print(f"Total Phase 4 BK hits: {total_phase4_hits}")
    print(f"Extra discoveries (no gene overlap): {total_extra}")
    print(f"Overall ratio (Phase4/GroundTruth): {total_phase4_hits/total_genes:.2f}x")
    print()

    # Families with issues
    print("=" * 80)
    print("PER-FAMILY SUMMARY")
    print("=" * 80)

    for _, row in summary_df.sort_values('ratio', ascending=False).iterrows():
        status = "✓" if row['recovery_rate'] == 1.0 else "⚠"
        ratio_str = f"{row['ratio']:.2f}x" if row['ratio'] != 1.0 else "1.00x"
        print(f"{status} {row['family']}: {row['recovered']}/{row['ground_truth_genes']} recovered "
              f"({100*row['recovery_rate']:.0f}%), P4 hits={row['phase4_hits']}, ratio={ratio_str}")
    print()

    # Show missed genes
    missed_genes = results_df[results_df['status'] == 'MISSED']
    if len(missed_genes) > 0:
        print("=" * 80)
        print("MISSED GENES (Phase 4 failed to find)")
        print("=" * 80)
        for _, row in missed_genes.iterrows():
            print(f"  {row['family']}/{row['locus_id']}/{row['gene_id']}: "
                  f"{row['chromosome']}:{row['gene_start']}-{row['gene_end']} ({row['gene_span_bp']}bp)")
        print()

    # Save results
    results_df.to_csv('phase4_gene_recovery.tsv', sep='\t', index=False)
    summary_df.to_csv('phase4_family_summary.tsv', sep='\t', index=False)

    if all_extra_hits:
        extra_df = pd.DataFrame(all_extra_hits)
        extra_df.to_csv('phase4_extra_discoveries.tsv', sep='\t', index=False)
        print(f"Saved extra discoveries to: phase4_extra_discoveries.tsv")

    print(f"Saved per-gene results to: phase4_gene_recovery.tsv")
    print(f"Saved family summary to: phase4_family_summary.tsv")

if __name__ == '__main__':
    main()
