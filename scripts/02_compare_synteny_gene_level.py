#!/usr/bin/env python3
print("DEBUG: Script starting...")  # Debug line to see if script runs at all
"""
Phase 2: Compare loci using gene-level synteny.

Determines which loci are:
- Orthologs: Different genomes, conserved synteny
- Paralogs: Same genome, different locations
- Tandem arrays: Adjacent copies in same genome

Usage:
    python 02_compare_synteny_gene_level.py <loci_dir>
"""

import sys
from pathlib import Path
import pandas as pd
from Bio import SeqIO
import numpy as np
from collections import defaultdict
import json
import subprocess

# Import universal synteny detection
from synteny_utils import SyntenyComparator, MIN_SYNTENY_MATCHES, SLIDING_WINDOW_SIZE, MIN_WINDOW_MATCHES, EVALUE_THRESHOLD

# Full landmark proteomes for proper synteny detection
LANDMARK_PROTEOMES = {
    'BK': Path('data/proteomes/GCF_010883055.1_Bkinseyi_NCBI.faa'),
    'LB': Path('data/proteomes/GCF_019393585.1_Lboulardi_NCBI.faa'),
    'TR': Path('data/proteomes/GCA_020615435.1.faa'),  # Telenomus remus (labeled as TR)
    'DR': Path('data/proteomes/GCA_030998225.1.faa')   # Diplolepis rosae (labeled as DR)
}


class Phase2SyntenyComparator(SyntenyComparator):
    """Compare loci based on flanking gene conservation."""

    def __init__(self):
        self.min_global = MIN_SYNTENY_MATCHES
        self.window_size = SLIDING_WINDOW_SIZE
        self.min_window = 3  # Lower threshold for divergent genomes - just need SOME clustering

    def compare_flanking(self, flanking_file1, flanking_file2, locus1_name, locus2_name):
        """Compare two sets of flanking proteins directly.

        Always compares flanking-to-flanking to avoid paralog contamination
        from full proteome searches.
        """

        # Compare flanking to flanking bidirectionally
        result1 = self._search_direction(flanking_file1, flanking_file2,
                                        f"{locus1_name.replace('/', '_')}_vs_{locus2_name.replace('/', '_')}")
        result2 = self._search_direction(flanking_file2, flanking_file1,
                                        f"{locus2_name.replace('/', '_')}_vs_{locus1_name.replace('/', '_')}")

        # Use the direction with higher match percentage
        if result1['percent'] >= result2['percent']:
            return result1
        else:
            # Return result2 but mark it as reverse search for clarity
            result2['search_direction'] = 'reverse'
            return result2

    def _search_direction(self, query_file, target_file, label):
        """Search query_file against target_file database."""

        # Create unique temp files to avoid race conditions in parallel array jobs
        import os
        pid = os.getpid()
        temp_db = Path(f"temp_db_{label}_{pid}")
        temp_output = Path(f"temp_compare_{label}_{pid}.tsv")

        try:
            # Check if this is a full proteome with existing database
            if str(target_file).endswith('.faa') and Path(str(target_file).replace('.faa', '.dmnd')).exists():
                # Use existing DIAMOND database
                temp_db = Path(str(target_file).replace('.faa', ''))
            else:
                # Make temporary database
                subprocess.run([
                    'diamond', 'makedb',
                    '--in', str(target_file),
                    '--db', str(temp_db)
                ], check=True, capture_output=True)

            # Search
            subprocess.run([
                'diamond', 'blastp',
                '--query', str(query_file),
                '--db', str(temp_db),
                '--out', str(temp_output),
                '--outfmt', '6',
                '--evalue', str(EVALUE_THRESHOLD),
                '--max-target-seqs', '1',
                '--threads', '1'
            ], check=True, capture_output=True)

            # Count unique query matches
            if temp_output.exists() and temp_output.stat().st_size > 0:
                hits = pd.read_csv(temp_output, sep='\t', header=None,
                                 names=['qseqid', 'sseqid', 'pident', 'length',
                                       'mismatch', 'gapopen', 'qstart', 'qend',
                                       'sstart', 'send', 'evalue', 'bitscore'])

                # Count unique queries that found a match
                unique_matches = hits['qseqid'].nunique()

                # Get total queries
                total_queries = len(list(SeqIO.parse(query_file, 'fasta')))

                return {
                    'matches': unique_matches,
                    'total': total_queries,
                    'percent': (unique_matches / total_queries * 100) if total_queries > 0 else 0,
                    'passes_global': unique_matches >= self.min_global,
                    'hits_df': hits
                }
            else:
                return {
                    'matches': 0,
                    'total': len(list(SeqIO.parse(query_file, 'fasta'))),
                    'percent': 0,
                    'passes_global': False,
                    'hits_df': pd.DataFrame()
                }

        finally:
            # Cleanup
            if temp_output.exists():
                temp_output.unlink()
            # Only delete database if it was a temp one we created
            if not str(target_file).endswith('.faa') and temp_db.with_suffix('.dmnd').exists():
                temp_db.with_suffix('.dmnd').unlink()

    def check_sliding_window(self, hits_df, query_order):
        """Check for concentrated synteny in sliding window."""

        if hits_df.empty:
            return False, 0, -1

        # Map queries to positions
        query_positions = {q: i for i, q in enumerate(query_order)}

        # Create binary array of matches
        match_array = np.zeros(len(query_order))
        for query in hits_df['qseqid'].unique():
            if query in query_positions:
                match_array[query_positions[query]] = 1

        # Slide window
        best_matches = 0
        best_position = -1

        for i in range(len(match_array) - self.window_size + 1):
            window = match_array[i:i + self.window_size]
            window_matches = int(np.sum(window))

            if window_matches > best_matches:
                best_matches = window_matches
                best_position = i

        passes_window = best_matches >= self.min_window

        return passes_window, best_matches, best_position

    def classify_relationship(self, comparison, same_genome, is_tandem1, is_tandem2, passes_window):
        """Classify the relationship between two loci.

        Requires BOTH global and window thresholds to confirm synteny:
        - Global: ≥9 matches (sufficient conservation)
        - Window: ≥3/10 matches (some clustering, not random scatter)
        """

        # Synteny requires both global AND window thresholds
        has_synteny = comparison['passes_global'] and passes_window

        if same_genome:
            if (is_tandem1 or is_tandem2) and has_synteny:
                return "tandem_array"  # Same genome, tandem, AND synteny
            elif has_synteny:
                return "duplicate"  # Same genome, high synteny
            else:
                return "paralog"  # Same genome, different location (no synteny)
        else:
            if has_synteny:
                return "ortholog"  # Different genome, conserved synteny
            else:
                return "unrelated"  # Different genome, no synteny


def main():
    if len(sys.argv) < 2:
        print("Usage: python 02_compare_synteny_gene_level.py <loci_dir>")
        sys.exit(1)

    loci_dir = Path(sys.argv[1])

    print("="*80)
    print("PHASE 2: SYNTENY-BASED LOCUS COMPARISON")
    print("="*80)
    print(f"\nInput directory: {loci_dir}")
    print(f"Synteny thresholds (BOTH required):")
    print(f"  Global: ≥{MIN_SYNTENY_MATCHES} matches")
    print(f"  Window: ≥3/{SLIDING_WINDOW_SIZE} matches (lowered for divergent genomes)")

    # Load loci information
    loci_file = loci_dir / "locus_definitions.tsv"
    if not loci_file.exists():
        print(f"ERROR: {loci_file} not found")
        print("Run Phase 1 first")
        sys.exit(1)

    loci_df = pd.read_csv(loci_file, sep='\t')
    print(f"\nLoaded {len(loci_df)} loci from {loci_df['genome'].nunique()} genomes")

    # Initialize comparator
    comparator = Phase2SyntenyComparator()

    # Compare all pairs
    print("\n[1] Comparing all locus pairs...")
    comparisons = []

    total_pairs = len(loci_df) * (len(loci_df) - 1) // 2
    pair_count = 0

    for i, locus1 in loci_df.iterrows():
        for j, locus2 in loci_df.iterrows():
            if j <= i:  # Skip self and avoid duplicate comparisons
                continue

            pair_count += 1
            print(f"  Comparing {pair_count}: {locus1['locus_id']} vs {locus2['locus_id']}...")

            flanking1 = Path(locus1['flanking_file'])
            flanking2 = Path(locus2['flanking_file'])

            if not flanking1.exists() or not flanking2.exists():
                continue

            # Compare flanking proteins
            result = comparator.compare_flanking(
                flanking1, flanking2,
                locus1['locus_id'], locus2['locus_id']
            )

            # Get query order for sliding window - depends on search direction
            if result.get('search_direction') == 'reverse':
                # For reverse search, hits_df has queries from flanking2
                query_order = [rec.id for rec in SeqIO.parse(flanking2, 'fasta')]
            else:
                # Normal search, queries from flanking1
                query_order = [rec.id for rec in SeqIO.parse(flanking1, 'fasta')]

            passes_window, window_matches, window_pos = comparator.check_sliding_window(
                result['hits_df'], query_order
            )

            # Classify relationship
            same_genome = locus1['genome'] == locus2['genome']
            relationship = comparator.classify_relationship(
                result, same_genome,
                locus1['is_tandem'], locus2['is_tandem'],
                passes_window
            )

            comparison = {
                'locus1': locus1['locus_id'],
                'locus2': locus2['locus_id'],
                'genome1': locus1['genome'],
                'genome2': locus2['genome'],
                'same_genome': same_genome,
                'matches': result['matches'],
                'total_flanking': result['total'],
                'percent': result['percent'],
                'passes_global': result['passes_global'],
                'passes_window': passes_window,
                'best_window': window_matches,
                'window_position': window_pos,
                'relationship': relationship
            }

            comparisons.append(comparison)

            # Show significant relationships
            if result['passes_global'] or passes_window:
                print(f"  {locus1['locus_id']} ↔ {locus2['locus_id']}:")
                print(f"    Matches: {result['matches']}/{result['total']} ({result['percent']:.1f}%)")
                print(f"    Window: {window_matches}/{SLIDING_WINDOW_SIZE}")
                print(f"    → {relationship}")

    # Save comparison results
    print(f"\n[2] Saving results...")
    comparisons_df = pd.DataFrame(comparisons)
    comparisons_df.to_csv(loci_dir / "synteny_comparisons.tsv", sep='\t', index=False)

    # Group loci by synteny
    print("\n[3] Grouping loci by synteny...")
    synteny_groups = []
    used_loci = set()

    for _, locus in loci_df.iterrows():
        if locus['locus_id'] in used_loci:
            continue

        # Start new group - use first locus as group name
        group = {
            'group_id': locus['locus_id'],  # Use first locus name as group ID
            'members': [locus['locus_id']],
            'genomes': [locus['genome']]
        }
        used_loci.add(locus['locus_id'])

        # Find all loci with synteny to this one
        # Include tandem_array - these are just paralogs in conserved synteny
        if len(comparisons_df) > 0:
            related = comparisons_df[
                ((comparisons_df['locus1'] == locus['locus_id']) |
                 (comparisons_df['locus2'] == locus['locus_id'])) &
                (comparisons_df['relationship'].isin(['ortholog', 'duplicate', 'tandem_array']))
            ]
        else:
            # No comparisons (only 1 locus), so no related loci
            related = pd.DataFrame()

        for _, rel in related.iterrows():
            other = rel['locus2'] if rel['locus1'] == locus['locus_id'] else rel['locus1']
            if other not in used_loci:
                group['members'].append(other)
                # Find genome for this locus
                other_genome = loci_df[loci_df['locus_id'] == other]['genome'].values[0]
                group['genomes'].append(other_genome)
                used_loci.add(other)

        synteny_groups.append(group)

    # Save synteny groups
    with open(loci_dir / "synteny_groups.json", 'w') as f:
        json.dump(synteny_groups, f, indent=2)

    # Summary
    print("\n" + "="*80)
    print("SUMMARY")
    print("-"*80)

    print(f"\nTotal comparisons: {len(comparisons_df)}")

    # Count relationships
    if len(comparisons_df) > 0:
        rel_counts = comparisons_df['relationship'].value_counts()
        print("\nRelationship types:")
        for rel, count in rel_counts.items():
            print(f"  {rel}: {count}")
    else:
        print("\nNo comparisons (only 1 locus found)")

    print(f"\nSynteny groups: {len(synteny_groups)}")
    for group in synteny_groups:
        genomes = set(group['genomes'])
        print(f"  {group['group_id']}: {len(group['members'])} loci in {genomes}")

    # Find conserved orthologs
    conserved = [g for g in synteny_groups if len(set(g['genomes'])) > 1]
    print(f"\nConserved ortholog groups: {len(conserved)}")
    for group in conserved:
        print(f"  {group['group_id']}: {set(group['genomes'])}")

    print(f"\nResults saved to {loci_dir}/")
    print("  - synteny_comparisons.tsv: All pairwise comparisons")
    print("  - synteny_groups.json: Loci grouped by synteny")


if __name__ == "__main__":
    main()