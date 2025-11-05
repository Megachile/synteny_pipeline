#!/usr/bin/env python3
"""
Universal synteny detection module - copied from Phase 2's working implementation.
Uses bidirectional search with full proteomes and variable flanking counts.
"""

import subprocess
import pandas as pd
import numpy as np
from pathlib import Path
from Bio import SeqIO

# Universal thresholds
MIN_SYNTENY_MATCHES = 9       # Global threshold
SLIDING_WINDOW_SIZE = 10      # Window size
MIN_WINDOW_MATCHES = 6        # Window threshold
EVALUE_THRESHOLD = 1e-3


class SyntenyComparator:
    """Universal synteny detection using proven Phase 2 logic."""

    def __init__(self):
        self.min_global = MIN_SYNTENY_MATCHES
        self.window_size = SLIDING_WINDOW_SIZE
        self.min_window = MIN_WINDOW_MATCHES

    def search_flanking_against_proteome(self, flanking_file, proteome_file, label):
        """
        Search flanking proteins against a full proteome.
        Returns match statistics.
        """
        # Put temp files in dedicated directory
        temp_dir = Path("temp_diamond_dbs")
        temp_dir.mkdir(exist_ok=True)

        temp_db = temp_dir / f"temp_db_{label}"
        temp_output = temp_dir / f"temp_compare_{label}.tsv"

        try:
            # Check if proteome has existing database
            if str(proteome_file).endswith('.faa') and Path(str(proteome_file).replace('.faa', '.dmnd')).exists():
                # Use existing DIAMOND database
                temp_db = Path(str(proteome_file).replace('.faa', ''))
            else:
                # Make temporary database
                subprocess.run([
                    'diamond', 'makedb',
                    '--in', str(proteome_file),
                    '--db', str(temp_db)
                ], check=True, capture_output=True)

            # Search
            subprocess.run([
                'diamond', 'blastp',
                '--query', str(flanking_file),
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
                total_queries = len(list(SeqIO.parse(flanking_file, 'fasta')))

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
                    'total': len(list(SeqIO.parse(flanking_file, 'fasta'))),
                    'percent': 0,
                    'passes_global': False,
                    'hits_df': pd.DataFrame()
                }

        finally:
            # Cleanup
            if temp_output.exists():
                temp_output.unlink()
            # Only delete database if it was a temp one we created
            if not str(proteome_file).endswith('.faa') and temp_db.with_suffix('.dmnd').exists():
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

    def evaluate_synteny(self, blast_hits, query_proteins):
        """
        Evaluate synteny from BLAST hits.

        Args:
            blast_hits: DataFrame with BLAST results (must have 'qseqid' column)
            query_proteins: List of query protein IDs

        Returns:
            dict with synteny metrics
        """
        # Count unique matches
        unique_matches = blast_hits['qseqid'].nunique() if not blast_hits.empty else 0
        total_queries = len(query_proteins)

        # Global threshold
        passes_global = unique_matches >= self.min_global
        global_percent = (unique_matches / total_queries * 100) if total_queries > 0 else 0

        # Window threshold
        passes_window, best_window, _ = self.check_sliding_window(blast_hits, query_proteins)

        return {
            'passes_either': passes_global or passes_window,
            'total_matches': unique_matches,
            'global_percent': global_percent,
            'passes_global': passes_global,
            'passes_window': passes_window,
            'best_window_matches': best_window
        }


# Keep old name for backward compatibility
SyntenyDetector = SyntenyComparator


def create_standard_detector():
    """Create detector with standard thresholds."""
    return SyntenyComparator()
