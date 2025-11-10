#!/usr/bin/env python3
# DEPRECATED: This configuration file is no longer used by the workflow.
# Thresholds based on percentages and sliding windows have been removed.
# The pipeline now uses count-based criteria (>=3 unique matches) and block ranking.
"""
Unified configuration for synteny pipeline - APPLIED TO CANONICAL SCRIPTS.

Key changes made to fix the pipeline:
1. Phase 1 (01_discover_paralogs_in_landmarks.py): Flanking excludes target (80 proteins)
2. Phase 2 (02_compare_flanking_and_deduplicate.py): Uses 9-gene threshold (0.1125)
3. Phase 3 (03_scan_proteomes_with_synteny.py): Uses 9-gene threshold
4. Matrix generation: Counts unique proteins, not total BLAST hits

Based on ferritin_MC102 analysis:
- Maximum cross-species flanking matches: 8 genes
- Threshold of 9 genes keeps all cross-species loci separate
- Enables detection of gene loss events (synteny preserved, target missing)
"""

# Flanking protein parameters (EXCLUDING target)
NUM_FLANKING = 40  # Number of proteins upstream AND downstream
TOTAL_FLANKING = 80  # Total flanking proteins (target NOT included)
MAX_PROTEIN_GAP = 10  # Maximum gap allowed in protein indices

# Synteny detection thresholds - CONSISTENT ACROSS ALL PHASES
MIN_PROTEINS_FOR_BLOCK = 9  # Minimum unique flanking matches to call synteny
SYNTENY_THRESHOLD_PERCENT = 0.1125  # 9/80 = 11.25%

# DIAMOND search parameters
EVALUE_THRESHOLD = 1e-3
MIN_PIDENT = 30.0  # Minimum percent identity for flanking protein match
MAX_TARGET_SEQS = 10

# Phase-specific settings
class Phase2Config:
    """Phase 2: Deduplication within landmark genomes"""
    FLANKING_PROTEINS = NUM_FLANKING
    MIN_SYNTENY_MATCHES = MIN_PROTEINS_FOR_BLOCK  # Use same as Phase 3!
    SYNTENY_THRESHOLD = SYNTENY_THRESHOLD_PERCENT

    # For same-genome comparisons (use Jaccard similarity)
    JACCARD_THRESHOLD = 0.5  # If proteins share >50% overlap, same region

class Phase3Config:
    """Phase 3: Synteny scanning across all genomes"""
    FLANKING_PROTEINS = NUM_FLANKING
    MIN_SYNTENY_MATCHES = MIN_PROTEINS_FOR_BLOCK
    SYNTENY_THRESHOLD = SYNTENY_THRESHOLD_PERCENT

    # For filtering hits
    MIN_BITSCORE = 50.0
    MAX_EVALUE = EVALUE_THRESHOLD

# Validation functions
def validate_synteny_block(num_matches, total_flanking=NUM_FLANKING*2):
    """Check if enough flanking matches to call synteny."""
    return num_matches >= MIN_PROTEINS_FOR_BLOCK

def calculate_synteny_score(num_matches, total_flanking=NUM_FLANKING*2):
    """Calculate synteny percentage score."""
    return (num_matches / total_flanking) * 100

def classify_locus_relationship(num_matches, same_genome=False):
    """Classify relationship between two loci based on flanking matches."""
    if same_genome and num_matches > NUM_FLANKING:
        # Same genome, most proteins match = same region (isoforms)
        return "isoform"
    elif num_matches >= MIN_PROTEINS_FOR_BLOCK:
        # Enough matches for synteny
        if same_genome:
            return "duplicate"  # Tandem duplication or assembly artifact
        else:
            return "ortholog"  # Same syntenic location, different species
    else:
        # Not enough matches
        if same_genome:
            return "paralog"  # Different loci in same genome
        else:
            return "unrelated"  # No syntenic relationship

# Report configuration
def print_config():
    """Print current configuration settings."""
    print("=" * 60)
    print("SYNTENY PIPELINE CONFIGURATION")
    print("=" * 60)
    print(f"Flanking proteins: ±{NUM_FLANKING} (total {NUM_FLANKING*2})")
    print(f"Min synteny matches: {MIN_PROTEINS_FOR_BLOCK} genes")
    print(f"Synteny threshold: {SYNTENY_THRESHOLD_PERCENT*100:.1f}%")
    print(f"DIAMOND e-value: {EVALUE_THRESHOLD}")
    print(f"Min percent identity: {MIN_PIDENT}%")
    print("=" * 60)
    print("Relationship classification:")
    print(f"  ≥{MIN_PROTEINS_FOR_BLOCK} matches, same genome → duplicate/isoform")
    print(f"  ≥{MIN_PROTEINS_FOR_BLOCK} matches, diff genome → ortholog")
    print(f"  <{MIN_PROTEINS_FOR_BLOCK} matches, same genome → paralog")
    print(f"  <{MIN_PROTEINS_FOR_BLOCK} matches, diff genome → unrelated")
    print("=" * 60)

if __name__ == "__main__":
    print_config()

    # Test with ferritin examples
    print("\nFerritin test cases:")
    print("-" * 40)

    test_cases = [
        ("BK vs BK (self)", 81, True),
        ("BK vs TR", 6, False),
        ("BK vs LB", 8, False),
        ("LB paralogs", 1, True),
        ("DR paralogs", 2, True),
    ]

    for desc, matches, same_genome in test_cases:
        relationship = classify_locus_relationship(matches, same_genome)
        score = calculate_synteny_score(matches)
        valid = validate_synteny_block(matches)

        print(f"{desc:15s}: {matches:2d} matches ({score:5.1f}%) → {relationship:10s} [synteny: {valid}]")
