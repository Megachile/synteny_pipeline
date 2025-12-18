"""
Helixer Pipeline Constants

Standardized column names, thresholds, and file naming conventions.
All scripts should import from this module to ensure consistency.
"""

# =============================================================================
# COLUMN NAMES - Standard names for TSV columns across all scripts
# =============================================================================

# Core identifiers
COL_GENOME = "genome"              # Genome accession (e.g., GCA_032370625.1)
COL_GENOME_ID = "genome_id"        # Alias for genome in some contexts
COL_TARGET_ID = "target_id"        # Target protein ID from Helixer
COL_TARGET_GENE_ID = "target_gene_id"  # Gene-level ID (may differ from protein)
COL_LOCUS_ID = "locus_id"          # Locus identifier (e.g., BK_LOC_chr1_1)

# Classification columns
COL_PLACEMENT = "placement"        # Classification result: 'synteny', 'unplaceable', 'novel'
COL_ASSIGNED_TO = "assigned_to"    # Which locus a target was assigned to

# Legacy column names (for backward compatibility during transition)
COL_CLASSIFICATION_LEGACY = "classification"  # Old name for placement
COL_LOCUS_ID_LEGACY = "locus_id"              # Old name for assigned_to in some contexts

# Flanking gene columns
COL_FLANKING_POS = "flanking_position"  # e.g., U1, U2, D1, D2
COL_FLANKING_GENE = "flanking_gene_id"  # Flanking gene protein ID

# Quality columns
COL_LENGTH_AA = "length_aa"        # Protein length in amino acids
COL_LENGTH_FLAG = "length_flag"    # QC flag: 'ok', 'short', 'long'
COL_BITSCORE = "bitscore"          # DIAMOND bitscore
COL_PIDENT = "pident"              # Percent identity
COL_EVALUE = "evalue"              # E-value

# SwissProt annotation columns
COL_SWISSPROT_ID = "swissprot_id"  # SwissProt accession
COL_SWISSPROT_NAME = "swissprot_name"  # SwissProt entry name
COL_SWISSPROT_DESC = "swissprot_desc"  # SwissProt description


# =============================================================================
# COLUMN NORMALIZATION - Helper function to standardize column names
# =============================================================================

def normalize_columns(df, inplace=True):
    """
    Normalize column names to standard conventions.

    Renames legacy column names to their standard equivalents.

    Args:
        df: pandas DataFrame to normalize
        inplace: If True, modify in place; if False, return a copy

    Returns:
        DataFrame with normalized column names
    """
    if not inplace:
        df = df.copy()

    rename_map = {}

    # classification -> placement
    if COL_CLASSIFICATION_LEGACY in df.columns and COL_PLACEMENT not in df.columns:
        rename_map[COL_CLASSIFICATION_LEGACY] = COL_PLACEMENT

    # When locus_id is used as assignment (in target classification context)
    # Note: locus_id remains locus_id in locus_definitions.tsv
    # Only rename to assigned_to when it represents a classification assignment

    if rename_map:
        df.rename(columns=rename_map, inplace=True)

    return df


def validate_required_columns(df, required_cols, source_name=""):
    """
    Validate that required columns exist in a DataFrame.

    Args:
        df: pandas DataFrame to validate
        required_cols: List of required column names
        source_name: Name of the source file for error messages

    Raises:
        ValueError: If any required columns are missing
    """
    missing = set(required_cols) - set(df.columns)
    if missing:
        source_str = f" in {source_name}" if source_name else ""
        raise ValueError(f"Missing required columns{source_str}: {sorted(missing)}")


# =============================================================================
# FILE NAMING - Standard file names for pipeline outputs
# =============================================================================

# Phase 2b outputs
FILE_PHASE2B_SYNTENY_BLOCKS = "phase2b_synteny_blocks.tsv"
FILE_PHASE2B_FLANKING_DETAILS = "phase2b_flanking_details.tsv"

# Legacy file names (for backward compatibility)
FILE_EMPTY_SYNTENY_BLOCKS_LEGACY = "empty_synteny_blocks.tsv"
FILE_EMPTY_BLOCK_DETAILS_LEGACY = "empty_block_details.tsv"

# Phase 5 outputs
FILE_SYNTENIC_TARGETS = "syntenic_targets.tsv"
FILE_UNPLACEABLE_TARGETS = "unplaceable_targets.tsv"

# Phase 5b outputs
FILE_NOVEL_LOCI_CLUSTERED = "novel_loci_clustered.tsv"

# Phase 6 outputs
FILE_LENGTH_QC = "length_qc.tsv"

# Phase 7 outputs
FILE_FLANKING_SWISSPROT = "flanking_swissprot.tsv"


# =============================================================================
# THRESHOLDS - Pipeline parameters
# =============================================================================

# Phase 5b: Novel loci validation
MIN_NOVEL_LOCI_GENOMES = 3  # Minimum genomes for valid novel locus

# Phase 6: Length QC
LENGTH_RATIO_SHORT = 0.8   # Ratio below this = "short"
LENGTH_RATIO_LONG = 1.3    # Ratio above this = "long"

# Phase 8b: Unplaceable classification
LOWQUALITY_BITSCORE_THRESHOLD = 100  # Below this = "lowquality"


# =============================================================================
# PLACEMENT VALUES - Standard classification values
# =============================================================================

PLACEMENT_SYNTENY = "synteny"
PLACEMENT_UNPLACEABLE = "unplaceable"
PLACEMENT_NOVEL = "novel"


# =============================================================================
# HELIXER ID HANDLING - Target ID format utilities
# =============================================================================

def strip_helixer_suffix(protein_id: str) -> str:
    """
    Remove .1 or .2 suffix from Helixer protein IDs.

    Helixer proteomes use suffixes like .1, .2 for isoforms.
    This function strips those for consistent matching.

    Args:
        protein_id: Protein ID potentially with suffix

    Returns:
        Protein ID without suffix
    """
    import re
    return re.sub(r'\.\d+$', '', protein_id)


def get_protein_variants(protein_id: str) -> list:
    """
    Get all possible variants of a protein ID for lookup.

    Returns the base ID plus common suffixes.

    Args:
        protein_id: Base protein ID

    Returns:
        List of possible ID variants to try
    """
    base = strip_helixer_suffix(protein_id)
    return [base, f"{base}.1", f"{base}.2"]


def lookup_protein(proteome: dict, protein_id: str) -> str:
    """
    Look up a protein sequence trying ID variants.

    Helixer proteomes use .1, .2 suffixes for isoforms.
    This function tries the base ID and common suffixes.

    Args:
        proteome: Dict of protein_id -> sequence
        protein_id: Target protein ID

    Returns:
        Protein sequence if found, else None
    """
    for variant in get_protein_variants(protein_id):
        if variant in proteome:
            return proteome[variant]
    return None
