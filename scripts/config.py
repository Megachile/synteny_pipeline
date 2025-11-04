#!/usr/bin/env python3
"""
Configuration file for synteny detection pipeline.
Paths are relative to support standalone deployment.
"""

from pathlib import Path
import os

# ==============================================================================
# INPUT CONFIGURATION
# ==============================================================================

# Base directory (parent of scripts/)
PIPELINE_DIR = Path(__file__).parent.parent

# Locus definitions file
LOCI_DEFINITIONS_FILE = Path(os.environ.get('SYNTENY_INPUT_FILE', 'inputs/mc6_locus_definitions.tsv'))

# Reference genome files (BK = GCF_010883055.1)
BK_GFF_FILE = PIPELINE_DIR / "data/reference/genomic.gff"
BK_PROTEINS_FILE = PIPELINE_DIR / "data/reference/protein.faa"

# Target genomes (RagTag assemblies)
RAGTAG_DB_DIR = PIPELINE_DIR / "data/ragtag_dbs"
RAGTAG_FASTA_DIR = PIPELINE_DIR / "data/ragtag_output"  # Genome FASTA files for Exonerate

# Species mapping and phylogenetic order
SPECIES_MAP_FILE = PIPELINE_DIR / "data/gca_to_species.tsv"

# SwissProt database for flanking protein annotation
SWISSPROT_DB = PIPELINE_DIR / "data/databases/swissprot.dmnd"
USE_DIAMOND = True  # Set to False to use traditional BLASTP instead

# ==============================================================================
# PROTEOME MODE (NEW)
# ==============================================================================

# Use proteomes instead of genomes
USE_PROTEOME_MODE = False  # If True, use DIAMOND blastp; if False, use tblastn (GENOME MODE for hybrid workflow)

# Proteome databases and files (only used if USE_PROTEOME_MODE = True)
PROTEOME_DB_DIR = PIPELINE_DIR / "data/proteome_dbs"  # DIAMOND databases
PROTEOME_FAA_DIR = PIPELINE_DIR / "data/proteomes"   # Protein FASTA files
PROTEOMES_DIR = PIPELINE_DIR / "data/proteomes"   # Alias for compatibility
PROTEOME_SUMMARY = PIPELINE_DIR / "data/final_proteome_summary.tsv"

# Proteome synteny parameters (protein index-based instead of genomic bp-based)
FLANKING_PROTEIN_COUNT = 20  # Number of proteins on each side to extract (±20 = 40 total query proteins)
NUM_FLANKING = 40  # Total flanking proteins on each side for extraction
SYNTENY_MAX_PROTEIN_GAP = 10  # Max gap in protein indices to be in same block
TARGET_CLUSTER_GAP_PROTEINS = 10  # Max gap in protein indices for target clustering

# ==============================================================================
# PIPELINE PARAMETERS
# ==============================================================================

# BLAST parameters
FLANKING_BLAST_EVALUE = "1e-5"
FLANKING_BLAST_MAX_TARGETS = 50  # Sufficient for synteny detection
MAX_FLANKING_FOR_SYNTENY = 15  # Limit flanking proteins used in Phase 3 (reduces from 100s to 15)
TARGET_BLAST_EVALUE = "1e-5"
TARGET_BLAST_MAX_TARGETS = 10000
SWISSPROT_BLAST_EVALUE = "1e-5"

# Synteny clustering
SYNTENY_MAX_GAP_KB = 500  # Maximum gap between proteins to be in same block
TARGET_CLUSTER_GAP_KB = 10  # Maximum gap between target hits to cluster

# Filtering
MIN_PROTEINS_FOR_BLOCK = 5  # Minimum query proteins matching in a block to call synteny

# ==============================================================================
# OUTPUT DIRECTORIES
# ==============================================================================

# Base output directory (override with SYNTENY_OUTPUTS_DIR environment variable)
OUTPUT_BASE = Path(os.environ.get('SYNTENY_OUTPUTS_DIR', 'outputs/mc6_phase3'))

# Step-specific outputs
STEP01_PROTEINS = OUTPUT_BASE / "01_extracted_proteins"
STEP02_SYNTENY = OUTPUT_BASE / "02_synteny_blocks"
STEP03_FILTERED = OUTPUT_BASE / "03_filtered_blocks"
STEP04_TARGETS = OUTPUT_BASE / "04_target_genes"
STEP05_CLASSIFIED = OUTPUT_BASE / "05_classified"
STEP06_SWISSPROT = OUTPUT_BASE / "06_swissprot"  # Old name, kept for compatibility
STEP07_SWISSPROT = OUTPUT_BASE / "07_swissprot"
STEP07_MATRICES = OUTPUT_BASE / "07_matrices"

# Final outputs
LOCUS_MATRICES_DIR = STEP07_MATRICES / "locus_specific"
SUMMARY_MATRICES_DIR = STEP07_MATRICES / "gene_type_summaries"

# ==============================================================================
# RUNTIME OPTIONS
# ==============================================================================

# Number of threads for BLAST
BLAST_THREADS = 16  # Use all available CPUs

# Verbose output
VERBOSE = True

# Clean intermediate files after completion
CLEAN_INTERMEDIATE = False

# ==============================================================================
# LOCUS DISCOVERY PARAMETERS (Step 00a)
# ==============================================================================

# BK DIAMOND database for self-discovery
BK_DIAMOND_DB = PIPELINE_DIR / "data/reference/bk_proteins"

# Discovery search parameters
DISCOVERY_EVALUE = 1e-10  # E-value threshold for DIAMOND search
DISCOVERY_PIDENT = 30.0   # Minimum percent identity for homologs
MAX_LOCUS_GAP_KB = 100    # Maximum gap (kb) between genes in same locus
