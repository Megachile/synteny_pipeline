#!/bin/bash
#SBATCH --job-name=synteny_pipeline
#SBATCH --output=logs/pipeline_%x_%j.out
#SBATCH --error=logs/pipeline_%x_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --partition=general

################################################################################
# SYNTENY PIPELINE MASTER COORDINATOR
#
# Usage: sbatch run_synteny_pipeline_master.slurm <gene_family> <output_dir>
#
# Example: sbatch run_synteny_pipeline_master.slurm ferritin ferritin_phase3_real
#
# This script runs the entire synteny detection pipeline end-to-end:
#   Phase 1-2: Discover paralogs in landmark genomes and deduplicate by synteny
#   Phase 3:   Detect synteny blocks genome-wide (tBLASTn + clustering)
#   Phase 4:   BLAST for target genes in synteny blocks
#   Phase 5:   Classify targets (syntenic/tandem/unplaceable)
#   Phase 6:   Extract full gene structures with Exonerate
#   Phase 7:   SwissProt annotation of flanking proteins
#   Phase 8:   Generate matrices (locus-specific + summary)
#
# Date: 2025-11-04
# Updated: Fixed Phase 7 signature, added Phase 3b filter step
################################################################################

set -e  # Exit on error
set -u  # Exit on undefined variable
set -o pipefail  # Exit on pipe failure

# Parse arguments
GENE_FAMILY=${1:-}
OUTPUT_DIR=${2:-}

if [ -z "$GENE_FAMILY" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "ERROR: Missing required arguments"
    echo "Usage: sbatch run_synteny_pipeline_master.slurm <gene_family> <output_dir>"
    echo "Example: sbatch run_synteny_pipeline_master.slurm ferritin ferritin_run1"
    exit 1
fi

# Change to project directory
cd /carc/scratch/projects/emartins/2016456/adam/synteny_scanning/hpc_deployment_package/hybrid_workflow

# Create directories
mkdir -p logs
mkdir -p outputs/${OUTPUT_DIR}

# Logging functions
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a logs/pipeline_${GENE_FAMILY}_progress.log
}

log_phase() {
    echo ""
    echo "================================================================================"
    echo "  $1"
    echo "================================================================================"
    echo ""
}

check_file() {
    if [ ! -f "$1" ]; then
        log "ERROR: Required file not found: $1"
        exit 1
    fi
}

log "================================================================================"
log "SYNTENY PIPELINE - ${GENE_FAMILY}"
log "================================================================================"
log "Job ID: $SLURM_JOB_ID"
log "Node: $(hostname)"
log "Output directory: outputs/${OUTPUT_DIR}/"
log "Start time: $(date)"
log ""

# Load environment
log "Loading environment..."
set +u  # Temporarily disable for bashrc
source ~/.bashrc || { log "ERROR: Failed to source ~/.bashrc"; exit 1; }
micromamba activate busco_env || { log "ERROR: Failed to activate busco_env"; exit 1; }
set -u
log "Environment activated: busco_env"

# Verify tools available
log "Verifying tools..."
python --version 2>&1 | head -1 || { log "ERROR: Python not found"; exit 1; }
tblastn -version 2>&1 | head -1 || { log "ERROR: tblastn not found"; exit 1; }
diamond --version 2>&1 | head -1 || { log "ERROR: diamond not found"; exit 1; }

# Load Exonerate module
log "Loading Exonerate module..."
module load exonerate/2.4.0-yl7q || { log "ERROR: Failed to load exonerate module"; exit 1; }
which exonerate > /dev/null || { log "ERROR: exonerate not found in PATH"; exit 1; }

log "Environment ready"
log ""

################################################################################
# PHASE 1: PARALOG DISCOVERY IN LANDMARK GENOMES
################################################################################
log_phase "PHASE 1: PARALOG DISCOVERY IN LANDMARK GENOMES"

PHASE12_DIR="outputs/${OUTPUT_DIR}/phase12_landmark"
mkdir -p ${PHASE12_DIR}

log "Input: inputs/${GENE_FAMILY}*_test.tsv (auto-detected)"
log "Output: ${PHASE12_DIR}/"

# Find input file (prefer *_test.tsv, skip *_locus_definitions.tsv which is Phase 2 output)
INPUT_FILE=$(ls inputs/${GENE_FAMILY}*_test.tsv 2>/dev/null | head -1)
if [ -z "$INPUT_FILE" ]; then
    # Fallback: try any .tsv file that's NOT locus_definitions
    INPUT_FILE=$(ls inputs/${GENE_FAMILY}*.tsv 2>/dev/null | grep -v "locus_definitions" | head -1)
fi
if [ -z "$INPUT_FILE" ]; then
    log "ERROR: No input file found matching inputs/${GENE_FAMILY}*_test.tsv"
    exit 1
fi
log "Using input: $INPUT_FILE"

# Extract LOC IDs from the target_genes column (skip header, get 2nd column)
LOC_IDS=$(tail -n +2 "$INPUT_FILE" | cut -f2 | tr -d '\r')
if [ -z "$LOC_IDS" ]; then
    log "ERROR: No LOC IDs found in $INPUT_FILE"
    exit 1
fi
log "LOC IDs: $LOC_IDS"

if python scripts/01_discover_paralogs_gene_level.py \
    --loc-ids "$LOC_IDS" \
    --output-dir ${PHASE12_DIR} \
    --gene-family ${GENE_FAMILY} \
    2>&1 | tee -a logs/phase1_${GENE_FAMILY}.log; then

    log "✓ Phase 1 complete"

    # Verify outputs
    check_file "${PHASE12_DIR}/locus_definitions.tsv"
    NUM_LOCI=$(tail -n +2 "${PHASE12_DIR}/locus_definitions.tsv" | wc -l)
    log "  Found ${NUM_LOCI} loci across landmark genomes"
else
    log "❌ Phase 1 failed - check logs/phase1_${GENE_FAMILY}.log"
    exit 1
fi

################################################################################
# PHASE 2: SYNTENY-BASED DEDUPLICATION
################################################################################
log_phase "PHASE 2: SYNTENY-BASED DEDUPLICATION"

log "Input: ${PHASE12_DIR}/locus_definitions.tsv"
log "Output: ${PHASE12_DIR}/synteny_comparisons.tsv"

if python scripts/02_compare_synteny_gene_level.py ${PHASE12_DIR} \
    2>&1 | tee -a logs/phase2_${GENE_FAMILY}.log; then

    log "✓ Phase 2 complete"

    # Verify outputs
    check_file "${PHASE12_DIR}/synteny_comparisons.tsv"
    check_file "${PHASE12_DIR}/synteny_groups.json"

    NUM_GROUPS=$(python -c "import json; print(len(json.load(open('${PHASE12_DIR}/synteny_groups.json'))))")
    log "  Created ${NUM_GROUPS} synteny groups after deduplication"
else
    log "❌ Phase 2 failed - check logs/phase2_${GENE_FAMILY}.log"
    exit 1
fi

################################################################################
# PHASE 3a: GENOME-WIDE SYNTENY DETECTION (tBLASTn + CLUSTERING)
################################################################################
log_phase "PHASE 3a: GENOME-WIDE SYNTENY DETECTION"

PHASE3_DIR="outputs/${OUTPUT_DIR}"
mkdir -p ${PHASE3_DIR}

log "Input: ${PHASE12_DIR}/locus_definitions.tsv"
log "Genome databases: data/ragtag_dbs/"
log "Output: ${PHASE3_DIR}/02_synteny_blocks/"

if python scripts/02_synteny_detection.py \
    --locus-defs ${PHASE12_DIR}/locus_definitions.tsv \
    --genome-db-dir data/ragtag_dbs \
    --output-dir ${PHASE3_DIR}/02_synteny_blocks \
    --evalue 1e-5 \
    --max-targets 50 \
    --max-gap-kb 500 \
    --threads 16 \
    2>&1 | tee -a logs/phase3a_${GENE_FAMILY}.log; then

    log "✓ Phase 3a complete"

    # Check for aggregated output
    if [ -f "${PHASE3_DIR}/02_synteny_blocks/all_synteny_blocks.tsv" ]; then
        NUM_BLOCKS=$(tail -n +2 "${PHASE3_DIR}/02_synteny_blocks/all_synteny_blocks.tsv" | wc -l)
        NUM_GENOMES=$(tail -n +2 "${PHASE3_DIR}/02_synteny_blocks/all_synteny_blocks.tsv" | cut -f2 | sort -u | wc -l)
        log "  Found ${NUM_BLOCKS} synteny blocks across ${NUM_GENOMES} genomes"
    fi
else
    log "❌ Phase 3a failed - check logs/phase3a_${GENE_FAMILY}.log"
    exit 1
fi

################################################################################
# PHASE 3b: FILTER TO BEST BLOCK PER GENOME
################################################################################
log_phase "PHASE 3b: FILTER TO BEST BLOCK PER GENOME"

log "Input: ${PHASE3_DIR}/02_synteny_blocks/all_synteny_blocks.tsv"
log "Output: ${PHASE3_DIR}/03_filtered_blocks/synteny_blocks_filtered.tsv"

# Need to aggregate blocks first if not already done
if [ ! -f "${PHASE3_DIR}/02_synteny_blocks/all_synteny_blocks.tsv" ]; then
    log "Aggregating synteny blocks from per-locus files..."

    # Create header
    echo -e "locus_id\tgenome\tblock_id\tscaffold\tstrand\tstart\tend\tspan_kb\tnum_target_proteins\tnum_query_matches" \
        > "${PHASE3_DIR}/02_synteny_blocks/all_synteny_blocks.tsv"

    # Concatenate all locus blocks
    for locus_file in ${PHASE3_DIR}/02_synteny_blocks/*_synteny_blocks.tsv; do
        if [ -f "$locus_file" ]; then
            tail -n +2 "$locus_file" >> "${PHASE3_DIR}/02_synteny_blocks/all_synteny_blocks.tsv"
        fi
    done

    log "  Aggregated blocks from per-locus files"
fi

if python scripts/03_filter_blocks.py \
    --input ${PHASE3_DIR}/02_synteny_blocks/all_synteny_blocks.tsv \
    --output-dir ${PHASE3_DIR}/03_filtered_blocks \
    --verbose \
    2>&1 | tee -a logs/phase3b_${GENE_FAMILY}.log; then

    log "✓ Phase 3b complete"

    check_file "${PHASE3_DIR}/03_filtered_blocks/synteny_blocks_filtered.tsv"
    NUM_FILTERED=$(tail -n +2 "${PHASE3_DIR}/03_filtered_blocks/synteny_blocks_filtered.tsv" | wc -l)
    log "  Filtered to ${NUM_FILTERED} best blocks (one per genome per locus)"
else
    log "❌ Phase 3b failed - check logs/phase3b_${GENE_FAMILY}.log"
    exit 1
fi

################################################################################
# PHASE 4: TARGET GENE DETECTION
################################################################################
log_phase "PHASE 4: TARGET GENE DETECTION"

log "Input: ${PHASE3_DIR}/03_filtered_blocks/synteny_blocks_filtered.tsv"
log "Output: ${PHASE3_DIR}/04_target_genes/"

if python scripts/04_blast_targets.py \
    --blocks ${PHASE3_DIR}/03_filtered_blocks/synteny_blocks_filtered.tsv \
    --locus-defs ${PHASE12_DIR}/locus_definitions.tsv \
    --genome-fasta-dir data/ragtag_output \
    --output-dir ${PHASE3_DIR}/04_target_genes \
    --threads 16 \
    2>&1 | tee -a logs/phase4_${GENE_FAMILY}.log; then

    log "✓ Phase 4 complete"

    if [ -f "${PHASE3_DIR}/04_target_genes/all_target_loci.tsv" ]; then
        NUM_TARGETS=$(tail -n +2 "${PHASE3_DIR}/04_target_genes/all_target_loci.tsv" | wc -l)
        NUM_TARGET_GENOMES=$(tail -n +2 "${PHASE3_DIR}/04_target_genes/all_target_loci.tsv" | cut -f1 | sort -u | wc -l)
        log "  Found ${NUM_TARGETS} target loci across ${NUM_TARGET_GENOMES} genomes"
    fi
else
    log "❌ Phase 4 failed - check logs/phase4_${GENE_FAMILY}.log"
    exit 1
fi

################################################################################
# PHASE 5: FUNCTIONAL CLASSIFICATION
################################################################################
log_phase "PHASE 5: FUNCTIONAL CLASSIFICATION"

log "Input: ${PHASE3_DIR}/04_target_genes/all_target_loci.tsv"
log "Output: ${PHASE3_DIR}/05_classified/"

if python scripts/05_classify_targets.py \
    --targets ${PHASE3_DIR}/04_target_genes/all_target_loci.tsv \
    --blocks ${PHASE3_DIR}/03_filtered_blocks/synteny_blocks_filtered.tsv \
    --output-dir ${PHASE3_DIR}/05_classified \
    2>&1 | tee -a logs/phase5_${GENE_FAMILY}.log; then

    log "✓ Phase 5 complete"

    if [ -f "${PHASE3_DIR}/05_classified/syntenic_targets.tsv" ]; then
        NUM_SYNTENIC=$(tail -n +2 "${PHASE3_DIR}/05_classified/syntenic_targets.tsv" | wc -l)
        log "  Classified ${NUM_SYNTENIC} syntenic targets"
    fi

    if [ -f "${PHASE3_DIR}/05_classified/unplaceable_targets.tsv" ]; then
        NUM_UNPLACEABLE=$(tail -n +2 "${PHASE3_DIR}/05_classified/unplaceable_targets.tsv" | wc -l)
        log "  Found ${NUM_UNPLACEABLE} unplaceable targets"
    fi
else
    log "❌ Phase 5 failed - check logs/phase5_${GENE_FAMILY}.log"
    exit 1
fi

################################################################################
# PHASE 6: EXTRACT FULL GENE STRUCTURES WITH EXONERATE
################################################################################
log_phase "PHASE 6: EXTRACT FULL GENE STRUCTURES (EXONERATE)"

log "Input: ${PHASE3_DIR}/05_classified/"
log "Output: ${PHASE3_DIR}/06_extracted_sequences/"

if python scripts/06_extract_sequences.py \
    --syntenic ${PHASE3_DIR}/05_classified/syntenic_targets.tsv \
    --unplaceable ${PHASE3_DIR}/05_classified/unplaceable_targets.tsv \
    --query-proteins ${PHASE3_DIR}/04_target_genes/combined_targets.faa \
    --genome-fasta-dir data/ragtag_output \
    --output-dir ${PHASE3_DIR}/06_extracted_sequences \
    --unplaceable-evalue 1e-10 \
    --threads 16 \
    --verbose \
    2>&1 | tee -a logs/phase6_${GENE_FAMILY}.log; then

    log "✓ Phase 6 complete"

    # Count extracted sequences
    NUM_EXTRACTED=$(find ${PHASE3_DIR}/06_extracted_sequences -name "*.faa" 2>/dev/null | wc -l)
    log "  Extracted ${NUM_EXTRACTED} gene structures"
else
    log "❌ Phase 6 failed - check logs/phase6_${GENE_FAMILY}.log"
    exit 1
fi

################################################################################
# PHASE 7: SWISSPROT ANNOTATION OF FLANKING PROTEINS
################################################################################
log_phase "PHASE 7: SWISSPROT ANNOTATION (SYNTENY BLOCK PROTEINS ONLY)"

log "Input: ${PHASE3_DIR}/02_synteny_blocks/, ${PHASE3_DIR}/03_filtered_blocks/synteny_blocks_filtered.tsv"
log "Database: data/reference/swissprot.dmnd"
log "Output: ${PHASE3_DIR}/07_swissprot_annotations/"

if python scripts/07_swissprot_annotation.py \
    --synteny-dir ${PHASE3_DIR}/02_synteny_blocks \
    --filtered-blocks ${PHASE3_DIR}/03_filtered_blocks/synteny_blocks_filtered.tsv \
    --swissprot-db data/reference/swissprot.dmnd \
    --output-dir ${PHASE3_DIR}/07_swissprot_annotations \
    --evalue 1e-5 \
    --threads 16 \
    2>&1 | tee -a logs/phase7_${GENE_FAMILY}.log; then

    log "✓ Phase 7 complete"

    if [ -f "${PHASE3_DIR}/07_swissprot_annotations/genome_specific_swissprot_annotations.tsv" ]; then
        NUM_ANNOTATIONS=$(tail -n +2 "${PHASE3_DIR}/07_swissprot_annotations/genome_specific_swissprot_annotations.tsv" | wc -l)
        log "  Generated ${NUM_ANNOTATIONS} protein annotations"
    fi
else
    log "❌ Phase 7 failed - check logs/phase7_${GENE_FAMILY}.log"
    exit 1
fi

################################################################################
# PHASE 8a: GENERATE LOCUS-SPECIFIC MATRICES
################################################################################
log_phase "PHASE 8a: GENERATE LOCUS-SPECIFIC MATRICES"

log "Input: Multiple sources (locus defs, synteny blocks, targets, SwissProt)"
log "Output: ${PHASE3_DIR}/08_matrices/"

# Copy locus_definitions.tsv to phase3 dir for easy access
cp ${PHASE12_DIR}/locus_definitions.tsv ${PHASE3_DIR}/locus_definitions.tsv

if python scripts/08a_generate_locus_matrices.py \
    --locus-defs ${PHASE3_DIR}/locus_definitions.tsv \
    --synteny-dir ${PHASE3_DIR}/02_synteny_blocks \
    --blocks ${PHASE3_DIR}/03_filtered_blocks/synteny_blocks_filtered.tsv \
    --targets ${PHASE3_DIR}/05_classified/all_targets_classified.tsv \
    --swissprot ${PHASE3_DIR}/07_swissprot_annotations/genome_specific_swissprot_annotations.tsv \
    --reference-proteins data/reference/protein.faa \
    --species-map data/gca_to_species.tsv \
    --output-dir ${PHASE3_DIR}/08_matrices \
    2>&1 | tee -a logs/phase8a_${GENE_FAMILY}.log; then

    log "✓ Phase 8a complete"

    NUM_MATRICES=$(find ${PHASE3_DIR}/08_matrices -name "*_genome_swissprot_matrix.tsv" 2>/dev/null | wc -l)
    log "  Created ${NUM_MATRICES} locus-specific matrices"
else
    log "❌ Phase 8a failed - check logs/phase8a_${GENE_FAMILY}.log"
    exit 1
fi

################################################################################
# PHASE 8b: GENERATE SUMMARY MATRICES
################################################################################
log_phase "PHASE 8b: GENERATE SUMMARY MATRICES"

log "Input: ${PHASE3_DIR}/08_matrices/*_genome_swissprot_matrix.tsv"
log "Output: ${PHASE3_DIR}/08_matrices/*_summary.tsv"

if python scripts/08b_generate_summary_matrices.py \
    --locus-matrices-dir ${PHASE3_DIR}/08_matrices \
    --output-dir ${PHASE3_DIR}/08_matrices \
    2>&1 | tee -a logs/phase8b_${GENE_FAMILY}.log; then

    log "✓ Phase 8b complete"

    # List summary matrices
    log "  Summary matrices:"
    find ${PHASE3_DIR}/08_matrices -name "*summary*.tsv" 2>/dev/null | while read matrix; do
        log "    - $(basename $matrix)"
    done
else
    log "❌ Phase 8b failed - check logs/phase8b_${GENE_FAMILY}.log"
    exit 1
fi

################################################################################
# FINAL SUMMARY
################################################################################
log ""
log "================================================================================"
log "PIPELINE COMPLETE - SUCCESS!"
log "================================================================================"
log ""
log "Gene family: ${GENE_FAMILY}"
log "Output directory: outputs/${OUTPUT_DIR}/"
log "Total runtime: $(($SECONDS / 3600))h $(($SECONDS % 3600 / 60))m $(($SECONDS % 60))s"
log ""
log "RESULTS SUMMARY:"
log "  Phase 1-2: ${NUM_DEFS} loci identified"
log "  Phase 3:   ${NUM_FILTERED} synteny blocks (best per genome)"
log "  Phase 4:   ${NUM_TARGET_GENOMES} genomes with target genes"
log "  Phase 5:   ${NUM_SYNTENIC} syntenic, ${NUM_UNPLACEABLE} unplaceable"
log "  Phase 6:   ${NUM_EXTRACTED} sequences extracted"
log "  Phase 7:   ${NUM_ANNOTATIONS} protein annotations"
log "  Phase 8:   ${NUM_MATRICES} locus matrices created"
log ""
log "KEY OUTPUT FILES:"
log "  - Synteny blocks: ${PHASE3_DIR}/03_filtered_blocks/synteny_blocks_filtered.tsv"
log "  - Target genes: ${PHASE3_DIR}/04_target_genes/all_target_loci.tsv"
log "  - Locus matrices: ${PHASE3_DIR}/08_matrices/*_genome_swissprot_matrix.tsv"
log "  - Summary matrices: ${PHASE3_DIR}/08_matrices/*_summary.tsv"
log ""
log "LOGS:"
log "  - Progress: logs/pipeline_${GENE_FAMILY}_progress.log"
log "  - Per-phase: logs/phase*_${GENE_FAMILY}.log"
log ""
log "================================================================================"

exit 0
