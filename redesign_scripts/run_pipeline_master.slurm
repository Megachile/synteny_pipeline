#!/bin/bash
#SBATCH --job-name=pipeline_v2_master
#SBATCH --output=logs/pipeline_v2_%x_%j.out
#SBATCH --error=logs/pipeline_v2_%x_%j.err
#SBATCH --time=48:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --partition=general

################################################################################
# REDESIGNED SYNTENY PIPELINE MASTER COORDINATOR
#
# Usage: sbatch redesign_scripts/run_pipeline_master.slurm <gene_family> <output_dir>
#
# Example: sbatch redesign_scripts/run_pipeline_master.slurm MC6 MC6_test_run
#
# This script orchestrates the redesigned pipeline end-to-end:
#   Phase 1:   Discover paralogs and build locus definitions (single job)
#   Phase 2:   Synteny detection per locus (array job - auto-sized)
#   Phase 3:   Aggregate and filter blocks (single job)
#   Phase 4:   Coverage-driven target detection (single job)
#   Phase 5:   Classify targets (single job)
#   Phase 6:   Extract sequences with Exonerate (array job - auto-sized)
#   Phase 7:   SwissProt annotation (array job - auto-sized)
#   Phase 8:   Generate matrices (single job)
#
# Date: 2025-11-15
################################################################################

set -e  # Exit on error
set -u  # Exit on undefined variable
set -o pipefail  # Exit on pipe failure

# Parse arguments
GENE_FAMILY=${1:-}
OUTPUT_DIR=${2:-${1}}  # Default output dir = gene_family name

if [ -z "$GENE_FAMILY" ]; then
    echo "ERROR: Missing required argument"
    echo "Usage: sbatch redesign_scripts/run_pipeline_master.slurm <gene_family> [output_dir]"
    echo "Example: sbatch redesign_scripts/run_pipeline_master.slurm MC6"
    exit 1
fi

# Change to project directory
cd /carc/scratch/projects/emartins/2016456/adam/synteny_scanning/hpc_deployment_package/hybrid_workflow

# Create directories
mkdir -p logs
mkdir -p outputs/${OUTPUT_DIR}

# Logging functions
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a logs/pipeline_v2_${OUTPUT_DIR}_progress.log
}

log_phase() {
    echo ""
    echo "================================================================================"
    echo "  $1"
    echo "================================================================================"
    echo ""
}

wait_for_job() {
    local job_id=$1
    local phase_name=$2
    log "Waiting for ${phase_name} (job ${job_id}) to complete..."

    while squeue -j ${job_id} 2>/dev/null | grep -q ${job_id}; do
        sleep 30
    done

    # Check if job succeeded (exit code 0)
    if sacct -j ${job_id} --format=State --noheader | grep -q "COMPLETED"; then
        log "✓ ${phase_name} completed successfully"
        return 0
    else
        log "❌ ${phase_name} failed - check logs"
        exit 1
    fi
}

log "================================================================================"
log "REDESIGNED PIPELINE - ${GENE_FAMILY} → ${OUTPUT_DIR}"
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

################################################################################
# PHASE 1: LOCUS DISCOVERY
################################################################################
log_phase "PHASE 1: LOCUS DISCOVERY"

INPUT_FILE=$(ls inputs/${GENE_FAMILY}*_test.tsv 2>/dev/null | head -1)
if [ -z "$INPUT_FILE" ]; then
    INPUT_FILE=$(ls inputs/${GENE_FAMILY}*.tsv 2>/dev/null | grep -v "locus_definitions" | head -1)
fi
if [ -z "$INPUT_FILE" ]; then
    log "ERROR: No input file found matching inputs/${GENE_FAMILY}*.tsv"
    exit 1
fi
log "Using input: $INPUT_FILE"

# Extract LOC IDs
LOC_IDS=$(tail -n +2 "$INPUT_FILE" | cut -f2 | tr '\n' ',' | sed 's/,$//')
if [ -z "$LOC_IDS" ]; then
    log "ERROR: No LOC IDs found in $INPUT_FILE"
    exit 1
fi
log "LOC IDs: $LOC_IDS"

PHASE1_DIR="outputs/${OUTPUT_DIR}/phase1_v2"
python redesign_scripts/01_phase1.py \
    --loc-ids "$LOC_IDS" \
    --gene-family "${GENE_FAMILY}" \
    --output-dir "${PHASE1_DIR}" \
    2>&1 | tee -a logs/phase1_v2_${OUTPUT_DIR}.log

if [ ! -f "${PHASE1_DIR}/locus_definitions.tsv" ]; then
    log "ERROR: Phase 1 failed - locus_definitions.tsv not found"
    exit 1
fi

NUM_LOCI=$(tail -n +2 "${PHASE1_DIR}/locus_definitions.tsv" | wc -l)
log "✓ Phase 1 complete: ${NUM_LOCI} loci identified"

################################################################################
# PHASE 1.5: BK-LB LOCUS UNIFICATION
################################################################################
log_phase "PHASE 1.5: BK-LB LOCUS UNIFICATION"

LB_PROTEOME_DB="data/proteomes/GCF_019393585.1_Lboulardi_NCBI.dmnd"
if [ -f "${LB_PROTEOME_DB}" ]; then
    log "Merging BK and LB loci based on flanking context..."
    python redesign_scripts/01b_merge_bk_lb_loci.py \
        --locus-defs "${PHASE1_DIR}/locus_definitions.tsv" \
        --lb-proteome-db "${LB_PROTEOME_DB}" \
        --diamond-threads 8 \
        2>&1 | tee -a logs/phase1.5_v2_${OUTPUT_DIR}.log

    NUM_LOCI_MERGED=$(tail -n +2 "${PHASE1_DIR}/locus_definitions.tsv" | wc -l)
    log "✓ Phase 1.5 complete: ${NUM_LOCI_MERGED} loci after BK-LB merge"
    NUM_LOCI=${NUM_LOCI_MERGED}  # Update for Phase 2 array sizing
else
    log "⚠ Skipping Phase 1.5: LB proteome DB not found at ${LB_PROTEOME_DB}"
fi

################################################################################
# PHASE 2: SYNTENY DETECTION (ARRAY)
################################################################################
log_phase "PHASE 2: SYNTENY DETECTION (ARRAY)"

ARRAY_SIZE=$((NUM_LOCI - 1))
log "Submitting array job for ${NUM_LOCI} loci (0-${ARRAY_SIZE})"

PHASE2_JOB=$(sbatch --parsable --array=0-${ARRAY_SIZE} \
    redesign_scripts/run_phase2_array.slurm ${OUTPUT_DIR})
log "Phase 2 job ID: ${PHASE2_JOB}"

wait_for_job ${PHASE2_JOB} "Phase 2"

################################################################################
# PHASE 3: AGGREGATE AND FILTER
################################################################################
log_phase "PHASE 3: AGGREGATE AND FILTER"

PHASE3_JOB=$(sbatch --parsable \
    redesign_scripts/run_phase3_single.slurm ${OUTPUT_DIR})
log "Phase 3 job ID: ${PHASE3_JOB}"

wait_for_job ${PHASE3_JOB} "Phase 3"

################################################################################
# PHASE 4: TARGET DETECTION (COVERAGE-DRIVEN)
################################################################################
log_phase "PHASE 4: TARGET DETECTION (COVERAGE-DRIVEN)"

PHASE4_JOB=$(sbatch --parsable \
    redesign_scripts/run_phase4_single.slurm ${OUTPUT_DIR})
log "Phase 4 job ID: ${PHASE4_JOB}"

wait_for_job ${PHASE4_JOB} "Phase 4"

################################################################################
# PHASE 5: CLASSIFICATION
################################################################################
log_phase "PHASE 5: CLASSIFICATION"

PHASE5_JOB=$(sbatch --parsable \
    redesign_scripts/run_phase5_single.slurm ${OUTPUT_DIR})
log "Phase 5 job ID: ${PHASE5_JOB}"

wait_for_job ${PHASE5_JOB} "Phase 5"

################################################################################
# PHASE 6: SEQUENCE EXTRACTION (ARRAY)
################################################################################
log_phase "PHASE 6: SEQUENCE EXTRACTION"

# Create family list for array job
echo "${OUTPUT_DIR}" > inputs/${OUTPUT_DIR}_list.txt

PHASE6_JOB=$(sbatch --parsable --array=1 \
    redesign_scripts/run_phase6_array.slurm inputs/${OUTPUT_DIR}_list.txt)
log "Phase 6 job ID: ${PHASE6_JOB}"

wait_for_job ${PHASE6_JOB} "Phase 6"

################################################################################
# PHASE 7: SWISSPROT ANNOTATION (ARRAY)
################################################################################
log_phase "PHASE 7: SWISSPROT ANNOTATION"

PHASE7_JOB=$(sbatch --parsable --array=1 \
    redesign_scripts/run_phase7_array.slurm inputs/${OUTPUT_DIR}_list.txt)
log "Phase 7 job ID: ${PHASE7_JOB}"

wait_for_job ${PHASE7_JOB} "Phase 7"

################################################################################
# PHASE 8: MATRIX GENERATION
################################################################################
log_phase "PHASE 8: MATRIX GENERATION"

PHASE8_JOB=$(sbatch --parsable \
    redesign_scripts/run_phase8_single.slurm ${OUTPUT_DIR})
log "Phase 8 job ID: ${PHASE8_JOB}"

wait_for_job ${PHASE8_JOB} "Phase 8"

################################################################################
# PIPELINE COMPLETE
################################################################################
log "================================================================================"
log "PIPELINE COMPLETE!"
log "================================================================================"
log "Gene family: ${GENE_FAMILY}"
log "Output directory: outputs/${OUTPUT_DIR}/"
log "End time: $(date)"
log "Total runtime: $SECONDS seconds"
log ""
log "Final outputs:"
log "  - Phase 1: ${PHASE1_DIR}/locus_definitions.tsv"
log "  - Phase 4: outputs/${OUTPUT_DIR}/phase4_v2/all_target_loci.tsv"
log "  - Phase 5: outputs/${OUTPUT_DIR}/phase5_v2/all_targets_classified.tsv"
log "  - Phase 8: outputs/${OUTPUT_DIR}/phase8_v2/*_matrix.tsv"
log ""
log "Check logs/pipeline_v2_${OUTPUT_DIR}_progress.log for full details"
