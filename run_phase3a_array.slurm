#!/bin/bash
#SBATCH --job-name=phase3a_array
#SBATCH --output=logs/phase3a_array_%A_%a.out
#SBATCH --error=logs/phase3a_array_%A_%a.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --partition=general

# Phase 3a Array Job - Run genome-wide synteny detection per locus in parallel
#
# Usage: sbatch --array=0-N run_phase3a_array.slurm <locus_definitions.tsv> <output_dir>
#   where N = number of loci - 1
#
# Example:
#   # First get locus count:
#   NLOCI=$(tail -n +2 outputs/MC18_BK_LB_batch/locus_definitions.tsv | wc -l)
#   # Then submit array:
#   sbatch --array=0-$((NLOCI-1)) run_phase3a_array.slurm \
#       outputs/MC18_BK_LB_batch/locus_definitions.tsv \
#       outputs/MC18_BK_LB_batch

set -euo pipefail

# Change to workflow directory
cd /carc/scratch/projects/emartins/2016456/adam/synteny_scanning/hpc_deployment_package/hybrid_workflow

# Check arguments
if [ $# -ne 2 ]; then
    echo "Usage: sbatch --array=0-N run_phase3a_array.slurm <locus_definitions.tsv> <output_dir>"
    exit 1
fi

LOCUS_DEFS="$1"
OUTPUT_DIR="$2"

# Load environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate busco_env

# Get list of loci from locus definitions (skip header)
if [ ! -f "$LOCUS_DEFS" ]; then
    echo "ERROR: Locus definitions file not found: $LOCUS_DEFS"
    exit 1
fi

mapfile -t LOCI < <(tail -n +2 "$LOCUS_DEFS" | cut -f1)

NLOCI=${#LOCI[@]}

if [ ${SLURM_ARRAY_TASK_ID} -ge ${NLOCI} ]; then
    echo "ERROR: Array task ID ${SLURM_ARRAY_TASK_ID} exceeds number of loci (${NLOCI})"
    exit 1
fi

LOCUS="${LOCI[$SLURM_ARRAY_TASK_ID]}"
OUTPUT_FILE="${OUTPUT_DIR}/02_synteny_blocks/${LOCUS}_synteny_blocks.tsv"

echo "================================================================================"
echo "PHASE 3a ARRAY TASK ${SLURM_ARRAY_TASK_ID}: Processing locus ${LOCUS}"
echo "================================================================================"
echo "Output directory: ${OUTPUT_DIR}/02_synteny_blocks"
echo "Total loci: ${NLOCI}"
echo "Job ID: ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "Node: $(hostname)"
echo ""

# Check if output already exists and is non-empty
if [ -f "$OUTPUT_FILE" ] && [ -s "$OUTPUT_FILE" ]; then
    echo "✓ Phase 3a already complete for locus ${LOCUS} (skipping)"
    echo "Output: ${OUTPUT_FILE}"
    exit 0
fi

# Run Phase 3a for this specific locus
python scripts/02_synteny_detection.py \
    --locus-defs "${LOCUS_DEFS}" \
    --locus "${LOCUS}" \
    --genome-db-dir data/ragtag_dbs \
    --output-dir "${OUTPUT_DIR}/02_synteny_blocks" \
    --evalue 1e-5 \
    --max-targets 50 \
    --max-gap-kb 500 \
    --threads 16

echo ""
echo "✓ Phase 3a complete for locus ${LOCUS}"
echo "Output: ${OUTPUT_FILE}"
