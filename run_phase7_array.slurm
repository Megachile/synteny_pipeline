#!/bin/bash
#SBATCH --job-name=phase7_array
#SBATCH --output=logs/phase7_array_%A_%a.out
#SBATCH --error=logs/phase7_array_%A_%a.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --partition=general

# Phase 7 Array Job - Run SwissProt annotation per locus in parallel
#
# Usage: sbatch --array=0-N run_phase7_array.slurm <gene_family_output_dir>
#   where N = number of loci - 1
#
# Example:
#   # First get locus count:
#   NLOCI=$(tail -n +2 outputs/ferritin_MC117_BK_LB_batch/03_filtered_blocks/synteny_blocks_filtered.tsv | cut -f1 | sort -u | wc -l)
#   # Then submit array:
#   sbatch --array=0-$((NLOCI-1)) run_phase7_array.slurm outputs/ferritin_MC117_BK_LB_batch

set -euo pipefail

# Change to workflow directory
cd /carc/scratch/projects/emartins/2016456/adam/synteny_scanning/hpc_deployment_package/hybrid_workflow

# Check arguments
if [ $# -ne 1 ]; then
    echo "Usage: sbatch --array=0-N run_phase7_array.slurm <gene_family_output_dir>"
    exit 1
fi

OUTPUT_DIR="$1"

# Load environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate busco_env

# Get list of unique loci from filtered blocks
FILTERED_BLOCKS="${OUTPUT_DIR}/03_filtered_blocks/synteny_blocks_filtered.tsv"

if [ ! -f "$FILTERED_BLOCKS" ]; then
    echo "ERROR: Filtered blocks file not found: $FILTERED_BLOCKS"
    exit 1
fi

# Extract unique loci (skip header)
mapfile -t LOCI < <(tail -n +2 "$FILTERED_BLOCKS" | cut -f1 | sort -u)

NLOCI=${#LOCI[@]}

if [ ${SLURM_ARRAY_TASK_ID} -ge ${NLOCI} ]; then
    echo "ERROR: Array task ID ${SLURM_ARRAY_TASK_ID} exceeds number of loci (${NLOCI})"
    exit 1
fi

LOCUS="${LOCI[$SLURM_ARRAY_TASK_ID]}"

echo "================================================================================"
echo "PHASE 7 ARRAY TASK ${SLURM_ARRAY_TASK_ID}: Processing locus ${LOCUS}"
echo "================================================================================"
echo "Output directory: ${OUTPUT_DIR}"
echo "Total loci: ${NLOCI}"
echo "Job ID: ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "Node: $(hostname)"
echo ""

# Run Phase 7 for this specific locus
python scripts/07_swissprot_annotation.py \
    --synteny-dir "${OUTPUT_DIR}/02_synteny_blocks" \
    --filtered-blocks "${FILTERED_BLOCKS}" \
    --swissprot-db data/reference/swissprot.dmnd \
    --output-dir "${OUTPUT_DIR}/07_swissprot_annotations" \
    --locus "${LOCUS}" \
    --evalue 1e-5 \
    --threads 8

echo ""
echo "✓ Phase 7 complete for locus ${LOCUS}"
echo "Output: ${OUTPUT_DIR}/07_swissprot_annotations/${LOCUS}_swissprot_annotations.tsv"
