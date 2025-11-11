#!/bin/bash
#SBATCH --job-name=phase2_synteny_v2
#SBATCH --array=0-0
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH --time=24:00:00
#SBATCH --partition=general
#SBATCH --output=logs/phase2_v2_%A_%a.out
#SBATCH --error=logs/phase2_v2_%A_%a.err

# Usage:
#   # After Phase 1 (redesign) writes outputs/<FAMILY>/phase1_v2/locus_definitions.tsv
#   FAMILY=ferritin_MC102
#   LOCUS_DEFS=outputs/${FAMILY}/phase1_v2/locus_definitions.tsv
#   N=$(( $(wc -l < ${LOCUS_DEFS}) - 1 ))
#   sbatch --array=0-$((N-1)) redesign_scripts/run_phase2_array.slurm ${FAMILY}
#
# Notes:
# - This script expects micromamba env providing tblastn (adjust as needed).
# - Each array task processes one locus (by row index in locus_definitions.tsv).

set -euo pipefail

cd /carc/scratch/projects/emartins/2016456/adam/synteny_scanning/hpc_deployment_package/hybrid_workflow

mkdir -p logs

# Activate environment (adjust if your cluster uses modules)
if command -v micromamba >/dev/null 2>&1; then
  eval "$(micromamba shell hook --shell bash)"
  micromamba activate busco_env || true
fi

if [ $# -lt 1 ]; then
  echo "Usage: $0 <FAMILY>" >&2
  exit 1
fi

FAMILY="$1"
PHASE1_DIR="outputs/${FAMILY}/phase1_v2"
PHASE2_DIR="outputs/${FAMILY}/phase2_synteny_v2"
LOCUS_DEFS="${PHASE1_DIR}/locus_definitions.tsv"

if [ ! -f "${LOCUS_DEFS}" ]; then
  echo "ERROR: locus_definitions.tsv not found: ${LOCUS_DEFS}" >&2
  exit 1
fi

mkdir -p "${PHASE2_DIR}"

# Extract locus_id for this array index
idx=${SLURM_ARRAY_TASK_ID:-0}
locus_id=$(awk -v i=$((idx+1)) 'BEGIN{FS="\t"} NR==1{next} NR==i+1{print $1; exit}' "${LOCUS_DEFS}" )
if [ -z "${locus_id}" ]; then
  echo "ERROR: Could not resolve locus_id for index ${idx}" >&2
  exit 1
fi

echo "Processing family: ${FAMILY} | locus: ${locus_id} | task ${SLURM_ARRAY_TASK_ID}"

python redesign_scripts/02_synteny_detection.py \
  --locus-defs "${LOCUS_DEFS}" \
  --genome-db-dir data/ragtag_dbs \
  --output-dir "${PHASE2_DIR}" \
  --evalue 1e-5 \
  --base-max-targets 1000 \
  --threads ${SLURM_CPUS_PER_TASK:-16} \
  --base-max-gap-kb 300 \
  --base-max-span-kb 800 \
  --base-merge-gap-kb 50 \
  --base-merge-span-kb 800 \
  --min-proteins 3 \
  --locus "${locus_id}"

echo "Done: ${FAMILY} ${locus_id}"

