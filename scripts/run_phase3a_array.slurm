#!/bin/bash
#SBATCH --job-name=phase3a_synteny
#SBATCH --output=logs/phase3a_%A_%a.out
#SBATCH --error=logs/phase3a_%A_%a.err
#SBATCH --time=6:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=16
#SBATCH --partition=general
#SBATCH --chdir=/carc/scratch/projects/emartins/2016456/adam/synteny_scanning/hpc_deployment_package/hybrid_workflow

# Usage: sbatch --array=0-(N-1) run_phase3a_array.slurm <locus_definitions.tsv> <family_root_dir>
set -euo pipefail

LOCUS_DEFS=${1:?"locus_definitions.tsv required"}
FAMILY_DIR=${2:?"family root dir required"}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${WORKFLOW_DIR}"

# Activate environment (match main organizer script)
eval "$(micromamba shell hook --shell bash)" || true
micromamba activate busco_env || true

mkdir -p logs
mkdir -p "${FAMILY_DIR}/02_synteny_blocks"

TASK_ID=${SLURM_ARRAY_TASK_ID:-0}

# Extract locus_id at this index (skip header)
LOCUS_ID=$(tail -n +2 "$LOCUS_DEFS" | awk -F '\t' -v idx=$((TASK_ID+1)) 'NR==idx{print $1}')
if [ -z "$LOCUS_ID" ]; then
  echo "ERROR: Could not resolve locus for index $TASK_ID from $LOCUS_DEFS" >&2
  exit 1
fi

echo "[Phase3a] Locus: $LOCUS_ID"

# Tunable params via env with sensible defaults
MAX_GAP_KB=${MAX_GAP_KB:-100}
MAX_SPAN_KB=${MAX_SPAN_KB:-300}
CLOSEST_EACH_SIDE=${CLOSEST_EACH_SIDE:-5}
# New: include internal anchors inside cluster span (I1..Ii)
CLOSEST_INTERNAL=${CLOSEST_INTERNAL:-5}
THREADS=${SLURM_CPUS_PER_TASK:-16}

python scripts/02_synteny_detection.py \
  --locus-defs "$LOCUS_DEFS" \
  --genome-db-dir data/ragtag_dbs \
  --output-dir "${FAMILY_DIR}/02_synteny_blocks" \
  --threads "$THREADS" \
  --max-gap-kb "$MAX_GAP_KB" \
  --max-span-kb "$MAX_SPAN_KB" \
  --closest-each-side "$CLOSEST_EACH_SIDE" \
  --closest-internal "$CLOSEST_INTERNAL" \
  --locus "$LOCUS_ID"
