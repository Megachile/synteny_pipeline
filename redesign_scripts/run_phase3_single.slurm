#!/bin/bash
#SBATCH --job-name=phase3_v2
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --partition=general
#SBATCH --output=logs/phase3_v2_%j.out
#SBATCH --error=logs/phase3_v2_%j.err

# Usage: sbatch redesign_scripts/run_phase3_single.slurm <FAMILY> [MIN_PROTEINS]

set -euo pipefail

cd /carc/scratch/projects/emartins/2016456/adam/synteny_scanning/hpc_deployment_package/hybrid_workflow

if [ $# -lt 1 ]; then
  echo "Usage: $0 <FAMILY> [MIN_PROTEINS]" >&2
  exit 1
fi

FAMILY="$1"
MIN_PROTEINS=${2:-3}

if command -v micromamba >/dev/null 2>&1; then
  eval "$(micromamba shell hook --shell bash)"
  micromamba activate busco_env || true
fi

SYNTENY_DIR=outputs/${FAMILY}/phase2_synteny_v2
OUT_DIR=outputs/${FAMILY}/phase3_v2

echo "Phase 3 (aggregate/filter) for ${FAMILY} with min_proteins=${MIN_PROTEINS}"
python redesign_scripts/03_aggregate_threshold_only.py \
  --synteny-dir "${SYNTENY_DIR}" \
  --output-dir "${OUT_DIR}" \
  --min-proteins ${MIN_PROTEINS}

echo "Done Phase 3 for ${FAMILY}"
