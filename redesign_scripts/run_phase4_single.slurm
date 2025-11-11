#!/bin/bash
#SBATCH --job-name=phase4_targets_v2
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH --time=06:00:00
#SBATCH --partition=general
#SBATCH --output=logs/phase4_v2_%j.out
#SBATCH --error=logs/phase4_v2_%j.err

# Single-family Phase 4 (targets) submission for redesigned flow
# Usage:
#   sbatch redesign_scripts/run_phase4_single.slurm <FAMILY>

set -euo pipefail

cd /carc/scratch/projects/emartins/2016456/adam/synteny_scanning/hpc_deployment_package/hybrid_workflow

mkdir -p logs

if [ $# -lt 1 ]; then
  echo "Usage: $0 <FAMILY>" >&2
  exit 1
fi

FAMILY="$1"

# Activate environment (adjust if your cluster uses modules)
if command -v micromamba >/dev/null 2>&1; then
  eval "$(micromamba shell hook --shell bash)"
  micromamba activate busco_env || true
fi

PHASE1_DIR="outputs/${FAMILY}/phase1_v2"
PHASE4_DIR="outputs/${FAMILY}/phase4_v2"
LOCUS_DEFS="${PHASE1_DIR}/locus_definitions.tsv"
DB_DIR="data/ragtag_dbs"

echo "Running Phase 4 (v2) target detection for family: ${FAMILY}"
echo "  Locus defs: ${LOCUS_DEFS}"
echo "  Genome DBs: ${DB_DIR}"
echo "  Output dir: ${PHASE4_DIR}"

if [ ! -f "${LOCUS_DEFS}" ]; then
  echo "ERROR: Missing locus_definitions.tsv: ${LOCUS_DEFS}" >&2
  exit 1
fi

mkdir -p "${PHASE4_DIR}"

python redesign_scripts/04_detect_targets.py \
  --locus-defs "${LOCUS_DEFS}" \
  --phase1-dir "${PHASE1_DIR}" \
  --genome-db-dir "${DB_DIR}" \
  --output-dir "${PHASE4_DIR}" \
  --evalue 1e-5 \
  --max-targets 10000 \
  --threads ${SLURM_CPUS_PER_TASK:-16} \
  --base-min-split-kb 40

RC=$?
if [ ${RC} -ne 0 ]; then
  echo "ERROR: Phase 4 (v2) failed for ${FAMILY} (rc=${RC})" >&2
  exit ${RC}
fi

echo "Phase 4 (v2) complete for family: ${FAMILY}"

