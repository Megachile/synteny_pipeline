#!/bin/bash
#SBATCH --job-name=qc_exonerate
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --partition=general
#SBATCH --output=logs/qc_exonerate_%j.out
#SBATCH --error=logs/qc_exonerate_%j.err

# Usage: sbatch redesign_scripts/run_qc_single.slurm <FAMILY>

set -euo pipefail

cd /carc/scratch/projects/emartins/2016456/adam/synteny_scanning/hpc_deployment_package/hybrid_workflow

# Activate environment (align with other phases)
if command -v micromamba >/dev/null 2>&1; then
  eval "$(micromamba shell hook --shell bash)"
  micromamba activate busco_env || true
fi

if [ $# -lt 1 ]; then
  echo "Usage: $0 <FAMILY>" >&2
  exit 1
fi

FAMILY="$1"

echo "QC: Phase 5 synteny-only checks for ${FAMILY}"
python scripts/qc_phase5_synteny.py --family "${FAMILY}"

echo "QC: Exonerate candidates vs selected blocks for ${FAMILY}"
python scripts/exonerate_qc.py --family "${FAMILY}"

echo "QC done for ${FAMILY}"
