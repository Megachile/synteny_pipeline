#!/bin/bash
#SBATCH --job-name=phase8_v2
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --time=02:00:00
#SBATCH --partition=general
#SBATCH --output=logs/phase8_v2_%j.out
#SBATCH --error=logs/phase8_v2_%j.err

# Usage: sbatch redesign_scripts/run_phase8_single.slurm <FAMILY>

set -euo pipefail

cd /carc/scratch/projects/emartins/2016456/adam/synteny_scanning/hpc_deployment_package/hybrid_workflow

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

echo "Phase 8a (locus matrices) for ${FAMILY}"
python redesign_scripts/08a_generate_locus_matrices.py --family "${FAMILY}" --species-map data/gca_to_species.tsv --reference-proteins data/reference/protein.faa

echo "Phase 8b (summary matrices) for ${FAMILY}"
python redesign_scripts/08b_generate_summary_matrices.py --family "${FAMILY}" --species-map data/gca_to_species.tsv

echo "Done Phase 8 for ${FAMILY}"

