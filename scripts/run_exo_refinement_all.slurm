#!/bin/bash
#SBATCH --job-name=synteny_06c6b_all
#SBATCH --output=logs/06c6b_all_%j.out
#SBATCH --error=logs/06c6b_all_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --partition=general

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WF_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "$WF_DIR"

mkdir -p logs

echo "Node: $(hostname)"
echo "Started: $(date)"

source ~/.bashrc
eval "$(micromamba shell hook --shell bash)"
micromamba activate busco_env

python scripts/run_exo_refinement_all.py

echo "Finished: $(date)"
