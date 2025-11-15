#!/bin/bash
#SBATCH --job-name=phase5_v2_single
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --partition=general
#SBATCH --output=logs/phase5_v2_%j.out
#SBATCH --error=logs/phase5_v2_%j.err

# Usage: sbatch redesign_scripts/run_phase5_single.slurm <FAMILY> [NEARBY_UPPER_KB]

set -euo pipefail

cd /carc/scratch/projects/emartins/2016456/adam/synteny_scanning/hpc_deployment_package/hybrid_workflow

if [ $# -lt 1 ]; then
  echo "Usage: $0 <FAMILY> [NEARBY_UPPER_KB]" >&2
  exit 1
fi

FAMILY="$1"
NEARBY=${2:-350}

if command -v micromamba >/dev/null 2>&1; then
  eval "$(micromamba shell hook --shell bash)"
  micromamba activate busco_env || true
fi

IN_T=outputs/${FAMILY}/phase4_v2/all_target_loci.tsv
IN_B=outputs/${FAMILY}/phase3_v2/all_synteny_blocks_ge3.tsv
OUT=outputs/${FAMILY}/phase5_v2
mkdir -p "$OUT"

echo "Phase 5 (test) for ${FAMILY} with nearby_upper_kb=${NEARBY}"
python redesign_scripts/05_classify_targets.py \
  --targets "$IN_T" \
  --blocks "$IN_B" \
  --output-dir "$OUT" \
  --pad-kb 20 \
  --nearby-upper-kb ${NEARBY} \
  --rescue-same-scaffold \
  --rescue-upper-kb 350 \
  --rescue-min-evalue 1e-50

# Create descriptive symlink for reference
SYMLINK=outputs/${FAMILY}/phase5_v2_nearby${NEARBY}
if [ ! -e "$SYMLINK" ]; then
  ln -s phase5_v2 "$SYMLINK"
fi

echo "Done Phase 5 test for ${FAMILY}"
