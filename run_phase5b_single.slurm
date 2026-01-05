#!/bin/bash
#SBATCH --job-name=helixer_p5b
#SBATCH --output=logs/helixer_p5b_%j.out
#SBATCH --error=logs/helixer_p5b_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8

# Phase 5b: Validate novel loci candidates
# Usage: sbatch run_phase5b_single.slurm FAMILY_NAME

cd /carc/scratch/projects/emartins/2016456/adam/synteny_scanning/hpc_deployment_package/hybrid_workflow

FAMILY=${1:-"CAP_MC28"}

echo "============================================================"
echo "HELIXER PIPELINE: Phase 5b for family: $FAMILY"
echo "  Date: $(date)"
echo "============================================================"

eval "$(/carc/scratch/projects/emartins/2016456/micromamba_root/bin/micromamba shell hook -s bash)"
micromamba activate busco_env

PHASE1_DIR="outputs/$FAMILY/phase1"
PHASE5_DIR="outputs/$FAMILY/phase5_helixer_filtered"
PHASE5B_OUT="outputs/$FAMILY/phase5b_helixer"
HELIXER_DIR="data/ragtag_output"

# Check for symlink and use actual directory
if [ -L "$PHASE5_DIR" ]; then
    PHASE5_DIR=$(readlink -f "$PHASE5_DIR")
fi

if [ ! -f "${PHASE5_DIR}/unplaceable_targets.tsv" ]; then
    echo "ERROR: Missing Phase 5 unplaceable targets"
    exit 1
fi

mkdir -p "$PHASE5B_OUT"

python helixer_pipeline/05b_validate_novel_loci.py \
    --family "$FAMILY" \
    --phase1-dir "$PHASE1_DIR" \
    --phase5-dir "$PHASE5_DIR" \
    --helixer-dir "$HELIXER_DIR" \
    --output-dir "$PHASE5B_OUT" \
    --threads 8

echo ""
echo "============================================================"
echo "COMPLETE: Phase 5b for $FAMILY"
echo "  Output: $PHASE5B_OUT"
echo "  Date: $(date)"
echo "============================================================"
