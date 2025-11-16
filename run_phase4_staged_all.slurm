#!/bin/bash
#SBATCH --job-name=phase4_staged
#SBATCH --array=0-32
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --time=01:00:00
#SBATCH --output=logs/phase4_staged_%A_%a.out
#SBATCH --error=logs/phase4_staged_%A_%a.err

cd /carc/scratch/projects/emartins/2016456/adam/synteny_scanning/hpc_deployment_package/hybrid_workflow

# Read family list
mapfile -t families < inputs/all_families.txt
family=${families[$SLURM_ARRAY_TASK_ID]}

echo "=========================================="
echo "Phase 4 (Staged Calibration): $family"
echo "Array task: $SLURM_ARRAY_TASK_ID"
echo "Started: $(date)"
echo "=========================================="

# Activate environment
eval "$(micromamba shell hook --shell=bash)"
micromamba activate busco_env

# Check inputs
LOCUS_DEFS="outputs/$family/phase1_v2/locus_definitions.tsv"
PHASE1_DIR="outputs/$family/phase1_v2"

if [ ! -f "$LOCUS_DEFS" ]; then
    echo "ERROR: Locus definitions not found: $LOCUS_DEFS"
    exit 1
fi

if [ ! -d "$PHASE1_DIR" ]; then
    echo "ERROR: Phase 1 directory not found: $PHASE1_DIR"
    exit 1
fi

# Run Phase 4 with staged calibration
OUT_DIR="outputs/$family/phase4_staged"
mkdir -p "$OUT_DIR"

python redesign_scripts/04_detect_targets.py \
    --locus-defs "$LOCUS_DEFS" \
    --phase1-dir "$PHASE1_DIR" \
    --genome-db-dir data/ragtag_dbs \
    --output-dir "$OUT_DIR"

RC=$?
if [ $RC -ne 0 ]; then
    echo "ERROR: Phase 4 failed for $family (rc=$RC)"
    exit $RC
fi

# Compare with previous Phase 4
echo ""
echo "=========================================="
echo "Comparison with Previous Phase 4"
echo "=========================================="

# Count BK genes from different versions
v3_bk=$(grep -c "Belonocnema_kinseyi" "outputs/$family/phase4_v3/all_target_loci.tsv" 2>/dev/null || echo 0)
staged_bk=$(grep -c "Belonocnema_kinseyi" "outputs/$family/phase4_staged/all_target_loci.tsv" 2>/dev/null || echo 0)

# Count Phase 1 BK ground truth
phase1_bk=$(ls -d outputs/$family/phase1_v2/BK_* 2>/dev/null | wc -l)

echo "Phase 1 ground truth: $phase1_bk BK loci"
echo "Phase 4 v3 (fixed 0.5): $v3_bk BK genes"
echo "Phase 4 staged: $staged_bk BK genes"

if [ $staged_bk -gt $v3_bk ]; then
    improvement=$((staged_bk - v3_bk))
    echo "✓ Improvement: +$improvement genes"
elif [ $staged_bk -lt $v3_bk ]; then
    regression=$((v3_bk - staged_bk))
    echo "⚠️  Regression: -$regression genes"
else
    echo "→ No change"
fi

echo ""
echo "Finished: $(date)"
