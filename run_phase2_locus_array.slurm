#!/bin/bash
#SBATCH --job-name=phase2_locus
#SBATCH --array=1-293
#SBATCH --cpus-per-task=16
#SBATCH --mem=20G
#SBATCH --time=03:00:00
#SBATCH --partition=general
#SBATCH --output=logs/phase2_locus_%A_%a.out
#SBATCH --error=logs/phase2_locus_%A_%a.err

cd /carc/scratch/projects/emartins/2016456/adam/synteny_scanning/hpc_deployment_package/hybrid_workflow

# Activate environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate busco_env

# Get locus info for this array task (line number = array task ID)
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" all_loci.txt)
family=$(echo "$line" | cut -f1)
locus_id=$(echo "$line" | cut -f2)
locus_defs=$(echo "$line" | cut -f3)

echo "Processing locus: $locus_id from family: $family (task $SLURM_ARRAY_TASK_ID)"

# Output directory for this family
phase2_output="outputs/${family}/phase2_synteny"

# Create output directory
mkdir -p "$phase2_output"

# Run synteny detection FOR THIS SPECIFIC LOCUS
python scripts/synteny_detection.py \
    --locus-defs "$locus_defs" \
    --genome-db-dir data/ragtag_dbs \
    --output-dir "$phase2_output" \
    --locus "$locus_id" \
    --evalue 1e-5 \
    --max-targets 50 \
    --threads 16

echo "Phase 2 complete for locus $locus_id"