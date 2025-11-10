#!/bin/bash
#SBATCH --job-name=phase2_synteny
#SBATCH --array=0-31
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH --time=24:00:00
#SBATCH --partition=general
#SBATCH --output=logs/phase2_v2_%A_%a.out
#SBATCH --error=logs/phase2_v2_%A_%a.err

cd /carc/scratch/projects/emartins/2016456/adam/synteny_scanning/hpc_deployment_package/hybrid_workflow

# Activate environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate busco_env

# Get list of families from outputs directory
families=($(ls -d outputs/*/phase1/locus_definitions.tsv | cut -d'/' -f2 | sort))

# Get family for this array task
family=${families[$SLURM_ARRAY_TASK_ID]}

echo "Processing family: $family (task $SLURM_ARRAY_TASK_ID)"

# Paths
locus_defs="outputs/${family}/phase1/locus_definitions.tsv"
phase2_output="outputs/${family}/phase2_synteny"

# Check if input exists
if [ ! -f "$locus_defs" ]; then
    echo "ERROR: Input file not found: $locus_defs"
    exit 1
fi

# Create output directory
mkdir -p "$phase2_output"

# Run synteny detection
python scripts/synteny_detection.py \
    --locus-defs "$locus_defs" \
    --genome-db-dir data/ragtag_dbs \
    --output-dir "$phase2_output" \
    --evalue 1e-5 \
    --max-targets 50 \
    --threads 16

echo "Phase 2 complete for $family"