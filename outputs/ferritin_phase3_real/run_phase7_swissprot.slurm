#!/bin/bash
#SBATCH --job-name=ferritin_phase7
#SBATCH --output=phase7_swissprot_%j.out
#SBATCH --error=phase7_swissprot_%j.err
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --partition=general

# CRITICAL: Change to project directory first
cd /carc/scratch/projects/emartins/2016456/adam/synteny_scanning/hpc_deployment_package/hybrid_workflow

# Activate environment (diamond is already available in PATH)
source ~/.bashrc
micromamba activate busco_env

# Run Phase 7: SwissProt annotation of flanking proteins
python scripts/07_swissprot_annotation.py \
    --synteny-dir outputs/ferritin_phase3_real/02_synteny_blocks \
    --swissprot-db data/reference/swissprot.dmnd \
    --output-dir outputs/ferritin_phase3_real/07_swissprot_annotations \
    --evalue 1e-5 \
    --threads 16

echo "Phase 7 SwissProt annotation complete!"
