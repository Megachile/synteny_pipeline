#!/bin/bash
#SBATCH --job-name=phase4_recluster
#SBATCH --output=logs/phase4_recluster_%A_%a.out
#SBATCH --error=logs/phase4_recluster_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --array=0-32

# Re-run Phase 4 clustering with locus-agnostic approach
# Reuses existing BLAST XML files from phase4_v2 - no tBLASTn needed

cd /carc/scratch/projects/emartins/2016456/adam/synteny_scanning/hpc_deployment_package/hybrid_workflow

# Family list (33 families)
FAMILIES=(
    "alaserpin_MC35"
    "ApoD_MC13"
    "beta_lectin_MC147"
    "CAP_MC28"
    "carboxypeptidase_MC26"
    "der_f_21_MC2"
    "endoglucanase_Z_MC6"
    "ester_hydrolase_MC19"
    "extensin-like_MC192"
    "ferritin_MC102"
    "Glutenin_MC21"
    "LRR_15_MC91"
    "LRR_4_MC10"
    "MC117"
    "MC174"
    "MC18"
    "MC211"
    "MC258"
    "MC3"
    "natterin-4-like_MC59"
    "neurexin_MC14"
    "OS-D_MC11"
    "pectin_lyase_PL"
    "peptidoglycan-recognition_protein_SA-like___peptidoglycan-recognition_protein_LC-like_MC125"
    "rhamnogalacturonate_lyase-like_MC47"
    "ribonuclease_MC33"
    "sprouty_MC112"
    "thioredoxin-2-like_MC108"
    "TIL_MC9"
    "venom_acid_phosphatase_Acph-1-like_transcript_variant_X1_MC63"
    "venom_R-like_MC1"
    "venom_serine_carboxypeptidase_MC203"
    "venom_serine_protease_34-like_MC175"
)

FAMILY=${FAMILIES[$SLURM_ARRAY_TASK_ID]}

echo "Re-clustering Phase 4 (locus-agnostic) for family: $FAMILY"
echo "  Array index: $SLURM_ARRAY_TASK_ID"
echo "  Output dir: outputs/$FAMILY/phase4_v2 (reusing existing BLAST XML)"

# Activate environment
eval "$(/carc/scratch/projects/emartins/2016456/micromamba_root/bin/micromamba shell hook -s bash)"
micromamba activate busco_env

# Run Phase 4 - will reuse existing XML files in phase4_v2/blast_xml/
python redesign_scripts/04_detect_targets.py \
    --locus-defs "outputs/$FAMILY/phase1_v2/locus_definitions.tsv" \
    --genome-db-dir data/ragtag_dbs \
    --output-dir "outputs/$FAMILY/phase4_v2" \
    --phase1-dir "outputs/$FAMILY/phase1_v2"

echo "Phase 4 re-clustering complete for family: $FAMILY"
