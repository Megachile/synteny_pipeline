#!/bin/bash
#SBATCH --job-name=phase4_adaptive_all
#SBATCH --output=phase4_adaptive_%A_%a.out
#SBATCH --error=phase4_adaptive_%A_%a.err
#SBATCH --time=01:30:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=16
#SBATCH --array=0-34

cd /carc/scratch/projects/emartins/2016456/adam/synteny_scanning/hpc_deployment_package/hybrid_workflow

# All 35 families
families=(
    "ApoD_MC13"
    "CAP_MC28"
    "Glutenin_MC21"
    "LRR_15_MC91"
    "LRR_4_MC10"
    "MC117"
    "MC174"
    "MC18"
    "MC211"
    "MC258"
    "MC3"
    "MC6"
    "MC6_pipeline_test"
    "OS-D_MC11"
    "TIL_MC9"
    "alaserpin_MC35"
    "carboxypeptidase_MC26"
    "der_f_21_MC2"
    "endoglucanase_Z_MC6"
    "ester_hydrolase_MC19"
    "extensin-like_MC192"
    "ferritin_MC102"
    "ferritin_test_fresh"
    "natterin-4-like_MC59"
    "neurexin_MC14"
    "pectin_lyase_PL"
    "peptidoglycan-recognition_protein_SA-like___peptidoglycan-recognition_protein_LC-like_MC125"
    "rhamnogalacturonate_lyase-like_MC47"
    "ribonuclease_MC33"
    "sprouty_MC112"
    "thioredoxin-2-like_MC108"
    "venom_R-like_MC1"
    "venom_acid_phosphatase_Acph-1-like_transcript_variant_X1_MC63"
    "venom_serine_carboxypeptidase_MC203"
    "venom_serine_protease_34-like_MC175"
)

FAMILY="${families[$SLURM_ARRAY_TASK_ID]}"

echo "=========================================="
echo "Phase 4 Adaptive Calibration: $FAMILY"
echo "Array task: $SLURM_ARRAY_TASK_ID"
echo "=========================================="
echo ""

# Activate environment
source /users/akranz8174/.bashrc
micromamba activate busco_env

# Copy existing BLAST XMLs from phase4_v3 (don't re-run BLAST!)
mkdir -p "outputs/$FAMILY/phase4_adaptive"
if [ -d "outputs/$FAMILY/phase4_v3/blast_xml" ]; then
    echo "Copying existing BLAST XMLs from phase4_v3..."
    cp -r "outputs/$FAMILY/phase4_v3/blast_xml" "outputs/$FAMILY/phase4_adaptive/"
    echo "  $(ls outputs/$FAMILY/phase4_adaptive/blast_xml/*.xml 2>/dev/null | wc -l) BLAST XMLs copied"
else
    echo "WARNING: No phase4_v3/blast_xml found for $FAMILY"
fi

echo ""

# Run Phase 4 with adaptive calibration
python redesign_scripts/04_detect_targets.py \
    --locus-defs "outputs/$FAMILY/phase1_v2/locus_definitions.tsv" \
    --phase1-dir "outputs/$FAMILY/phase1_v2" \
    --genome-db-dir "data/ragtag_dbs" \
    --output-dir "outputs/$FAMILY/phase4_adaptive" \
    --threads 16

echo ""
echo "=========================================="
echo "RESULTS FOR $FAMILY"
echo "=========================================="
echo "Phase 1 BK targets: $(for dir in outputs/$FAMILY/phase1_v2/BK_*; do [ -d "$dir" ] && grep -c "^>" $dir/*_targets.faa 2>/dev/null; done | awk '{sum+=$1} END {print sum}')"
echo "Phase 4 v3 (thresh 0.2): $(grep -c "Belonocnema_kinseyi" "outputs/$FAMILY/phase4_v3/all_target_loci.tsv" 2>/dev/null || echo 0)"
echo "Phase 4 adaptive: $(grep -c "Belonocnema_kinseyi" "outputs/$FAMILY/phase4_adaptive/all_target_loci.tsv" 2>/dev/null || echo 0)"
