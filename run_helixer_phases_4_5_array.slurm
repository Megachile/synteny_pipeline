#!/bin/bash
#SBATCH --job-name=helixer_p4_5
#SBATCH --output=logs/helixer_p4_5_%A_%a.out
#SBATCH --error=logs/helixer_p4_5_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --array=0-32

# Helixer Pipeline: Phase 4 (DIAMOND) + Phase 4b (flanking) + Phase 5 (classify)
# Runs with 30% query coverage filter to remove partial domain matches

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

echo "============================================================"
echo "HELIXER PIPELINE: Phase 4 + 4b + 5 for family: $FAMILY"
echo "  Array index: $SLURM_ARRAY_TASK_ID"
echo "  Date: $(date)"
echo "============================================================"

# Activate environment
eval "$(/carc/scratch/projects/emartins/2016456/micromamba_root/bin/micromamba shell hook -s bash)"
micromamba activate busco_env

# Paths
PHASE1_DIR="outputs/$FAMILY/phase1"
PHASE4_OUT="outputs/$FAMILY/phase4_helixer_filtered"
PHASE4B_OUT="outputs/$FAMILY/phase4b_helixer_filtered"
PHASE5_OUT="outputs/$FAMILY/phase5_helixer_filtered"
HELIXER_DIR="data/ragtag_output"
GENOME_LIST="helixer_genomes.tsv"

# Check inputs
if [ ! -f "${PHASE1_DIR}/locus_definitions.tsv" ]; then
    echo "ERROR: Missing Phase 1 locus definitions: ${PHASE1_DIR}/locus_definitions.tsv"
    exit 1
fi

# =========================================
# Phase 4: DIAMOND against Helixer proteomes
# =========================================
echo ""
echo "=== PHASE 4: DIAMOND blastp against Helixer proteomes ==="
echo "  Output: $PHASE4_OUT"
echo "  Min query coverage: 30% (filtering partial domain matches)"

python helixer_pipeline/04_detect_targets_helixer.py \
    --family "$FAMILY" \
    --phase1-dir "$PHASE1_DIR" \
    --helixer-dir "$HELIXER_DIR" \
    --genome-list "$GENOME_LIST" \
    --output-dir "$PHASE4_OUT" \
    --threads 8 \
    --min-query-coverage 0.3

if [ $? -ne 0 ]; then
    echo "ERROR: Phase 4 failed for $FAMILY"
    exit 1
fi

# =========================================
# Phase 4b: Extract flanking genes from GFF3
# =========================================
echo ""
echo "=== PHASE 4b: Extract flanking genes from GFF3 ==="
echo "  Output: $PHASE4B_OUT"

python helixer_pipeline/04b_extract_flanking_helixer.py \
    --family "$FAMILY" \
    --phase4-output "${PHASE4_OUT}/all_target_loci.tsv" \
    --helixer-dir "$HELIXER_DIR" \
    --output-dir "$PHASE4B_OUT" \
    --n-flanking 10

if [ $? -ne 0 ]; then
    echo "ERROR: Phase 4b failed for $FAMILY"
    exit 1
fi

# =========================================
# Phase 5: Classify by flanking gene synteny
# =========================================
echo ""
echo "=== PHASE 5: Classify targets by flanking gene synteny ==="
echo "  Output: $PHASE5_OUT"

python helixer_pipeline/05_classify_targets_helixer.py \
    --family "$FAMILY" \
    --phase4-dir "$PHASE4_OUT" \
    --phase4b-dir "$PHASE4B_OUT" \
    --phase1-dir "$PHASE1_DIR" \
    --output-dir "$PHASE5_OUT" \
    --threads 8

if [ $? -ne 0 ]; then
    echo "ERROR: Phase 5 failed for $FAMILY"
    exit 1
fi

echo ""
echo "============================================================"
echo "COMPLETE: Helixer pipeline for $FAMILY"
echo "  Phase 4: $PHASE4_OUT"
echo "  Phase 4b: $PHASE4B_OUT"
echo "  Phase 5: $PHASE5_OUT"
echo "  Date: $(date)"
echo "============================================================"
