#!/bin/bash
#SBATCH --job-name=helixer_p6_8
#SBATCH --output=logs/helixer_p6_8_%A_%a.out
#SBATCH --error=logs/helixer_p6_8_%A_%a.err
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --array=0-32

# Helixer Pipeline: Phase 6 + Phase 7 (SwissProt) + Phase 8a/8b
#
# This produces outputs matching the original tblastn pipeline format:
#   Phase 6:  outputs/{family}/phase6_helixer_filtered/{genome}/{locus}/gene1_protein.fasta
#   Phase 7:  outputs/{family}/phase7_helixer_filtered/flanking_matches_annotated.tsv
#   Phase 8a: outputs/{family}/phase8_helixer_filtered/{locus}_genome_swissprot_matrix.tsv
#   Phase 8b: outputs/{family}/phase8_helixer_filtered/{family}_summary_matrix.tsv

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
echo "HELIXER PIPELINE: Phase 6 + 8a + 8b for family: $FAMILY"
echo "  Array index: $SLURM_ARRAY_TASK_ID"
echo "  Date: $(date)"
echo "============================================================"

# Activate environment
eval "$(/carc/scratch/projects/emartins/2016456/micromamba_root/bin/micromamba shell hook -s bash)"
micromamba activate busco_env

# Paths
PHASE1_DIR="outputs/$FAMILY/phase1_v2"
PHASE4_DIR="outputs/$FAMILY/phase4_helixer_filtered"
PHASE4B_DIR="outputs/$FAMILY/phase4b_helixer_filtered"
PHASE5_DIR="outputs/$FAMILY/phase5_helixer_filtered"
PHASE6_OUT="outputs/$FAMILY/phase6_helixer_filtered"
PHASE7_OUT="outputs/$FAMILY/phase7_helixer_filtered"
PHASE8_OUT="outputs/$FAMILY/phase8_helixer_filtered"
HELIXER_DIR="data/ragtag_output"
SPECIES_MAP="data/gca_to_species.tsv"
SWISSPROT_DB="/carc/scratch/projects/emartins/2016456/adam/databases/uniprot_sprot.dmnd"

# Check if Phase 5 outputs exist (skip families with no targets)
if [ ! -f "${PHASE5_DIR}/syntenic_targets.tsv" ] && [ ! -f "${PHASE5_DIR}/unplaceable_targets.tsv" ]; then
    echo "SKIP: No Phase 5 outputs for $FAMILY (no Helixer targets)"
    exit 0
fi

# =========================================
# Phase 6: Extract and organize sequences
# =========================================
echo ""
echo "=== PHASE 6: Extract and organize target sequences ==="
echo "  Output: $PHASE6_OUT"

python helixer_pipeline/06_extract_sequences_helixer.py \
    --family "$FAMILY" \
    --phase4-dir "$PHASE4_DIR" \
    --phase5-dir "$PHASE5_DIR" \
    --helixer-dir "$HELIXER_DIR" \
    --output-dir "$PHASE6_OUT"

if [ $? -ne 0 ]; then
    echo "ERROR: Phase 6 failed for $FAMILY"
    exit 1
fi

# =========================================
# Phase 7: SwissProt annotation for flanking genes
# =========================================
echo ""
echo "=== PHASE 7: Annotate flanking genes with SwissProt ==="
echo "  Output: $PHASE7_OUT"

python helixer_pipeline/07_annotate_flanking_swissprot.py \
    --family "$FAMILY" \
    --phase4b-dir "$PHASE4B_DIR" \
    --phase5-dir "$PHASE5_DIR" \
    --swissprot-db "$SWISSPROT_DB" \
    --output-dir "$PHASE7_OUT" \
    --threads 4

if [ $? -ne 0 ]; then
    echo "ERROR: Phase 7 failed for $FAMILY"
    exit 1
fi

# =========================================
# Phase 8a: Generate locus matrices
# =========================================
echo ""
echo "=== PHASE 8a: Generate locus-specific matrices ==="
echo "  Output: $PHASE8_OUT"

python helixer_pipeline/08a_generate_locus_matrices_helixer.py \
    --family "$FAMILY" \
    --phase1-dir "$PHASE1_DIR" \
    --phase5-dir "$PHASE5_DIR" \
    --phase6-dir "$PHASE6_OUT" \
    --phase7-dir "$PHASE7_OUT" \
    --species-map "$SPECIES_MAP" \
    --output-dir "$PHASE8_OUT"

if [ $? -ne 0 ]; then
    echo "ERROR: Phase 8a failed for $FAMILY"
    exit 1
fi

# =========================================
# Phase 8b: Generate summary matrix
# =========================================
echo ""
echo "=== PHASE 8b: Generate summary matrix ==="
echo "  Output: $PHASE8_OUT"

python helixer_pipeline/08b_generate_summary_matrices_helixer.py \
    --family "$FAMILY" \
    --phase1-dir "$PHASE1_DIR" \
    --phase5-dir "$PHASE5_DIR" \
    --phase6-dir "$PHASE6_OUT" \
    --phase8a-dir "$PHASE8_OUT" \
    --species-map "$SPECIES_MAP" \
    --output-dir "$PHASE8_OUT"

if [ $? -ne 0 ]; then
    echo "ERROR: Phase 8b failed for $FAMILY"
    exit 1
fi

echo ""
echo "============================================================"
echo "COMPLETE: Helixer pipeline Phase 6-8 for $FAMILY"
echo "  Phase 6: $PHASE6_OUT"
echo "  Phase 7: $PHASE7_OUT"
echo "  Phase 8: $PHASE8_OUT"
echo "  Date: $(date)"
echo "============================================================"
