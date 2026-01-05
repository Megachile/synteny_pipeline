#!/bin/bash
#SBATCH --job-name=helixer_phase1
#SBATCH --output=logs/helixer_phase1_%A_%a.out
#SBATCH --error=logs/helixer_phase1_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --array=0-32

################################################################################
# Helixer Pipeline: Phase 1 - Locus Discovery
#
# Discovers BK/LB reference loci for each gene family using DIAMOND against
# the landmark genomes (BK and LB). This is the prerequisite step before
# running Phase 4-8 which search Helixer-annotated genomes.
#
# Phase 1: Discover paralogs → build locus definitions
# Phase 1.5: Merge BK-LB orthologs into unified locus space
################################################################################

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
echo "HELIXER PIPELINE: Phase 1 (Locus Discovery) for family: $FAMILY"
echo "  Array index: $SLURM_ARRAY_TASK_ID"
echo "  Date: $(date)"
echo "============================================================"

# Activate environment
eval "$(/carc/scratch/projects/emartins/2016456/micromamba_root/bin/micromamba shell hook -s bash)"
micromamba activate busco_env

# Paths
INPUT_FILE="inputs/${FAMILY}_test.tsv"
if [ ! -f "$INPUT_FILE" ]; then
    INPUT_FILE=$(ls inputs/${FAMILY}*.tsv 2>/dev/null | grep -v "locus_definitions" | head -1)
fi

if [ -z "$INPUT_FILE" ] || [ ! -f "$INPUT_FILE" ]; then
    echo "ERROR: No input file found for $FAMILY"
    exit 1
fi

echo "Using input file: $INPUT_FILE"

# Extract LOC IDs from input file
LOC_IDS=$(tail -n +2 "$INPUT_FILE" | cut -f2 | tr '\n' ',' | sed 's/,$//')
if [ -z "$LOC_IDS" ]; then
    echo "ERROR: No LOC IDs found in $INPUT_FILE"
    exit 1
fi
echo "LOC IDs: $LOC_IDS"

# Output directories
PHASE1_DIR="outputs/$FAMILY/phase1"
GENOME_GFF_DIR="data/reference"
LB_PROTEOME_DB="data/proteomes/GCF_019393585.1_Lboulardi_NCBI.dmnd"

mkdir -p "$PHASE1_DIR"
mkdir -p logs

# =========================================
# Phase 1: Locus Discovery
# =========================================
echo ""
echo "=== PHASE 1: Locus Discovery ==="
echo "  Output: $PHASE1_DIR"

python helixer_pipeline/01_phase1.py \
    --loc-ids "$LOC_IDS" \
    --gene-family "$FAMILY" \
    --output-dir "$PHASE1_DIR" \
    --genome-gff-dir "$GENOME_GFF_DIR"

if [ $? -ne 0 ]; then
    echo "ERROR: Phase 1 failed for $FAMILY"
    exit 1
fi

if [ ! -f "${PHASE1_DIR}/locus_definitions.tsv" ]; then
    echo "ERROR: Phase 1 failed - locus_definitions.tsv not found"
    exit 1
fi

NUM_LOCI=$(tail -n +2 "${PHASE1_DIR}/locus_definitions.tsv" | wc -l)
echo "Phase 1 complete: ${NUM_LOCI} loci identified"

# =========================================
# Phase 1.5: BK-LB Locus Unification
# =========================================
if [ -f "$LB_PROTEOME_DB" ]; then
    echo ""
    echo "=== PHASE 1.5: BK-LB Locus Unification ==="

    python helixer_pipeline/01b_merge_bk_lb_loci.py \
        --locus-defs "${PHASE1_DIR}/locus_definitions.tsv" \
        --lb-proteome-db "$LB_PROTEOME_DB" \
        --diamond-threads 4

    if [ $? -ne 0 ]; then
        echo "WARNING: Phase 1.5 failed, continuing with original loci"
    else
        echo "Phase 1.5 complete: BK-LB loci unified"
    fi
else
    echo "SKIP Phase 1.5: LB proteome database not found"
fi

echo ""
echo "============================================================"
echo "COMPLETE: Phase 1 for $FAMILY"
echo "  Output: $PHASE1_DIR"
echo "  Loci: ${NUM_LOCI}"
echo "  Date: $(date)"
echo "============================================================"
