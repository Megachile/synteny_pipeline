#!/bin/bash
#SBATCH --job-name=helixer_full
#SBATCH --output=logs/helixer_full_%A_%a.out
#SBATCH --error=logs/helixer_full_%A_%a.err
#SBATCH --time=03:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --array=0-32

################################################################################
# Helixer Pipeline: Complete single-job execution
#
# Runs ALL phases for one family:
#   Phase 1:   Locus discovery
#   Phase 1.5: BK-LB unification
#   Phase 2b:  Flanking-first synteny detection
#   Phase 4:   Target detection (DIAMOND)
#   Phase 4b:  Extract flanking genes
#   Phase 5:   Classify targets
#   Phase 5b:  Validate novel loci
#   Phase 6:   Extract sequences
#   Phase 7:   SwissProt annotation
#   Phase 8a:  Locus matrices
#   Phase 8b:  Summary matrix
#   Phase 9:   Novel loci detection
################################################################################

set -e  # Exit on error

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
echo "HELIXER PIPELINE: FULL RUN for family: $FAMILY"
echo "  Array index: $SLURM_ARRAY_TASK_ID"
echo "  Date: $(date)"
echo "============================================================"

# Activate environment
eval "$(/carc/scratch/projects/emartins/2016456/micromamba_root/bin/micromamba shell hook -s bash)"
micromamba activate busco_env

# Common paths
HELIXER_DIR="data/ragtag_output"
GENOME_LIST="helixer_genomes.tsv"
GENOME_GFF_DIR="data/reference"
SPECIES_MAP="data/gca_to_species.tsv"
SWISSPROT_DB="/carc/scratch/projects/emartins/2016456/adam/databases/uniprot_sprot.dmnd"
LB_PROTEOME_DB="data/proteomes/GCF_019393585.1_Lboulardi_NCBI.dmnd"
BK_PROTEOME_DB="data/proteomes/GCF_010883055.1_bkinseyi.dmnd"

# Output directories
PHASE1_DIR="outputs/$FAMILY/phase1"
PHASE2B_DIR="outputs/$FAMILY/phase2b_helixer"
PHASE4_DIR="outputs/$FAMILY/phase4_helixer_filtered"
PHASE4B_DIR="outputs/$FAMILY/phase4b_helixer_filtered"
PHASE5_DIR="outputs/$FAMILY/phase5_helixer_filtered"
PHASE5B_DIR="outputs/$FAMILY/phase5b_helixer"
PHASE6_DIR="outputs/$FAMILY/phase6_helixer_filtered"
PHASE7_DIR="outputs/$FAMILY/phase7_helixer_filtered"
PHASE8_DIR="outputs/$FAMILY/phase8_helixer_filtered"
PHASE9_DIR="outputs/$FAMILY/phase9_novel_loci"

mkdir -p logs

################################################################################
# PHASE 1: Locus Discovery
################################################################################
echo ""
echo "=== PHASE 1: Locus Discovery ==="

INPUT_FILE="inputs/${FAMILY}_test.tsv"
if [ ! -f "$INPUT_FILE" ]; then
    INPUT_FILE=$(ls inputs/${FAMILY}*.tsv 2>/dev/null | grep -v "locus_definitions" | head -1)
fi

if [ -z "$INPUT_FILE" ] || [ ! -f "$INPUT_FILE" ]; then
    echo "ERROR: No input file found for $FAMILY"
    exit 1
fi

LOC_IDS=$(tail -n +2 "$INPUT_FILE" | cut -f2 | tr '\n' ',' | sed 's/,$//')
echo "  Input: $INPUT_FILE"
echo "  LOC IDs: $LOC_IDS"

python helixer_pipeline/01_phase1.py \
    --loc-ids "$LOC_IDS" \
    --gene-family "$FAMILY" \
    --output-dir "$PHASE1_DIR" \
    --genome-gff-dir "$GENOME_GFF_DIR"

################################################################################
# PHASE 1.5: BK-LB Unification
################################################################################
echo ""
echo "=== PHASE 1.5: BK-LB Unification ==="

if [ -f "$LB_PROTEOME_DB" ]; then
    python helixer_pipeline/01b_merge_bk_lb_loci.py \
        --locus-defs "${PHASE1_DIR}/locus_definitions.tsv" \
        --lb-proteome-db "$LB_PROTEOME_DB" \
        --diamond-threads 8 || echo "  WARNING: Phase 1.5 failed, continuing"
else
    echo "  SKIP: LB proteome not found"
fi

################################################################################
# PHASE 4: Target Detection
################################################################################
echo ""
echo "=== PHASE 4: Target Detection (DIAMOND) ==="

python helixer_pipeline/04_detect_targets_helixer.py \
    --family "$FAMILY" \
    --phase1-dir "$PHASE1_DIR" \
    --helixer-dir "$HELIXER_DIR" \
    --genome-list "$GENOME_LIST" \
    --output-dir "$PHASE4_DIR" \
    --threads 8 \
    --min-query-coverage 0.3

################################################################################
# PHASE 4b: Extract Flanking
################################################################################
echo ""
echo "=== PHASE 4b: Extract Flanking Genes ==="

python helixer_pipeline/04b_extract_flanking_helixer.py \
    --family "$FAMILY" \
    --phase4-output "${PHASE4_DIR}/all_target_loci.tsv" \
    --helixer-dir "$HELIXER_DIR" \
    --output-dir "$PHASE4B_DIR" \
    --n-flanking 10

################################################################################
# PHASE 2b: Flanking-First Synteny (can run after Phase 4)
################################################################################
echo ""
echo "=== PHASE 2b: Flanking-First Synteny Detection ==="

python helixer_pipeline/02b_detect_empty_blocks_helixer.py \
    --family "$FAMILY" \
    --phase1-dir "$PHASE1_DIR" \
    --phase4-dir "$PHASE4_DIR" \
    --helixer-dir "$HELIXER_DIR" \
    --output-dir "$PHASE2B_DIR" \
    --threads 8 || echo "  WARNING: Phase 2b failed, continuing"

################################################################################
# PHASE 5: Classify Targets
################################################################################
echo ""
echo "=== PHASE 5: Classify Targets ==="

python helixer_pipeline/05_classify_targets_helixer.py \
    --family "$FAMILY" \
    --phase4-dir "$PHASE4_DIR" \
    --phase4b-dir "$PHASE4B_DIR" \
    --phase1-dir "$PHASE1_DIR" \
    --output-dir "$PHASE5_DIR" \
    --threads 8

# Check if we have any targets
if [ ! -f "${PHASE5_DIR}/syntenic_targets.tsv" ] && [ ! -f "${PHASE5_DIR}/unplaceable_targets.tsv" ]; then
    echo "SKIP: No targets found for $FAMILY"
    exit 0
fi

################################################################################
# PHASE 5b: Validate Novel Loci
################################################################################
echo ""
echo "=== PHASE 5b: Validate Novel Loci ==="

if [ -f "${PHASE5_DIR}/unplaceable_targets.tsv" ]; then
    python helixer_pipeline/05b_validate_novel_loci.py \
        --family "$FAMILY" \
        --phase5-dir "$PHASE5_DIR" \
        --phase4b-dir "$PHASE4B_DIR" \
        --helixer-dir "$HELIXER_DIR" \
        --bk-proteome "$BK_PROTEOME_DB" \
        --output-dir "$PHASE5B_DIR" \
        --threads 8 || echo "  WARNING: Phase 5b failed, continuing"
else
    echo "  SKIP: No unplaceable targets"
fi

################################################################################
# PHASE 6: Extract Sequences
################################################################################
echo ""
echo "=== PHASE 6: Extract Sequences ==="

python helixer_pipeline/06_extract_sequences_helixer.py \
    --family "$FAMILY" \
    --phase1-dir "$PHASE1_DIR" \
    --phase4-dir "$PHASE4_DIR" \
    --phase5-dir "$PHASE5_DIR" \
    --helixer-dir "$HELIXER_DIR" \
    --output-dir "$PHASE6_DIR"

################################################################################
# PHASE 7: SwissProt Annotation
################################################################################
echo ""
echo "=== PHASE 7: SwissProt Annotation ==="

PHASE2B_ARG=""
if [ -d "$PHASE2B_DIR" ]; then
    PHASE2B_ARG="--phase2b-dir $PHASE2B_DIR"
fi

python helixer_pipeline/07_annotate_flanking_swissprot.py \
    --family "$FAMILY" \
    --phase4b-dir "$PHASE4B_DIR" \
    --phase5-dir "$PHASE5_DIR" \
    $PHASE2B_ARG \
    --swissprot-db "$SWISSPROT_DB" \
    --output-dir "$PHASE7_DIR" \
    --threads 8

# Also annotate BK flanking proteins
BK_SWISSPROT_DIR="${PHASE1_DIR}/bk_swissprot"
BK_SWISSPROT_FILE="${BK_SWISSPROT_DIR}/bk_flanking_swissprot.tsv"
if [ ! -f "$BK_SWISSPROT_FILE" ]; then
    mkdir -p "$BK_SWISSPROT_DIR"
    cat ${PHASE1_DIR}/*_flanking.faa > /tmp/bk_flanking_${FAMILY}.faa 2>/dev/null || true
    cat ${PHASE1_DIR}/*_flanking_dedup.faa >> /tmp/bk_flanking_${FAMILY}.faa 2>/dev/null || true
    if [ -s /tmp/bk_flanking_${FAMILY}.faa ]; then
        diamond blastp \
            --query /tmp/bk_flanking_${FAMILY}.faa \
            --db ${SWISSPROT_DB%.dmnd} \
            --out "$BK_SWISSPROT_FILE" \
            --outfmt 6 qseqid sseqid pident length evalue bitscore stitle \
            --evalue 1e-5 --max-target-seqs 1 --threads 8 --quiet
        rm -f /tmp/bk_flanking_${FAMILY}.faa
    fi
fi

################################################################################
# PHASE 8a: Locus Matrices
################################################################################
echo ""
echo "=== PHASE 8a: Locus Matrices ==="

BK_SWISSPROT_ARG=""
if [ -f "$BK_SWISSPROT_FILE" ]; then
    BK_SWISSPROT_ARG="--bk-swissprot $BK_SWISSPROT_FILE"
fi

PHASE5B_ARG=""
if [ -d "$PHASE5B_DIR" ] && [ -f "$PHASE5B_DIR/novel_loci_clustered.tsv" ]; then
    PHASE5B_ARG="--phase5b-dir $PHASE5B_DIR"
fi

python helixer_pipeline/08a_generate_locus_matrices_helixer.py \
    --family "$FAMILY" \
    --phase1-dir "$PHASE1_DIR" \
    --phase5-dir "$PHASE5_DIR" \
    $PHASE5B_ARG \
    --phase2b-dir "$PHASE2B_DIR" \
    --phase6-dir "$PHASE6_DIR" \
    --phase7-dir "$PHASE7_DIR" \
    --species-map "$SPECIES_MAP" \
    --output-dir "$PHASE8_DIR" \
    $BK_SWISSPROT_ARG || echo "  WARNING: Phase 8a partial failure"

################################################################################
# PHASE 8b: Summary Matrix
################################################################################
echo ""
echo "=== PHASE 8b: Summary Matrix ==="

python helixer_pipeline/08b_generate_summary_matrices_helixer.py \
    --family "$FAMILY" \
    --phase1-dir "$PHASE1_DIR" \
    --phase2b-dir "$PHASE2B_DIR" \
    --phase5-dir "$PHASE5_DIR" \
    $PHASE5B_ARG \
    --phase6-dir "$PHASE6_DIR" \
    --species-map "$SPECIES_MAP" \
    --output-dir "$PHASE8_DIR"

################################################################################
# PHASE 9: Novel Loci Detection
################################################################################
echo ""
echo "=== PHASE 9: Novel Loci Detection ==="

if [ -d "$PHASE5B_DIR" ]; then
    python helixer_pipeline/09_detect_novel_loci_helixer.py \
        --family "$FAMILY" \
        --phase5b-dir "$PHASE5B_DIR" \
        --phase4b-dir "$PHASE4B_DIR" \
        --helixer-dir "$HELIXER_DIR" \
        --output-dir "$PHASE9_DIR" || echo "  WARNING: Phase 9 failed"
else
    echo "  SKIP: No Phase 5b outputs"
fi

################################################################################
# DONE
################################################################################
echo ""
echo "============================================================"
echo "COMPLETE: Full pipeline for $FAMILY"
echo "  Phase 1:  $PHASE1_DIR"
echo "  Phase 2b: $PHASE2B_DIR"
echo "  Phase 4:  $PHASE4_DIR"
echo "  Phase 5:  $PHASE5_DIR"
echo "  Phase 5b: $PHASE5B_DIR"
echo "  Phase 6:  $PHASE6_DIR"
echo "  Phase 7:  $PHASE7_DIR"
echo "  Phase 8:  $PHASE8_DIR"
echo "  Phase 9:  $PHASE9_DIR"
echo "  Date: $(date)"
echo "============================================================"
