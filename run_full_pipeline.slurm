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
#   Step 01: Locus discovery (BK GFF parsing or coordinates file)
#   Step 02: BK-LB unification (DIAMOND merge)
#   Step 03: Target detection (DIAMOND vs all proteomes)
#   Step 04: Extract flanking genes
#   Step 05: Synteny block detection (flanking-first)
#   Step 06: Classify targets (syntenic vs novel)
#   Step 07: Validate novel loci (synteny-aware clustering)
#   Step 08: Extract sequences (proteins from targets + NR annotations)
#   Step 09: Annotate flanking (parse NR annotations from headers)
#   Step 10: Generate matrices (locus + summary)
################################################################################

set -e  # Exit on error

cd /carc/scratch/projects/emartins/2016456/adam/analysis/helixer_pipeline

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

python 01_discover_loci.py \
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
    python 02_merge_reference_loci.py \
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

python 03_detect_targets.py \
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

python 04_extract_flanking.py \
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

python 05_detect_synteny_blocks.py \
    --family "$FAMILY" \
    --phase1-dir "$PHASE1_DIR" \
    --phase4-targets "${PHASE4_DIR}/all_target_loci.tsv" \
    --helixer-dir "$HELIXER_DIR" \
    --output-dir "$PHASE2B_DIR" \
    --threads 8 || echo "  WARNING: Phase 2b failed, continuing"

################################################################################
# PHASE 5: Classify Targets
################################################################################
echo ""
echo "=== PHASE 5: Classify Targets ==="

python 06_classify_targets.py \
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
    python 07_validate_novel_loci.py \
        --family "$FAMILY" \
        --phase5-dir "$PHASE5_DIR" \
        --phase1-dir "$PHASE1_DIR" \
        --helixer-dir "$HELIXER_DIR" \
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

python 08_extract_sequences.py \
    --family "$FAMILY" \
    --phase1-dir "$PHASE1_DIR" \
    --phase4-dir "$PHASE4_DIR" \
    --phase5-dir "$PHASE5_DIR" \
    --helixer-dir "$HELIXER_DIR" \
    --output-dir "$PHASE6_DIR"

################################################################################
# PHASE 7: Annotate Flanking (NR annotations from Helixer FASTA headers)
################################################################################
echo ""
echo "=== PHASE 7: Annotate Flanking (NR) ==="

PHASE2B_ARG=""
if [ -d "$PHASE2B_DIR" ]; then
    PHASE2B_ARG="--phase2b-dir $PHASE2B_DIR"
fi

python 09_annotate_flanking.py \
    --family "$FAMILY" \
    --phase4b-dir "$PHASE4B_DIR" \
    --phase5-dir "$PHASE5_DIR" \
    --helixer-dir "$HELIXER_DIR" \
    $PHASE2B_ARG \
    --output-dir "$PHASE7_DIR"

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

python 10a_generate_locus_matrices.py \
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

python 10b_generate_summary_matrices.py \
    --family "$FAMILY" \
    --phase1-dir "$PHASE1_DIR" \
    --phase2b-dir "$PHASE2B_DIR" \
    --phase5-dir "$PHASE5_DIR" \
    $PHASE5B_ARG \
    --phase6-dir "$PHASE6_DIR" \
    --species-map "$SPECIES_MAP" \
    --output-dir "$PHASE8_DIR"

################################################################################
# PHASE 9: Novel Loci Detection (DEPRECATED 2026-01-04)
# - Script moved to deprecated/ folder
# - Input/output mismatch with Phase 5b prevents execution
# - Phase 5b's novel_loci_clustered.tsv provides sufficient novel locus info
################################################################################
# Skipping deprecated Phase 9
echo ""
echo "=== PHASE 9: Skipped (deprecated) ==="

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
echo "  Date: $(date)"
echo "============================================================"
