#!/bin/bash
#SBATCH --job-name=test_helixer
#SBATCH --output=logs/test_helixer_%j.out
#SBATCH --error=logs/test_helixer_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4

################################################################################
# Test script for refactored Helixer Pipeline
# Tests one family end-to-end with fresh outputs
################################################################################

set -e  # Exit on error

cd /carc/scratch/projects/emartins/2016456/adam/synteny_scanning/hpc_deployment_package/hybrid_workflow

FAMILY="sprouty_MC112"

echo "============================================================"
echo "TESTING REFACTORED HELIXER PIPELINE"
echo "  Family: $FAMILY"
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

# Output directories - use _test suffix for this run
PHASE1_DIR="outputs/$FAMILY/phase1_test"
PHASE2B_DIR="outputs/$FAMILY/phase2b_test"
PHASE4_DIR="outputs/$FAMILY/phase4_test"
PHASE4B_DIR="outputs/$FAMILY/phase4b_test"
PHASE5_DIR="outputs/$FAMILY/phase5_test"
PHASE5B_DIR="outputs/$FAMILY/phase5b_test"
PHASE6_DIR="outputs/$FAMILY/phase6_test"
PHASE7_DIR="outputs/$FAMILY/phase7_test"
PHASE8_DIR="outputs/$FAMILY/phase8_test"

mkdir -p logs

# Clean any previous test outputs
rm -rf outputs/$FAMILY/*_test

################################################################################
# PHASE 1: Locus Discovery
################################################################################
echo ""
echo "=== PHASE 1: Locus Discovery ==="

INPUT_FILE="inputs/${FAMILY}_test.tsv"
LOC_IDS=$(tail -n +2 "$INPUT_FILE" | cut -f2 | tr '\n' ',' | sed 's/,$//')
echo "  Input: $INPUT_FILE"
echo "  LOC IDs: $LOC_IDS"

python helixer_pipeline/01_phase1.py \
    --loc-ids "$LOC_IDS" \
    --gene-family "$FAMILY" \
    --output-dir "$PHASE1_DIR" \
    --genome-gff-dir "$GENOME_GFF_DIR"

echo "  Phase 1 complete: $(ls $PHASE1_DIR | wc -l) files"

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
    --threads 4 \
    --min-query-coverage 0.3

echo "  Phase 4 complete"

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

echo "  Phase 4b complete"

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
    --threads 4

echo "  Phase 5 complete"

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

python helixer_pipeline/05b_validate_novel_loci.py \
    --family "$FAMILY" \
    --phase5-dir "$PHASE5_DIR" \
    --helixer-dir "$HELIXER_DIR" \
    --phase1-dir "$PHASE1_DIR" \
    --output-dir "$PHASE5B_DIR" \
    --threads 4

echo "  Phase 5b complete"

################################################################################
# PHASE 6: Extract Sequences (REFACTORED)
################################################################################
echo ""
echo "=== PHASE 6: Extract Sequences (REFACTORED) ==="

python helixer_pipeline/06_extract_sequences_helixer.py \
    --family "$FAMILY" \
    --phase1-dir "$PHASE1_DIR" \
    --phase4-dir "$PHASE4_DIR" \
    --phase5-dir "$PHASE5_DIR" \
    --phase5b-dir "$PHASE5B_DIR" \
    --helixer-dir "$HELIXER_DIR" \
    --output-dir "$PHASE6_DIR"

echo "  Phase 6 complete"

################################################################################
# PHASE 7: SwissProt Annotation (REFACTORED)
################################################################################
echo ""
echo "=== PHASE 7: SwissProt Annotation (REFACTORED) ==="

python helixer_pipeline/07_annotate_flanking_swissprot.py \
    --family "$FAMILY" \
    --phase4b-dir "$PHASE4B_DIR" \
    --phase5-dir "$PHASE5_DIR" \
    --phase5b-dir "$PHASE5B_DIR" \
    --swissprot-db "$SWISSPROT_DB" \
    --output-dir "$PHASE7_DIR" \
    --threads 4

echo "  Phase 7 complete"

################################################################################
# PHASE 8a: Locus Matrices (REFACTORED)
################################################################################
echo ""
echo "=== PHASE 8a: Locus Matrices (REFACTORED) ==="

python helixer_pipeline/08a_generate_locus_matrices_helixer.py \
    --family "$FAMILY" \
    --phase1-dir "$PHASE1_DIR" \
    --phase5-dir "$PHASE5_DIR" \
    --phase5b-dir "$PHASE5B_DIR" \
    --phase6-dir "$PHASE6_DIR" \
    --phase7-dir "$PHASE7_DIR" \
    --species-map "$SPECIES_MAP" \
    --output-dir "$PHASE8_DIR"

echo "  Phase 8a complete"

################################################################################
# PHASE 8b: Summary Matrix (REFACTORED)
################################################################################
echo ""
echo "=== PHASE 8b: Summary Matrix (REFACTORED) ==="

python helixer_pipeline/08b_generate_summary_matrices_helixer.py \
    --family "$FAMILY" \
    --phase1-dir "$PHASE1_DIR" \
    --phase5-dir "$PHASE5_DIR" \
    --phase5b-dir "$PHASE5B_DIR" \
    --phase6-dir "$PHASE6_DIR" \
    --species-map "$SPECIES_MAP" \
    --output-dir "$PHASE8_DIR"

echo "  Phase 8b complete"

################################################################################
# VERIFY OUTPUTS
################################################################################
echo ""
echo "============================================================"
echo "OUTPUT VERIFICATION"
echo "============================================================"

echo ""
echo "Phase 1:"
ls -la $PHASE1_DIR/*.tsv 2>/dev/null || echo "  No TSV files"

echo ""
echo "Phase 4:"
ls -la $PHASE4_DIR/*.tsv 2>/dev/null || echo "  No TSV files"

echo ""
echo "Phase 5:"
ls -la $PHASE5_DIR/*.tsv 2>/dev/null || echo "  No TSV files"

echo ""
echo "Phase 6:"
ls -la $PHASE6_DIR/*.tsv 2>/dev/null || echo "  No TSV files"
ls -la $PHASE6_DIR/*.faa 2>/dev/null || echo "  No FAA files"

echo ""
echo "Phase 7:"
ls -la $PHASE7_DIR/*.tsv 2>/dev/null || echo "  No TSV files"

echo ""
echo "Phase 8:"
ls -la $PHASE8_DIR/*.tsv 2>/dev/null || echo "  No TSV files"

echo ""
echo "============================================================"
echo "TEST COMPLETE: $FAMILY"
echo "  Date: $(date)"
echo "============================================================"
