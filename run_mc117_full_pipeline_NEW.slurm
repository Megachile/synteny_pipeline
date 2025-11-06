#!/bin/bash
#SBATCH --job-name=mc117_new_pipeline
#SBATCH --output=logs/mc117_pipeline_NEW_%j.out
#SBATCH --error=logs/mc117_pipeline_NEW_%j.err
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --partition=general

# CRITICAL: Change to working directory
cd /carc/scratch/projects/emartins/2016456/adam/synteny_scanning/hpc_deployment_package/hybrid_workflow

# Activate environment
source ~/.bashrc
eval "$(micromamba shell hook --shell bash)"
micromamba activate busco_env

# Set variables
FAMILY="outputs/MC117_BK_LB_batch"
GENOME_DB="../data/ragtag_dbs"

echo "======================================================================"
echo "MC117 FULL PIPELINE TEST WITH NEW PHASE 2 CLUSTERING"
echo "======================================================================"
echo ""
echo "Settings:"
echo "  - Phase 2: --max-gap-kb 100 (no --max-span-kb, using MAX_BLOCK_SPAN_KB=300 constant)"
echo "  - Output directories: _NEW suffix to avoid overwriting"
echo "  - All 27 BK chr5 loci"
echo ""

# Step 02: Synteny detection with new clustering
echo "[1/4] Running Phase 2: Synteny detection..."
python scripts/02_synteny_detection.py \
    --locus-defs ${FAMILY}/locus_definitions.tsv \
    --genome-db-dir ${GENOME_DB} \
    --output-dir ${FAMILY}/02_synteny_blocks_NEW \
    --max-gap-kb 100 \
    --threads 16

if [ $? -ne 0 ]; then
    echo "ERROR: Phase 2 failed"
    exit 1
fi

# Step 02b: Aggregate synteny blocks
echo ""
echo "[2/4] Running Phase 2b: Aggregate synteny blocks..."
python scripts/02b_aggregate_synteny_blocks.py \
    --input-dir ${FAMILY}/02_synteny_blocks_NEW \
    --output ${FAMILY}/02_synteny_blocks_NEW/all_synteny_blocks.tsv

if [ $? -ne 0 ]; then
    echo "ERROR: Phase 2b failed"
    exit 1
fi

# Step 03: Filter best blocks
echo ""
echo "[3/4] Running Phase 3: Filter best blocks..."
python scripts/03_filter_blocks.py \
    --input ${FAMILY}/02_synteny_blocks_NEW/all_synteny_blocks.tsv \
    --output-dir ${FAMILY}/03_filtered_blocks_NEW \
    --verbose

if [ $? -ne 0 ]; then
    echo "ERROR: Phase 3 failed"
    exit 1
fi

# Step 05: Classify targets (using existing Phase 4 targets)
echo ""
echo "[4/4] Running Phase 5: Classify targets..."
python scripts/05_classify_targets.py \
    --targets ${FAMILY}/04_target_genes/all_target_loci.tsv \
    --blocks ${FAMILY}/03_filtered_blocks_NEW/synteny_blocks_filtered.tsv \
    --output-dir ${FAMILY}/05_classified_NEW

if [ $? -ne 0 ]; then
    echo "ERROR: Phase 5 failed"
    exit 1
fi

echo ""
echo "======================================================================"
echo "PIPELINE COMPLETE"
echo "======================================================================"
echo ""
echo "Results:"
echo "  Phase 2: ${FAMILY}/02_synteny_blocks_NEW/"
echo "  Phase 3: ${FAMILY}/03_filtered_blocks_NEW/"
echo "  Phase 5: ${FAMILY}/05_classified_NEW/"
echo ""

# Compare BK target recovery
echo "Comparing BK target recovery (GCA_010883055.1):"
echo ""
echo "OLD Phase 5:"
if [ -f "${FAMILY}/05_classified/syntenic_targets.tsv" ]; then
    OLD_COUNT=$(grep "GCA_010883055" ${FAMILY}/05_classified/syntenic_targets.tsv 2>/dev/null | wc -l)
    echo "  BK syntenic targets: ${OLD_COUNT}"
fi
echo ""
echo "NEW Phase 5:"
if [ -f "${FAMILY}/05_classified_NEW/syntenic_targets.tsv" ]; then
    NEW_COUNT=$(grep "GCA_010883055" ${FAMILY}/05_classified_NEW/syntenic_targets.tsv 2>/dev/null | wc -l)
    echo "  BK syntenic targets: ${NEW_COUNT}"
fi
echo ""

# Count blocks per locus for BK
echo "BK blocks per locus comparison:"
echo "OLD (should have many mega-blocks):"
grep "GCA_010883055" ${FAMILY}/02_synteny_blocks/all_synteny_blocks.tsv 2>/dev/null | \
    awk '{print $1}' | sort | uniq -c | head -10
echo ""
echo "NEW (should have more, smaller blocks):"
grep "GCA_010883055" ${FAMILY}/02_synteny_blocks_NEW/all_synteny_blocks.tsv 2>/dev/null | \
    awk '{print $1}' | sort | uniq -c | head -10
