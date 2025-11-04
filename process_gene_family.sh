#!/bin/bash
# Wrapper script to process one gene family through Phases 1-5
# Usage: process_gene_family.sh <family_name> <loc_list>

set -e  # Exit on error

FAMILY_NAME=$1
LOC_LIST=$2

if [ -z "$FAMILY_NAME" ] || [ -z "$LOC_LIST" ]; then
    echo "ERROR: Missing arguments"
    echo "Usage: $0 <family_name> <loc_list>"
    exit 1
fi

echo "================================================================================"
echo "PROCESSING GENE FAMILY: $FAMILY_NAME"
echo "================================================================================"
echo "LOCs: $LOC_LIST"
echo "Started: $(date)"
echo ""

# Set up directories
BASE_DIR="/carc/scratch/projects/emartins/2016456/adam/synteny_scanning/hpc_deployment_package/hybrid_workflow"
OUTPUT_BASE="$BASE_DIR/outputs/batch_run"
FAMILY_OUTPUT="$OUTPUT_BASE/$FAMILY_NAME"

mkdir -p "$FAMILY_OUTPUT"

cd "$BASE_DIR"

# Load environment
source ~/.bashrc
micromamba activate busco_env

# Set gene family for downstream scripts
export GENE_FAMILY="$FAMILY_NAME"

echo "=== PHASE 1: Paralog Discovery ==="
python scripts/01_discover_paralogs_gene_level.py "$LOC_LIST" "$FAMILY_OUTPUT/01_paralogs" "$FAMILY_NAME"
if [ $? -ne 0 ]; then
    echo "ERROR: Phase 1 failed for $FAMILY_NAME"
    exit 1
fi
echo ""

# Set environment variables for subsequent phases
export SYNTENY_INPUT_FILE="$FAMILY_OUTPUT/01_paralogs/unique_loci.tsv"
export SYNTENY_OUTPUTS_DIR="$FAMILY_OUTPUT"

echo "=== PHASE 3: Synteny Detection ==="
python scripts/04_blast_targets.py
if [ $? -ne 0 ]; then
    echo "ERROR: Phase 3 failed for $FAMILY_NAME"
    exit 1
fi
echo ""

echo "=== PHASE 4: Target Gene Detection ==="
# Phase 4 runs as part of the same script (04_blast_targets.py handles both)
echo "Completed as part of Phase 3"
echo ""

echo "=== PHASE 5: Target Classification ==="
python scripts/05_classify_targets.py
if [ $? -ne 0 ]; then
    echo "ERROR: Phase 5 failed for $FAMILY_NAME"
    exit 1
fi
echo ""

echo "================================================================================"
echo "COMPLETED: $FAMILY_NAME"
echo "Finished: $(date)"
echo "Output directory: $FAMILY_OUTPUT"
echo "================================================================================"
