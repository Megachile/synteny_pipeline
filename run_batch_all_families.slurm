#!/bin/bash
#SBATCH --job-name=synteny_batch
#SBATCH --output=logs/batch_array_%A_%a.out
#SBATCH --error=logs/batch_array_%A_%a.err
#SBATCH --time=48:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --partition=general
#SBATCH --array=0-31

################################################################################
# BATCH ARRAY: RUN FULL PIPELINE FOR ALL GENE FAMILIES
#
# Reads from hpc_input_loci.tsv and runs the complete pipeline (Phases 1-8)
# for each gene family
#
# Usage: sbatch run_batch_all_families.slurm
################################################################################

set -e
set -u
set -o pipefail

# Change to workflow directory
cd /carc/scratch/projects/emartins/2016456/adam/synteny_scanning/hpc_deployment_package/hybrid_workflow

# Create logs directory
mkdir -p logs

# Read gene families from input file (skip header and empty lines)
# Use mapfile to preserve spaces in gene family names
mapfile -t GENE_FAMILIES < <(tail -n +2 hpc_input_loci.tsv | grep -v "^$" | cut -f1)

# Get this array task's gene family
GENE_FAMILY="${GENE_FAMILIES[$SLURM_ARRAY_TASK_ID]}"
TARGET_GENES=$(tail -n +2 hpc_input_loci.tsv | awk -F'\t' -v gf="$GENE_FAMILY" '$1 == gf {print $2}')

if [ -z "$GENE_FAMILY" ] || [ -z "$TARGET_GENES" ]; then
    echo "ERROR: Could not find gene family at index $SLURM_ARRAY_TASK_ID"
    exit 1
fi

echo "================================================================================"
echo "ARRAY TASK $SLURM_ARRAY_TASK_ID: Processing $GENE_FAMILY"
echo "================================================================================"
echo "Target genes: $TARGET_GENES"
echo "Job ID: ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "Node: $(hostname)"
echo ""

# Create inputs directory
mkdir -p inputs

# Create input file for this gene family
INPUT_FILE="inputs/${GENE_FAMILY}_test.tsv"
echo -e "gene_family\ttarget_genes" > "$INPUT_FILE"
echo -e "${GENE_FAMILY}\t${TARGET_GENES}" >> "$INPUT_FILE"

echo "Created input file: $INPUT_FILE"

# Set output directory name (clean version of gene family name)
OUTPUT_DIR="${GENE_FAMILY}_BK_LB_batch"

# Load environment
echo "Loading environment..."
set +e  # Disable exit on error for environment setup
set +u  # Disable unbound variable check
source ~/.bashrc 2>&1 | grep -v "cannot access" | grep -v "No such file" || true
eval "$(micromamba shell hook --shell bash)" 2>/dev/null || true
micromamba activate busco_env 2>&1 | grep -v "already activated" || true
module load exonerate/2.4.0-yl7q 2>/dev/null || true
set -u  # Re-enable unbound variable check
set -e  # Re-enable exit on error

echo "Environment ready"
echo ""

################################################################################
# RUN COMPLETE PIPELINE
################################################################################

echo "================================================================================"
echo "STARTING FULL PIPELINE (PHASES 1-8)"
echo "================================================================================"

# Phase 1: Discover paralogs
PHASE12_DIR="outputs/${OUTPUT_DIR}/phase12_landmark"
mkdir -p ${PHASE12_DIR}

echo ""
echo "--- PHASE 1: PARALOG DISCOVERY ---"
if python scripts/01_discover_paralogs_gene_level.py \
    --loc-ids "$TARGET_GENES" \
    --output-dir ${PHASE12_DIR} \
    --gene-family ${GENE_FAMILY} \
    2>&1 | tee logs/phase1_${GENE_FAMILY}_${SLURM_ARRAY_TASK_ID}.log; then
    echo "✓ Phase 1 complete"
else
    echo "❌ Phase 1 failed"
    exit 1
fi

# Phase 2: Synteny comparison and deduplication
echo ""
echo "--- PHASE 2: SYNTENY COMPARISON & DEDUPLICATION ---"
if python scripts/02_compare_synteny_gene_level.py ${PHASE12_DIR} \
    2>&1 | tee logs/phase2_${GENE_FAMILY}_${SLURM_ARRAY_TASK_ID}.log; then
    echo "✓ Phase 2 complete"
else
    echo "❌ Phase 2 failed"
    exit 1
fi

# Phase 3a: Genome-wide synteny detection
PHASE3_DIR="outputs/${OUTPUT_DIR}"
mkdir -p ${PHASE3_DIR}

echo ""
echo "--- PHASE 3a: GENOME-WIDE SYNTENY DETECTION ---"
if python scripts/02_synteny_detection.py \
    --locus-defs ${PHASE12_DIR}/locus_definitions.tsv \
    --genome-db-dir data/ragtag_dbs \
    --output-dir ${PHASE3_DIR}/02_synteny_blocks \
    --evalue 1e-5 \
    --max-targets 50 \
    --max-gap-kb 500 \
    --threads 16 \
    2>&1 | tee logs/phase3a_${GENE_FAMILY}_${SLURM_ARRAY_TASK_ID}.log; then
    echo "✓ Phase 3a complete"
else
    echo "❌ Phase 3a failed"
    exit 1
fi

# Phase 3b: Filter to best block per genome
echo ""
echo "--- PHASE 3b: FILTER TO BEST BLOCKS ---"

# Aggregate blocks if needed
if [ ! -f "${PHASE3_DIR}/02_synteny_blocks/all_synteny_blocks.tsv" ]; then
    echo "Aggregating synteny blocks..."
    echo -e "locus_id\tgenome\tblock_id\tscaffold\tstrand\tstart\tend\tspan_kb\tnum_target_proteins\tnum_query_matches" \
        > "${PHASE3_DIR}/02_synteny_blocks/all_synteny_blocks.tsv"

    for locus_file in ${PHASE3_DIR}/02_synteny_blocks/*_synteny_blocks.tsv; do
        if [ -f "$locus_file" ]; then
            tail -n +2 "$locus_file" >> "${PHASE3_DIR}/02_synteny_blocks/all_synteny_blocks.tsv"
        fi
    done
fi

if python scripts/03_filter_blocks.py \
    --input ${PHASE3_DIR}/02_synteny_blocks/all_synteny_blocks.tsv \
    --output-dir ${PHASE3_DIR}/03_filtered_blocks \
    --verbose \
    2>&1 | tee logs/phase3b_${GENE_FAMILY}_${SLURM_ARRAY_TASK_ID}.log; then
    echo "✓ Phase 3b complete"
else
    echo "❌ Phase 3b failed"
    exit 1
fi

# Phase 3c: Extract target proteins from proteomes
echo ""
echo "--- PHASE 3c: EXTRACT TARGET PROTEINS ---"
if python scripts/03c_extract_target_proteins.py \
    --locus-defs ${PHASE12_DIR}/locus_definitions.tsv \
    --proteome-dir data/proteomes \
    --genome-dir data/genomes \
    --output-dir ${PHASE3_DIR}/01_extracted_proteins \
    2>&1 | tee logs/phase3c_${GENE_FAMILY}_${SLURM_ARRAY_TASK_ID}.log; then
    echo "✓ Phase 3c complete"
else
    echo "❌ Phase 3c failed"
    exit 1
fi

# Phase 4: Target gene detection
echo ""
echo "--- PHASE 4: TARGET GENE DETECTION ---"
if python scripts/04_blast_targets.py \
    --locus-defs ${PHASE12_DIR}/locus_definitions.tsv \
    --query-proteins ${PHASE3_DIR}/01_extracted_proteins \
    --genome-db-dir data/ragtag_dbs \
    --output-dir ${PHASE3_DIR}/04_target_genes \
    --threads 16 \
    2>&1 | tee logs/phase4_${GENE_FAMILY}_${SLURM_ARRAY_TASK_ID}.log; then
    echo "✓ Phase 4 complete"
else
    echo "❌ Phase 4 failed"
    exit 1
fi

# Phase 5: Functional classification
echo ""
echo "--- PHASE 5: FUNCTIONAL CLASSIFICATION ---"
if python scripts/05_classify_targets.py \
    --targets ${PHASE3_DIR}/04_target_genes/all_target_loci.tsv \
    --blocks ${PHASE3_DIR}/03_filtered_blocks/synteny_blocks_filtered.tsv \
    --output-dir ${PHASE3_DIR}/05_classified \
    2>&1 | tee logs/phase5_${GENE_FAMILY}_${SLURM_ARRAY_TASK_ID}.log; then
    echo "✓ Phase 5 complete"
else
    echo "❌ Phase 5 failed"
    exit 1
fi

# Phase 6: Extract sequences with Exonerate
echo ""
echo "--- PHASE 6: EXTRACT GENE STRUCTURES ---"
if python scripts/06_extract_sequences.py \
    --syntenic ${PHASE3_DIR}/05_classified/syntenic_targets.tsv \
    --unplaceable ${PHASE3_DIR}/05_classified/unplaceable_targets.tsv \
    --query-proteins ${PHASE3_DIR}/04_target_genes/combined_targets.faa \
    --genome-fasta-dir data/ragtag_output \
    --output-dir ${PHASE3_DIR}/06_extracted_sequences \
    --unplaceable-evalue 1e-10 \
    --threads 16 \
    --verbose \
    2>&1 | tee logs/phase6_${GENE_FAMILY}_${SLURM_ARRAY_TASK_ID}.log; then
    echo "✓ Phase 6 complete"
else
    echo "❌ Phase 6 failed"
    exit 1
fi

# Phase 7: SwissProt annotation
echo ""
echo "--- PHASE 7: SWISSPROT ANNOTATION ---"
if python scripts/07_swissprot_annotation.py \
    --synteny-dir ${PHASE3_DIR}/02_synteny_blocks \
    --filtered-blocks ${PHASE3_DIR}/03_filtered_blocks/synteny_blocks_filtered.tsv \
    --swissprot-db data/reference/swissprot.dmnd \
    --output-dir ${PHASE3_DIR}/07_swissprot_annotations \
    --evalue 1e-5 \
    --threads 16 \
    2>&1 | tee logs/phase7_${GENE_FAMILY}_${SLURM_ARRAY_TASK_ID}.log; then
    echo "✓ Phase 7 complete"
else
    echo "❌ Phase 7 failed"
    exit 1
fi

# Phase 8a: Generate locus matrices
echo ""
echo "--- PHASE 8a: GENERATE LOCUS MATRICES ---"
cp ${PHASE12_DIR}/locus_definitions.tsv ${PHASE3_DIR}/locus_definitions.tsv

if python scripts/08a_generate_locus_matrices.py \
    --locus-defs ${PHASE3_DIR}/locus_definitions.tsv \
    --synteny-dir ${PHASE3_DIR}/02_synteny_blocks \
    --blocks ${PHASE3_DIR}/03_filtered_blocks/synteny_blocks_filtered.tsv \
    --targets ${PHASE3_DIR}/05_classified/all_targets_classified.tsv \
    --swissprot ${PHASE3_DIR}/07_swissprot_annotations/genome_specific_swissprot_annotations.tsv \
    --reference-proteins data/reference/protein.faa \
    --species-map data/gca_to_species.tsv \
    --output-dir ${PHASE3_DIR}/08_matrices \
    2>&1 | tee logs/phase8a_${GENE_FAMILY}_${SLURM_ARRAY_TASK_ID}.log; then
    echo "✓ Phase 8a complete"
else
    echo "❌ Phase 8a failed"
    exit 1
fi

# Phase 8b: Generate summary matrices
echo ""
echo "--- PHASE 8b: GENERATE SUMMARY MATRICES ---"
if python scripts/08b_generate_summary_matrices.py \
    --locus-matrices-dir ${PHASE3_DIR}/08_matrices \
    --output-dir ${PHASE3_DIR}/08_matrices \
    2>&1 | tee logs/phase8b_${GENE_FAMILY}_${SLURM_ARRAY_TASK_ID}.log; then
    echo "✓ Phase 8b complete"
else
    echo "❌ Phase 8b failed"
    exit 1
fi

################################################################################
# FINAL SUMMARY
################################################################################

echo ""
echo "================================================================================"
echo "PIPELINE COMPLETE - ${GENE_FAMILY}"
echo "================================================================================"
echo ""
echo "Output directory: outputs/${OUTPUT_DIR}/"
echo "All phases (1-8) completed successfully"
echo ""
echo "Runtime: $(($SECONDS / 3600))h $(($SECONDS % 3600 / 60))m $(($SECONDS % 60))s"
echo ""

exit 0
