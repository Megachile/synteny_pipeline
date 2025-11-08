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

# CRITICAL: Change to workflow directory using absolute path
# SLURM jobs start in the directory where sbatch is called, NOT where the script lives
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

# Sanitize gene family name for filesystem use (remove spaces, pipes, commas, special chars)
CLEAN_GENE_FAMILY=$(echo "$GENE_FAMILY" | sed 's/ /_/g; s/|/_/g; s/,//g; s/[^a-zA-Z0-9_-]/_/g')

echo "================================================================================"
echo "ARRAY TASK $SLURM_ARRAY_TASK_ID: Processing $GENE_FAMILY"
echo "================================================================================"
echo "Sanitized name: $CLEAN_GENE_FAMILY"
echo "Target genes: $TARGET_GENES"
echo "Job ID: ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "Node: $(hostname)"
echo ""

# Create inputs directory
mkdir -p inputs

# Create input file for this gene family (use sanitized name for filename)
INPUT_FILE="inputs/${CLEAN_GENE_FAMILY}_test.tsv"
echo -e "gene_family\ttarget_genes" > "$INPUT_FILE"
echo -e "${GENE_FAMILY}\t${TARGET_GENES}" >> "$INPUT_FILE"

echo "Created input file: $INPUT_FILE"

# Set output directory name (sanitized version of gene family name)
OUTPUT_DIR="${CLEAN_GENE_FAMILY}_BK_LB_batch"

# Load environment
echo "Loading environment..."

# Initialize micromamba
eval "$(micromamba shell hook --shell bash)"
micromamba activate busco_env

# Load Exonerate module
module load exonerate/2.4.0-yl7q 2>/dev/null || echo "  Exonerate module not loaded (may already be available)"

echo "Environment ready"
echo "  Python: $(which python)"
echo "  Conda env: $CONDA_DEFAULT_ENV"
echo ""

################################################################################
# RUN COMPLETE PIPELINE
################################################################################

echo "================================================================================"
echo "STARTING FULL PIPELINE (PHASES 1-8)"
echo "================================================================================"

# Phase 1: Discover paralogs
PHASE12_DIR="outputs/${OUTPUT_DIR}/phase12_landmark"
mkdir -p "${PHASE12_DIR}"

echo ""
echo "--- PHASE 1: PARALOG DISCOVERY ---"
if [ -f "${PHASE12_DIR}/locus_definitions.tsv" ]; then
    echo "✓ Phase 1 already complete (skipping)"
else
    if python scripts/01_discover_paralogs_gene_level.py \
        --loc-ids "$TARGET_GENES" \
        --output-dir "${PHASE12_DIR}" \
        --gene-family "${GENE_FAMILY}" \
        2>&1 | tee "logs/phase1_${CLEAN_GENE_FAMILY}_${SLURM_ARRAY_TASK_ID}.log"; then
        echo "✓ Phase 1 complete"
    else
        echo "❌ Phase 1 failed"
        exit 1
    fi
fi

# Phase 2: Synteny comparison and deduplication
echo ""
echo "--- PHASE 2: SYNTENY COMPARISON & DEDUPLICATION ---"
if [ -f "${PHASE12_DIR}/synteny_groups.json" ]; then
    echo "✓ Phase 2 already complete (skipping)"
else
    if python scripts/02_compare_synteny_gene_level.py "${PHASE12_DIR}" \
        2>&1 | tee "logs/phase2_${CLEAN_GENE_FAMILY}_${SLURM_ARRAY_TASK_ID}.log"; then
        echo "✓ Phase 2 complete"
    else
        echo "❌ Phase 2 failed"
        exit 1
    fi
fi

# Phase 3a: Genome-wide synteny detection (parallelized per locus)
PHASE3_DIR="outputs/${OUTPUT_DIR}"
mkdir -p "${PHASE3_DIR}"
mkdir -p "${PHASE3_DIR}/02_synteny_blocks"

echo ""
echo "--- PHASE 3a: GENOME-WIDE SYNTENY DETECTION ---"

SYNTENY_BLOCKS_DIR="${PHASE3_DIR}/02_synteny_blocks"
ALL_BLOCKS="${SYNTENY_BLOCKS_DIR}/all_synteny_blocks.tsv"

# Get locus definitions
LOCUS_DEFS="${PHASE12_DIR}/locus_definitions.tsv"
if [ ! -f "$LOCUS_DEFS" ]; then
    echo "❌ Phase 3a failed - locus definitions not found"
    exit 1
fi

NLOCI=$(tail -n +2 "$LOCUS_DEFS" | wc -l)

if [ $NLOCI -eq 0 ]; then
    echo "❌ Phase 3a failed - no loci found in locus definitions"
    exit 1
fi

# Check which loci already have synteny files
shopt -s nullglob
existing_locus_files=("${SYNTENY_BLOCKS_DIR}"/*_synteny_blocks.tsv)
existing_count=0
for f in "${existing_locus_files[@]}"; do
    if [[ ! "$f" =~ all_synteny_blocks.tsv$ ]]; then
        existing_count=$((existing_count + 1))
    fi
done

if [ $existing_count -eq $NLOCI ]; then
    echo "✓ Phase 3a already complete for all ${NLOCI} loci (skipping array job)"
    # Just aggregate if needed
    if [ ! -f "$ALL_BLOCKS" ] || [ $(wc -l < "$ALL_BLOCKS") -le 1 ]; then
        echo "  Aggregating existing results..."
        echo -e "locus_id\tgenome\tblock_id\tscaffold\tstrand\tstart\tend\tspan_kb\tnum_target_proteins\tnum_query_matches" > "$ALL_BLOCKS"
        for f in "${existing_locus_files[@]}"; do
            if [[ ! "$f" =~ all_synteny_blocks.tsv$ ]]; then
                tail -n +2 "$f" >> "$ALL_BLOCKS"
            fi
        done
        TOTAL_BLOCKS=$(tail -n +2 "$ALL_BLOCKS" | wc -l)
        echo "  Aggregated ${TOTAL_BLOCKS} total blocks"
    fi
else
    echo "  Detected ${NLOCI} total loci, ${existing_count} already complete"
    echo "  Submitting Phase 3a array job for all loci (0-$((NLOCI-1)))..."

    PHASE3A_JOB=$(sbatch --parsable \
        --array=0-$((NLOCI-1)) \
        --job-name="phase3a_${GENE_FAMILY}" \
        --output="logs/phase3a_${CLEAN_GENE_FAMILY}_%A_%a.out" \
        --error="logs/phase3a_${CLEAN_GENE_FAMILY}_%A_%a.err" \
        --cpus-per-task=16 \
        --mem=32G \
        --time=4:00:00 \
        run_phase3a_array.slurm "${LOCUS_DEFS}" "${PHASE3_DIR}")

    echo "  Phase 3a array job submitted: ${PHASE3A_JOB}"
    echo "  Waiting for completion..."

    # Wait for array job to complete
    while true; do
        if squeue -j ${PHASE3A_JOB} 2>/dev/null | grep -q "${PHASE3A_JOB}"; then
            sleep 30
        else
            break
        fi
    done

    # Check for failures
    FAILED_TASKS=0
    if compgen -G "logs/phase3a_${CLEAN_GENE_FAMILY}_${PHASE3A_JOB}_*.err" > /dev/null; then
        FAILED_TASKS=$(grep -l "Phase 3a failed\|ERROR" logs/phase3a_${CLEAN_GENE_FAMILY}_${PHASE3A_JOB}_*.err 2>/dev/null | wc -l || true)
    fi

    if [ $FAILED_TASKS -gt 0 ]; then
        echo "❌ Phase 3a failed - ${FAILED_TASKS} loci failed"
        exit 1
    fi

    # Aggregate results from all loci
    echo "  Aggregating results from ${NLOCI} loci..."

    # Find per-locus output files
    shopt -s nullglob
    locus_files=("${SYNTENY_BLOCKS_DIR}"/*_synteny_blocks.tsv)
    filtered_files=()
    for f in "${locus_files[@]}"; do
        if [[ ! "$f" =~ all_synteny_blocks.tsv$ ]]; then
            filtered_files+=("$f")
        fi
    done

    if [ ${#filtered_files[@]} -eq 0 ]; then
        echo "❌ Phase 3a failed - no locus output files found"
        exit 1
    fi

    echo "  Found ${#filtered_files[@]} locus files"

    # Create header
    echo -e "locus_id\tgenome\tblock_id\tscaffold\tstrand\tstart\tend\tspan_kb\tnum_target_proteins\tnum_query_matches" \
        > "$ALL_BLOCKS"

    # Append all per-locus files
    for locus_file in "${filtered_files[@]}"; do
        tail -n +2 "$locus_file" >> "$ALL_BLOCKS"
    done

    TOTAL_BLOCKS=$(tail -n +2 "$ALL_BLOCKS" | wc -l)
    echo "✓ Phase 3a complete"
    echo "  Aggregated output: ${ALL_BLOCKS} (${TOTAL_BLOCKS} total blocks)"
fi

# Phase 3b: Filter to best block per genome
echo ""
echo "--- PHASE 3b: FILTER TO BEST BLOCKS ---"

if [ -f "${PHASE3_DIR}/03_filtered_blocks/synteny_blocks_filtered.tsv" ]; then
    echo "✓ Phase 3b already complete (skipping)"
else
    # Aggregate blocks if needed
    if [ ! -f "${PHASE3_DIR}/02_synteny_blocks/all_synteny_blocks.tsv" ]; then
        echo "Aggregating synteny blocks..."
        echo -e "locus_id\tgenome\tblock_id\tscaffold\tstrand\tstart\tend\tspan_kb\tnum_target_proteins\tnum_query_matches" \
            > "${PHASE3_DIR}/02_synteny_blocks/all_synteny_blocks.tsv"

        for locus_file in "${PHASE3_DIR}"/02_synteny_blocks/*_synteny_blocks.tsv; do
            if [ -f "$locus_file" ]; then
                tail -n +2 "$locus_file" >> "${PHASE3_DIR}/02_synteny_blocks/all_synteny_blocks.tsv"
            fi
        done
    fi

    if python scripts/03_filter_blocks.py \
        --input "${PHASE3_DIR}/02_synteny_blocks/all_synteny_blocks.tsv" \
        --output-dir "${PHASE3_DIR}/03_filtered_blocks" \
        --verbose \
        2>&1 | tee "logs/phase3b_${CLEAN_GENE_FAMILY}_${SLURM_ARRAY_TASK_ID}.log"; then
        echo "✓ Phase 3b complete"
    else
        echo "❌ Phase 3b failed"
        exit 1
    fi
fi

# Phase 3c: Extract target proteins from proteomes
echo ""
echo "--- PHASE 3c: EXTRACT TARGET PROTEINS ---"
# Check if extracted proteins exist (they're in subdirectories per locus)
if [ -d "${PHASE3_DIR}/01_extracted_proteins" ] && [ -n "$(find ${PHASE3_DIR}/01_extracted_proteins -name '*_targets.faa' -type f 2>/dev/null)" ]; then
    echo "✓ Phase 3c already complete (skipping)"
else
    if python scripts/03c_extract_target_proteins.py \
        --locus-defs "${PHASE12_DIR}/locus_definitions.tsv" \
        --proteome-dir data/proteomes \
        --genome-dir data/genomes \
        --output-dir "${PHASE3_DIR}/01_extracted_proteins" \
        2>&1 | tee "logs/phase3c_${CLEAN_GENE_FAMILY}_${SLURM_ARRAY_TASK_ID}.log"; then
        echo "✓ Phase 3c complete"
    else
        echo "❌ Phase 3c failed"
        exit 1
    fi
fi

# Phase 4: Target gene detection
echo ""
echo "--- PHASE 4: TARGET GENE DETECTION ---"
if [ -f "${PHASE3_DIR}/04_target_genes/all_target_loci.tsv" ]; then
    echo "✓ Phase 4 already complete (skipping)"
else
    if python scripts/04_blast_targets.py \
        --locus-defs "${PHASE12_DIR}/locus_definitions.tsv" \
        --query-proteins "${PHASE3_DIR}/01_extracted_proteins" \
        --genome-db-dir data/ragtag_dbs \
        --output-dir "${PHASE3_DIR}/04_target_genes" \
        --gene-family "${GENE_FAMILY}" \
        --threads 16 \
        2>&1 | tee "logs/phase4_${CLEAN_GENE_FAMILY}_${SLURM_ARRAY_TASK_ID}.log"; then
        echo "✓ Phase 4 complete"
    else
        echo "❌ Phase 4 failed"
        exit 1
    fi
fi

# Phase 5: Functional classification
echo ""
echo "--- PHASE 5: FUNCTIONAL CLASSIFICATION ---"
if [ -f "${PHASE3_DIR}/05_classified/all_targets_classified.tsv" ]; then
    echo "✓ Phase 5 already complete (skipping)"
else
    if python scripts/05_classify_targets.py \
        --targets "${PHASE3_DIR}/04_target_genes/all_target_loci.tsv" \
        --blocks "${PHASE3_DIR}/03_filtered_blocks/synteny_blocks_filtered.tsv" \
        --output-dir "${PHASE3_DIR}/05_classified" \
        2>&1 | tee "logs/phase5_${CLEAN_GENE_FAMILY}_${SLURM_ARRAY_TASK_ID}.log"; then
        # Verify outputs exist and are non-empty
        if [ ! -s "${PHASE3_DIR}/05_classified/all_targets_classified.tsv" ]; then
            echo "❌ Phase 5 produced no classifications (all_targets_classified.tsv missing or empty)"
            exit 1
        fi
        echo "✓ Phase 5 complete"
    else
        echo "❌ Phase 5 failed"
        exit 1
    fi
fi

# Phase 6: Extract sequences with Exonerate
echo ""
echo "--- PHASE 6: EXTRACT GENE STRUCTURES ---"
if [ -f "${PHASE3_DIR}/06_extracted_sequences/all_extracted_genes.faa" ]; then
    echo "✓ Phase 6 already complete (skipping)"
else
    if python scripts/06_extract_sequences.py \
        --syntenic "${PHASE3_DIR}/05_classified/syntenic_targets.tsv" \
        --unplaceable "${PHASE3_DIR}/05_classified/unplaceable_targets.tsv" \
        --query-proteins "${PHASE3_DIR}/04_target_genes/combined_targets.faa" \
        --genome-fasta-dir data/ragtag_output \
        --output-dir "${PHASE3_DIR}/06_extracted_sequences" \
        --unplaceable-evalue 1e-10 \
        --verbose \
        2>&1 | tee "logs/phase6_${CLEAN_GENE_FAMILY}_${SLURM_ARRAY_TASK_ID}.log"; then
        # Verify aggregated file was created
        if [ ! -f "${PHASE3_DIR}/06_extracted_sequences/all_extracted_genes.faa" ]; then
            echo "❌ Phase 6 failed - aggregated file not created"
            exit 1
        fi
        echo "✓ Phase 6 complete"
    else
        echo "❌ Phase 6 failed"
        exit 1
    fi
fi

# Phase 7: SwissProt annotation (parallelized per locus)
echo ""
echo "--- PHASE 7: SWISSPROT ANNOTATION ---"

ANNOT_DIR="${PHASE3_DIR}/07_swissprot_annotations"
FINAL_OUTPUT="${ANNOT_DIR}/genome_specific_swissprot_annotations.tsv"

if [ -f "$FINAL_OUTPUT" ]; then
    echo "✓ Phase 7 already complete (skipping)"
else
    # Count loci from filtered blocks
    FILTERED_BLOCKS="${PHASE3_DIR}/03_filtered_blocks/synteny_blocks_filtered.tsv"
    if [ ! -f "$FILTERED_BLOCKS" ]; then
        echo "❌ Phase 7 failed - filtered blocks not found"
        exit 1
    fi

    NLOCI=$(tail -n +2 "$FILTERED_BLOCKS" | cut -f1 | sort -u | wc -l)
    echo "  Detected ${NLOCI} loci to annotate"

    if [ $NLOCI -eq 0 ]; then
        echo "❌ Phase 7 failed - no loci found in filtered blocks"
        exit 1
    fi

    # Submit Phase 7 array job (one task per locus)
    echo "  Submitting Phase 7 array job (0-$((NLOCI-1)))..."
    PHASE7_JOB=$(sbatch --parsable \
        --array=0-$((NLOCI-1)) \
        --job-name="phase7_${GENE_FAMILY}" \
        --output="logs/phase7_${CLEAN_GENE_FAMILY}_%A_%a.out" \
        --error="logs/phase7_${CLEAN_GENE_FAMILY}_%A_%a.err" \
        --cpus-per-task=8 \
        --mem=16G \
        --time=2:00:00 \
        run_phase7_array.slurm "${PHASE3_DIR}")

    echo "  Phase 7 array job submitted: ${PHASE7_JOB}"
    echo "  Waiting for completion..."

    # Wait for array job to complete
    while true; do
        # Check if job is still in queue (running or pending)
        if squeue -j ${PHASE7_JOB} 2>/dev/null | grep -q "${PHASE7_JOB}"; then
            sleep 30
        else
            # Job not in queue - either completed or failed
            break
        fi
    done

    # Check if any tasks failed
    FAILED_TASKS=0
    if compgen -G "logs/phase7_${CLEAN_GENE_FAMILY}_${PHASE7_JOB}_*.err" > /dev/null; then
        # grep returns 1 when no matches found, so use || true to prevent pipefail
        FAILED_TASKS=$(grep -l "Phase 7 failed" logs/phase7_${CLEAN_GENE_FAMILY}_${PHASE7_JOB}_*.err 2>/dev/null | wc -l || true)
    fi
    if [ $FAILED_TASKS -gt 0 ]; then
        echo "❌ Phase 7 failed - ${FAILED_TASKS} array tasks failed"
        exit 1
    fi

    # Aggregate per-locus results into single file
    echo "  Aggregating results from ${NLOCI} loci..."
    echo "  Starting aggregation" >&2

    # Create header from first locus file
    shopt -s nullglob
    locus_files=(${ANNOT_DIR}/*_swissprot_annotations.tsv)
    if [ ${#locus_files[@]} -eq 0 ]; then
        echo "❌ Phase 7 failed - no locus annotation files found"
        exit 1
    fi

    # Get first non-genome_specific file for header
    for first_file in "${locus_files[@]}"; do
        if [[ ! "$first_file" =~ "genome_specific" ]]; then
            head -1 "$first_file" > "$FINAL_OUTPUT"
            break
        fi
    done

    # Concatenate all per-locus results (skip headers)
    for locus_file in "${locus_files[@]}"; do
        if [[ ! "$locus_file" =~ "genome_specific" ]]; then
            tail -n +2 "$locus_file" >> "$FINAL_OUTPUT"
        fi
    done

    echo "✓ Phase 7 complete (${NLOCI} loci processed in parallel)"
    echo "  Aggregated output: ${FINAL_OUTPUT}"
fi

# Phase 8a: Generate locus matrices
echo ""
echo "--- PHASE 8a: GENERATE LOCUS MATRICES ---"

if [ -d "${PHASE3_DIR}/08_matrices" ] && [ -n "$(ls -A ${PHASE3_DIR}/08_matrices/*_locus_matrix.tsv 2>/dev/null)" ]; then
    echo "✓ Phase 8a already complete (skipping)"
else
    cp "${PHASE12_DIR}/locus_definitions.tsv" "${PHASE3_DIR}/locus_definitions.tsv"

    if python scripts/08a_generate_locus_matrices.py \
        --locus-defs "${PHASE3_DIR}/locus_definitions.tsv" \
        --synteny-dir "${PHASE3_DIR}/02_synteny_blocks" \
        --blocks "${PHASE3_DIR}/03_filtered_blocks/synteny_blocks_filtered.tsv" \
        --targets "${PHASE3_DIR}/05_classified/all_targets_classified.tsv" \
        --swissprot "${PHASE3_DIR}/07_swissprot_annotations/genome_specific_swissprot_annotations.tsv" \
        --reference-proteins data/reference/protein.faa \
        --extracted-seqs "${PHASE3_DIR}/06_extracted_sequences" \
        --species-map data/gca_to_species.tsv \
        --output-dir "${PHASE3_DIR}/08_matrices" \
        2>&1 | tee "logs/phase8a_${CLEAN_GENE_FAMILY}_${SLURM_ARRAY_TASK_ID}.log"; then
        echo "✓ Phase 8a complete"
    else
        echo "❌ Phase 8a failed"
        exit 1
    fi
fi

# Phase 8b: Generate summary matrices
echo ""
echo "--- PHASE 8b: GENERATE SUMMARY MATRICES ---"

if [ -f "${PHASE3_DIR}/08_matrices/summary_matrix.tsv" ]; then
    echo "✓ Phase 8b already complete (skipping)"
else
    if python scripts/08b_generate_summary_matrices.py \
        --locus-defs "${PHASE3_DIR}/locus_definitions.tsv" \
        --blocks "${PHASE3_DIR}/03_filtered_blocks/synteny_blocks_filtered.tsv" \
        --targets "${PHASE3_DIR}/05_classified/all_targets_classified.tsv" \
        --species-map data/gca_to_species.tsv \
        --extracted-seqs "${PHASE3_DIR}/06_extracted_sequences" \
        --output-dir "${PHASE3_DIR}/08_matrices" \
        2>&1 | tee "logs/phase8b_${CLEAN_GENE_FAMILY}_${SLURM_ARRAY_TASK_ID}.log"; then
        echo "✓ Phase 8b complete"
    else
        echo "❌ Phase 8b failed"
        exit 1
    fi
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
