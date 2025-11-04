#!/bin/bash
#SBATCH --job-name=synteny_batch
#SBATCH --array=0-31
#SBATCH --output=logs/batch_%A_%a.out
#SBATCH --error=logs/batch_%A_%a.err
#SBATCH --time=8:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --partition=general

# Array job to process all 32 gene families through Phases 1-5
# Each array task processes one gene family in parallel

cd /carc/scratch/projects/emartins/2016456/adam/synteny_scanning/hpc_deployment_package/hybrid_workflow

echo "=========================================="
echo "Batch Synteny Array Job"
echo "Job ID: $SLURM_ARRAY_JOB_ID"
echo "Task ID: $SLURM_ARRAY_TASK_ID"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo "=========================================="

# Read gene family data from input file
INPUT_FILE="hpc_input_loci.tsv"

# Skip header, get line number SLURM_ARRAY_TASK_ID + 2 (header is line 1, tasks start at 0)
LINE_NUM=$((SLURM_ARRAY_TASK_ID + 2))

# Read the specific line for this task
FAMILY_DATA=$(sed -n "${LINE_NUM}p" "$INPUT_FILE")

if [ -z "$FAMILY_DATA" ]; then
    echo "ERROR: Could not read line $LINE_NUM from $INPUT_FILE"
    exit 1
fi

# Parse TSV columns: gene_family, target_genes, num_locs
FAMILY_NAME=$(echo "$FAMILY_DATA" | cut -f1)
LOC_LIST=$(echo "$FAMILY_DATA" | cut -f2)
NUM_LOCS=$(echo "$FAMILY_DATA" | cut -f3)

echo ""
echo "Processing family: $FAMILY_NAME"
echo "Number of LOCs: $NUM_LOCS"
echo "LOC list: $LOC_LIST"
echo ""

# Run the processing script
bash process_gene_family.sh "$FAMILY_NAME" "$LOC_LIST"

EXIT_CODE=$?

echo ""
echo "=========================================="
echo "Task $SLURM_ARRAY_TASK_ID completed with exit code: $EXIT_CODE"
echo "Finished: $(date)"
echo "=========================================="

exit $EXIT_CODE
