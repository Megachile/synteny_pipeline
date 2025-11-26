#!/bin/bash
#SBATCH --job-name=phase6_extract_v2
#SBATCH --array=1-1
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --time=04:00:00
#SBATCH --partition=general
#SBATCH --output=logs/phase6_extract_v2_%A_%a.out
#SBATCH --error=logs/phase6_extract_v2_%A_%a.err

# Change to working directory
cd /carc/scratch/projects/emartins/2016456/adam/synteny_scanning/hpc_deployment_package/hybrid_workflow

# Activate environment (adjust if your cluster uses modules)
if command -v micromamba >/dev/null 2>&1; then
  eval "$(micromamba shell hook --shell bash)"
  micromamba activate busco_env || true
fi

# Load exonerate module - critical for Phase 6 extraction!
module load exonerate/2.4.0-yl7q || echo "WARNING: Could not load exonerate module"

# Verify exonerate availability
if ! command -v exonerate &> /dev/null; then
    echo "ERROR: exonerate not found in PATH"
    exit 1
fi

# Family list handling (one family per array task)
# Default list file: families.txt, or pass a custom file path as $1
if [ -n "$1" ]; then
    FAMILY_FILE="$1"
else
    FAMILY_FILE="families.txt"
fi

FAMILY=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$FAMILY_FILE")
if [ -z "$FAMILY" ]; then
    echo "ERROR: No family found for task ${SLURM_ARRAY_TASK_ID} in $FAMILY_FILE"
    exit 1
fi

echo "Processing Phase 6 (v2) extraction for family: $FAMILY (task ${SLURM_ARRAY_TASK_ID})"
echo "=============================================================="

# Use canonical ragtag_output directory (Python script defaults to this)
# But BK needs NCBI reference (NC_* chromosome names) - create minimal override
BK_NCBI="/carc/scratch/projects/emartins/2016456/adam/genomes/belonocnema_kinseyi/GCF_010883055.1_B_treatae_v1_genomic.fna"
BK_OVERRIDE_DIR="outputs/${FAMILY}/phase6_bk_override"

if [ -f "${BK_NCBI}" ]; then
    # Create BK override symlink (pipeline uses NC_* coords, ragtag has CM_* names)
    rm -rf "${BK_OVERRIDE_DIR}"
    mkdir -p "${BK_OVERRIDE_DIR}/Belonocnema_kinseyi_GCF"
    ln -s "${BK_NCBI}" "${BK_OVERRIDE_DIR}/Belonocnema_kinseyi_GCF/ragtag.scaffold.fasta"
    echo "BK genome: using NCBI NC_* FASTA override"
fi

# Inputs from redesigned phases
SYNTENIC_TSV="outputs/${FAMILY}/phase5_v2/syntenic_targets.tsv"
UNPLACEABLE_TSV="outputs/${FAMILY}/phase5_v2/unplaceable_targets.tsv"
COMBINED_TARGETS="outputs/${FAMILY}/phase4_v2/combined_targets.faa"
OUT_DIR="outputs/${FAMILY}/phase6_extracted_v2"

# Validate inputs
if [ ! -f "${SYNTENIC_TSV}" ]; then
    echo "ERROR: Missing syntenic targets: ${SYNTENIC_TSV}"
    exit 1
fi
if [ ! -f "${COMBINED_TARGETS}" ]; then
    echo "ERROR: Missing combined targets FASTA: ${COMBINED_TARGETS}"
    exit 1
fi

# Rotate old outputs
echo "Preparing output directory..."
rm -rf "${OUT_DIR}_OLD"
mv "${OUT_DIR}" "${OUT_DIR}_OLD" 2>/dev/null || true
mkdir -p "${OUT_DIR}"

# Run extraction (uses default genome-fasta-dir, with BK override if available)
echo "Running Phase 6 extraction (exonerate)..."
CMD=( python redesign_scripts/06_extract_sequences.py \
    --syntenic "${SYNTENIC_TSV}" \
    --query-proteins "${COMBINED_TARGETS}" \
    --output-dir "${OUT_DIR}" \
    --verbose )

# Add BK override directory if it was created
if [ -d "${BK_OVERRIDE_DIR}" ]; then
    CMD+=( --genome-override-dir "${BK_OVERRIDE_DIR}" )
fi

# Include strong unplaceables when available; extractor will filter by --unplaceable-evalue
if [ -f "${UNPLACEABLE_TSV}" ]; then
  CMD+=( --unplaceable "${UNPLACEABLE_TSV}" )
fi

"${CMD[@]}"

RC=$?
if [ ${RC} -ne 0 ]; then
    echo "ERROR: Phase 6 extraction failed for ${FAMILY} (rc=${RC})"
    exit ${RC}
fi

echo ""
echo "Phase 6 (v2) complete for family: ${FAMILY}"

# Clean up BK override directory
rm -rf "${BK_OVERRIDE_DIR}"

# Summarize outputs
CDS_COUNT=$(find "${OUT_DIR}" -name "*_gene1_cds.fasta" 2>/dev/null | wc -l)
PROT_COUNT=$(find "${OUT_DIR}" -name "*_gene1_protein.fasta" 2>/dev/null | wc -l)

echo "Extracted: ${CDS_COUNT} CDS files | ${PROT_COUNT} protein files"
if [ "${CDS_COUNT}" -gt 0 ]; then
    echo "Sample CDS files:"
    find "${OUT_DIR}" -name "*_gene1_cds.fasta" 2>/dev/null | head -3
fi

# Aggregate proteins into one FASTA
AGG="${OUT_DIR}/all_extracted_genes.faa"
> "${AGG}"
for pf in $(find "${OUT_DIR}" -name "*_gene1_protein.fasta" 2>/dev/null); do
    cat "$pf" >> "${AGG}"
done
TOTAL_SEQS=$(grep -c "^>" "${AGG}" 2>/dev/null || echo 0)
echo "Total extracted protein sequences: ${TOTAL_SEQS}"
