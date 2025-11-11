#!/bin/bash
#SBATCH --job-name=phase7_swissprot_v2
#SBATCH --array=1-1
#SBATCH --cpus-per-task=16
#SBATCH --mem=20G
#SBATCH --time=04:00:00
#SBATCH --partition=general
#SBATCH --output=logs/phase7_v2_%A_%a.out
#SBATCH --error=logs/phase7_v2_%A_%a.err

# Change to working directory
cd /carc/scratch/projects/emartins/2016456/adam/synteny_scanning/hpc_deployment_package/hybrid_workflow

# Activate environment (same as working script)
source ~/.bashrc
micromamba activate busco_env

# Family list handling (one family per array task)
if [ -n "$1" ]; then
  FAMILY_FILE="$1"
else
  FAMILY_FILE="families.txt"
fi

FAMILY=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$FAMILY_FILE")
if [ -z "$FAMILY" ]; then
  echo "ERROR: No family found for task ${SLURM_ARRAY_TASK_ID} in $FAMILY_FILE" >&2
  exit 1
fi

echo "Processing Phase 7 (v2) SwissProt annotation for family: ${FAMILY} (task ${SLURM_ARRAY_TASK_ID})"
echo "==========================================================================="

PHASE1_DIR="outputs/${FAMILY}/phase1_v2"
PHASE2_DIR="outputs/${FAMILY}/phase2_synteny_v2"
PHASE3_DIR="outputs/${FAMILY}/phase3_v2"
OUT_DIR="outputs/${FAMILY}/phase7_v2"
mkdir -p "${OUT_DIR}"

LOCUS_DEFS="${PHASE1_DIR}/locus_definitions.tsv"
if [ ! -f "${LOCUS_DEFS}" ]; then
  echo "ERROR: Missing locus_definitions.tsv: ${LOCUS_DEFS}" >&2
  exit 1
fi

# SwissProt DB resolution: prefer in-repo DB, fallback to shared path
SWISSPROT_DB="data/databases/swissprot.dmnd"
if [ ! -f "${SWISSPROT_DB}" ]; then
  ALT_DB="/carc/scratch/projects/emartins/2016456/adam/databases/uniprot_sprot.dmnd"
  if [ -f "${ALT_DB}" ]; then
    SWISSPROT_DB="${ALT_DB}"
  else
    echo "WARNING: SwissProt DB not found at ${SWISSPROT_DB} nor ${ALT_DB}. The step will be skipped." >&2
    exit 0
  fi
fi

# Determine filtered blocks file for this family (best-per-genome required by annotate script)
FILTERED="${PHASE5_DIR:-outputs/${FAMILY}/phase5_v2}/selected_synteny_blocks.tsv"
if [ ! -f "${FILTERED}" ]; then
  # Fallback to traditional best-per-genome from Phase 3 or derive one if needed
  FILTERED="${PHASE3_DIR}/synteny_blocks_filtered.tsv"
fi
if [ ! -f "${FILTERED}" ]; then
  THR="${PHASE3_DIR}/all_synteny_blocks_ge3.tsv"
  if [ -f "${THR}" ]; then
    echo "Creating best-per-genome filtered blocks from threshold-only file..."
    python - "$THR" "$FILTERED" <<'PY'
import sys, pandas as pd
thr = sys.argv[1]; out = sys.argv[2]
df = pd.read_csv(thr, sep='\t')
# choose count column
cnt = None
for c in ('num_query_matches','num_target_proteins','num_proteins'):
    if c in df.columns:
        cnt = c; break
if cnt is None:
    cnt = 'num_query_matches'; df[cnt] = 0
# best per (locus, genome)
idx = df.groupby(['locus_id','genome'])[cnt].idxmax()
best = df.loc[idx].copy()
best.to_csv(out, sep='\t', index=False)
print(f"Wrote best-per-genome: {out} ({len(best)})")
PY
  else
    COMB="${PHASE2_DIR}/combined_synteny_blocks.tsv"
    if [ -f "${COMB}" ]; then
      echo "Creating best-per-genome filtered blocks from combined synteny blocks..."
      python - "$COMB" "$FILTERED" <<'PY'
import sys, pandas as pd
inp = sys.argv[1]; out = sys.argv[2]
df = pd.read_csv(inp, sep='\t')
cnt = None
for c in ('num_query_matches','num_target_proteins','num_proteins'):
    if c in df.columns:
        cnt = c; break
if cnt is None:
    cnt = 'num_query_matches'; df[cnt] = 0
idx = df.groupby(['locus_id','genome'])[cnt].idxmax()
best = df.loc[idx].copy()
best.to_csv(out, sep='\t', index=False)
print(f"Wrote best-per-genome: {out} ({len(best)})")
PY
    else
      echo "ERROR: No blocks file found to derive best-per-genome (missing ${THR} and ${COMB})." >&2
      exit 1
    fi
  fi
fi

# Iterate loci in this family and annotate
LOCI=$(awk 'NR>1{print $1}' "${LOCUS_DEFS}")
for L in ${LOCI}; do
  echo "Annotating locus: ${L}"
  python scripts/annotate_swissprot.py \
    --synteny-dir "${PHASE2_DIR}" \
    --filtered-blocks "${FILTERED}" \
    --swissprot-db "${SWISSPROT_DB}" \
    --output-dir "${OUT_DIR}" \
    --locus "${L}" \
    --evalue 1e-5 \
    --threads ${SLURM_CPUS_PER_TASK:-16}
done

echo "Phase 7 (v2) complete for family: ${FAMILY}"

