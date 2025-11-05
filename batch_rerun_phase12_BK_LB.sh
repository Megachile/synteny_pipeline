#!/bin/bash
#SBATCH --job-name=batch_ph12_BK_LB
#SBATCH --output=logs/batch_ph12_%A_%a.out
#SBATCH --error=logs/batch_ph12_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --mem=8GB
#SBATCH --cpus-per-task=4
#SBATCH --array=0-31  # 32 gene families from batch_run

cd /carc/scratch/projects/emartins/2016456/adam/synteny_scanning/hpc_deployment_package/hybrid_workflow

# Activate environment
source ~/.bashrc
micromamba activate busco_env

# Gene families from batch_run (in alphabetical order)
families=(
    "alaserpin_MC35"
    "ApoD_MC13"
    "CAP_MC28"
    "carboxypeptidase_MC26"
    "der f 21_MC2"
    "endoglucanase Z_MC6"
    "ester hydrolase_MC19"
    "extensin-like_MC192"
    "ferritin_MC102"
    "Glutenin_MC21"
    "LRR 15_MC91"
    "LRR 4_MC10"
    "MC117"
    "MC174"
    "MC18"
    "MC211"
    "MC258"
    "MC3"
    "natterin-4-like_MC59"
    "neurexin_MC14"
    "OS-D_MC11"
    "pectin lyase_PL"
    "peptidoglycan-recognition protein SA-like | peptidoglycan-recognition protein LC-like_MC125"
    "rhamnogalacturonate lyase-like_MC47"
    "ribonuclease_MC33"
    "sprouty_MC112"
    "thioredoxin-2-like_MC108"
    "TIL_MC9"
    "venom acid phosphatase Acph-1-like, transcript variant X1_MC63"
    "venom R-like_MC1"
    "venom serine carboxypeptidase_MC203"
    "venom serine protease 34-like_MC175"
)

family="${families[$SLURM_ARRAY_TASK_ID]}"

echo "========================================="
echo "Processing gene family: $family"
echo "Array task ID: $SLURM_ARRAY_TASK_ID"
echo "========================================="

# Read LOC IDs from gene_family_loc_map.tsv
loc_map_file="gene_family_loc_map.tsv"

if [ ! -f "$loc_map_file" ]; then
    echo "ERROR: LOC map file not found: $loc_map_file"
    echo "Run extract_loc_ids_from_batch.py first"
    exit 1
fi

# Extract LOC IDs for this family (skip header, find matching family)
loc_ids=$(awk -F'\t' -v fam="$family" '$1 == fam {print $2}' "$loc_map_file")

if [ -z "$loc_ids" ]; then
    echo "ERROR: No LOC IDs found for family: $family"
    exit 1
fi

echo "LOC IDs: $loc_ids"

# Create output directory
output_dir="outputs/${family}_BK_LB_rerun"
mkdir -p "$output_dir"

# Phase 1: Discover paralogs (BK+LB only)
# Using genome-specific distance caps: BK: 1500 kb, LB: 400 kb, MAX_FLANKING_GENES: 50
echo ""
echo "[Phase 1] Discovering paralogs..."
python scripts/01_discover_paralogs_gene_level.py \
    --loc-ids "$loc_ids" \
    --output-dir "$output_dir/phase12_landmark" \
    --gene-family "$family" \
    --evalue 1e-3 \
    --min-pident 40

if [ $? -ne 0 ]; then
    echo "ERROR: Phase 1 failed for $family"
    exit 1
fi

# Phase 2: Compare synteny
echo ""
echo "[Phase 2] Comparing locus synteny..."
python scripts/02_compare_synteny_gene_level.py "$output_dir/phase12_landmark"

if [ $? -ne 0 ]; then
    echo "ERROR: Phase 2 failed for $family"
    exit 1
fi

echo ""
echo "========================================="
echo "COMPLETED: $family"
echo "========================================="
