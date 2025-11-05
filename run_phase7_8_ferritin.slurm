#!/bin/bash
#SBATCH --job-name=ferritin_ph7_8
#SBATCH --output=ferritin_ph7_8_%j.out
#SBATCH --error=ferritin_ph7_8_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G

cd /carc/scratch/projects/emartins/2016456/adam/synteny_scanning/hpc_deployment_package/hybrid_workflow

# Activate environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate busco_env

echo "===== PHASE 7: SwissProt Annotation ====="
python3 scripts/07_swissprot_annotation.py \
  --extracted-dir outputs/ferritin_phase3_real/06_extracted_sequences \
  --swissprot-db data/databases/uniprot_sprot.fasta \
  --output-dir outputs/ferritin_phase3_real/07_swissprot_annotations \
  --threads 8

echo ""
echo "===== PHASE 8a: Generate Locus Matrices ====="
python3 scripts/08a_generate_locus_matrices.py \
  --locus-defs outputs/ferritin_phase3_real/locus_definitions.tsv \
  --synteny-dir outputs/ferritin_phase3_real/02_synteny_blocks \
  --blocks outputs/ferritin_phase3_real/03_filtered_blocks/synteny_blocks_filtered.tsv \
  --targets outputs/ferritin_phase3_real/05_classified/all_targets_classified.tsv \
  --swissprot outputs/ferritin_phase3_real/07_swissprot_annotations/genome_specific_swissprot_annotations.tsv \
  --extracted-seqs outputs/ferritin_phase3_real/04_target_genes/extracted_sequences \
  --species-map data/gca_to_species.tsv \
  --output-dir outputs/ferritin_phase3_real/08_matrices

echo ""
echo "===== PHASE 8b: Generate Summary Matrices ====="
python3 scripts/08b_generate_summary_matrices.py \
  --locus-matrices-dir outputs/ferritin_phase3_real/08_matrices \
  --output-dir outputs/ferritin_phase3_real/08_matrices

echo ""
echo "===== FINAL SUMMARY: BK Genome Results ====="
echo ""
echo "Matrix files created:"
ls -lh outputs/ferritin_phase3_real/08_matrices/*.tsv

echo ""
echo "Checking summary matrix for BK (GCA_010883055.1):"
grep "GCA_010883055.1" outputs/ferritin_phase3_real/08_matrices/*summary*.tsv 2>/dev/null || echo "Not found in summary"
