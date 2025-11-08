#!/bin/bash
# Surgical script: Apply Phase 7 grading to families that completed Phase 6
# Then regenerate Phase 9 matrices using graded targets

set -euo pipefail

cd /carc/scratch/projects/emartins/2016456/adam/synteny_scanning/hpc_deployment_package/hybrid_workflow

# Load environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate busco_env

echo "================================================================================"
echo "SURGICAL APPLICATION OF GRADING SYSTEM"
echo "================================================================================"
echo ""

# Families that completed Phase 6
FAMILIES=(
    "MC117"
    "ApoD_MC13"
    "Glutenin_MC21"
    "ester_hydrolase_MC19"
    "TIL_MC9"
    "ribonuclease_MC33"
    "OS-D_MC11"
    "alaserpin_MC35"
    "endoglucanase_Z_MC6"
    "MC211"
    "rhamnogalacturonate_lyase-like_MC47"
    "peptidoglycan-recognition_protein_SA-like___peptidoglycan-recognition_protein_LC-like_MC125"
    "natterin-4-like_MC59"
    "ferritin_MC102"
    "thioredoxin-2-like_MC108"
    "sprouty_MC112"
    "MC174"
    "extensin-like_MC192"
    "venom_serine_carboxypeptidase_MC203"
    "MC258"
    "pectin_lyase_PL"
)

echo "Applying grading to ${#FAMILIES[@]} families..."
echo ""

SUCCESS_COUNT=0
FAIL_COUNT=0

for family in "${FAMILIES[@]}"; do
    OUTPUT_DIR="outputs/${family}_BK_LB_batch"
    GRADED_FILE="${OUTPUT_DIR}/05_classified/syntenic_targets_graded.tsv"

    echo "--- Processing $family ---"

    # Check if Phase 6 completed
    if [ ! -f "${OUTPUT_DIR}/06_extracted_sequences/all_extracted_genes.faa" ]; then
        echo "  ⚠ Phase 6 not complete, skipping..."
        ((FAIL_COUNT++))
        continue
    fi

    # Phase 7: Grade matches
    if [ -f "$GRADED_FILE" ]; then
        echo "  ✓ Phase 7 already complete (graded file exists)"
    else
        echo "  Running Phase 7: Grade match quality..."
        if python scripts/07_grade_matches.py \
            --syntenic "${OUTPUT_DIR}/05_classified/syntenic_targets.tsv" \
            --extracted-dir "${OUTPUT_DIR}/06_extracted_sequences" \
            --query-proteins "${OUTPUT_DIR}/04_target_genes/combined_targets.faa" \
            2>&1 | tee "logs/grade_surgical_${family}.log"; then

            if [ -f "$GRADED_FILE" ]; then
                echo "  ✓ Phase 7 complete"

                # Count grades
                INTACT=$(grep -c "intact" "$GRADED_FILE" || echo 0)
                DEGRADED_FRAG=$(grep -c "degraded_fragment" "$GRADED_FILE" || echo 0)
                DEGRADED_PSEUDO=$(grep -c "degraded_pseudogene" "$GRADED_FILE" || echo 0)
                DEGRADED_NO_CDS=$(grep -c "degraded_no_cds" "$GRADED_FILE" || echo 0)
                BAD=$(grep -c "bad_match" "$GRADED_FILE" || echo 0)

                echo "    Grades: intact=$INTACT, degraded_fragment=$DEGRADED_FRAG, pseudogene=$DEGRADED_PSEUDO, no_cds=$DEGRADED_NO_CDS, bad=$BAD"
            else
                echo "  ❌ Phase 7 failed - graded file not created"
                ((FAIL_COUNT++))
                continue
            fi
        else
            echo "  ❌ Phase 7 failed"
            ((FAIL_COUNT++))
            continue
        fi
    fi

    # Phase 9a: Regenerate locus matrices with graded targets
    echo "  Regenerating Phase 9a matrices with graded targets..."

    # Back up old matrices
    TIMESTAMP=$(date +%Y%m%d_%H%M%S)
    if [ -d "${OUTPUT_DIR}/08_matrices" ]; then
        mv "${OUTPUT_DIR}/08_matrices" "${OUTPUT_DIR}/08_matrices_OLD_${TIMESTAMP}"
    fi

    # Check if Phase 8 (SwissProt) completed
    if [ ! -f "${OUTPUT_DIR}/07_swissprot_annotations/genome_specific_swissprot_annotations.tsv" ]; then
        echo "  ⚠ Phase 8 (SwissProt) not complete, skipping matrix regeneration..."
        ((FAIL_COUNT++))
        continue
    fi

    # Copy locus definitions
    cp "${OUTPUT_DIR}/phase12_landmark/locus_definitions.tsv" "${OUTPUT_DIR}/locus_definitions.tsv"

    # Run Phase 9a with grading
    if python scripts/08a_generate_locus_matrices.py \
        --locus-defs "${OUTPUT_DIR}/locus_definitions.tsv" \
        --synteny-dir "${OUTPUT_DIR}/02_synteny_blocks" \
        --blocks "${OUTPUT_DIR}/03_filtered_blocks/synteny_blocks_filtered.tsv" \
        --targets "${OUTPUT_DIR}/05_classified/all_targets_classified.tsv" \
        --graded-syntenic "$GRADED_FILE" \
        --keep-grades "intact,degraded_fragment" \
        --swissprot "${OUTPUT_DIR}/07_swissprot_annotations/genome_specific_swissprot_annotations.tsv" \
        --reference-proteins data/reference/protein.faa \
        --extracted-seqs "${OUTPUT_DIR}/06_extracted_sequences" \
        --species-map data/gca_to_species.tsv \
        --output-dir "${OUTPUT_DIR}/08_matrices" \
        2>&1 | tee "logs/phase9a_surgical_${family}.log"; then

        echo "  ✓ Phase 9a complete (with grading)"
    else
        echo "  ❌ Phase 9a failed"
        ((FAIL_COUNT++))
        continue
    fi

    # Phase 9b: Generate summary matrices
    if python scripts/08b_generate_summary_matrices.py \
        --locus-defs "${OUTPUT_DIR}/locus_definitions.tsv" \
        --blocks "${OUTPUT_DIR}/03_filtered_blocks/synteny_blocks_filtered.tsv" \
        --targets "${OUTPUT_DIR}/05_classified/all_targets_classified.tsv" \
        --species-map data/gca_to_species.tsv \
        --extracted-seqs "${OUTPUT_DIR}/06_extracted_sequences" \
        --output-dir "${OUTPUT_DIR}/08_matrices" \
        2>&1 | tee "logs/phase9b_surgical_${family}.log"; then

        echo "  ✓ Phase 9b complete"
        ((SUCCESS_COUNT++))
    else
        echo "  ❌ Phase 9b failed"
        ((FAIL_COUNT++))
    fi

    echo ""
done

echo "================================================================================"
echo "SURGICAL GRADING COMPLETE"
echo "================================================================================"
echo "  Successful: $SUCCESS_COUNT families"
echo "  Failed: $FAIL_COUNT families"
echo "================================================================================"
