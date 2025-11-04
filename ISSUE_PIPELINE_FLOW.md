# Issue: Pipeline Flow Confusion - Multiple "Step 01" Scripts

## Problem

The pipeline has TWO different scripts both claiming to be "Step 01":

1. **`01_discover_paralogs_gene_level.py`** - NEW approach
   - Takes LOC IDs as input
   - Detects tandem duplications
   - Uses DIAMOND to find paralogs in landmark genomes
   - Outputs: `unique_loci.tsv`, flanking protein FAA files

2. **`01_extract_proteins.py`** - OLD approach
   - Takes locus definitions with protein IDs
   - Extracts specific proteins from reference
   - Outputs: `LOCUS_ID/LOCUS_ID_targets.faa`

## Impact

- Batch wrapper (`process_gene_family.sh`) calls `01_discover_paralogs_gene_level.py`
- But `04_blast_targets.py` expects output from `01_extract_proteins.py`
- This causes batch runs to fail when Phase 4 can't find expected files

## Working Example

The ferritin test works because it:
1. Uses a manually created `locus_definitions.tsv`
2. Ran `01_extract_proteins.py` (the OLD script)
3. Created directory structure: `01_extracted_proteins/BK_chr2_a/BK_chr2_a_targets.faa`

## Failed Example

Batch runs fail because:
1. Wrapper calls `01_discover_paralogs_gene_level.py` (the NEW script)
2. Creates: `01_paralogs/unique_loci.tsv` and `*_flanking.faa` files
3. Phase 4 looks for: `01_extracted_proteins/LOCUS_ID/LOCUS_ID_targets.faa`
4. Files don't exist → Phase 4 fails

## Questions

1. **Which script is the canonical Step 01?**
2. **Do we need both approaches or should we consolidate?**
3. **What's the complete, correct pipeline flow for batch processing?**

## Next Steps

1. Manually trace through MC6 to understand correct flow
2. Document the actual working pipeline sequence
3. Either:
   - Update wrapper to use correct script sequence, OR
   - Modify scripts to have compatible inputs/outputs
4. Remove or rename deprecated scripts to avoid confusion
