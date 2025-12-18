# Helixer Pipeline Session Status - 2025-12-02

## Session Summary

**Goal**: Implement proper tracking and filtering of novel loci discovered through Phase 5b validation.

**Completed**:
1. Added clustering logic to `05b_validate_novel_loci.py` to deduplicate novel loci
2. Fixed species name lookup (GCA IDs → species names)
3. Added CM→chr mapping for BK chromosomes (CM021338-CM021347 → chr1-chr10)
4. Added CM→chr mapping for T. remus native assembly (CM036336-CM036345 → chr1-chr10)
5. Added locus suffix logic (a, b, c...) for multiple loci on same chromosome
6. Exempt T. remus from chromosome comparison penalty in Phase 2b
7. Added `--min-total-genomes` filter to remove spurious single-genome hits
8. Deprecated `09_detect_novel_loci_helixer.py` (superseded by 05b)
9. **Phase 6**: Added `--phase5b-dir` argument to extract novel loci targets
10. **Phase 7**: Added `--phase5b-dir` argument to annotate novel loci flanking genes
11. **Phase 8a**: Added novel loci matrix generation (`{family}_NOVEL_{locus}_matrix.tsv`)

## Test Results (CAP_MC28)

### Before Filter
```
39 candidates → 7 unique novel loci

NOVEL_1: Alloxysta_arcuata_chr5_a     (1 member, 12 empty) - total: 13
NOVEL_2: Alloxysta_arcuata_chr8_a     (1 member, 12 empty) - total: 13
NOVEL_3: Synergus_coniferae_chr2_a    (11 members, 5 empty) - total: 16 ← strongest
NOVEL_4: Callaspidia_notata_chr5_a    (1 member, 6 empty) - total: 7
NOVEL_5: Telenomus_remus_chr8_a       (1 member, 0 empty) - total: 1 ← spurious
NOVEL_6: Telenomus_remus_chr8_b       (1 member, 0 empty) - total: 1 ← spurious
NOVEL_7: Telenomus_remus_chr8_c       (1 member, 0 empty) - total: 1 ← spurious
```

### After Filter (--min-total-genomes 3)
```
NOVEL_1: Alloxysta_arcuata_chr5_a     ✓ PASS (block in 13 genomes)
NOVEL_2: Alloxysta_arcuata_chr8_a     ✓ PASS (block in 13 genomes)
NOVEL_3: Synergus_coniferae_chr2_a    ✓ PASS (block in 16 genomes)
NOVEL_4: Callaspidia_notata_chr5_a    ✓ PASS (block in 7 genomes)
NOVEL_5-7: Telenomus_remus_chr8       ✗ FILTERED (block only in 1 genome)
```

**Interpretation**:
- NOVEL_1,2,4: Synteny block conserved across many genomes, but target rare = likely novel insertion
- NOVEL_3: Synteny block AND target both conserved = strong novel locus
- NOVEL_5-7: Block only in T. remus = spurious/outgroup-specific, filtered out

## Git Commits Today

```
b5158f2  Phase 8a: Add novel loci matrix generation
ba7e8f5  Phase 7: Add novel loci flanking gene SwissProt annotation
ccbfbe5  Phase 6: Add novel loci sequence extraction from Phase 5b
ff11833  Phase 5b: Add min_total_genomes filter for novel loci
2a7bced  Phase 5b: Add locus suffix (a, b, c...) for multiple loci on same chromosome
76989b7  Phase 5b: Novel loci clustering + T. remus chromosome exemption
```

## Key Files Modified

### 05b_validate_novel_loci.py
- `cluster_novel_loci()`: Deduplicate via mutual flanking synteny
- `load_chromosome_mapping()`: CM accessions → chr# conversion
- Locus naming: `Species_name_chr#_suffix` format
- `--min-total-genomes` filter (default: 3)

### 02b_detect_empty_blocks_helixer.py
- `NATIVE_CHROMOSOME_GENOMES` set for non-RagTag'd genomes
- T. remus exempt from chromosome comparison penalty
- New tier 2b: `native_chr_best` → MEDIUM confidence

### 06_extract_sequences_helixer.py
- Added `--phase5b-dir` argument
- Extract one target per member genome per novel locus (best bitscore)
- Output: `NOVEL_{locus_name}/` directories with protein FASTAs

### 07_annotate_flanking_swissprot.py
- Added `--phase5b-dir` argument
- Load `NOVEL_*_flanking.faa` files and include in DIAMOND run
- Output: `novel_loci_flanking_annotated.tsv`

### 08a_generate_locus_matrices_helixer.py
- Added `--phase5b-dir` argument
- New `create_novel_locus_matrix()` function
- Shows MEMBER, EMPTY_SYNTENY, NO_MATCH status per genome
- Output: `{family}_NOVEL_{locus_name}_matrix.tsv`

### run_helixer_phases_6_8_array.slurm
- Updated to pass `--phase5b-dir` to Phases 6, 7, 8a when available

### Deprecated
- `09_detect_novel_loci_helixer.py` → moved to deprecated/

## T. remus Chromosome Handling

**Problem**: T. remus (GCA_020615435.1) has native chromosome-level assembly, NOT RagTag'd against BK. Its chromosome numbers are NOT syntenic with BK.

**Solution**:
- Phase 2b: Exempt from chr comparison, get MEDIUM confidence (not penalized)
- Phase 5b: Use native chr numbers for naming, but filter by total genome count

## Novel Loci Filter Logic

A genuine novel locus should have:
1. **Conserved synteny block** across multiple genomes (even if target is rare)
2. Block in only 1 genome = spurious hit or outgroup-specific

Filter: `n_total_genomes >= 3` (members + empty)

## Chromosome Mapping Reference

### BK (Belonocnema kinseyi)
| RefSeq | GenBank | Chr |
|--------|---------|-----|
| NC_046657.1 | CM021338.1 | chr1 |
| NC_046658.1 | CM021339.1 | chr2 |
| NC_046659.1 | CM021340.1 | chr3 |
| NC_046660.1 | CM021341.1 | chr4 |
| NC_046661.1 | CM021342.1 | chr5 |
| NC_046662.1 | CM021343.1 | chr6 |
| NC_046663.1 | CM021344.1 | chr7 |
| NC_046664.1 | CM021345.1 | chr8 |
| NC_046665.1 | CM021346.1 | chr9 |
| NC_046666.1 | CM021347.1 | chr10 |

### T. remus (GCA_020615435.1) - Native Assembly
| GenBank | Chr |
|---------|-----|
| CM036336.1 | chr1 |
| CM036337.1 | chr2 |
| CM036338.1 | chr3 |
| CM036339.1 | chr4 |
| CM036340.1 | chr5 |
| CM036341.1 | chr6 |
| CM036342.1 | chr7 |
| CM036343.1 | chr8 |
| CM036344.1 | chr9 |
| CM036345.1 | chr10 |

## Next Steps

1. ~~**Integrate novel loci into Phases 6, 7, 8a**~~: ✓ Complete
2. **Make helixer pipeline standalone**: Copy Phase 1 dependencies into helixer_pipeline/
3. **Test on more families**: Validate filtering works across different gene families
4. **Run full pipeline**: Submit array job for all 33 families with novel loci support
