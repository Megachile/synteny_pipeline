# Proposed Pipeline Structure

## What the Pipeline Actually Does

### Current Flow (with confusing numbering)

```
INPUT: List of BK LOC IDs for a gene family

Phase 1:  Define loci from BK reference (LOC → target + 20 flanking genes)
Phase 1b: Add LB orthologs to locus definitions
    ↓
    ├──────────────────────────────────┐
    ↓                                  ↓
Phase 4:  DIAMOND targets vs all      Phase 2b: DIAMOND BK flanking vs all
          Helixer proteomes                     Helixer proteomes
    ↓                                           (finds synteny blocks)
Phase 4b: Extract flanking around              ↓
          each target hit                       │
    ↓                                           │
Phase 5:  Classify targets by        ←──────────┘
          flanking synteny                (2b concordance check)
    ↓
Phase 5b: Cluster unplaceables into novel loci
    ↓
Phase 6:  Extract protein sequences
    ↓
Phase 7:  Annotate flanking genes (from FASTA headers)
    ↓
Phase 8a: Per-locus matrices (genome × flanking presence)
Phase 8b: Summary matrix (genome × all loci)

OUTPUT: Classified targets + matrices showing locus presence across genomes
```

---

## Proposed Renaming

| Old Name | What It Does | Proposed Name |
|----------|--------------|---------------|
| `01_phase1.py` | Parse BK GFF, extract target + flanking for each LOC ID | `01_define_loci.py` |
| `01b_merge_bk_lb_loci.py` | BLAST BK targets vs LB to find orthologs | `02_add_outgroup.py` |
| `04_detect_targets_helixer.py` | DIAMOND BK targets vs all proteomes | `03_find_targets.py` |
| `04b_extract_flanking_helixer.py` | Get flanking genes around each hit | `04_get_flanking.py` |
| `02b_detect_empty_blocks_helixer.py` | DIAMOND BK flanking vs all, cluster into blocks | `05_find_synteny_blocks.py` |
| `05_classify_targets_helixer.py` | Score targets by flanking match to BK loci | `06_classify_targets.py` |
| `05b_validate_novel_loci.py` | Cluster unplaceables with shared flanking | `07_find_novel_loci.py` |
| `06_extract_sequences_helixer.py` | Pull protein sequences for all targets | `08_extract_sequences.py` |
| `07_annotate_flanking.py` | Parse NR annotations from headers | `09_annotate_flanking.py` |
| `08a_generate_locus_matrices_helixer.py` | Per-locus presence matrix | `10a_locus_matrices.py` |
| `08b_generate_summary_matrices_helixer.py` | Family summary matrix | `10b_summary_matrix.py` |

---

## Plain English Pipeline Description

### Step 1: Define Reference Loci
From user-provided LOC IDs, look up each gene in the BK reference genome. For each:
- Extract the target protein sequence
- Extract 10 upstream + 10 downstream flanking gene proteins
- Record chromosome location and gene boundaries

### Step 2: Add Outgroup Reference
BLAST BK targets against Leptopilina boulardi (parasitoid outgroup) to find orthologs. This helps with loci that might be BK-specific vs broadly conserved.

### Step 3: Find Target Homologs
DIAMOND BLAST the BK target proteins against every Helixer-annotated proteome. Collect hits above threshold as candidate targets.

### Step 4: Get Flanking Context
For each candidate target hit, extract the N genes upstream and downstream from the Helixer GFF. These flanking genes will be used for synteny scoring.

### Step 5: Find Synteny Blocks
DIAMOND BLAST the BK flanking proteins against all Helixer proteomes. Cluster co-located hits into "synteny blocks" - regions that share the BK flanking pattern. Blocks without a target hit are "empty blocks" (gene loss or annotation gap).

### Step 6: Classify Targets
For each target from Step 3:
- BLAST its flanking genes against BK proteome
- Count how many match each locus's flanking set
- Assign to best-matching locus if above threshold
- Mark as "unplaceable" if no good match

### Step 7: Discover Novel Loci
Take unplaceable targets that have good BLAST scores but don't match any BK locus. Check if their flanking pattern appears in multiple genomes. If so, cluster them as a "novel locus" not present in the BK reference.

### Step 8: Extract Sequences and Gene Structure
Pull sequences and structural information from Helixer outputs for all classified targets.

**Outputs per target:**

| File | Content | Use Case |
|------|---------|----------|
| `*_protein.fasta` | Amino acid sequence | Phylogenetics, domain analysis |
| `*_cds.fasta` | Spliced coding nucleotide sequence | dN/dS selection analysis, codon usage |
| `*_structure.tsv` | Intron positions and lengths | Gene structure evolution |

**Structure TSV columns:**
```
target_id  genome  scaffold  strand  gene_start  gene_end  n_exons  exon_lengths  n_introns  intron_lengths  intron_phases  total_cds_length
```

Where:
- `exon_lengths` is semicolon-delimited (e.g., "245;187;412")
- `intron_lengths` is semicolon-delimited (e.g., "1523;892;2104")
- `intron_phases` indicates reading frame interruption (0/1/2 for each intron)

**Implementation:**
1. Parse Helixer GFF3 for each target gene to get exon/CDS coordinates
2. Extract CDS regions from genome assembly, splice together, handle strand
3. Record intron positions (gaps between exons) with lengths
4. Calculate intron phases from CDS frame

**Also extract:**
- NR annotation from each target's FASTA header for identity confirmation
- Add `helixer_annotation` column to `length_qc.tsv`

This sets up downstream positive selection analysis (dN/dS requires CDS alignments) and gene structure comparisons.

### Step 9: Annotate Flanking
Parse NR functional annotations from the Helixer proteome FASTA headers. Add to flanking gene tables for matrix display.

### Step 10: Generate Matrices
Create output tables:
- **Locus matrices**: One per locus, showing which flanking genes are present in each genome
- **Summary matrix**: One per family, showing all loci × all genomes with target presence/absence

---

## Dependency Graph (for parallel SLURM)

```
01_define_loci
      │
      ├────────────────┐
      ↓                ↓
02_add_outgroup   03_find_targets
                       │
                       ↓
                  04_get_flanking
                       │
      ┌────────────────┼────────────────┐
      ↓                ↓                ↓
05_find_synteny   06_classify     (parallel)
      │                │
      │                ↓
      │         07_find_novel
      │                │
      └───────┬────────┘
              ↓
        08_extract_sequences
              │
              ↓
        09_annotate_flanking
              │
              ↓
        10_generate_matrices
```

**Parallelizable steps:**
- Steps 3-5 could potentially run in parallel (all depend only on step 1/2)
- The bottleneck is step 6 which needs both step 4 and step 5 outputs

---

## Expected Runtime (per family)

| Step | Operation | Est. Time |
|------|-----------|-----------|
| 01 | GFF parsing, FASTA extraction | <5s |
| 02 | Single DIAMOND search | 10-30s |
| 03 | DIAMOND vs 33 proteomes | 1-2 min |
| 04 | GFF parsing per genome | <30s |
| 05 | DIAMOND vs 33 proteomes + clustering | 1-2 min |
| 06 | DIAMOND + scoring | 30-60s |
| 07 | DIAMOND + clustering | 30-60s |
| 08 | FASTA I/O | <10s |
| 09 | Header parsing | <5s |
| 10 | Pandas operations | <10s |

**Total: ~4-6 minutes per family**

With 33 families and 30-minute SLURM allocation, array jobs should complete with margin.
