# Synteny Scanning Hybrid Workflow - DETAILED BLUEPRINT

## Overview

This pipeline detects conserved gene families across 75 wasp genomes using synteny-based comparative genomics. The approach: (1) find paralogs of target genes in 4 landmark genomes (BK, LB, TR, DR), (2) extract flanking proteins around each locus, (3) use flanking proteins to detect syntenic regions in 75 target genomes, (4) find the actual target gene within each syntenic region, (5) extract complete gene structures with Exonerate, (6) annotate for functional context, (7) generate presence/absence matrices.

**Key insight:** Synteny (conserved gene order) helps identify true orthologs vs. distant paralogs or horizontal transfers.

---

## Phase 1: Paralog Discovery & Locus Definition

**Script:** `scripts/01_discover_paralogs_gene_level.py`

**Input:** Single LOC ID (e.g., `LOC117167432`) or protein ID from B. kinseyi

**Purpose:** Find all paralogs of a query gene across 4 landmark genomes and extract their flanking proteins for synteny detection

### Detailed Workflow:

#### Step 1: Build Genome Mappers (lines 350-354)
For each of 4 landmark genomes (BK, LB, TR, DR):

**Creates a `GenomicMapper` object that:**
1. Parses GFF files to build protein→gene and gene→position mappings
2. Handles two annotation formats:
   - **NCBI format (BK, LB):** Uses LOC IDs (e.g., `LOC117167432`), protein IDs (e.g., `XP_033211301.1`)
     - Extracts from GFF: `Name=LOC117167432` for genes, `protein_id=XP_...;gene=LOC...` for CDS
   - **BRAKER3 format (TR, DR):** Uses gene IDs (e.g., `g10440`), transcript IDs (e.g., `g10440.t1`)
     - Extracts from GFF3: `ID=g10440;` for genes, `ID=g10440.t1;Parent=g10440;` for mRNA
3. Builds data structures:
   - `protein_to_gene`: Map each protein/transcript ID to its parent gene
   - `gene_to_proteins`: Map each gene to all its isoforms (handles alternative splicing)
   - `gene_positions`: Map each gene to (chromosome, genomic_position)
   - `chr_genes`: For each chromosome, list of (position, gene_id) sorted by position

**Why this matters:** Multiple protein isoforms (e.g., XP_033211301.1, XP_033211301.2) come from the same LOC. We need to search at GENE level, not protein level, to avoid counting isoforms as paralogs.

#### Step 2: Find Query Protein in BK (lines 356-380)
1. If input is LOC ID (starts with `LOC`):
   - Look up all proteins for this gene in `gene_to_proteins`
   - Pick first protein found in BK proteome as representative
2. If input is protein ID:
   - Look up directly in BK proteome
   - Map back to gene ID using `protein_to_gene`
3. Extract full protein sequence from BK proteome FASTA

**Output:** Query protein sequence + its gene ID

#### Step 3: Find Paralogs in All Landmark Genomes (lines 382-392)

**For each genome (BK, LB, TR, DR):**

1. **Run DIAMOND blastp:**
   - Query: BK protein sequence (from Step 2)
   - Database: Landmark genome proteome (`.dmnd` file)
   - Parameters: E-value ≤ 1e-3, max 100 targets, 2 threads
   - Output: TSV with 12 columns (qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore)

2. **Filter hits:**
   - Keep only hits with ≥40% protein identity (`MIN_PIDENT`)
   - This avoids distant homologs that aren't true paralogs

3. **Map protein hits to genes:**
   - For each hit protein ID (sseqid), look up parent gene using `mapper.get_gene_for_protein()`
   - Handles isoforms: if `XP_123.1` and `XP_123.2` both hit, they collapse to single gene `LOC123`

4. **Deduplicate to unique genes:**
   - Result: List of unique gene IDs that are paralogs of the query

**Example output:**
```
BK: Found 45 protein hits (≥40% pident) → 12 unique genes
LB: Found 38 protein hits (≥40% pident) → 8 unique genes
TR: Found 52 protein hits (≥40% pident) → 15 unique genes
DR: Found 41 protein hits (≥40% pident) → 11 unique genes
```

#### Step 4: Process Each Unique Locus (lines 394-451)

**For each paralog gene found:**

**A. Extract flanking genes (lines 403-408):**
- Get chromosome and position of target gene
- Find all genes on same chromosome within 1000 kb (1 Mb) upstream/downstream
- Uses `get_flanking_genes(target_gene, flanking_distance_kb=1000)`
- **For NCBI genomes:** Uses actual genomic positions from GFF
- **For BRAKER3 genomes:** Uses protein order as proxy (no scaffold coordinates in some annotations)

**B. Detect tandem duplications (lines 410-411):**
- Check if any flanking genes are ALSO in the paralog list
- Tandem = when multiple paralogs sit adjacent on same chromosome
- Example: `LOC117167432` and `LOC117167433` sitting next to each other
- Flags: `tandem_count`, `is_tandem` (boolean)

**C. Get chromosome info (lines 413-417):**
- Extract chromosome ID from gene_positions
- Used for systematic locus naming

**D. Create systematic locus name (lines 419-421):**
- Uses `name_locus(genome_code, chromosome, locus_number)`
- Format: `<genome>_<chr>_<letter>`
- Examples:
  - BK chromosome 2, first locus → `BK_chr2_a`
  - BK chromosome 2, second locus → `BK_chr2_b`
  - LB scaffold NW_026137864.1 → `LB_scf7864_a` (uses last 4 digits)
  - BRAKER3 scaffolds → `TR_scf123_a`

**E. Extract flanking proteins (lines 423-433):**
1. For each flanking gene, get representative protein:
   - If gene has multiple isoforms, pick first alphabetically (for consistency)
   - Example: gene `LOC123` has proteins `XP_123.1`, `XP_123.2` → picks `XP_123.1`
2. Extract protein sequence from proteome FASTA
3. Save to `<locus_name>_flanking.faa`
   - Header format: `>XP_033211300.1|LOC117167431` (protein_id|gene_id)

**Why flanking proteins matter:** These define the conserved genomic neighborhood. If we find these same proteins in similar order in another genome, we've found a syntenic region where the target gene likely exists.

**F. Create locus entry (lines 435-451):**
Dictionary with:
- `locus_id`: Systematic name (e.g., `BK_chr3_a`)
- `genome`: Genome code (BK, LB, TR, or DR)
- `chromosome`: Chromosome/scaffold ID
- `target_gene`: Gene ID of the paralog at this locus
- `flanking_count`: Number of flanking genes found
- `tandem_count`: Number of tandem duplicates detected
- `is_tandem`: Boolean flag
- `flanking_file`: Path to `*_flanking.faa` file

#### Step 5: Save Results (lines 453-470)

**A. unique_loci.tsv:**
```tsv
locus_id	genome	chromosome	target_gene	flanking_count	tandem_count	is_tandem	flanking_file
BK_chr2_a	BK	NC_046658.1	LOC117167432	156	0	False	outputs/.../BK_chr2_a_flanking.faa
BK_chr3_a	BK	NC_046659.1	LOC117169185	142	0	False	outputs/.../BK_chr3_a_flanking.faa
LB_scf7864_a	LB	NW_026137864.1	LOC12345	98	1	True	outputs/.../LB_scf7864_a_flanking.faa
```

**B. paralog_details.json:**
Full DIAMOND results for each genome including all protein hits before gene-level deduplication

**C. *_flanking.faa files:**
One FASTA file per locus with all flanking proteins (typically 50-200 proteins per locus)

### Key Parameters:
- `FLANKING_DISTANCE_KB = 1000` (1 Mb window)
- `EVALUE_THRESHOLD = 1e-3` (DIAMOND search)
- `MIN_PIDENT = 40.0` (minimum 40% identity)
- `MIN_SYNTENY_MATCHES = 9` (not used in Phase 1, but defined for downstream)

### Outputs:
```
outputs/<gene_family>_phase12/loc_<LOC>/
├── unique_loci.tsv                 # Summary table
├── paralog_details.json            # Full BLAST results
├── BK_chr2_a_flanking.faa         # Flanking proteins
├── BK_chr3_a_flanking.faa
├── LB_scf7864_a_flanking.faa
└── ... (one .faa per locus)
```

### Notes:
- This phase is **optional** if you already have locus definitions
- Can be run per-LOC in parallel for batch processing
- For ferritin: 1 LOC (LOC117167432) → 9 loci across 4 genomes
- The `unique_loci.tsv` becomes the `locus_definitions.tsv` input for Phase 3+

---

## Phase 2: Synteny-Based Deduplication

**Status:** DEPRECATED - now handled in Phase 1

**Original purpose:** Merge duplicate loci detected across landmarks by comparing their synteny

**Why deprecated:** Phase 1 now handles deduplication at gene level using GFF mappings, making this step redundant

---

## Phase 3: Genome-Wide Synteny Detection

**Script:** `scripts/02_synteny_detection.py` *(numbering mismatch!)*

**Input:** `locus_definitions.tsv` (from Phase 1 or pre-built)

**Purpose:** Use flanking proteins to find syntenic regions across 75 target genomes

### Detailed Workflow:

#### Step 1: Load Locus Definitions (lines 344-356)
- Read `locus_definitions.tsv` from `SYNTENY_INPUT_FILE` env variable
- If `--locus <locus_id>` arg provided, filter to that specific locus (for SLURM array jobs)
- **Locus definition format:**
  ```tsv
  locus_id	genome	chromosome	flanking_file	gene_family
  BK_chr2_a	BK	NC_046658.1	.../BK_chr2_a_flanking.faa	ferritin_MC102
  ```

#### Step 2: Find Genome Databases (lines 358-364)
- Scan `data/ragtag_dbs/` for `.nhr` files (BLAST database headers)
- Each `.nhr` file indicates a BLAST-formatted nucleotide database
- Example filenames: `GCA_010883055.1_ASM1088305v1_genomic.nhr`
- Extract database paths (without extension): `data/ragtag_dbs/GCA_010883055.1_ASM1088305v1_genomic`
- Result: Dict of `{genome_name: database_path}` for 75 genomes

#### Step 3: Process Each Locus (lines 366-374)

**For each locus, calls `process_locus()` which does:**

**A. Filter Flanking Proteins (lines 199-239, called in process_locus):**

**Why filtering?** Some loci have 100+ flanking proteins. Searching all of them is slow and unnecessary.

1. Read first `MAX_FLANKING_FOR_SYNTENY` proteins (default: 15) from `*_flanking.faa`
2. Save to temporary filtered file: `<locus_id>_flanking_filtered.faa`
3. Use this smaller set for tBLASTn searches

**Example:**
```
Input: BK_chr3_a_flanking.faa (156 proteins)
Output: BK_chr3_a_flanking_filtered.faa (15 proteins)
```

**B. Run tBLASTn for Each Genome (lines 265-285):**

**For each of 75 target genomes:**

1. **Check for existing results:**
   - Look for `blast_xml/<genome_name>.xml`
   - If exists, skip BLAST and parse existing file
   - This allows resuming interrupted runs

2. **Run tBLASTn if needed:**
   ```bash
   tblastn \
     -query <locus>_flanking_filtered.faa \
     -db data/ragtag_dbs/<genome>.fasta \
     -outfmt 5 \             # XML output
     -evalue 1e-5 \
     -max_target_seqs 50 \
     -num_threads 16
   ```

3. **What tBLASTn does:**
   - Takes protein queries (flanking proteins)
   - Searches nucleotide database (genome)
   - Translates genome in all 6 reading frames (+1, +2, +3, -1, -2, -3)
   - Finds protein similarity matches (HSPs = High-Scoring Segment Pairs)
   - Returns: scaffold, coordinates, strand, alignment, percent identity

**C. Parse BLAST XML (lines 38-94, parse_blast_xml function):**

For each HSP in XML:
1. Extract query protein ID (e.g., `XP_033211300.1|LOC117167431`)
2. Extract hit details:
   - Scaffold ID (e.g., `CM021340.1_RagTag`)
   - Coordinates: start, end (nucleotide positions)
   - Strand: `+` or `-`
   - E-value, bitscore
   - Alignment length, identity count
   - Calculate percent identity: `(identity / align_length) × 100`
3. **Extract HSP sequences:**
   - `hit_protein_seq`: Translated protein sequence from genome (this is what we'll annotate later!)
   - `query_protein_seq`: Query protein sequence (for comparison)

**Output:** List of hit dictionaries:
```python
{
  'qseqid': 'XP_033211300.1|LOC117167431',
  'sseqid': 'CM021340.1_RagTag',
  'scaffold_desc': 'Diachasmimorpha longicaudata chromosome 3',
  'strand': '+',
  'coord_start': 12450000,
  'coord_end': 12450687,
  'evalue': 1.2e-45,
  'bitscore': 185.2,
  'pident': 65.3,
  'length': 229,
  'hit_protein_seq': 'MSLVRQNFHD...',  # Translated from genome
  'query_protein_seq': 'MSLVRQNFHD...'
}
```

**D. Save HSP Sequences (lines 161-197, save_hit_sequences function):**

**Why save these?** These are the actual protein sequences found in each genome - we'll annotate them with SwissProt in Phase 7.

1. For each HSP, remove alignment gaps from protein sequence
2. Create descriptive FASTA header:
   ```
   >GCA_123456.1_XP_033211300.1|LOC117167431_CM021340_12450000-12450687_hit1
   ```
3. Save all HSPs to: `<locus>_output/hit_sequences/<genome>_hit_proteins.fasta`

**E. Cluster Hits into Synteny Blocks (lines 96-157, cluster_into_blocks function):**

**Goal:** Group nearby BLAST hits into "synteny blocks" representing conserved genomic regions.

**Algorithm:**

1. **Group by scaffold and strand:**
   - Hits on `CM021340.1 (+)` separate from `CM021340.1 (-)`
   - Hits on `CM021340.1` separate from `CM021341.1`
   - Key: `(scaffold_id, strand)`

2. **For each (scaffold, strand) group:**

   a. **Sort hits by genomic position** (coord_start)

   b. **Cluster by proximity:**
      - Start with first hit as a "block"
      - For each subsequent hit:
        - Calculate gap from end of current block to start of this hit
        - If gap ≤ 500 kb (`SYNTENY_MAX_GAP_KB`): add to current block
        - If gap > 500 kb: save current block, start new block

   c. **Filter blocks:**
      - Count unique target proteins in block (not query proteins!)
      - Unique = distinct (scaffold, start, end) positions
      - Keep block only if ≥ 5 unique target proteins (`MIN_PROTEINS_FOR_BLOCK`)
      - This filters out spurious matches

   d. **Create block record:**
      ```python
      {
        'block_id': 'block_00001',
        'scaffold': 'CM021340.1_RagTag',
        'strand': '+',
        'start': 12400000,      # Leftmost hit
        'end': 12650000,        # Rightmost hit
        'span_kb': 250.0,       # (end - start) / 1000
        'num_target_proteins': 12,      # Unique proteins found
        'num_query_matches': 11,        # Unique query proteins that matched
        'hits': [<list of all HSP dicts in this block>],
        'genome': 'GCA_010883055.1',
        'locus_id': 'BK_chr3_a'
      }
      ```

**Why this matters:** A synteny block with 11/15 query proteins matching means we found a region with conserved gene order - likely the orthologous locus!

**F. Save Per-Locus Outputs (lines 308-327):**

**For each locus:**

1. `<locus>/flanking_blast_all.tsv`: All BLAST hits (without full sequences)
   - Columns: qseqid, sseqid, scaffold_desc, strand, coord_start, coord_end, evalue, bitscore, pident, length, genome, block_id

2. `<locus>/hit_sequences/<genome>_hit_proteins.fasta`: HSP protein sequences for annotation

#### Step 4: Filter to Best Block Per Genome (lines 376-383)

**Problem:** Some genomes may have multiple synteny blocks for same locus (duplications, paralogs)

**Solution:**
- Group blocks by (locus_id, genome)
- For each group, keep block with highest `num_query_matches`
- This selects the "best" synteny block per genome

**Example:**
```
Before filtering:
  Locus BK_chr3_a, Genome GCA_123456.1:
    - block_00001: 11 query matches
    - block_00002: 7 query matches
    - block_00003: 9 query matches

After filtering:
  Locus BK_chr3_a, Genome GCA_123456.1:
    - block_00001: 11 query matches  ← KEPT
```

#### Step 5: Save Combined Results (lines 385-417)

**A. `<locus_id>_synteny_blocks.tsv`** (for single locus array jobs):
```tsv
locus_id	genome	block_id	scaffold	strand	start	end	span_kb	num_target_proteins	num_query_matches
BK_chr3_a	GCA_010883055.1	block_00001	CM021340.1_RagTag	+	12400000	12650000	250.0	12	11
BK_chr3_a	GCA_011057455.1	block_00012	RICB01_scaffold_5	-	8950000	9120000	170.0	8	9
```

This file is used by Phase 4 to know where to search for target genes.

### Key Parameters:
- `FLANKING_BLAST_EVALUE = 1e-5` (tBLASTn sensitivity)
- `MAX_FLANKING_FOR_SYNTENY = 15` (proteins to use for search)
- `SYNTENY_MAX_GAP_KB = 500` (max gap between proteins in same block)
- `MIN_PROTEINS_FOR_BLOCK = 5` (min unique proteins for valid block)
- `BLAST_THREADS = 16` (parallelization per BLAST run)

### Outputs:
```
outputs/<gene_family>/02_synteny_blocks/
├── <locus_id>/
│   ├── flanking_blast_all.tsv
│   ├── hit_sequences/
│   │   ├── GCA_010883055.1_hit_proteins.fasta
│   │   ├── GCA_011057455.1_hit_proteins.fasta
│   │   └── ... (75 genomes)
│   └── blast_xml/
│       ├── GCA_010883055.1_ASM1088305v1_genomic.xml
│       └── ... (75 genomes)
└── <locus_id>_synteny_blocks.tsv
```

### Parallelization:
- Can run per-locus with `--locus <locus_id>` flag
- SLURM array: `#SBATCH --array=0-8` for 9 loci
- Each locus processes all 75 genomes independently

### Notes:
- For ferritin (9 loci, 75 genomes): ~30 minutes total
- Most time spent in tBLASTn searches
- XML files are large but allow resuming interrupted runs

---

## Phase 4: Target Gene Detection & Clustering

**Script:** `scripts/04_blast_targets.py`

**Input:** Synteny blocks from Phase 3 + BK reference proteins

**Purpose:** Find the actual target gene (the ortholog of the original query) within each synteny block

### Detailed Workflow:

#### Background: Why Do We Need This Phase?

Phase 3 found synteny blocks using **flanking proteins**. These blocks tell us "there's conserved gene order here", but they don't tell us which specific gene in that region is the actual target gene we care about.

**Example:**
- Query: BK ferritin gene `LOC117169185`
- Phase 3 found: Synteny block on scaffold CM021340.1, positions 12,400,000-12,650,000 (250 kb region)
- Question: Where exactly is the ferritin gene in this 250 kb region?

**Answer:** Search with the actual target protein (ferritin) to find its precise location.

#### Step 1: Load Locus Definitions (lines 185-191)
Same as Phase 3 - reads `locus_definitions.tsv`

#### Step 2: Find Unique Target Proteins (lines 195-225)

**Problem:** Multiple loci may share the same target protein sequence (paralogs that are identical/near-identical at protein level).

**Solution:** Create combined query file with unique sequences to avoid redundant BLAST searches.

**Algorithm:**

1. **For each locus:**
   - Look for `01_extracted_proteins/<locus_id>/<locus_id>_targets.faa`
   - This file contains the query protein from BK reference
   - Read protein sequence

2. **Deduplicate by sequence:**
   - Use protein sequence as dictionary key
   - If same sequence seen before, add this locus to its loci list
   - Track which protein ID goes with which loci

3. **Create combined query file:**
   - Save all unique proteins to `combined_targets.faa`
   - One BLAST search per unique protein, assigned to multiple loci if needed

**Example output:**
```
Found 3 unique target proteins for 9 loci
  Protein 1 (XP_033211301.1): used by 5 loci [BK_chr2_a, LB_scf7864_a, ...]
  Protein 2 (XP_033211302.1): used by 3 loci [BK_chr3_a, TR_CM036338_a, ...]
  Protein 3 (XP_033211303.1): used by 1 locus [DR_CM021344_b]
```

#### Step 3: Find Genome Databases (lines 228-232)
Same as Phase 3 - scans for 75 genome BLAST databases

#### Step 4: Run tBLASTn Per Genome (lines 235-264)

**Key difference from Phase 3:** Now searching with TARGET protein, not flanking proteins.

**For each of 75 genomes:**

1. **Check for existing results** (allows resuming)
2. **Run tBLASTn:**
   ```bash
   tblastn \
     -query combined_targets.faa \    # Target proteins (e.g., ferritin)
     -db <genome>.fasta \
     -outfmt 5 \                       # XML
     -evalue 1e-5 \
     -max_target_seqs 10000 \          # Higher than Phase 3!
     -num_threads 16
   ```

**Why max_target_seqs = 10000?**
- Target gene may have many paralogs across genome
- Need to find ALL copies to cluster them properly
- Phase 3 only needed top 50 hits per flanking protein

3. **Parse XML:** Same function as Phase 3 (`parse_blast_xml`)

**Result per genome:** List of all HSPs for target protein

#### Step 5: Cluster Hits by Reading Frame (lines 266-309, cluster_into_loci function)

**Critical fix (Nov 3):** Hits in different reading frames are now clustered separately!

**Why this matters:**
- tBLASTn translates genome in 6 frames: +1, +2, +3, -1, -2, -3
- Overlapping genes in different frames are DIFFERENT genes
- Example: Ferritin in +1 frame vs. unrelated gene in +2 frame at same locus

**Algorithm:**

1. **Extract reading frame from coordinates:**
   ```python
   frame = (start - 1) % 3 + 1  # Results in 1, 2, or 3
   if strand == '-':
       frame = -frame  # Results in -1, -2, or -3
   ```

2. **Group hits by (scaffold, strand, frame):**
   - Before fix: grouped by (scaffold, strand) only
   - After fix: separate clusters for different frames

3. **For each (scaffold, strand, frame) group:**

   a. **Sort hits by position**

   b. **Cluster by proximity:**
      - Gap threshold: 10 kb (`TARGET_CLUSTER_GAP_KB`)
      - Much smaller than Phase 3's 500 kb
      - Target gene hits should be tightly clustered (exons of same gene)

   c. **Filter clusters:**
      - Keep if ≥ 1 hit (`min_hits=1`)
      - Less stringent than Phase 3's 5-protein requirement

   d. **Create target cluster:**
      ```python
      {
        'locus_id': 'BK_chr3_a',
        'genome': 'GCA_010883055.1',
        'gene_family': 'ferritin_MC102',
        'scaffold': 'CM021340.1_RagTag',
        'strand': '+',
        'frame': 1,                    # NEW: reading frame
        'start': 12891127,             # Leftmost HSP
        'end': 12891312,               # Rightmost HSP
        'num_hsps': 3,                 # Number of HSPs (likely exons)
        'max_bitscore': 185.2,
        'avg_pident': 87.3,
        'target_id': 'BK_chr3_a_CM021340.1_RagTag_002',  # Unique cluster ID
        'hits': [<HSP list>]
      }
      ```

**Example - before and after frame-aware clustering:**

**Before (WRONG):**
```
Scaffold CM021340.1, strand +, position ~12,890,000:
  - HSPs in frame +1 (ferritin exons)
  - HSPs in frame +2 (different gene overlapping ferritin)
  → Clustered together as 1 target ❌
```

**After (CORRECT):**
```
Scaffold CM021340.1, strand +, frame +1, position ~12,890,000:
  - HSPs in frame +1 (ferritin exons)
  → Target cluster 1 ✓

Scaffold CM021340.1, strand +, frame +2, position ~12,890,000:
  - HSPs in frame +2 (different gene)
  → Target cluster 2 ✓
```

#### Step 6: Assign Clusters to Synteny Blocks (lines 311-337)

**Goal:** Determine which clusters fall within (or near) synteny blocks from Phase 3.

**For each locus:**

1. **Load synteny blocks** for this locus:
   - Read `02_synteny_blocks/<locus>_synteny_blocks.tsv`
   - Get list of blocks with scaffold, start, end for each genome

2. **Get target clusters** for this locus (from Step 5)

3. **For each target cluster:**

   a. **Find nearest synteny block:**
      - Check if cluster scaffold matches any block scaffold in this genome
      - Calculate distance from cluster midpoint to block
      - Keep block with smallest distance

   b. **Calculate proximity:**
      ```python
      cluster_mid = (cluster['start'] + cluster['end']) / 2
      if cluster_mid >= block['start'] and cluster_mid <= block['end']:
          distance_kb = 0  # Inside block
      else:
          # Distance from cluster to nearest block edge
          distance_kb = min(
              abs(cluster_mid - block['start']),
              abs(cluster_mid - block['end'])
          ) / 1000
      ```

   c. **Assign block_id and distance:**
      ```python
      cluster['block_id'] = block['block_id']
      cluster['distance_to_block_kb'] = distance_kb
      ```

4. **Create unique target_id:**
   - Format: `<locus>_<scaffold>_<number>`
   - Example: `BK_chr3_a_CM021340.1_RagTag_001`
   - Number increments for multiple clusters on same scaffold

**Result:** Each target cluster now knows which synteny block it belongs to (or is nearest to).

#### Step 7: Save Outputs (lines 339-365)

**A. Per-genome outputs:** `04_target_genes/<genome>/`
- `target_hits.tsv`: All raw tBLASTn HSPs
- `target_clusters.tsv`: Clustered targets with block assignments

**B. Combined output:** `04_target_genes/all_targets.tsv`
```tsv
locus_id	genome	target_id	scaffold	strand	frame	start	end	num_hsps	max_bitscore	block_id	distance_to_block_kb
BK_chr3_a	GCA_010883055.1	BK_chr3_a_CM021340.1_RagTag_001	CM021340.1_RagTag	+	1	12891127	12891312	3	185.2	block_00001	0.0
BK_chr3_a	GCA_010883055.1	BK_chr3_a_CM021340.1_RagTag_002	CM021340.1_RagTag	+	2	12891100	12891290	2	92.4	block_00001	0.0
```

Note the two separate clusters in frames +1 and +2!

#### Step 8: Run Phase 5 Classification (embedded in Phase 4, lines 343-345)

Phase 4 automatically calls `scripts/05_classify_targets.py` as a subprocess.

**What Phase 5 does:**

1. **Load target clusters** from `all_targets.tsv`
2. **Load synteny blocks** from `02_synteny_blocks/*.tsv`
3. **For each target:**
   - Check `distance_to_block_kb`
   - If ≤ 500 kb: classify as **"syntenic"** (near conserved region)
   - If > 500 kb: classify as **"unplaceable"** (isolated hit, possible horizontal transfer or distant paralog)
4. **Save classified_targets.tsv** with new column: `classification`
5. **Save syntenic_targets.tsv**: Filtered to syntenic targets only

**Output:** `05_classified/syntenic_targets.tsv` used by Phase 6

### Key Parameters:
- `TARGET_BLAST_EVALUE = 1e-5` (same as Phase 3)
- `TARGET_CLUSTER_GAP_KB = 10` (much smaller than Phase 3's 500 kb)
- `TARGET_BLAST_MAX_TARGETS = 10000` (finds all paralogs)
- `max_distance_kb = 500` (for syntenic classification in Phase 5)

### Outputs:
```
outputs/<gene_family>/04_target_genes/
├── combined_targets.faa
├── blast_xml/
│   ├── GCA_010883055.1.xml
│   └── ...
├── GCA_010883055.1/
│   ├── target_hits.tsv
│   └── target_clusters.tsv
├── all_targets.tsv
└── (Phase 5 creates 05_classified/ automatically)
```

### Notes:
- **Optimized:** Uses combined query file → 1 BLAST per genome instead of 1 per locus per genome
- **Frame-aware:** Critical bug fix prevents merging genes in different reading frames
- For ferritin: ~1 hour for 75 genomes (vs. many hours with old per-locus approach)

---

## Phase 5: Classify Targets (Integrated into Phase 4)

**Script:** `scripts/05_classify_targets.py`

**Input:** `04_target_genes/all_targets.tsv` + `02_synteny_blocks/*.tsv`

**Purpose:** Filter targets to keep only those near synteny blocks (removes spurious distant hits)

### Classification Logic:

```python
def check_proximity(target, block, max_distance_kb=500):
    # Check if on same scaffold
    if target['scaffold'] != block['scaffold']:
        return False, float('inf')

    # Calculate distance from target to block
    target_mid = (target['start'] + target['end']) / 2

    if target_mid >= block['start'] and target_mid <= block['end']:
        return True, 0  # Inside block

    # Distance to nearest edge
    dist_to_start = abs(target_mid - block['start'])
    dist_to_end = abs(target_mid - block['end'])
    distance_kb = min(dist_to_start, dist_to_end) / 1000

    is_near = (distance_kb <= max_distance_kb)
    return is_near, distance_kb
```

**Classifications:**
- **syntenic:** Target within 500 kb of a synteny block → likely true ortholog
- **unplaceable:** Target >500 kb from any block → possible horizontal transfer, distant paralog, or assembly error

**Output:** Only **syntenic** targets proceed to Phase 6.

### Why This Matters:

Without this filter, we'd try to extract genes from:
- Distant paralogs (duplicated and diverged)
- Horizontal gene transfers
- Assembly errors
- Random BLAST hits

Syntenic targets are the most reliable orthologs.

---

## Phase 6: Extract Gene Structures with Exonerate

**Script:** `scripts/06_extract_sequences.py`
**Modules:** `scripts/extract_with_exonerate.py`, `scripts/exonerate_extract.py`

**Input:** `05_classified/syntenic_targets.tsv` + RagTag genome FASTA files

**Purpose:** Extract complete gene structures (all exons, introns, CDS, protein) using Exonerate protein2genome alignment

### Background: Why Exonerate?

**tBLASTn (Phases 3-4) limitations:**
- Finds HSPs (high-scoring segment pairs) = local alignments
- Doesn't predict gene structure (exon boundaries, splice sites)
- Doesn't handle large introns well
- For ferritin: tBLASTn found only middle exon (76aa), missed exons 8-13kb away

**Exonerate solution:**
- Protein-to-genome alignment tool designed for gene structure prediction
- Respects splice donor/acceptor sites
- Can bridge large introns
- Outputs GFF with gene, exon, CDS features
- Returns complete gene models

### Detailed Workflow:

#### Step 1: Load Syntenic Targets (lines 19-29, in main())
- Read `05_classified/syntenic_targets.tsv`
- Filter to keep only rows with classification == "syntenic"
- Group by (genome, scaffold, target_id)

#### Step 2: For Each Target Cluster (lines 31-onwards, main() calls extract_block_genes())

**Calls:** `extract_with_exonerate.extract_block_genes(block, query_protein_file, genome_fasta, output_dir)`

**Inputs per target:**
- **block:** Dict with scaffold, start, end, strand, target_id
- **query_protein_file:** Path to BK reference protein FASTA (e.g., ferritin XP_033211301.1)
- **genome_fasta:** Path to genome assembly FASTA
- **output_dir:** Where to save extracted sequences

#### Step 3: Adaptive Windowing Strategy (lines 48-113, extract_with_exonerate.py)

**Problem:** We don't know how large the introns are. Ferritin spans 22 kb total.

**Solution:** Try progressively larger search windows until we find complete gene.

**Window sizes (flanking on each side of HSP):**
```python
flank_sizes = [0, 1000, 5000, 10000, 15000, 20000, 50000, 100000, 200000]  # bp
```

**For each window size:**

**A. Extract genomic region (lines 56-62):**
```python
region_file = exonerate_extract.extract_genomic_region(
    genome_fasta=genome_fasta,
    scaffold=block['scaffold'],
    start=block['start'],    # HSP start from Phase 4
    end=block['end'],        # HSP end from Phase 4
    flank=flank            # Add this much on each side
)
```

**What this does:**
1. Parse genome FASTA to find scaffold
2. Extract subsequence: `[start - flank : end + flank]`
3. Save to temp file: `temp_<scaffold>_<start>_<end>.fasta`
4. Return path to temp file

**Example:**
- HSP coordinates: 12,891,127 - 12,891,312 (middle exon, 186 bp)
- Window 0 kb: Extract 12,891,127 - 12,891,312 (186 bp total)
- Window 1 kb: Extract 12,890,127 - 12,892,312 (2.2 kb total)
- Window 20 kb: Extract 12,871,127 - 12,911,312 (40 kb total)

**B. Run Exonerate (lines 68-74):**
```bash
exonerate \
  --model protein2genome \
  --query <BK_ferritin_protein>.faa \
  --target <temp_genomic_region>.fasta \
  --showtargetgff yes \
  --showvulgar yes \
  --showquerygff no \
  --showcigar no \
  --ryo "# Alignment %qi vs %ti\n..."
```

**What Exonerate does:**
1. Aligns query protein to genomic DNA
2. Predicts exon boundaries based on:
   - Splice donor sites (GT)
   - Splice acceptor sites (AG)
   - Reading frame maintenance
   - Protein similarity
3. Outputs GFF with gene structure

**Example Exonerate output:**
```gff
CM021340.1_RagTag	exonerate	gene	8373	21744	127	+	.	gene_id 1 ; ...
CM021340.1_RagTag	exonerate	exon	8373	8470	.	+	.	identity 39.29 ; ...
CM021340.1_RagTag	exonerate	exon	16742	16966	.	+	.	identity 39.29 ; ...
CM021340.1_RagTag	exonerate	exon	21558	21744	.	+	.	identity 39.29 ; ...
```

This shows 3 exons spanning 13.4 kb!

**C. Parse Exonerate GFF (lines 76-85):**

Uses `exonerate_extract.parse_exonerate_gff()` to extract:
- gene features (with gene_id attribute)
- exon features
- CDS features
- Each feature has: type, start, end, strand, attributes

**D. Calculate coverage (lines 87-90):**
```python
cds_features = [f for f in features if f['type'].lower() in ['cds', 'coding_exon', 'exon']]
total_cds_length = sum(f['end'] - f['start'] + 1 for f in cds_features)
coverage = total_cds_length / (query_length * 3)  # Divide by 3: protein aa → nucleotides
```

**Example:**
- Query protein: 174 aa → 522 nt expected
- CDS found: 99 + 225 + 186 = 510 nt
- Coverage: 510 / 522 = 97.7% ✓

**E. Check completeness (lines 92-112):**
```python
is_complete = (coverage >= 0.90)  # 90% threshold

if coverage > best_coverage:
    # This is better than previous attempts
    best_features = features
    best_region_file = region_file
    best_flank = flank
    best_coverage = coverage

if is_complete:
    print(f"Complete gene at {flank}bp flanking (coverage: {coverage:.1%})")
    break  # Stop expanding window
```

**Adaptive strategy in action:**
```
Window 0kb (186 bp):   Coverage 35.6%  → keep searching
Window 1kb (2.2 kb):   Coverage 35.6%  → keep searching
Window 5kb (10 kb):    Coverage 63.2%  → keep searching
Window 20kb (40 kb):   Coverage 97.7%  → COMPLETE! Stop.
```

#### Step 4: Group Exons into Genes (lines 133-185, extract_with_exonerate.py)

**Bug fix (Nov 3):** Previously treated each exon as separate gene!

**Correct approach:**

1. **Find gene features:**
   ```python
   gene_features = [f for f in features if f['type'].lower() == 'gene']
   ```

2. **For each gene feature:**
   - Extract gene_id from attributes: `gene_id 1 ; ...` → gene_id = "1"
   - Get gene span: start, end, strand

3. **Collect CDS/exon features within gene span:**
   ```python
   for feat in features:
       if feat['type'].lower() in ['cds', 'coding_exon']:
           if feat['start'] >= gene_start and feat['end'] <= gene_end:
               genes_by_id[gene_id]['cds_features'].append(feat)
   ```

4. **If no explicit gene features:**
   - Treat all CDS/exons as one gene
   - Some Exonerate versions don't output gene wrapper

**Result:** Dictionary mapping gene_id → {gene_feature, cds_features, exon_features, start, end, strand}

#### Step 5: Extract Sequences for Each Gene (lines 187-273)

**For each gene:**

**A. Extract CDS sequence (lines 196-201):**

Calls `exonerate_extract.extract_cds_sequence(region_file, cds_features, gene_id)`:

1. **Load temp genomic region file** (NOT full genome - coordinates are relative!)
2. **For each CDS feature:**
   - Extract subsequence using GFF coordinates (1-based inclusive)
   - Convert to 0-based Python indexing: `seq[start-1:end]`
3. **Concatenate CDS sequences** in genomic order
4. **If minus strand:** Reverse complement the concatenated sequence

**Example:**
```python
CDS 1: positions 8373-8470   → extract 98 bp
CDS 2: positions 16742-16966  → extract 225 bp
CDS 3: positions 21558-21744  → extract 187 bp
Concatenate: 510 bp total
```

**B. Extract full genomic sequence (lines 207-221):**

1. Extract gene span from temp region file: `seq[gene_start-1:gene_end]`
2. This includes all introns
3. If minus strand: reverse complement

**C. Translate to protein (lines 240-249):**
```python
from Bio.Seq import Seq
protein_seq = str(Seq(cds_seq).translate(to_stop=True))
```

BioPython's translate():
- Translates nucleotides → amino acids
- Stops at first stop codon (`to_stop=True`)
- Handles frame automatically (expects in-frame CDS)

**D. Classify functional status (line 204):**

Uses `exonerate_extract.classify_gene_status(cds_features, cds_seq)`:

Checks for:
- Has start codon (ATG)?
- Has stop codon (TAA, TAG, TGA)?
- In-frame stop codons (premature stops)?
- Length compared to query

**Classifications:**
- **functional:** Has start, has stop, no internal stops, reasonable length
- **partial:** Missing start or stop (truncated at assembly edge)
- **pseudogene:** Internal stop codons (frameshifts, mutations)

**E. Save sequences (lines 224-268):**

**Three files per gene:**

1. **CDS FASTA** (`<block>_gene<id>_cds.fasta`):
   ```fasta
   >gene1 CM021340.1_RagTag:12891127-12891312 strand:+ status:functional exons:3
   ATGAGCCTGGTGAGGCAGAACTTCCACGACGAGTGCGAGACGGC...
   ```

2. **Genomic FASTA** (`<block>_gene<id>_genomic.fasta`):
   ```fasta
   >gene1 CM021340.1_RagTag:12891127-12891312 type:genomic_with_introns length:13371bp
   ATGAGCCTGGTGAGGCAGAACTTCCACGACGAGTGCGAGACGGC...GTAAG...exon1...GTAAG
   T...intron1...AGGTAAG...exon2...GTAAGT...intron2...AGAT...exon3...TGATAA
   ```
   (Shows full gene with introns)

3. **Protein FASTA** (`<block>_gene<id>_protein.fasta`):
   ```fasta
   >gene1 CM021340.1_RagTag:12891127-12891312 strand:+ length:170aa
   MSLVRQNFHDECETALNIQINLELYASYVYLSMAYHFDRSDVALQGLSKYFKEASEEERKHAMK...
   ```

#### Step 6: Clean Up (line 271-272)
Delete temp genomic region file (only keep extracted sequences)

### Key Parameters:
- **Window sizes:** 0, 1kb, 5kb, 10kb, 15kb, 20kb, 50kb, 100kb, 200kb flanking
- **Completeness threshold:** 90% coverage
- **Exonerate model:** protein2genome (splice-aware)

### Outputs:
```
outputs/<gene_family>/04_target_genes/extracted_sequences/
└── <genome>/
    └── <target_id>/
        ├── <target>_exonerate_flank<size>.txt    # Exonerate alignment details
        ├── <target>_gene1_cds.fasta              # CDS sequence
        ├── <target>_gene1_genomic.fasta          # Gene with introns
        └── <target>_gene1_protein.fasta          # Translated protein
```

### Performance:
- For ferritin (736 targets, 75 genomes): ~42 hours
- Most genes complete at 0-20kb window (fast)
- Large intron genes need 50-200kb window (slower)
- Adaptive strategy optimizes speed vs completeness

### Notes:
- **Critical for multi-exon genes:** Without adaptive windowing, would miss distant exons
- **Handles partial genes:** Reports coverage % for incomplete extractions
- **Three output types:** Enables downstream CDS analysis, genomic analysis, and protein annotation

---

## Phase 7: SwissProt Annotation

**Script:** `scripts/07_swissprot_annotation.py`

**Input:** HSP protein sequences from Phase 3 (in `02_synteny_blocks/*/hit_sequences/`)

**Purpose:** Annotate all flanking protein hits with SwissProt for functional context

### Detailed Workflow:

#### Background

In Phase 3, we saved the **translated protein sequences** from tBLASTn HSPs. These are the actual proteins found in each genome that matched our flanking queries. We want to know what these proteins are (functionally) to understand the genomic context.

**Example:**
- Query flanking protein: BK `XP_033211300.1` (unknown function)
- HSP from genome GCA_123: Translated sequence `MKVLPGRS...`
- Question: What is this protein?
- Answer: Search SwissProt → finds match to "Cytochrome P450 3A4"

This gives us functional labels for synteny validation.

#### Step 1: Collect All HSP Sequences (lines 30-60, in main())

1. **Scan for hit_sequences directories:**
   ```python
   hit_dirs = list(config.STEP02_SYNTENY.glob("*/hit_sequences"))
   ```

2. **Read all FASTA files:**
   ```python
   for fasta_file in hit_dir.glob("*.fasta"):
       for record in SeqIO.parse(fasta_file, 'fasta'):
           all_seqs[record.id] = record
   ```

3. **Deduplicate:**
   - Use sequence ID as key
   - Avoids re-annotating same protein multiple times
   - IDs are unique: `<genome>_<query>_<scaffold>_<coords>_hit<N>`

4. **Save combined file:**
   ```python
   combined_faa = config.STEP07_SWISSPROT / "all_hits.faa"
   SeqIO.write(all_seqs.values(), combined_faa, 'fasta')
   ```

**Output:** One FASTA with all unique HSP sequences (for ferritin: 293,478 sequences)

#### Step 2: Run DIAMOND Search (lines 62-85)

**Using DIAMOND instead of BLASTP:**
- DIAMOND is 100-10,000x faster than BLASTP
- Designed for large-scale protein searches
- Comparable sensitivity for close homologs

**Command:**
```bash
diamond blastp \
  --query all_hits.faa \                 # All HSP sequences
  --db data/databases/swissprot.dmnd \    # SwissProt database
  --out swissprot_results.tsv \
  --outfmt 6 \                           # Tabular format
  --evalue 1e-5 \
  --max-target-seqs 1 \                  # Best hit only
  --threads 16
```

**Output format (tab-separated):**
```
qseqid	sseqid	pident	length	mismatch	gapopen	qstart	qend	sstart	send	evalue	bitscore
```

#### Step 3: Parse and Filter Results (lines 87-120)

1. **Read DIAMOND output:**
   ```python
   results = pd.read_csv(diamond_out, sep='\t', header=None, names=[...])
   ```

2. **For each query, keep best hit:**
   - Group by qseqid
   - Sort by bitscore (descending)
   - Keep first row per group
   ```python
   best_hits = results.sort_values('bitscore', ascending=False).groupby('qseqid').first()
   ```

3. **Extract SwissProt annotations:**
   - SwissProt IDs have format: `sp|P12345|PROTEIN_NAME_SPECIES`
   - Parse to extract protein name and species
   ```python
   # Example: sp|P00789|CYPE_HUMAN → CYPE (Cytochrome P450)
   sseqid = 'sp|P00789|CYPE_HUMAN'
   protein_name = sseqid.split('|')[2].split('_')[0]  # "CYPE"
   species = sseqid.split('|')[2].split('_')[1]       # "HUMAN"
   ```

4. **Create annotation table:**
   ```python
   annotations = {
       'query_id': query,
       'swissprot_id': sseqid,
       'protein_name': protein_name,
       'species': species,
       'pident': pident,
       'evalue': evalue,
       'bitscore': bitscore
   }
   ```

#### Step 4: Save Annotations (lines 122-140)

**Two output files:**

1. **swissprot_results.tsv:** Raw DIAMOND output (all fields)

2. **annotations.tsv:** Clean annotation table
   ```tsv
   query_id	swissprot_id	protein_name	species	pident	evalue	bitscore
   GCA_010883055.1_XP_033211300_..._hit1	sp|P00789|CYPE_HUMAN	CYPE	HUMAN	67.3	1.2e-45	182.5
   ```

### Key Parameters:
- `SWISSPROT_BLAST_EVALUE = 1e-5`
- `USE_DIAMOND = True` (much faster than BLASTP)
- `--max-target-seqs 1` (best hit only per query)
- `BLAST_THREADS = 16`

### Outputs:
```
outputs/<gene_family>/07_swissprot/
├── all_hits.faa               # Combined HSP sequences
├── swissprot_results.tsv      # Raw DIAMOND output
└── annotations.tsv            # Clean annotations
```

### Performance:
- For ferritin: 293,478 proteins annotated in <1 hour
- 100% match rate (every protein got a SwissProt hit at E≤1e-5)
- DIAMOND is critical for speed at this scale

### Notes:
- Annotations used in Phase 8 matrices for functional context
- Helps validate synteny (do flanking genes make sense functionally?)
- Can reveal horizontal gene transfers (bacterial genes in insects)

---

## Phase 8a: Generate Locus-Specific Matrices

**Script:** `scripts/08a_generate_locus_matrices.py`

**Input:** All previous phase outputs

**Purpose:** Create detailed presence/absence matrix for each locus showing all flanking proteins and target gene status

### Matrix Structure:

**Rows:** 75 genomes (sorted phylogenetically)
**Columns:**
- Flanking proteins (Upstream1, Upstream2, ..., Downstream1, Downstream2, ...)
- Target gene
- Metadata (genome, species, assembly)

**Cell values:**
- **For flanking proteins:**
  - `1` if protein found in genome's synteny block
  - `0` if not found
  - Annotation from SwissProt in adjacent column
- **For target gene:**
  - Length in amino acids if present
  - Functional status (functional/partial/pseudogene)
  - `ABSENT` if no synteny block in this genome

### Example Matrix:

```tsv
Genome	Species	U1_ferritin	U1_annot	U2_enolase	U2_annot	Target_ferritin	Target_status
GCA_010883055.1	Bkinseyi	1	FERRITIN	1	ENOLASE	174	functional
GCA_011057455.1	Lboulardi	1	FERRITIN	1	ENOLASE	171	functional
GCA_015476485.1	Pnear	0	-	1	ENOLASE	ABSENT	-
```

### Key Features:
- Phylogenetic ordering reveals evolutionary patterns
- SwissProt annotations show conserved functional neighborhoods
- Target status distinguishes functional vs. degraded genes

---

## Phase 8b: Generate Summary Matrices

**Script:** `scripts/08b_generate_summary_matrices.py`

**Input:** All locus matrices from Phase 8a

**Purpose:** Aggregate locus-level data into gene family summaries

### Summary Structure:

**Rows:** 75 genomes
**Columns:**
- Functional_count: # of functional genes
- Partial_count: # of partial genes
- Pseudogene_count: # of pseudogenes
- Absent_count: # of loci with no synteny
- Total_loci: # of loci searched
- Functional_rate: Functional / Total

### Example Summary:

```tsv
Genome	Species	Functional	Partial	Pseudogene	Absent	Total	Functional_rate
GCA_010883055.1	Bkinseyi	9	0	0	0	9	100%
GCA_011057455.1	Lboulardi	7	1	0	1	9	77.8%
GCA_015476485.1	Pnear	3	2	1	3	9	33.3%
```

### Applications:
- Identify genome-wide gene gain/loss patterns
- Detect lineage-specific duplications
- Quality check: High pseudogene rate may indicate assembly issues

---

## Key Configuration (config.py)

### File Paths:
```python
LOCI_DEFINITIONS_FILE = os.environ.get('SYNTENY_INPUT_FILE', 'inputs/locus_definitions.tsv')
OUTPUT_BASE = os.environ.get('SYNTENY_OUTPUTS_DIR', 'outputs/default')

BK_GFF_FILE = "data/reference/genomic.gff"
BK_PROTEINS_FILE = "data/reference/protein.faa"

RAGTAG_DB_DIR = "data/ragtag_dbs/"           # tBLASTn databases
RAGTAG_FASTA_DIR = "data/ragtag_output/"     # Genome FASTA files
SWISSPROT_DB = "data/databases/swissprot.dmnd"
```

### BLAST Parameters:
```python
FLANKING_BLAST_EVALUE = "1e-5"
TARGET_BLAST_EVALUE = "1e-5"
SWISSPROT_BLAST_EVALUE = "1e-5"
BLAST_THREADS = 16

FLANKING_BLAST_MAX_TARGETS = 50    # Phase 3
TARGET_BLAST_MAX_TARGETS = 10000   # Phase 4 (needs to find all paralogs)
```

### Synteny Parameters:
```python
MAX_FLANKING_FOR_SYNTENY = 15      # Use top 15 flanking proteins (speed optimization)
SYNTENY_MAX_GAP_KB = 500           # Max gap between proteins in synteny block
MIN_PROTEINS_FOR_BLOCK = 5         # Min proteins required for valid block

TARGET_CLUSTER_GAP_KB = 10         # Max gap for clustering target HSPs
```

### Output Directories:
```python
STEP02_SYNTENY = OUTPUT_BASE / "02_synteny_blocks"
STEP04_TARGETS = OUTPUT_BASE / "04_target_genes"
STEP05_CLASSIFIED = OUTPUT_BASE / "05_classified"
STEP07_SWISSPROT = OUTPUT_BASE / "07_swissprot"
STEP07_MATRICES = OUTPUT_BASE / "07_matrices"
```

---

## Complete Output Structure

```
outputs/<gene_family>_<timestamp>/
├── 02_synteny_blocks/
│   ├── <locus_id>/
│   │   ├── <locus>_flanking_filtered.faa         # Top 15 flanking proteins used
│   │   ├── flanking_blast_all.tsv                # All tBLASTn hits
│   │   ├── hit_sequences/
│   │   │   ├── GCA_010883055.1_hit_proteins.fasta  # HSP sequences for SwissProt
│   │   │   ├── GCA_011057455.1_hit_proteins.fasta
│   │   │   └── ... (75 genomes)
│   │   └── blast_xml/
│   │       ├── GCA_010883055.1.xml               # Raw BLAST output
│   │       └── ...
│   └── <locus_id>_synteny_blocks.tsv             # Filtered blocks summary
│
├── 04_target_genes/
│   ├── combined_targets.faa                      # Unique target queries
│   ├── blast_xml/
│   │   └── GCA_*.xml                            # Target BLAST results
│   ├── GCA_*/
│   │   ├── target_hits.tsv                       # Raw target HSPs
│   │   └── target_clusters.tsv                   # Clustered by frame
│   ├── all_targets.tsv                           # All clusters combined
│   └── extracted_sequences/
│       └── GCA_*/
│           └── <target_id>/
│               ├── <target>_exonerate_flank<kb>.txt  # Alignment details
│               ├── <target>_gene1_cds.fasta           # CDS
│               ├── <target>_gene1_genomic.fasta       # With introns
│               └── <target>_gene1_protein.fasta       # Translated
│
├── 05_classified/
│   ├── classified_targets.tsv                    # All with classification
│   └── syntenic_targets.tsv                      # Filtered to syntenic only
│
├── 07_swissprot/
│   ├── all_hits.faa                              # All HSP sequences
│   ├── swissprot_results.tsv                     # DIAMOND results
│   └── annotations.tsv                           # Best hit per protein
│
└── 07_matrices/
    ├── locus_specific/
    │   └── <locus>_matrix.tsv                    # Detailed matrix per locus
    └── gene_type_summaries/
        └── <gene_family>_summary.tsv             # Aggregate summary
```

---

## Execution Modes

### Single Gene Family (Manual):
```bash
# Set environment
export SYNTENY_INPUT_FILE="outputs/ferritin_phase12/locus_definitions.tsv"
export SYNTENY_OUTPUTS_DIR="outputs/ferritin_results"

# Run phases sequentially
sbatch run_phase3.slurm      # Synteny detection
sbatch run_phase4.slurm      # Target finding (includes Phase 5)
sbatch run_phase6.slurm      # Sequence extraction
sbatch run_phase7.slurm      # SwissProt annotation
python scripts/08a_generate_locus_matrices.py
python scripts/08b_generate_summary_matrices.py
```

### Per-Locus Parallelization (Array Jobs):
```bash
# Phase 3 array (9 loci)
#SBATCH --array=0-8
LOCI=(BK_chr2_a BK_chr3_a LB_scf7864_a ...)
python scripts/02_synteny_detection.py --locus ${LOCI[$SLURM_ARRAY_TASK_ID]}
```

### Batch Processing (32 Gene Families):
**TO BE IMPLEMENTED** - would use array jobs per gene family

---

## Performance Benchmarks (Ferritin Test Case)

**Specs:** 9 loci, 75 genomes, 736 target clusters

| Phase | Time | Parallelization | Bottleneck |
|-------|------|----------------|------------|
| 1 (Optional) | ~1 hour | Per-LOC | DIAMOND searches |
| 3 | ~30 min | Per-locus array | tBLASTn searches |
| 4 | ~1 hour | None (optimized) | tBLASTn searches |
| 5 | ~5 min | None | I/O bound |
| 6 | ~42 hours | Could parallelize | Exonerate + adaptive windows |
| 7 | ~1 hour | None | DIAMOND search |
| 8 | <10 min | None | Matrix generation |

**Total:** ~45 hours for single gene family
**With array jobs:** Could reduce to ~10 hours (Phase 3 + 6 parallelized)

---

## TODO / Future Improvements

1. **Batch mode for 32 gene families** - Array job infrastructure
2. **Parallelize Phase 6** - Per-genome or per-target arrays
3. **Validate window size limits** - Check if 200kb is sufficient
4. **Add QC metrics** - Coverage plots, synteny conservation scores
5. **Optimize Exonerate** - Tune parameters for speed vs accuracy
6. **Add visualization** - Synteny plots, phylogenetic trees with presence/absence
