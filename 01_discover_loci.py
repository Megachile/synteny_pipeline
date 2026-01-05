#!/usr/bin/env python3
"""
Phase 1: Locus discovery + locus envelope metadata + dedup flanking preference.

Part of the Helixer proteome pipeline (canonical synteny-scanning pipeline).

TWO MODES OF OPERATION:

1. Legacy mode (--loc-ids): Takes LOC IDs, parses BK GFF to find coordinates
   python 01_phase1.py \
       --loc-ids LOC117167432,LOC117167433 \
       --gene-family ferritin_MC102 \
       --output-dir outputs/ferritin_MC102/phase1

2. Coordinates mode (--coordinates-file): Uses pre-computed genomic coordinates
   python 01_phase1.py \
       --coordinates-file top50_effector_genomic_coordinates.tsv \
       --gene-family ferritin_MC102 \
       --output-dir outputs/ferritin_MC102/phase1

   Expected columns in coordinates file:
   - subcluster_id: Used as locus_id
   - gene_id: LOC ID (e.g., LOC117175425)
   - chromosome: RefSeq chromosome ID (e.g., NC_046662.1)
   - start, end: Genomic coordinates
   - strand: + or -

This wraps the discovery implementation to produce flanking/targets, then
post-processes locus_definitions.tsv to:
  - Prefer deduplicated flanking files ("*_flanking_dedup.faa")
  - Add expected chromosome from locus_id (e.g., BK_chr2_x → chr2)
  - Add locus envelope coordinates from debug flanking files when available
  - Compute a per-locus scale factor (locus_scale) based on span vs cohort median
"""

from __future__ import annotations

from pathlib import Path
import sys
import argparse
import pandas as pd


def _import_legacy_scripts():
    """Import phase1_discovery_impl from the same directory as this script."""
    here = Path(__file__).resolve().parent
    # phase1_discovery_impl.py is now co-located in helixer_pipeline/
    sys.path.insert(0, str(here))
    return here


def run_legacy_discovery(loc_ids: str, gene_family: str, out_dir: Path) -> None:
    """Invoke existing phase1_discovery_impl in-process, writing to out_dir.

    This mirrors scripts/01_phase1_combined.py's approach.
    """
    _import_legacy_scripts()
    import phase1_discovery_impl as _disc  # type: ignore
    old_argv = sys.argv[:]
    try:
        sys.argv = [
            _disc.__file__,
            '--loc-ids', loc_ids,
            '--output-dir', str(out_dir),
            '--gene-family', gene_family,
        ]
        if hasattr(_disc, 'main'):
            _disc.main()
        else:
            raise RuntimeError('phase1_discovery_impl has no main()')
    finally:
        sys.argv = old_argv


def _expected_chr_from_locus_id(locus_id: str) -> str:
    # Examples: BK_chr2_a → chr2; LB_scf7864_a → scf7864
    try:
        parts = locus_id.split('_')
        if len(parts) >= 3 and parts[1].startswith('chr'):
            return parts[1]
        if len(parts) >= 3 and parts[1].startswith('scf'):
            return parts[1]
        return ''
    except Exception:
        return ''


def _load_gff_gene_positions(gff_path: Path) -> dict[str, tuple[str, int, int]]:
    """Load gene positions from GFF. Returns {gene_id: (chromosome, start, end)}."""
    positions = {}
    try:
        with open(gff_path) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 9 or parts[2] != 'gene':
                    continue
                chrom = parts[0]
                start = int(parts[3])
                end = int(parts[4])
                attrs = parts[8]
                # Extract gene ID from attributes (e.g., ID=gene-LOC117173775)
                for attr in attrs.split(';'):
                    if attr.startswith('ID=gene-'):
                        gene_id = attr.replace('ID=gene-', '')
                        positions[gene_id] = (chrom, start, end)
                        break
    except Exception as e:
        print(f"Warning: Failed to load GFF {gff_path}: {e}")
    return positions


def _calculate_flanking_spans_from_fasta(phase1_dir: Path, genome_gff_dir: Path) -> pd.DataFrame:
    """
    Calculate flanking gene spans directly from flanking FASTA files and GFF.

    For each locus, reads the flanking FASTA, extracts gene IDs from headers,
    looks up their coordinates in the GFF, and returns min/max as the span.
    """
    from Bio import SeqIO

    rows = []

    # Find BK GFF for coordinate lookup - try multiple locations
    bk_gff = None
    candidates = [
        genome_gff_dir / 'belonocnema_kinseyi' / 'GCF_010883055.1_B_treatae_v1_genomic.gff',
        Path('/carc/scratch/projects/emartins/2016456/adam/genomes/belonocnema_kinseyi/GCF_010883055.1_B_treatae_v1_genomic.gff'),
        Path('/carc/scratch/projects/emartins/2016456/adam/genomes/belonocnema_kinseyi/bkinseyi_genome.gff'),
    ]
    for candidate in candidates:
        if candidate.exists():
            bk_gff = candidate
            break

    # Also try glob search if still not found
    if bk_gff is None and genome_gff_dir.exists():
        for candidate in genome_gff_dir.glob('**/GCF_010883055*.gff'):
            bk_gff = candidate
            break
        if bk_gff is None:
            for candidate in genome_gff_dir.glob('**/bkinseyi*.gff'):
                bk_gff = candidate
                break

    if bk_gff is None or not bk_gff.exists():
        print(f"Warning: BK GFF not found (tried {len(candidates)} locations)")
        return pd.DataFrame(columns=['locus_id', 'flanking_span_start', 'flanking_span_end', 'flanking_span_kb'])

    print(f"  Loading BK gene positions from {bk_gff.name}...")
    gene_positions = _load_gff_gene_positions(bk_gff)
    print(f"  Loaded positions for {len(gene_positions)} genes")

    # Find all flanking FASTA files
    flanking_files = list(phase1_dir.glob('*_flanking*.faa'))
    if not flanking_files:
        print(f"  No flanking FASTA files found in {phase1_dir}")
        return pd.DataFrame(columns=['locus_id', 'flanking_span_start', 'flanking_span_end', 'flanking_span_kb'])

    for fasta_file in flanking_files:
        # Extract locus_id from filename (e.g., BK_chr6_a_flanking.faa -> BK_chr6_a)
        locus_id = fasta_file.stem.replace('_flanking_dedup', '').replace('_flanking', '')

        try:
            coords = []
            for record in SeqIO.parse(fasta_file, 'fasta'):
                # Header format: >XP_033221061.1|LOC117175465 U1 XP_033221061.1 description
                # Extract LOC ID
                header_parts = record.id.split('|')
                if len(header_parts) >= 2:
                    loc_id = header_parts[1].split()[0]  # Get LOC117175465
                    if loc_id in gene_positions:
                        chrom, start, end = gene_positions[loc_id]
                        coords.append((start, end))

            if coords:
                min_coord = min(c[0] for c in coords)
                max_coord = max(c[1] for c in coords)
                span_kb = (max_coord - min_coord) / 1000.0

                rows.append({
                    'locus_id': locus_id,
                    'flanking_span_start': min_coord,
                    'flanking_span_end': max_coord,
                    'flanking_span_kb': span_kb
                })
                print(f"    {locus_id}: {len(coords)} genes, span={span_kb:.1f} kb ({min_coord:,}-{max_coord:,})")
        except Exception as e:
            print(f"  Warning: Failed to process {fasta_file.name}: {e}")
            continue

    if not rows:
        return pd.DataFrame(columns=['locus_id', 'flanking_span_start', 'flanking_span_end', 'flanking_span_kb'])

    return pd.DataFrame(rows)


def _calculate_flanking_spans(phase1_dir: Path, genome_gff_dir: Path) -> pd.DataFrame:
    """Calculate flanking gene spans from debug_flanking files and GFF."""
    # Genome ID mapping for common abbreviations
    genome_mappings = {
        'BK': 'kinseyi',
        'LB': 'boulardi',
        'AF': 'quercusfoliatus',
        'CB': 'quercusbatatoides',
        'DC': 'cinerosa',
        'DL': 'quercuslanigerum',
        'NH': 'howertoni'
    }

    rows = []

    for dbg in phase1_dir.glob('debug_flanking_*.tsv'):
        try:
            # Extract genome from filename (e.g., debug_flanking_BK.tsv → BK)
            genome_id = dbg.stem.replace('debug_flanking_', '')

            # Find corresponding GFF file (try direct match first, then mapping)
            gff_candidates = list(genome_gff_dir.glob(f'*{genome_id}*.gff'))
            if not gff_candidates and genome_id in genome_mappings:
                # Try with species name
                species_hint = genome_mappings[genome_id]
                gff_candidates = list(genome_gff_dir.glob(f'*{species_hint}*.gff'))

            if not gff_candidates:
                print(f"Warning: No GFF found for genome {genome_id}")
                continue

            gff_path = gff_candidates[0]
            print(f"  Using GFF: {gff_path.name} for genome {genome_id}")
            gene_positions = _load_gff_gene_positions(gff_path)

            # Read debug_flanking file
            df = pd.read_csv(dbg, sep='\t')
            if not {'gene', 'chromosome', 'upstream_genes', 'downstream_genes'} <= set(df.columns):
                continue

            for _, row in df.iterrows():
                target_gene = row['gene']
                upstream_str = str(row['upstream_genes']) if pd.notna(row['upstream_genes']) else ''
                downstream_str = str(row['downstream_genes']) if pd.notna(row['downstream_genes']) else ''

                upstream_list = [g.strip() for g in upstream_str.split(',') if g.strip()]
                downstream_list = [g.strip() for g in downstream_str.split(',') if g.strip()]

                # Get positions
                positions_all = []

                # Add upstream genes
                for gene_id in upstream_list:
                    if gene_id in gene_positions:
                        positions_all.append(gene_positions[gene_id])

                # Add downstream genes
                for gene_id in downstream_list:
                    if gene_id in gene_positions:
                        positions_all.append(gene_positions[gene_id])

                if not positions_all:
                    continue

                # Calculate span from first to last flanking gene
                all_starts = [pos[1] for pos in positions_all]
                all_ends = [pos[2] for pos in positions_all]

                flanking_span_start = min(all_starts)
                flanking_span_end = max(all_ends)

                # Find locus_id (need to match target_gene to locus_id)
                # For now, use gene ID as proxy; will match properly in postprocess
                rows.append({
                    'target_gene': target_gene,
                    'flanking_span_start': flanking_span_start,
                    'flanking_span_end': flanking_span_end
                })

        except Exception as e:
            print(f"Warning: Failed to process {dbg}: {e}")
            continue

    if not rows:
        return pd.DataFrame(columns=['target_gene', 'flanking_span_start', 'flanking_span_end'])

    return pd.DataFrame(rows)


def _read_debug_boundaries(phase1_dir: Path) -> pd.DataFrame:
    rows = []
    for dbg in phase1_dir.glob('debug_flanking_*.tsv'):
        try:
            df = pd.read_csv(dbg, sep='\t')
            # Expect columns: boundary_start, boundary_end, locus_id
            if {'boundary_start', 'boundary_end', 'locus_id'} <= set(df.columns):
                rows.append(df[['locus_id', 'boundary_start', 'boundary_end']])
        except Exception:
            continue
    if not rows:
        return pd.DataFrame(columns=['locus_id', 'boundary_start', 'boundary_end'])
    merged = pd.concat(rows, ignore_index=True)
    # Some files may duplicate; keep first per locus_id
    merged = merged.drop_duplicates(subset=['locus_id'], keep='first')
    return merged


def _compute_locus_scale(df: pd.DataFrame) -> pd.Series:
    # Compute per-genome median span (kb), then ratio per locus, clamped to [0.5, 2.0]
    if 'genome' not in df.columns or 'locus_span_kb' not in df.columns:
        return pd.Series([1.0] * len(df))
    med_by_genome = df.groupby('genome')['locus_span_kb'].median().to_dict()
    def scale(row):
        base = med_by_genome.get(row['genome'], max(1.0, df['locus_span_kb'].median()))
        if base <= 0:
            base = 1.0
        r = float(row['locus_span_kb']) / float(base)
        return max(0.5, min(2.0, r))
    return df.apply(scale, axis=1)


def _find_gff_file(genome_gff_dir: Path, prefix: str) -> Path | None:
    """Find GFF file for a genome prefix (BK or LB) with flexible path matching."""
    # Define search patterns for each genome
    search_patterns = {
        'BK': [
            # NCBI reference GFF (has LOC IDs) - prioritize this!
            'Belonocnema_kinseyi*/GCF_010883055*.gff',
            '**/GCF_010883055*_genomic.gff',
            '**/GCF_010883055*.gff',
            # Direct file in root
            'genomic.gff',
            # Subdirectory patterns
            'belonocnema_kinseyi/*.gff',
            'belonocnema_kinseyi/**/*.gff',
            '**/bkinseyi*.gff',
            # Helixer GFF (no LOC IDs) - only as fallback
            'Belonocnema_kinseyi*/*_helixer.gff3',
        ],
        'LB': [
            # Subdirectory patterns
            'Leptopilina_boulardi/**/*.gff',
            'leptopilina_boulardi/**/*.gff',
            '**/GCF_015476425*.gff',
            '**/GCF_019393585*.gff',  # Alternative LB assembly
            '**/GCF_019393585*.gff3',
            '**/GCA_019393585*_helixer.gff3',
        ],
    }

    patterns = search_patterns.get(prefix, [])
    for pattern in patterns:
        matches = list(genome_gff_dir.glob(pattern))
        if matches:
            return matches[0]

    return None


def extract_cluster_gene_coordinates(phase1_dir: Path, genome_gff_dir: Path) -> None:
    """Extract and save individual gene coordinates for each cluster member from GFF files."""
    locus_defs = phase1_dir / 'locus_definitions.tsv'
    if not locus_defs.exists():
        return

    df = pd.read_csv(locus_defs, sep='\t')
    all_gene_coords = []

    # Cache GFF positions by genome prefix to avoid re-parsing
    gff_cache: dict[str, dict] = {}

    for _, row in df.iterrows():
        locus_id = row['locus_id']
        cluster_members = row.get('cluster_members', '')

        if not cluster_members:
            continue

        # Parse cluster members
        gene_ids = [g.strip() for g in str(cluster_members).split(',')]

        # Determine genome prefix
        prefix = None
        if locus_id.startswith('BK'):
            prefix = 'BK'
        elif locus_id.startswith('LB'):
            prefix = 'LB'

        if prefix is None:
            continue

        # Load GFF positions (cached)
        if prefix not in gff_cache:
            gff_file = _find_gff_file(genome_gff_dir, prefix)
            if gff_file and gff_file.exists():
                print(f"  Loading {prefix} GFF: {gff_file}")
                gff_cache[prefix] = _load_gff_gene_positions(gff_file)
            else:
                print(f"  WARNING: No GFF file found for {prefix} in {genome_gff_dir}")
                gff_cache[prefix] = {}

        positions = gff_cache.get(prefix, {})
        for gene_id in gene_ids:
            if gene_id in positions:
                chrom, start, end = positions[gene_id]
                all_gene_coords.append({
                    'locus_id': locus_id,
                    'gene_id': gene_id,
                    'chromosome': chrom,
                    'start': start,
                    'end': end,
                    'length': end - start + 1
                })

    if all_gene_coords:
        output_file = phase1_dir / 'cluster_gene_coordinates.tsv'
        coords_df = pd.DataFrame(all_gene_coords)
        coords_df.to_csv(output_file, sep='\t', index=False)
        print(f"  Saved {len(coords_df)} gene coordinates to {output_file}")


def postprocess_phase1(phase1_dir: Path, genome_gff_dir: Path | None = None) -> Path:
    locus_defs = phase1_dir / 'locus_definitions.tsv'
    if not locus_defs.exists():
        raise FileNotFoundError(f'Expected locus_definitions.tsv at {locus_defs}')

    df = pd.read_csv(locus_defs, sep='\t')

    # Prefer deduplicated flanking file
    def pick_dedup(path_str: str) -> str:
        p = Path(path_str)
        dedup = p.with_name(f"{p.stem}_dedup{p.suffix}")
        return str(dedup if dedup.exists() else p)

    if 'flanking_file' in df.columns:
        df['flanking_file'] = df['flanking_file'].astype(str).apply(pick_dedup)

    # Add expected chromosome from locus_id
    df['expected_chromosome'] = df['locus_id'].astype(str).apply(_expected_chr_from_locus_id)

    # Calculate flanking gene spans from FASTA files + GFF (preferred method)
    if genome_gff_dir and genome_gff_dir.exists():
        print("\n[Calculating flanking gene spans from FASTA files...]")
        flanking_spans = _calculate_flanking_spans_from_fasta(phase1_dir, genome_gff_dir)

        if not flanking_spans.empty:
            # Merge on locus_id (new method returns locus_id directly)
            df = df.merge(flanking_spans, on='locus_id', how='left', suffixes=('', '_new'))
            # Use new values, keeping any existing non-null values
            for col in ['flanking_span_start', 'flanking_span_end', 'flanking_span_kb']:
                if f'{col}_new' in df.columns:
                    df[col] = df[f'{col}_new'].fillna(df.get(col, pd.NA))
                    df.drop(columns=[f'{col}_new'], inplace=True, errors='ignore')
            print(f"Added flanking spans for {df['flanking_span_kb'].notna().sum()} loci")
        else:
            # Fallback to old debug file method
            flanking_spans = _calculate_flanking_spans(phase1_dir, genome_gff_dir)
            if not flanking_spans.empty and 'target_gene' in df.columns:
                df = df.merge(flanking_spans, on='target_gene', how='left')
                if {'flanking_span_start', 'flanking_span_end'} <= set(df.columns):
                    df['flanking_span_kb'] = ((df['flanking_span_end'].fillna(0).astype(float) -
                                              df['flanking_span_start'].fillna(0).astype(float)).clip(lower=0)) / 1000.0
                    print(f"Added flanking spans for {df['flanking_span_kb'].notna().sum()} loci (from debug files)")

    # Join locus envelope boundaries from debug files (legacy format)
    dbg = _read_debug_boundaries(phase1_dir)
    if not dbg.empty:
        df = df.merge(dbg, on='locus_id', how='left')
        df.rename(columns={'boundary_start': 'locus_span_start', 'boundary_end': 'locus_span_end'}, inplace=True)
    else:
        # Fallback to target extents if no flanking spans calculated
        if 'locus_span_start' not in df.columns and {'target_start', 'target_end'} <= set(df.columns):
            df['locus_span_start'] = df['target_start']
            df['locus_span_end'] = df['target_end']

    # Compute span kb if possible (for target gene span)
    if {'locus_span_start', 'locus_span_end'} <= set(df.columns):
        try:
            df['locus_span_kb'] = ((df['locus_span_end'].astype(float) - df['locus_span_start'].astype(float)).clip(lower=0)) / 1000.0
        except Exception:
            if 'locus_span_kb' not in df.columns:
                df['locus_span_kb'] = 0.0

    # Compute locus_scale
    df['locus_scale'] = _compute_locus_scale(df)

    # Write updated locus_definitions.tsv (in-place)
    df.to_csv(locus_defs, sep='\t', index=False)
    return locus_defs


def run_coordinates_mode(
    coords_file: Path,
    gene_family: str,
    out_dir: Path,
    genome_gff_dir: Path,
    n_flanking: int = 20
) -> None:
    """
    Create Phase 1 outputs directly from pre-computed coordinates file.

    This bypasses GFF LOC lookup - coordinates are already provided.
    Still needs GFF for extracting flanking gene sequences.
    """
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    print(f"\n[COORDINATES MODE] Reading {coords_file}")
    coords_df = pd.read_csv(coords_file, sep='\t')

    # Validate required columns
    required = {'subcluster_id', 'gene_id', 'chromosome', 'start', 'end'}
    missing = required - set(coords_df.columns)
    if missing:
        raise ValueError(f"Coordinates file missing columns: {missing}")

    # Filter to only the target gene family
    coords_df = coords_df[coords_df['subcluster_id'] == gene_family]
    if coords_df.empty:
        raise ValueError(f"No coordinates found for gene family: {gene_family}")
    print(f"  Filtered to {len(coords_df)} genes for family {gene_family}")

    # Group by subcluster_id (each subcluster = one locus)
    locus_groups = coords_df.groupby('subcluster_id')
    print(f"  Found {len(locus_groups)} loci (subclusters)")

    # Load BK GFF for flanking gene extraction
    bk_gff = _find_gff_file(genome_gff_dir, 'BK')
    if not bk_gff:
        raise FileNotFoundError(f"BK GFF not found in {genome_gff_dir}")
    print(f"  Loading BK GFF: {bk_gff}")
    gene_positions = _load_gff_gene_positions(bk_gff)
    print(f"  Loaded {len(gene_positions)} gene positions")

    # Load BK proteome for flanking sequences
    bk_proteome_path = None
    proteome_candidates = [
        genome_gff_dir / 'belonocnema_kinseyi' / 'GCF_010883055.1_B_treatae_v1_protein.faa',
        Path('/carc/scratch/projects/emartins/2016456/adam/genomes/belonocnema_kinseyi/GCF_010883055.1_B_treatae_v1_protein.faa'),
    ]
    for candidate in proteome_candidates:
        if candidate.exists():
            bk_proteome_path = candidate
            break

    # Also try glob
    if not bk_proteome_path:
        for candidate in genome_gff_dir.glob('**/GCF_010883055*protein*.faa'):
            bk_proteome_path = candidate
            break
        if not bk_proteome_path:
            for candidate in genome_gff_dir.glob('**/*kinseyi*protein*.faa'):
                bk_proteome_path = candidate
                break

    # Load proteome sequences
    proteome_seqs = {}
    if bk_proteome_path and bk_proteome_path.exists():
        print(f"  Loading BK proteome: {bk_proteome_path}")
        for record in SeqIO.parse(bk_proteome_path, 'fasta'):
            # Extract LOC ID from header (e.g., XP_033221023.1 -> find LOC in description)
            # Header format varies: "XP_033221023.1 protein description [Belonocnema kinseyi]"
            protein_id = record.id.split()[0]
            proteome_seqs[protein_id] = record
        print(f"  Loaded {len(proteome_seqs)} protein sequences")
    else:
        print(f"  WARNING: BK proteome not found, flanking FASTAs will be empty")

    # Build chromosome -> sorted gene list for flanking extraction
    chrom_genes = {}
    for gene_id, (chrom, start, end) in gene_positions.items():
        if chrom not in chrom_genes:
            chrom_genes[chrom] = []
        chrom_genes[chrom].append((start, end, gene_id))

    for chrom in chrom_genes:
        chrom_genes[chrom].sort(key=lambda x: x[0])

    # Chromosome name mapping
    chrom_map = {
        'NC_046657.1': 'chr1', 'NC_046658.1': 'chr2', 'NC_046659.1': 'chr3',
        'NC_046660.1': 'chr4', 'NC_046661.1': 'chr5', 'NC_046662.1': 'chr6',
        'NC_046663.1': 'chr7', 'NC_046664.1': 'chr8', 'NC_046665.1': 'chr9',
        'NC_046666.1': 'chr10', 'NC_046667.1': 'chr11', 'NC_046668.1': 'chr12',
        'NC_046669.1': 'chr13', 'NC_046670.1': 'chr14',
    }

    # Cluster genes into physical tandem loci (within 50kb of each other)
    # This matches the logic in phase1_discovery_impl.detect_input_tandem_clusters
    MAX_TANDEM_DISTANCE_KB = 50

    def cluster_genes_into_loci(genes_df):
        """Cluster genes within MAX_TANDEM_DISTANCE_KB into tandem loci."""
        # Sort by chromosome and position
        sorted_df = genes_df.sort_values(['chromosome', 'start'])

        clusters = []
        current_cluster = []

        for _, row in sorted_df.iterrows():
            if not current_cluster:
                current_cluster.append(row)
            else:
                prev_row = current_cluster[-1]
                # Check if same chromosome and within distance
                if (row['chromosome'] == prev_row['chromosome'] and
                    abs(row['start'] - prev_row['start']) < MAX_TANDEM_DISTANCE_KB * 1000):
                    current_cluster.append(row)
                else:
                    clusters.append(pd.DataFrame(current_cluster))
                    current_cluster = [row]

        if current_cluster:
            clusters.append(pd.DataFrame(current_cluster))

        return clusters

    # Process each family's genes into tandem loci
    locus_rows = []
    gene_coord_rows = []

    for subcluster_id, group in locus_groups:
        # Cluster this family's genes into physical tandem loci
        tandem_clusters = cluster_genes_into_loci(group)
        print(f"  {subcluster_id}: {len(group)} genes → {len(tandem_clusters)} tandem loci")

        for cluster_idx, cluster_df in enumerate(tandem_clusters):
            # Get locus envelope from this tandem cluster
            locus_start = cluster_df['start'].min()
            locus_end = cluster_df['end'].max()
            chrom = cluster_df['chromosome'].iloc[0]
            strand = cluster_df['strand'].iloc[0] if 'strand' in cluster_df.columns else '+'

            # Get cluster members
            cluster_members = ','.join(cluster_df['gene_id'].tolist())
            protein_ids = cluster_df['protein_id'].tolist() if 'protein_id' in cluster_df.columns else []

            # Create locus_id with chromosome and letter suffix
            chr_name = chrom_map.get(chrom, chrom.replace('NC_', 'unk'))
            # Use letter suffix for multiple loci (a, b, c, ...)
            suffix = chr(ord('a') + cluster_idx) if len(tandem_clusters) > 1 else ''
            display_locus_id = f"BK_{chr_name}_{subcluster_id}{suffix}"

            # Find flanking genes on this chromosome
            flanking_genes = []
            if chrom in chrom_genes:
                genes_on_chrom = chrom_genes[chrom]

                # Find genes that OVERLAP the locus (not strict containment)
                # Gene coordinates from GFF include UTRs, input coords are CDS-only
                locus_gene_indices = []
                for i, (gstart, gend, gid) in enumerate(genes_on_chrom):
                    # Check for any overlap: gene overlaps locus if it doesn't end before or start after
                    if gend >= locus_start and gstart <= locus_end:
                        locus_gene_indices.append(i)

                if locus_gene_indices:
                    first_idx = min(locus_gene_indices)
                    last_idx = max(locus_gene_indices)

                    # Get N upstream and N downstream
                    upstream_start = max(0, first_idx - n_flanking)
                    downstream_end = min(len(genes_on_chrom), last_idx + n_flanking + 1)

                    for i in range(upstream_start, first_idx):
                        flanking_genes.append(('U', genes_on_chrom[i][2]))
                    for i in range(last_idx + 1, downstream_end):
                        flanking_genes.append(('D', genes_on_chrom[i][2]))

            # Write flanking FASTA
            flanking_file = out_dir / f"{display_locus_id}_flanking.faa"
            flanking_records = []
            for direction, flank_gene_id in flanking_genes:
                # Find protein ID for this gene (from proteome)
                # Try to match by LOC ID in proteome headers
                for prot_id, record in proteome_seqs.items():
                    if flank_gene_id in record.description:
                        new_record = SeqRecord(
                            record.seq,
                            id=f"{prot_id}|{flank_gene_id}",
                            description=f"{direction} {prot_id} {record.description}"
                        )
                        flanking_records.append(new_record)
                        break

            if flanking_records:
                SeqIO.write(flanking_records, flanking_file, 'fasta')
            else:
                # Write empty file
                flanking_file.touch()

            # Record locus definition
            locus_rows.append({
                'locus_id': display_locus_id,
                'subcluster_id': subcluster_id,
                'gene_family': gene_family,  # For downstream compatibility
                'genome': 'BK',
                'chromosome': chrom,
                'strand': strand,
                'locus_span_start': locus_start,
                'locus_span_end': locus_end,
                'locus_span_kb': (locus_end - locus_start) / 1000.0,
                'n_genes': len(cluster_df),
                'cluster_members': cluster_members,
                'flanking_file': str(flanking_file),
                'n_flanking': len(flanking_genes),
            })

            # Record individual gene coordinates
            for _, gene_row in cluster_df.iterrows():
                gene_coord_rows.append({
                    'locus_id': display_locus_id,
                    'gene_id': gene_row['gene_id'],
                    'protein_id': gene_row.get('protein_id', ''),
                    'chromosome': gene_row['chromosome'],
                    'start': gene_row['start'],
                    'end': gene_row['end'],
                    'strand': gene_row.get('strand', '+'),
                    'length': gene_row['end'] - gene_row['start'] + 1,
                })

    # Write outputs
    locus_defs_df = pd.DataFrame(locus_rows)
    locus_defs_df.to_csv(out_dir / 'locus_definitions.tsv', sep='\t', index=False)
    print(f"  Wrote {len(locus_defs_df)} loci to locus_definitions.tsv")

    gene_coords_df = pd.DataFrame(gene_coord_rows)
    gene_coords_df.to_csv(out_dir / 'cluster_gene_coordinates.tsv', sep='\t', index=False)
    print(f"  Wrote {len(gene_coords_df)} gene coordinates to cluster_gene_coordinates.tsv")

    # Create query_proteins.faa from cluster member protein IDs
    query_records = []
    if 'protein_id' in coords_df.columns and proteome_seqs:
        # Get unique protein IDs from this family only
        family_df = coords_df[coords_df['subcluster_id'] == gene_family]
        protein_ids_to_extract = family_df['protein_id'].dropna().unique().tolist()

        for prot_id in protein_ids_to_extract:
            if prot_id in proteome_seqs:
                record = proteome_seqs[prot_id]
                # Get the gene_id for this protein
                matching_row = family_df[family_df['protein_id'] == prot_id].iloc[0]
                gene_id = matching_row.get('gene_id', '')
                new_record = SeqRecord(
                    record.seq,
                    id=prot_id,
                    description=f"{gene_id} {record.description}"
                )
                query_records.append(new_record)

    query_faa = out_dir / 'query_proteins.faa'
    if query_records:
        SeqIO.write(query_records, query_faa, 'fasta')
        print(f"  Wrote {len(query_records)} query proteins to query_proteins.faa")
    else:
        # Try to extract from gene coordinates if protein_id column is available
        print(f"  WARNING: No query proteins extracted (no matching protein IDs in proteome)")
        query_faa.touch()


def parse_args():
    p = argparse.ArgumentParser(
        description='Phase 1: locus discovery + envelope + scaling',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Legacy mode (LOC ID lookup):
  python 01_phase1.py --loc-ids LOC117167432,LOC117167433 --gene-family MC102 --output-dir out/phase1

  # Coordinates mode (pre-computed):
  python 01_phase1.py --coordinates-file coords.tsv --gene-family MC102 --output-dir out/phase1
"""
    )

    # Input mode (mutually exclusive)
    input_group = p.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--loc-ids', help='Comma-separated LOC IDs (legacy mode)')
    input_group.add_argument('--coordinates-file', type=Path,
                             help='TSV with pre-computed genomic coordinates (coordinates mode)')

    p.add_argument('--gene-family', required=True, help='Gene family label')
    p.add_argument('--output-dir', required=True, type=Path, help='Output directory for phase1 outputs')
    p.add_argument('--genome-gff-dir', type=Path, default=Path('data/reference'),
                   help='Directory containing genome GFF files')
    p.add_argument('--n-flanking', type=int, default=20,
                   help='Number of flanking genes per side (default: 20)')
    return p.parse_args()


def main():
    args = parse_args()
    out_dir = args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    print('=' * 80)
    print('PHASE 1: LOCUS DISCOVERY')
    print('=' * 80)
    print(f'Output dir: {out_dir}')
    print(f'Gene family: {args.gene_family}')

    # Determine mode
    if args.coordinates_file:
        # COORDINATES MODE: Use pre-computed coordinates
        print(f'Mode: COORDINATES (using {args.coordinates_file.name})')
        print('=' * 80)

        run_coordinates_mode(
            coords_file=args.coordinates_file,
            gene_family=args.gene_family,
            out_dir=out_dir,
            genome_gff_dir=args.genome_gff_dir,
            n_flanking=args.n_flanking
        )

    else:
        # LEGACY MODE: LOC ID lookup
        print(f'Mode: LEGACY (LOC ID lookup)')
        print('=' * 80)

        # 1) Run legacy discovery into the requested directory
        print('\n[1] Running legacy discovery...')
        run_legacy_discovery(args.loc_ids, args.gene_family, out_dir)

        # 2) Post-process locus definitions
        print('\n[2] Post-processing locus definitions...')
        locus_defs = postprocess_phase1(out_dir, genome_gff_dir=args.genome_gff_dir)
        print(f'  Updated: {locus_defs}')

        # 3) Extract individual gene coordinates (REQUIRED for Phase 4 calibration)
        print('\n[3] Extracting cluster gene coordinates from GFF...')
        extract_cluster_gene_coordinates(out_dir, args.genome_gff_dir)

    # Validate outputs (both modes)
    gene_coords_file = out_dir / 'cluster_gene_coordinates.tsv'
    locus_defs_file = out_dir / 'locus_definitions.tsv'

    if not gene_coords_file.exists():
        print(f'ERROR: Failed to create {gene_coords_file}')
        print('       This file is REQUIRED for Phase 4 factorial calibration.')
        raise FileNotFoundError(f'Missing required output: {gene_coords_file}')

    if not locus_defs_file.exists():
        print(f'ERROR: Failed to create {locus_defs_file}')
        raise FileNotFoundError(f'Missing required output: {locus_defs_file}')

    coords_df = pd.read_csv(gene_coords_file, sep='\t')
    if len(coords_df) == 0:
        print(f'ERROR: {gene_coords_file} is empty')
        print('       No gene coordinates could be extracted.')
        raise ValueError(f'Empty gene coordinates file: {gene_coords_file}')

    print(f'\n[VALIDATION] Extracted {len(coords_df)} gene coordinates')

    # Summary
    print('\n[SUMMARY]')
    try:
        df = pd.read_csv(locus_defs_file, sep='\t')
        cols = [c for c in ['locus_id', 'genome', 'chromosome', 'locus_span_kb', 'n_genes', 'n_flanking'] if c in df.columns]
        if cols:
            print(df[cols].head(10).to_string(index=False))
            if len(df) > 10:
                print(f'  ... and {len(df) - 10} more loci')
    except Exception as e:
        print(f'  Could not display summary: {e}')

    print('\nDone.')


if __name__ == '__main__':
    main()

