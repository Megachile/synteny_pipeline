#!/usr/bin/env python3
"""
Phase 1: Locus discovery + locus envelope metadata + dedup flanking preference.

Part of the Helixer proteome pipeline (canonical synteny-scanning pipeline).

This wraps the discovery implementation to produce flanking/targets, then
post-processes locus_definitions.tsv to:
  - Prefer deduplicated flanking files ("*_flanking_dedup.faa")
  - Add expected chromosome from locus_id (e.g., BK_chr2_x → chr2)
  - Add locus envelope coordinates from debug flanking files when available
  - Compute a per-locus scale factor (locus_scale) based on span vs cohort median

Usage:
    python 01_phase1.py \
        --loc-ids LOC117167432,LOC117167433 \
        --gene-family ferritin_MC102 \
        --output-dir outputs/ferritin_MC102/phase1
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
            # Direct file in root
            'genomic.gff',
            # Subdirectory patterns
            'belonocnema_kinseyi/*.gff',
            'belonocnema_kinseyi/**/*.gff',
            '**/GCF_010883055*.gff',
            '**/bkinseyi*.gff',
        ],
        'LB': [
            # Subdirectory patterns
            'Leptopilina_boulardi/**/*.gff',
            'leptopilina_boulardi/**/*.gff',
            '**/GCF_015476425*.gff',
            '**/GCF_019393585*.gff',  # Alternative LB assembly
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


def parse_args():
    p = argparse.ArgumentParser(description='Phase 1 (redesign): locus discovery + envelope + scaling')
    p.add_argument('--loc-ids', required=True, help='Comma-separated LOC IDs (e.g., LOC...,LOC...)')
    p.add_argument('--gene-family', required=True, help='Gene family label')
    p.add_argument('--output-dir', required=True, type=Path, help='Output directory for phase1 outputs')
    p.add_argument('--genome-gff-dir', required=True, type=Path,
                   help='Directory containing genome GFF files (REQUIRED for per-gene coordinate extraction)')
    return p.parse_args()


def main():
    args = parse_args()
    out_dir = args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    print('=' * 80)
    print('PHASE 1 (REDESIGN): DISCOVERY + ENVELOPE + DEDUP PREFERENCE')
    print('=' * 80)
    print(f'Output dir: {out_dir}')

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

    # Validate that gene coordinates were extracted
    gene_coords_file = out_dir / 'cluster_gene_coordinates.tsv'
    if not gene_coords_file.exists():
        print(f'ERROR: Failed to create {gene_coords_file}')
        print('       This file is REQUIRED for Phase 4 factorial calibration.')
        print('       Check that GFF files exist in --genome-gff-dir')
        raise FileNotFoundError(f'Missing required output: {gene_coords_file}')

    coords_df = pd.read_csv(gene_coords_file, sep='\t')
    if len(coords_df) == 0:
        print(f'ERROR: {gene_coords_file} is empty')
        print('       No gene coordinates could be extracted from GFF files.')
        raise ValueError(f'Empty gene coordinates file: {gene_coords_file}')

    print(f'  ✓ Extracted {len(coords_df)} gene coordinates')

    # 4) Summary
    print('\n[4] Summary (preview):')
    try:
        df = pd.read_csv(locus_defs, sep='\t')
        cols = [c for c in ['locus_id','genome','expected_chromosome','locus_span_kb','locus_scale'] if c in df.columns]
        if cols:
            print(df[cols].to_string(index=False))
    except Exception:
        pass

    print('\nDone.')


if __name__ == '__main__':
    main()

