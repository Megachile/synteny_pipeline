#!/usr/bin/env python3
"""
Phase 8a (redesign wrapper): Generate locus-specific matrices using existing generator.

Resolves redesigned outputs to the expected inputs and calls scripts/generate_locus_matrices.py.
Also combines per-locus SwissProt annotation TSVs into a single file when needed.
"""

from __future__ import annotations

from pathlib import Path
import argparse
import subprocess
import sys


def combine_swissprot(swissprot_dir: Path) -> Path | None:
    """Combine per-locus SwissProt TSVs into a single file expected by Phase 8a."""
    if not swissprot_dir.exists():
        return None
    tsvs = sorted(swissprot_dir.glob('*_swissprot_annotations.tsv'))
    if not tsvs:
        return None
    out = swissprot_dir / 'genome_specific_swissprot_annotations.tsv'
    # Simple concat with header from first file
    with open(out, 'w') as w:
        for i, p in enumerate(tsvs):
            with open(p, 'r') as r:
                if i == 0:
                    w.write(r.read())
                else:
                    r.readline()  # skip header
                    w.write(r.read())
    return out


def parse_args():
    p = argparse.ArgumentParser(description='Phase 8a (redesign): locus matrices wrapper')
    p.add_argument('--family', required=True, help='Gene family name (folder under outputs)')
    p.add_argument('--reference-proteins', type=Path, default=Path('data/reference/protein.faa'))
    p.add_argument('--species-map', type=Path, default=Path('data/gca_to_species.tsv'))
    p.add_argument('--out-dir', type=Path, help='Override output directory (default outputs/<family>/phase8a_v2)')
    return p.parse_args()


def main():
    args = parse_args()

    base = Path('outputs') / args.family
    phase1 = base / 'phase1_v2'
    phase2 = base / 'phase2_synteny_v2'
    phase5 = base / 'phase5_v2'
    phase6 = base / 'phase6_extracted_v2'
    phase7 = base / 'phase7_v2'
    # Write Phase 8a outputs into the unified Phase 8 folder by default
    out_dir = args.out_dir if args.out_dir else (base / 'phase8_v2')
    out_dir.mkdir(parents=True, exist_ok=True)

    locus_defs = phase1 / 'locus_definitions.tsv'
    synteny_dir = phase2
    blocks = phase5 / 'selected_synteny_blocks.tsv'  # preferred
    if not blocks.exists():
        # Fallback to Phase 3_v2 filtered if present
        p3 = base / 'phase3_v2'
        cand = p3 / 'synteny_blocks_filtered.tsv'
        blocks = cand if cand.exists() else blocks

    targets = phase5 / 'all_targets_classified.tsv'
    swissprot = combine_swissprot(phase7)

    if not locus_defs.exists():
        print(f'ERROR: Missing locus_definitions.tsv: {locus_defs}', file=sys.stderr)
        sys.exit(1)
    if not synteny_dir.exists():
        print(f'ERROR: Missing synteny directory: {synteny_dir}', file=sys.stderr)
        sys.exit(1)
    if not targets.exists():
        print(f'ERROR: Missing all_targets_classified.tsv: {targets}', file=sys.stderr)
        sys.exit(1)

    cmd = [
        sys.executable,
        str(Path(__file__).resolve().parent.parent / 'scripts' / 'generate_locus_matrices.py'),
        '--locus-defs', str(locus_defs),
        '--synteny-dir', str(synteny_dir),
        '--blocks', str(blocks),
        '--targets', str(targets),
        '--species-map', str(args.species_map),
        '--output-dir', str(out_dir),
    ]
    if swissprot and swissprot.exists():
        cmd.extend(['--swissprot', str(swissprot)])
    if args.reference_proteins and args.reference_proteins.exists():
        cmd.extend(['--reference-proteins', str(args.reference_proteins)])
    if phase6.exists():
        cmd.extend(['--extracted-seqs', str(phase6)])

    print('Running:', ' '.join(cmd))
    subprocess.run(cmd, check=True)
    print('Phase 8a complete:', out_dir)


if __name__ == '__main__':
    main()
