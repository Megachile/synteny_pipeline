#!/usr/bin/env python3
"""
Run classification sweeps over existing 03 filtered blocks and 04 targets without re-running BLAST.

Variants tested per family:
  - overlap + pad_kb in {0,10,20}
  - allow_nearby + proximity_kb in {50,100,200} with pad_kb {10,20}
  - dynamic_nearby + nearby_frac_span in {0.3,0.5}, nearby_upper_kb in {200,500}, pad_kb {10,20}

Writes results to outputs/<FAM>/05_classified_tune_* and a summary TSV.
"""

from pathlib import Path
import argparse
import subprocess
import pandas as pd


def classify(fam_dir: Path, variant_name: str, args_list: list) -> int:
    blocks = fam_dir / '03_filtered_blocks' / 'synteny_blocks_filtered.tsv'
    targets = fam_dir / '04_target_genes' / 'all_target_loci.tsv'
    outdir = fam_dir / f'05_classified_{variant_name}'
    outdir.mkdir(exist_ok=True, parents=True)

    cmd = [
        'python', 'scripts/05_classify_targets.py',
        '--targets', str(targets),
        '--blocks', str(blocks),
        '--output-dir', str(outdir)
    ] + args_list
    print('Running:', ' '.join(cmd))
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        print('ERROR:', res.stderr[-800:])
        return 0
    # Count syntenic
    syn = outdir / 'syntenic_targets.tsv'
    if syn.exists():
        try:
            n = max(0, sum(1 for _ in open(syn)) - 1)
        except Exception:
            n = 0
        return n
    return 0


def sweep_family(fam_dir: Path) -> list:
    results = []
    # A) overlap only (pad)
    for pad in (0, 10, 20):
        name = f'overlap_pad{pad}'
        n = classify(fam_dir, name, ['--pad-kb', str(pad)])
        results.append({'family': fam_dir.name, 'variant': name, 'syntenic': n})

    # B) nearby fixed
    for prox in (50, 100, 200):
        for pad in (10, 20):
            name = f'nearby{prox}_pad{pad}'
            n = classify(fam_dir, name, ['--allow-nearby', '--proximity-kb', str(prox), '--pad-kb', str(pad)])
            results.append({'family': fam_dir.name, 'variant': name, 'syntenic': n})

    # C) dynamic nearby
    for frac in (0.3, 0.5):
        for upper in (200, 500):
            for pad in (10, 20):
                name = f'dyn_frac{frac}_cap{upper}_pad{pad}'
                n = classify(
                    fam_dir,
                    name,
                    ['--allow-nearby', '--dynamic-nearby', '--nearby-frac-span', str(frac), '--nearby-upper-kb', str(upper), '--pad-kb', str(pad)]
                )
                results.append({'family': fam_dir.name, 'variant': name, 'syntenic': n})

    return results


def find_families(outputs_root: Path):
    families = []
    for p in outputs_root.iterdir():
        if not p.is_dir():
            continue
        if (p / '04_target_genes' / 'all_target_loci.tsv').exists() and (p / '03_filtered_blocks' / 'synteny_blocks_filtered.tsv').exists():
            families.append(p)
    return sorted(families)


def main():
    ap = argparse.ArgumentParser(description='Sweep classification parameters across families using existing blocks and targets.')
    ap.add_argument('--family-dir', type=Path, help='Run only for this family (outputs/<FAM>)')
    ap.add_argument('--outputs-root', type=Path, default=Path('outputs'))
    args = ap.parse_args()

    fams = [args.family_dir] if args.family_dir else find_families(args.outputs_root)
    print(f'Found {len(fams)} families')

    all_rows = []
    for fam in fams:
        print(f'=== {fam.name} ===')
        rows = sweep_family(fam)
        all_rows.extend(rows)

    if all_rows:
        df = pd.DataFrame(all_rows)
        # Write per-family summaries to avoid array races
        if len(fams) == 1:
            out = fams[0] / 'reclassify_sweep_summary.tsv'
            df.to_csv(out, sep='\t', index=False)
            print(f'Summary: {out}')
        else:
            df.to_csv('reclassify_sweep_summary.tsv', sep='\t', index=False)
            print('Summary: reclassify_sweep_summary.tsv')


if __name__ == '__main__':
    main()
