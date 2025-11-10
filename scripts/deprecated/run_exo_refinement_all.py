#!/usr/bin/env python3
"""
Run 06c (Exonerate-guided block re-selection) + 6b (Exonerate refinement)
across all families finished at phase 8, then summarize BK recovery.

Outputs:
  - outputs/<FAMILY>/03_filtered_blocks_exo/synteny_blocks_filtered.tsv
  - outputs/<FAMILY>/05_classified_exonerate_refined_exo/refined_targets.tsv
  - outputs/refined_bk_recovery_summary_exo.tsv
"""

from pathlib import Path
import subprocess
import json
import pandas as pd
import sys


def main():
    root = Path(__file__).resolve().parent.parent
    out_root = root / 'outputs'
    comp = json.load(open(out_root / 'pipeline_completion_analysis.json'))
    finished = [f['name'] for f in comp['families'] if f.get('stopped_at') == 'phase8']

    rows = []
    BK = 'GCA_010883055.1'

    for fam in finished:
        fam_dir = out_root / fam
        cand = fam_dir / '02_synteny_blocks' / 'all_synteny_blocks.tsv'
        exo_ref = fam_dir / '05_classified_exonerate_refined' / 'refined_targets.tsv'

        # Ensure a 6b baseline exists
        if not exo_ref.exists():
            blocks0 = fam_dir / '03_filtered_blocks' / 'synteny_blocks_filtered.tsv'
            if blocks0.exists():
                subprocess.run([
                    'python', str(root/'scripts/06b_refine_targets_with_exonerate.py'),
                    '--family-dir', str(fam_dir),
                    '--blocks', str(blocks0),
                    '--output-dir', str(fam_dir/'05_classified_exonerate_refined')
                ], check=True)
            else:
                print(f"Skipping {fam} (no initial 03_filtered_blocks)")
                continue

        if not cand.exists() or not exo_ref.exists():
            print(f"Skipping {fam} (missing candidates or refined targets)")
            continue

        # 06c reselect blocks
        out_blocks = fam_dir / '03_filtered_blocks_exo'
        out_blocks.mkdir(parents=True, exist_ok=True)
        subprocess.run([
            'python', str(root/'scripts/06c_reselect_blocks_with_exonerate.py'),
            '--family-dir', str(fam_dir),
            '--candidates', str(cand),
            '--exonerate-refined', str(exo_ref),
            '--output-dir', str(out_blocks)
        ], check=True)

        # 6b using reselected blocks
        out_dir = fam_dir / '05_classified_exonerate_refined_exo'
        out_dir.mkdir(parents=True, exist_ok=True)
        blocks = out_blocks / 'synteny_blocks_filtered.tsv'
        subprocess.run([
            'python', str(root/'scripts/06b_refine_targets_with_exonerate.py'),
            '--family-dir', str(fam_dir),
            '--blocks', str(blocks),
            '--output-dir', str(out_dir)
        ], check=True)

        # Summarize BK
        ref = out_dir / 'refined_targets.tsv'
        if not ref.exists():
            print(f"No refined targets for {fam}")
            continue
        df = pd.read_csv(ref, sep='\t')
        dfbk = df[(df['genome'] == BK) & (df['placement'] == 'synteny')]
        observed = dfbk['locus_id'].nunique()

        ldf_path = fam_dir / 'locus_definitions.tsv'
        expected = 0
        if ldf_path.exists():
            ldf = pd.read_csv(ldf_path, sep='\t')
            ldf_bk = ldf[ldf['genome'] == 'BK']
            expected = ldf_bk.apply(
                lambda r: int(r['input_cluster_size'])
                if (bool(r['is_tandem']) and pd.notna(r['input_cluster_size']))
                else 1,
                axis=1
            ).sum()

        rows.append({'family': fam, 'bk_expected': int(expected), 'bk_refined_inblock_exo': int(observed)})

    out = out_root / 'refined_bk_recovery_summary_exo.tsv'
    pd.DataFrame(rows).to_csv(out, sep='\t', index=False)
    print('Wrote', out)


if __name__ == '__main__':
    sys.exit(main())

