#!/usr/bin/env python3
"""
Sync locus_definitions.tsv with actual synteny block directories after Phase 3 filtering.

PROBLEM: Phase 1/12 creates initial locus_definitions.tsv, but Phase 2-3 may:
- Create different locus IDs due to tandem clustering/re-ranking
- Remove loci during filtering
- Change locus IDs during re-ranking

This script updates locus_definitions.tsv to match reality after Phase 3.

USAGE: Called automatically after Phase 3 in batch pipeline, or manually:
    python scripts/sync_locus_definitions.py --output-dir outputs/FAMILY_batch
"""

from pathlib import Path
import pandas as pd
import argparse
import sys


def sync_locus_definitions(output_dir):
    """
    Sync locus_definitions.tsv with actual synteny block directories.

    Strategy:
    1. Read current locus_definitions.tsv (from Phase 1/12)
    2. Get actual locus IDs from 02_synteny_blocks/ directories
    3. Match old->new by chromosome base (e.g., BK_chr1_a -> BK_chr1_b)
    4. Update locus_definitions.tsv with correct IDs
    5. Remove loci that don't exist anymore
    """

    output_dir = Path(output_dir)
    locus_defs_file = output_dir / "locus_definitions.tsv"
    phase12_locus_defs = output_dir / "phase12_landmark" / "locus_definitions.tsv"
    synteny_dir = output_dir / "02_synteny_blocks"
    phase12_dir = output_dir / "phase12_landmark"

    print(f"Syncing locus definitions for: {output_dir.name}")

    # Check prerequisites
    if not locus_defs_file.exists() and not phase12_locus_defs.exists():
        print(f"  ⚠ No locus_definitions.tsv found")
        return False

    # Use Phase 12 definitions if root doesn't exist yet
    if not locus_defs_file.exists() and phase12_locus_defs.exists():
        locus_defs_file = phase12_locus_defs

    if not synteny_dir.exists():
        print(f"  ⚠ No synteny_blocks directory (Phase 2 not complete)")
        return False

    # Get actual locus IDs from directories
    locus_dirs = [d.name for d in synteny_dir.iterdir()
                  if d.is_dir() and (d.name.startswith('BK_') or d.name.startswith('LB_'))]

    if not locus_dirs:
        print(f"  ⚠ No locus directories found in synteny_blocks")
        return False

    locus_ids_actual = set(locus_dirs)

    # Read current definitions
    try:
        df = pd.read_csv(locus_defs_file, sep='\t')
        locus_ids_defs = set(df['locus_id'].values)
    except Exception as e:
        print(f"  ❌ Error reading locus_definitions.tsv: {e}")
        return False

    # Check if sync needed
    if locus_ids_defs == locus_ids_actual:
        print(f"  ✓ Already synced ({len(locus_ids_actual)} loci)")

        # Still copy to root if it's only in phase12
        if locus_defs_file == phase12_locus_defs:
            output_locus_defs = output_dir / "locus_definitions.tsv"
            df.to_csv(output_locus_defs, sep='\t', index=False)
            print(f"  ✓ Copied to {output_locus_defs}")

        return True

    print(f"  Found mismatch - syncing...")
    print(f"    Definitions: {len(locus_ids_defs)} loci")
    print(f"    Actual:      {len(locus_ids_actual)} loci")

    # Create mapping from old to new locus IDs
    # Strategy: Match by chromosome/scaffold base name (e.g., BK_chr1)
    old_to_new = {}
    unmatched_actual = set(locus_ids_actual)

    for old_id in locus_ids_defs:
        # Extract base (e.g., "BK_chr1" from "BK_chr1_a")
        parts = old_id.split('_')
        if len(parts) < 3:
            continue

        old_base = f"{parts[0]}_{parts[1]}"  # e.g., "BK_chr1" or "LB_scf7713"

        # Find matching new ID with same base
        for new_id in locus_ids_actual:
            new_parts = new_id.split('_')
            if len(new_parts) < 3:
                continue

            new_base = f"{new_parts[0]}_{new_parts[1]}"

            if old_base == new_base and new_id in unmatched_actual:
                old_to_new[old_id] = new_id
                unmatched_actual.remove(new_id)
                break

    print(f"    Mapped: {len(old_to_new)} loci")

    # Handle loci that exist but aren't in current definitions
    # Try to add them from phase12 definitions
    if unmatched_actual and phase12_locus_defs.exists() and locus_defs_file != phase12_locus_defs:
        try:
            df_phase12 = pd.read_csv(phase12_locus_defs, sep='\t')

            rows_to_add = []
            for locus in unmatched_actual:
                if locus in df_phase12['locus_id'].values:
                    row = df_phase12[df_phase12['locus_id'] == locus].iloc[0]
                    rows_to_add.append(row)
                    print(f"    Added from phase12: {locus}")

            if rows_to_add:
                df = pd.concat([df, pd.DataFrame(rows_to_add)], ignore_index=True)
                unmatched_actual -= set([r['locus_id'] for r in rows_to_add])

        except Exception as e:
            print(f"    ⚠ Could not read phase12 definitions: {e}")

    if unmatched_actual:
        print(f"    ⚠ {len(unmatched_actual)} loci exist but no metadata:")
        for locus in sorted(unmatched_actual)[:5]:
            print(f"      - {locus}")

    # Update dataframe
    changes_made = False

    # Update existing loci with new IDs
    for old_id, new_id in old_to_new.items():
        if old_id != new_id:
            df.loc[df['locus_id'] == old_id, 'locus_id'] = new_id

            # Update flanking_file path
            old_path = str(df.loc[df['locus_id'] == new_id, 'flanking_file'].iloc[0])
            new_path = old_path.replace(f"{old_id}_flanking.faa", f"{new_id}_flanking.faa")
            df.loc[df['locus_id'] == new_id, 'flanking_file'] = new_path

            print(f"    Updated: {old_id} -> {new_id}")
            changes_made = True

    # Remove loci that don't exist anymore
    loci_to_remove = locus_ids_defs - set(old_to_new.keys())
    if loci_to_remove:
        df = df[df['locus_id'].isin(locus_ids_actual)]
        print(f"    Removed {len(loci_to_remove)} non-existent loci")
        changes_made = True

    # Copy flanking files with new names if needed
    if phase12_dir.exists():
        for old_id, new_id in old_to_new.items():
            if old_id == new_id:
                continue

            old_file = phase12_dir / f"{old_id}_flanking.faa"
            new_file = phase12_dir / f"{new_id}_flanking.faa"

            if old_file.exists() and not new_file.exists():
                import shutil
                shutil.copy(old_file, new_file)
                print(f"    Copied: {old_id}_flanking.faa -> {new_id}_flanking.faa")

    # Save updated definitions to root directory
    output_locus_defs = output_dir / "locus_definitions.tsv"
    df.to_csv(output_locus_defs, sep='\t', index=False)

    if changes_made:
        print(f"  ✓ Updated {output_locus_defs} ({len(df)} loci)")
    else:
        print(f"  ✓ Synced {output_locus_defs} ({len(df)} loci)")

    return True


def main():
    parser = argparse.ArgumentParser(
        description="Sync locus_definitions.tsv with actual synteny blocks after Phase 3"
    )
    parser.add_argument('--output-dir', required=True, type=Path,
                       help="Family output directory (e.g., outputs/FAMILY_batch)")

    args = parser.parse_args()

    success = sync_locus_definitions(args.output_dir)
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
