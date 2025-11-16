#!/usr/bin/env python3
"""Compare Phase 4 adaptive calibration results vs fixed threshold 0.2."""

from pathlib import Path

def count_phase1_bk_targets(family):
    """Count BK targets from Phase 1 per-locus target files."""
    phase1_dir = Path(f"outputs/{family}/phase1_v2")
    if not phase1_dir.exists():
        return 0

    total = 0
    for locus_dir in phase1_dir.glob("BK_*"):
        if locus_dir.is_dir():
            for faa_file in locus_dir.glob("*_targets.faa"):
                with open(faa_file) as f:
                    total += sum(1 for line in f if line.startswith('>'))
    return total

def count_phase4_bk_targets(family, version="v3"):
    """Count BK targets from Phase 4 output."""
    targets_file = Path(f"outputs/{family}/phase4_{version}/all_target_loci.tsv")
    if not targets_file.exists():
        return None

    with open(targets_file) as f:
        lines = f.readlines()
        if len(lines) <= 1:
            return 0
        return sum(1 for line in lines[1:] if 'Belonocnema_kinseyi' in line)

def main():
    # Get all families with Phase 1
    families = []
    for phase1_dir in sorted(Path("outputs").glob("*/phase1_v2")):
        family = phase1_dir.parent.name
        families.append(family)

    print("=" * 120)
    print("BK GENE RECOVERY: Phase 4 Adaptive vs Fixed Threshold 0.2")
    print("=" * 120)
    print()

    improved = []
    worsened = []
    unchanged = []
    still_missing = []
    still_extra = []
    perfect = []

    for family in families:
        phase1_count = count_phase1_bk_targets(family)
        v3_count = count_phase4_bk_targets(family, "v3")
        adaptive_count = count_phase4_bk_targets(family, "adaptive")

        if adaptive_count is None:
            continue

        # Compare to Phase 1
        v3_diff = abs(v3_count - phase1_count) if v3_count is not None else 999
        adaptive_diff = abs(adaptive_count - phase1_count)

        # Determine improvement
        if adaptive_count == phase1_count:
            status = "✓ PERFECT"
            perfect.append((family, phase1_count, v3_count, adaptive_count))
        elif adaptive_diff < v3_diff:
            improvement = v3_diff - adaptive_diff
            status = f"📈 IMPROVED (±{adaptive_diff}, was ±{v3_diff})"
            improved.append((family, phase1_count, v3_count, adaptive_count, improvement))
        elif adaptive_diff > v3_diff:
            regression = adaptive_diff - v3_diff
            status = f"📉 WORSE (±{adaptive_diff}, was ±{v3_diff})"
            worsened.append((family, phase1_count, v3_count, adaptive_count, regression))
        else:
            status = f"→ UNCHANGED (±{adaptive_diff})"
            unchanged.append((family, phase1_count, v3_count, adaptive_count))

        # Track if still has issues
        if adaptive_count < phase1_count:
            missing = phase1_count - adaptive_count
            still_missing.append((family, phase1_count, adaptive_count, missing))
        elif adaptive_count > phase1_count:
            extra = adaptive_count - phase1_count
            still_extra.append((family, phase1_count, adaptive_count, extra))

        print(f"{family:65} P1:{phase1_count:3}  v3:{v3_count:3}  adaptive:{adaptive_count:3}  {status}")

    print()
    print("=" * 120)
    print("SUMMARY")
    print("=" * 120)
    print(f"Perfect recovery: {len(perfect)}/{len(families)} families")
    print(f"Improved from v3: {len(improved)} families")
    print(f"Worse than v3: {len(worsened)} families")
    print(f"Unchanged: {len(unchanged)} families")
    print()
    print(f"Still missing BK genes: {len(still_missing)} families")
    print(f"Still extra BK genes: {len(still_extra)} families")

    if improved:
        print()
        print("📈 FAMILIES IMPROVED BY ADAPTIVE CALIBRATION:")
        for family, p1, v3, adaptive, improvement in improved:
            print(f"  {family}: {v3} → {adaptive} (closer to {p1} by {improvement})")

    if worsened:
        print()
        print("📉 FAMILIES WORSE WITH ADAPTIVE CALIBRATION:")
        for family, p1, v3, adaptive, regression in worsened:
            print(f"  {family}: {v3} → {adaptive} (farther from {p1} by {regression})")

    if still_missing:
        print()
        print("⚠️  FAMILIES STILL MISSING BK GENES:")
        for family, p1, adaptive, missing in still_missing:
            print(f"  {family}: {p1} → {adaptive} (missing {missing})")

if __name__ == "__main__":
    main()
