#!/usr/bin/env python3
"""Create a list of all loci with their family and locus_definitions file paths."""

from pathlib import Path
import pandas as pd

# Output file for locus list
output_file = Path("all_loci.txt")

# Find all locus_definitions files
locus_defs_files = sorted(Path("outputs").glob("*/phase1/locus_definitions.tsv"))

all_loci = []

for locus_def_file in locus_defs_files:
    family = locus_def_file.parent.parent.name

    # Read locus definitions
    df = pd.read_csv(locus_def_file, sep='\t')

    for _, row in df.iterrows():
        locus_id = row['locus_id']
        all_loci.append(f"{family}\t{locus_id}\t{locus_def_file}")

# Write to file
with open(output_file, 'w') as f:
    for line in all_loci:
        f.write(line + '\n')

print(f"Created list of {len(all_loci)} loci in {output_file}")