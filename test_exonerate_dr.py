#!/usr/bin/env python3
"""
Quick test: Run exonerate on DR loci to see what we extract
"""

import subprocess
from pathlib import Path
from Bio import SeqIO

# DR genome
genome_fasta = "data/ragtag_output/GCA_030998225.1/ragtag.scaffold.fasta"

# Query protein (BK ferritin from Phase 1)
query_protein = "outputs/ferritin_phase3_real/01_extracted_proteins/DR_CM021340.1_a/DR_CM021340.1_a_targets.faa"

# Test cases
test_cases = [
    {
        "name": "DR_CM021340.1_a (NO TARGET)",
        "scaffold": "CM021340.1_RagTag",
        "start": 43822130,
        "end": 44409502,
        "span_kb": 587.4
    },
    {
        "name": "DR_CM021340.1_a_BLAST_HIT (600kb away)",
        "scaffold": "CM021340.1_RagTag",
        "start": 45041637 - 10000,  # Add flanking for exonerate
        "end": 45045491 + 10000,
        "span_kb": 24  # ~24kb with flanking
    },
    {
        "name": "DR_CM021344.1_a (HAS TARGET)",
        "scaffold": "CM021344.1_RagTag",
        "start": 23165057,
        "end": 37027159,
        "span_kb": 13862  # Very large block
    }
]

print("=" * 80)
print("EXONERATE TEST: DR Loci")
print("=" * 80)
print(f"\nQuery: {query_protein}")
print(f"Genome: {genome_fasta}\n")

for test in test_cases:
    print(f"\n{'=' * 80}")
    print(f"TEST: {test['name']}")
    print(f"Scaffold: {test['scaffold']}")
    print(f"Coordinates: {test['start']:,} - {test['end']:,}")
    print(f"Span: {test['span_kb']:.1f} kb")
    print(f"{'=' * 80}\n")

    # Extract region using BioPython
    region_file = f"test_region_{test['name'].replace(' ', '_').replace('(', '').replace(')', '')}.fasta"

    print(f"Extracting region with BioPython...")
    try:
        # Parse genome FASTA and find the scaffold
        found = False
        for record in SeqIO.parse(genome_fasta, "fasta"):
            if record.id == test['scaffold']:
                # Extract the subsequence (1-based to 0-based indexing)
                subseq = record.seq[test['start']-1:test['end']]

                # Write to file
                with open(region_file, 'w') as f:
                    f.write(f">{test['scaffold']}:{test['start']}-{test['end']}\n")
                    f.write(str(subseq) + "\n")

                found = True
                break

        if not found:
            print(f"  ERROR: Scaffold {test['scaffold']} not found in genome")
            continue

        # Check if region file exists and has content
        region_path = Path(region_file)
        if not region_path.exists():
            print(f"  ERROR: Region file not created")
            continue

        region_size = region_path.stat().st_size
        print(f"  Region extracted: {region_size:,} bytes ({len(subseq):,} bp)")

    except Exception as e:
        print(f"  ERROR extracting region: {e}")
        continue

    # Run exonerate
    exonerate_output = f"test_exonerate_{test['name'].replace(' ', '_').replace('(', '').replace(')', '')}.txt"

    exonerate_cmd = [
        'exonerate',
        '--model', 'protein2genome',
        '--query', query_protein,
        '--target', region_file,
        '--showtargetgff', 'yes',
        '--showvulgar', 'no',
        '--showalignment', 'no',
        '--ryo', '%qi\\t%ti\\t%qab\\t%qae\\t%tab\\t%tae\\t%s\\t%C\\n'
    ]

    print(f"\nRunning exonerate...")
    result = subprocess.run(exonerate_cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"  ERROR running exonerate: {result.stderr}")
        continue

    # Parse results
    with open(exonerate_output, 'w') as f:
        f.write(result.stdout)

    # Count hits
    lines = [line for line in result.stdout.split('\n') if line and not line.startswith('#') and not line.startswith('--') and 'Command line' not in line]
    gene_lines = [line for line in lines if 'gene' in line.lower() or 'cds' in line.lower()]

    print(f"\nResults:")
    print(f"  Total output lines: {len(lines)}")
    print(f"  Gene/CDS features: {len(gene_lines)}")

    if gene_lines:
        print(f"\n  First few features:")
        for line in gene_lines[:5]:
            print(f"    {line}")
    else:
        print(f"  NO GENE FEATURES FOUND")

    print(f"\n  Full output saved to: {exonerate_output}")

print(f"\n{'=' * 80}")
print("TEST COMPLETE")
print("=" * 80)
