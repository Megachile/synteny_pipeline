#!/usr/bin/env python3
"""
Helixer QC Inventory Script
Collects comprehensive metrics for all genomes in ragtag_output
Implements Phases 1-3 of HELIXER_QC_PLAN.md
"""

import os
import sys
import subprocess
from pathlib import Path
from collections import defaultdict
import json

# Base paths - use data_paths module for consistency
from data_paths import get_data_dir, list_genomes, get_genome_dir
RAGTAG_BASE = get_data_dir()

def get_genome_dirs():
    """Get all genome directories (excluding 'excluded')"""
    return sorted([d for d in RAGTAG_BASE.iterdir()
                   if d.is_dir() and d.name != 'excluded' and not d.name.endswith('.py')])

def get_fai_stats(genome_dir):
    """Get assembly stats from .fai file"""
    fai_files = list(genome_dir.glob("*.fai")) + list(genome_dir.glob("ragtag.scaffold.fasta.fai"))
    if not fai_files:
        return None

    fai_file = fai_files[0]
    scaffolds = []
    cm_scaffolds = []

    with open(fai_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            name = parts[0]
            length = int(parts[1])
            scaffolds.append((name, length))
            if name.startswith('CM'):
                cm_scaffolds.append((name, length))

    if not scaffolds:
        return None

    # Calculate N50
    total_length = sum(length for _, length in scaffolds)
    sorted_lengths = sorted([length for _, length in scaffolds], reverse=True)
    cumsum = 0
    n50 = 0
    for length in sorted_lengths:
        cumsum += length
        if cumsum >= total_length / 2:
            n50 = length
            break

    cm_total = sum(length for _, length in cm_scaffolds)

    return {
        'total_length_bp': total_length,
        'total_length_mb': total_length / 1_000_000,
        'scaffold_count': len(scaffolds),
        'n50': n50,
        'n50_kb': n50 / 1000,
        'cm_scaffold_count': len(cm_scaffolds),
        'cm_total_bp': cm_total,
        'cm_percent': (cm_total / total_length * 100) if total_length > 0 else 0,
        'largest_scaffold': max(length for _, length in scaffolds),
        'smallest_scaffold': min(length for _, length in scaffolds),
        'fai_file': str(fai_file)
    }

def get_helixer_gff_stats(genome_dir):
    """Get stats from Helixer GFF output"""
    genome_name = genome_dir.name
    # Look for GFF inside the genome directory
    gff_files = list(genome_dir.glob(f"{genome_name}_helixer.gff3")) + \
                list(genome_dir.glob(f"{genome_name}_helixer.gff")) + \
                list(genome_dir.glob("*_helixer.gff3")) + \
                list(genome_dir.glob("*_helixer.gff"))

    if not gff_files:
        return None

    gff_file = gff_files[0]

    gene_count = 0
    mrna_count = 0
    scaffolds_with_genes = set()
    cm_gene_count = 0
    last_scaffold = None

    with open(gff_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            scaffold = parts[0]
            feature_type = parts[2]
            last_scaffold = scaffold

            if feature_type == 'gene':
                gene_count += 1
                scaffolds_with_genes.add(scaffold)
                if scaffold.startswith('CM'):
                    cm_gene_count += 1
            elif feature_type == 'mRNA':
                mrna_count += 1

    return {
        'gff_exists': True,
        'gene_count': gene_count,
        'mrna_count': mrna_count,
        'scaffolds_with_genes': len(scaffolds_with_genes),
        'cm_gene_count': cm_gene_count,
        'cm_gene_percent': (cm_gene_count / gene_count * 100) if gene_count > 0 else 0,
        'last_scaffold_in_gff': last_scaffold,
        'gff_file': str(gff_file),
        'gff_size_mb': gff_file.stat().st_size / 1_000_000
    }

def get_protein_stats(genome_dir):
    """Get stats from extracted proteins"""
    genome_name = genome_dir.name
    # Look for protein file inside genome directory
    protein_files = list(genome_dir.glob(f"{genome_name}_helixer_proteins.faa")) + \
                    list(genome_dir.glob("*_helixer_proteins.faa")) + \
                    list(genome_dir.glob("*_proteins.faa"))

    if not protein_files:
        return None

    protein_file = protein_files[0]

    protein_count = 0
    total_length = 0
    lengths = []

    current_length = 0
    with open(protein_file) as f:
        for line in f:
            if line.startswith('>'):
                if current_length > 0:
                    lengths.append(current_length)
                    total_length += current_length
                protein_count += 1
                current_length = 0
            else:
                current_length += len(line.strip())
        if current_length > 0:
            lengths.append(current_length)
            total_length += current_length

    avg_length = total_length / protein_count if protein_count > 0 else 0

    return {
        'protein_exists': True,
        'protein_count': protein_count,
        'total_aa': total_length,
        'avg_protein_length': avg_length,
        'median_protein_length': sorted(lengths)[len(lengths)//2] if lengths else 0,
        'protein_file': str(protein_file),
        'protein_file_mb': protein_file.stat().st_size / 1_000_000
    }

def get_busco_stats(genome_dir):
    """Get BUSCO results if they exist"""
    genome_name = genome_dir.name
    # Look for BUSCO results inside genome directory
    busco_dirs = list(genome_dir.glob(f"busco_{genome_name}")) + \
                 list(genome_dir.glob(f"{genome_name}_busco")) + \
                 list(genome_dir.glob("busco_*"))

    summaries = []
    for busco_dir in busco_dirs:
        if busco_dir.is_dir():
            summaries.extend(list(busco_dir.glob("short_summary*.txt")))

    if not summaries:
        return None

    summary = summaries[0]

    with open(summary) as f:
        content = f.read()

    # Parse BUSCO results
    import re

    complete_match = re.search(r'C:(\d+\.?\d*)%', content)
    single_match = re.search(r'S:(\d+\.?\d*)%', content)
    dup_match = re.search(r'D:(\d+\.?\d*)%', content)
    frag_match = re.search(r'F:(\d+\.?\d*)%', content)
    missing_match = re.search(r'M:(\d+\.?\d*)%', content)

    return {
        'busco_exists': True,
        'complete': float(complete_match.group(1)) if complete_match else 0,
        'single': float(single_match.group(1)) if single_match else 0,
        'duplicated': float(dup_match.group(1)) if dup_match else 0,
        'fragmented': float(frag_match.group(1)) if frag_match else 0,
        'missing': float(missing_match.group(1)) if missing_match else 0,
        'summary_file': str(summary)
    }

def check_helixer_completeness(genome_name, assembly_stats, gff_stats):
    """Determine if Helixer ran to completion"""
    if not gff_stats:
        return {'complete': False, 'reason': 'No GFF file'}

    if not assembly_stats:
        return {'complete': None, 'reason': 'Cannot verify (no FAI)'}

    # Check if last scaffold in GFF matches something reasonable
    # For complete runs, we'd expect genes on most scaffolds
    scaffold_coverage = gff_stats['scaffolds_with_genes'] / assembly_stats['scaffold_count'] * 100

    # Heuristics for completeness
    if gff_stats['gene_count'] < 1000:
        return {'complete': False, 'reason': f'Very low gene count ({gff_stats["gene_count"]})', 'scaffold_coverage': scaffold_coverage}

    if scaffold_coverage < 1:
        return {'complete': False, 'reason': f'Very low scaffold coverage ({scaffold_coverage:.2f}%)', 'scaffold_coverage': scaffold_coverage}

    return {'complete': True, 'reason': 'Appears complete', 'scaffold_coverage': scaffold_coverage}

def assess_assembly_quality(stats):
    """Assess if assembly is suitable for Helixer"""
    if not stats:
        return {'quality': 'unknown', 'issues': ['No FAI file found']}

    issues = []

    if stats['scaffold_count'] > 1_000_000:
        issues.append(f"Extremely fragmented ({stats['scaffold_count']:,} scaffolds)")
    elif stats['scaffold_count'] > 500_000:
        issues.append(f"Very fragmented ({stats['scaffold_count']:,} scaffolds)")

    if stats['n50'] < 1000:
        issues.append(f"Very low N50 ({stats['n50']:,} bp)")
    elif stats['n50'] < 10000:
        issues.append(f"Low N50 ({stats['n50']:,} bp)")

    if stats['cm_percent'] < 10:
        issues.append(f"CM scaffolds only {stats['cm_percent']:.1f}% of genome")

    if stats['total_length_mb'] > 1500:
        issues.append(f"Very large genome ({stats['total_length_mb']:.0f} Mb)")

    if not issues:
        quality = 'good'
    elif len(issues) == 1 and 'fragmented' not in issues[0].lower():
        quality = 'acceptable'
    else:
        quality = 'poor'

    return {'quality': quality, 'issues': issues}

def main():
    print("=" * 80)
    print("HELIXER QC INVENTORY")
    print("=" * 80)

    genome_dirs = get_genome_dirs()
    print(f"\nFound {len(genome_dirs)} genome directories\n")

    # Collect all data
    all_data = []

    for genome_dir in genome_dirs:
        genome_name = genome_dir.name
        print(f"Processing {genome_name}...")

        entry = {'genome': genome_name}

        # Assembly stats
        assembly_stats = get_fai_stats(genome_dir)
        entry['assembly'] = assembly_stats

        # Assembly quality assessment
        entry['assembly_quality'] = assess_assembly_quality(assembly_stats)

        # Helixer GFF stats
        gff_stats = get_helixer_gff_stats(genome_dir)
        entry['helixer_gff'] = gff_stats

        # Protein stats
        protein_stats = get_protein_stats(genome_dir)
        entry['proteins'] = protein_stats

        # BUSCO stats
        busco_stats = get_busco_stats(genome_dir)
        entry['busco'] = busco_stats

        # Completeness check
        entry['completeness'] = check_helixer_completeness(genome_name, assembly_stats, gff_stats)

        all_data.append(entry)

    # Generate summary TSV
    print("\n" + "=" * 80)
    print("GENERATING SUMMARY")
    print("=" * 80)

    tsv_file = Path("/carc/scratch/projects/emartins/2016456/adam/synteny_scanning/helixer_pipeline/helixer_qc_summary.tsv")

    with open(tsv_file, 'w') as f:
        # Header
        headers = [
            'genome',
            'assembly_size_mb', 'scaffold_count', 'n50_kb', 'cm_scaffolds', 'cm_percent',
            'assembly_quality', 'assembly_issues',
            'helixer_gff_exists', 'gene_count', 'mrna_count', 'scaffolds_with_genes', 'cm_gene_percent',
            'helixer_complete', 'completion_reason',
            'protein_exists', 'protein_count', 'avg_protein_length',
            'busco_exists', 'busco_complete', 'busco_single', 'busco_dup', 'busco_frag', 'busco_missing',
            'status_category'
        ]
        f.write('\t'.join(headers) + '\n')

        for entry in all_data:
            asm = entry['assembly'] or {}
            gff = entry['helixer_gff'] or {}
            prot = entry['proteins'] or {}
            busco = entry['busco'] or {}
            quality = entry['assembly_quality']
            complete = entry['completeness']

            # Determine status category
            busco_complete = busco.get('complete', 0)
            if quality['quality'] == 'poor':
                category = 'ASSEMBLY_ISSUE'
            elif not gff:
                category = 'NO_HELIXER'
            elif not complete.get('complete', False):
                category = 'INCOMPLETE_HELIXER'
            elif busco_complete < 30:
                category = 'LOW_BUSCO'
            elif busco_complete < 60:
                category = 'ACCEPTABLE'
            else:
                category = 'GOOD'

            row = [
                entry['genome'],
                f"{asm.get('total_length_mb', ''):.1f}" if asm.get('total_length_mb') else '',
                str(asm.get('scaffold_count', '')),
                f"{asm.get('n50_kb', ''):.1f}" if asm.get('n50_kb') else '',
                str(asm.get('cm_scaffold_count', '')),
                f"{asm.get('cm_percent', ''):.1f}" if asm.get('cm_percent') else '',
                quality['quality'],
                '; '.join(quality['issues']) if quality['issues'] else '',
                'Y' if gff else 'N',
                str(gff.get('gene_count', '')),
                str(gff.get('mrna_count', '')),
                str(gff.get('scaffolds_with_genes', '')),
                f"{gff.get('cm_gene_percent', ''):.1f}" if gff.get('cm_gene_percent') is not None else '',
                'Y' if complete.get('complete') else 'N',
                complete.get('reason', ''),
                'Y' if prot else 'N',
                str(prot.get('protein_count', '')),
                f"{prot.get('avg_protein_length', ''):.1f}" if prot.get('avg_protein_length') else '',
                'Y' if busco else 'N',
                f"{busco.get('complete', '')}" if busco.get('complete') is not None else '',
                f"{busco.get('single', '')}" if busco.get('single') is not None else '',
                f"{busco.get('duplicated', '')}" if busco.get('duplicated') is not None else '',
                f"{busco.get('fragmented', '')}" if busco.get('fragmented') is not None else '',
                f"{busco.get('missing', '')}" if busco.get('missing') is not None else '',
                category
            ]
            f.write('\t'.join(row) + '\n')

    print(f"\nSummary written to: {tsv_file}")

    # Print category counts
    print("\n" + "=" * 80)
    print("STATUS SUMMARY")
    print("=" * 80)

    categories = defaultdict(list)
    for entry in all_data:
        asm = entry['assembly'] or {}
        gff = entry['helixer_gff'] or {}
        busco = entry['busco'] or {}
        quality = entry['assembly_quality']
        complete = entry['completeness']
        busco_complete = busco.get('complete', 0)

        if quality['quality'] == 'poor':
            categories['ASSEMBLY_ISSUE'].append(entry['genome'])
        elif not gff:
            categories['NO_HELIXER'].append(entry['genome'])
        elif not complete.get('complete', False):
            categories['INCOMPLETE_HELIXER'].append(entry['genome'])
        elif busco_complete < 30:
            categories['LOW_BUSCO'].append(entry['genome'])
        elif busco_complete < 60:
            categories['ACCEPTABLE'].append(entry['genome'])
        else:
            categories['GOOD'].append(entry['genome'])

    for cat in ['GOOD', 'ACCEPTABLE', 'LOW_BUSCO', 'INCOMPLETE_HELIXER', 'NO_HELIXER', 'ASSEMBLY_ISSUE']:
        count = len(categories[cat])
        print(f"\n{cat}: {count}")
        if count <= 10:
            for g in categories[cat]:
                print(f"  - {g}")
        else:
            for g in categories[cat][:5]:
                print(f"  - {g}")
            print(f"  ... and {count - 5} more")

    # Save full JSON for detailed analysis
    json_file = Path("/carc/scratch/projects/emartins/2016456/adam/synteny_scanning/helixer_pipeline/helixer_qc_full.json")
    with open(json_file, 'w') as f:
        json.dump(all_data, f, indent=2, default=str)
    print(f"\nFull data saved to: {json_file}")

    # List genomes needing BUSCO
    need_busco = [entry['genome'] for entry in all_data
                  if entry['proteins'] and not entry['busco']]

    if need_busco:
        print(f"\n" + "=" * 80)
        print(f"GENOMES NEEDING BUSCO ({len(need_busco)})")
        print("=" * 80)
        for g in need_busco:
            print(f"  - {g}")

        # Save list for batch BUSCO
        busco_list = Path("/carc/scratch/projects/emartins/2016456/adam/synteny_scanning/helixer_pipeline/need_busco.txt")
        with open(busco_list, 'w') as f:
            for g in need_busco:
                f.write(f"{g}\n")
        print(f"\nList saved to: {busco_list}")

if __name__ == "__main__":
    main()
