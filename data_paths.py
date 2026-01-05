#!/usr/bin/env python3
"""
Data Path Configuration for Helixer Pipeline

This module handles the per-species directory structure in pancynipoid_helixer.
All data access should go through this module for consistency.

Directory structure:
    data/ -> /carc/scratch/projects/emartins/2016456/adam/genomes/pancynipoid_helixer/
        {species_or_accession}/
            {name}_helixer.gff3
            {name}_helixer_proteins.faa
            {name}_helixer.dmnd
            ragtag.scaffold.fasta
            ragtag.scaffold.fasta.fai
"""

from pathlib import Path
from typing import Dict, Optional, List
import os

# =============================================================================
# BASE PATHS
# =============================================================================

# Pipeline root directory (where scripts live)
PIPELINE_DIR = Path(__file__).resolve().parent

# Data directory (symlink to pancynipoid_helixer)
DATA_DIR = PIPELINE_DIR / "data"

# Alternative: direct path if symlink not available
PANCYNIPOID_DIR = Path("/carc/scratch/projects/emartins/2016456/adam/genomes/pancynipoid_helixer")


def get_data_dir() -> Path:
    """Get the data directory, preferring symlink if available."""
    if DATA_DIR.exists():
        return DATA_DIR
    return PANCYNIPOID_DIR


# =============================================================================
# METADATA FILES
# =============================================================================

def get_species_mapping_file() -> Path:
    """Get path to gca_to_species.tsv (species names, families, BUSCO)."""
    return get_data_dir() / "gca_to_species.tsv"


def get_genome_quality_file() -> Path:
    """Get path to genome_quality_full.tsv (assembly metrics)."""
    return get_data_dir() / "genome_quality_full.tsv"


GCA_TO_SPECIES_FILE = get_species_mapping_file()
GENOME_QUALITY_FILE = get_genome_quality_file()


# =============================================================================
# GENOME DISCOVERY
# =============================================================================

def list_genomes() -> List[str]:
    """List all available genome directories."""
    data_dir = get_data_dir()
    if not data_dir.exists():
        return []
    return sorted([
        d.name for d in data_dir.iterdir()
        if d.is_dir() and not d.name.startswith('.')
    ])


def get_genome_dir(genome_id: str) -> Optional[Path]:
    """Get the directory for a specific genome."""
    data_dir = get_data_dir()
    genome_dir = data_dir / genome_id
    if genome_dir.exists():
        return genome_dir
    return None


# =============================================================================
# FILE ACCESSORS
# =============================================================================

def get_proteome_path(genome_id: str) -> Optional[Path]:
    """Get path to Helixer proteome FASTA for a genome."""
    genome_dir = get_genome_dir(genome_id)
    if genome_dir is None:
        return None
    # Try common naming patterns
    patterns = [
        f"{genome_id}_helixer_proteins.faa",
        f"{genome_id}_proteins.faa",
        f"{genome_id}.faa",
    ]
    for pattern in patterns:
        path = genome_dir / pattern
        if path.exists():
            return path
    return None


def get_gff_path(genome_id: str) -> Optional[Path]:
    """Get path to Helixer GFF3 annotation for a genome."""
    genome_dir = get_genome_dir(genome_id)
    if genome_dir is None:
        return None
    patterns = [
        f"{genome_id}_helixer.gff3",
        f"{genome_id}.gff3",
        f"{genome_id}.gff",
    ]
    for pattern in patterns:
        path = genome_dir / pattern
        if path.exists():
            return path
    return None


def get_diamond_db_path(genome_id: str) -> Optional[Path]:
    """Get path to DIAMOND database for a genome."""
    genome_dir = get_genome_dir(genome_id)
    if genome_dir is None:
        return None
    patterns = [
        f"{genome_id}_helixer.dmnd",
        f"{genome_id}.dmnd",
    ]
    for pattern in patterns:
        path = genome_dir / pattern
        if path.exists():
            return path
    return None


def get_assembly_path(genome_id: str) -> Optional[Path]:
    """Get path to scaffolded genome assembly for a genome."""
    genome_dir = get_genome_dir(genome_id)
    if genome_dir is None:
        return None
    patterns = [
        "ragtag.scaffold.fasta",
        f"{genome_id}.fasta",
        f"{genome_id}.fna",
    ]
    for pattern in patterns:
        path = genome_dir / pattern
        if path.exists():
            return path
    return None


# =============================================================================
# GENOME INFO
# =============================================================================

class GenomeInfo:
    """Container for genome paths and metadata."""

    def __init__(self, genome_id: str):
        self.genome_id = genome_id
        self.genome_dir = get_genome_dir(genome_id)
        self.proteome = get_proteome_path(genome_id)
        self.gff = get_gff_path(genome_id)
        self.diamond_db = get_diamond_db_path(genome_id)
        self.assembly = get_assembly_path(genome_id)

    @property
    def is_valid(self) -> bool:
        """Check if minimum required files exist."""
        return self.proteome is not None and self.gff is not None

    def to_dict(self) -> Dict[str, str]:
        """Return paths as dict (for compatibility with LANDMARKS format)."""
        return {
            'proteome': str(self.proteome) if self.proteome else None,
            'gff': str(self.gff) if self.gff else None,
            'db': str(self.diamond_db) if self.diamond_db else None,
            'assembly': str(self.assembly) if self.assembly else None,
        }

    def __repr__(self):
        return f"GenomeInfo({self.genome_id}, valid={self.is_valid})"


def get_genome_info(genome_id: str) -> GenomeInfo:
    """Get GenomeInfo for a specific genome."""
    return GenomeInfo(genome_id)


def get_all_genome_info() -> Dict[str, GenomeInfo]:
    """Get GenomeInfo for all available genomes."""
    return {gid: GenomeInfo(gid) for gid in list_genomes()}


# =============================================================================
# LANDMARK GENOMES
# =============================================================================

# Primary landmark genomes for synteny analysis
# These are well-annotated genomes used as anchors
LANDMARK_IDS = {
    'BK': 'Belonocnema_kinseyi_GCF',  # Primary cynipid reference
    'LB': 'GCA_019393585.1',           # Leptopilina boulardi (parasitoid outgroup)
}

def get_landmarks() -> Dict[str, GenomeInfo]:
    """Get GenomeInfo for landmark genomes."""
    return {code: GenomeInfo(gid) for code, gid in LANDMARK_IDS.items()}


def get_landmark_dict() -> Dict[str, Dict[str, str]]:
    """Get landmarks in legacy dict format for backward compatibility."""
    landmarks = {}
    for code, gid in LANDMARK_IDS.items():
        info = GenomeInfo(gid)
        if info.is_valid:
            landmarks[code] = info.to_dict()
    return landmarks


# =============================================================================
# VALIDATION
# =============================================================================

def validate_data_directory() -> Dict[str, any]:
    """Validate data directory and return status report."""
    data_dir = get_data_dir()
    report = {
        'data_dir': str(data_dir),
        'exists': data_dir.exists(),
        'genomes': [],
        'valid_count': 0,
        'invalid_count': 0,
        'landmarks_ok': True,
    }

    if not data_dir.exists():
        return report

    for gid in list_genomes():
        info = GenomeInfo(gid)
        genome_status = {
            'id': gid,
            'valid': info.is_valid,
            'has_proteome': info.proteome is not None,
            'has_gff': info.gff is not None,
            'has_diamond': info.diamond_db is not None,
        }
        report['genomes'].append(genome_status)
        if info.is_valid:
            report['valid_count'] += 1
        else:
            report['invalid_count'] += 1

    # Check landmarks
    for code, gid in LANDMARK_IDS.items():
        info = GenomeInfo(gid)
        if not info.is_valid:
            report['landmarks_ok'] = False

    return report


# =============================================================================
# CLI
# =============================================================================

if __name__ == "__main__":
    import json

    print("Data Path Configuration")
    print("=" * 60)

    report = validate_data_directory()
    print(f"Data directory: {report['data_dir']}")
    print(f"Exists: {report['exists']}")
    print(f"Valid genomes: {report['valid_count']}")
    print(f"Invalid genomes: {report['invalid_count']}")
    print(f"Landmarks OK: {report['landmarks_ok']}")

    print("\nLandmarks:")
    for code, info in get_landmarks().items():
        print(f"  {code}: {info}")

    print("\nFirst 5 genomes:")
    for gid in list_genomes()[:5]:
        info = GenomeInfo(gid)
        print(f"  {gid}: valid={info.is_valid}")
