"""
SNP Analyser CLI Package

This package provides tools for analyzing DNA sequences and their mutations,
including amino acid translation effects and reporting capabilities.
"""

from .analyser import Analyser, get_codon_table
from .cds_annotation import (
    CdsAnnotation,
    TranslationException,
    genetic_code_name,
    read_cds_annotation,
)
from .protein_alignment import ProteinAlignment, align_translations, group_residues
from .reporter import Reporter, SubReport, ReportFormatter
from .fasta_reader import FastaReader
from .data_manager import DataManager

__all__ = [
    'Analyser',
    'get_codon_table',
    'CdsAnnotation',
    'TranslationException',
    'genetic_code_name',
    'read_cds_annotation',
    'ProteinAlignment',
    'align_translations',
    'group_residues',
    'Reporter', 
    'SubReport',
    'ReportFormatter',
    "FastaReader",
    "DataManager",
]

__version__ = "0.1.0"
