"""
SNP Analyser CLI Package

This package provides tools for analyzing DNA sequences and their mutations,
including amino acid translation effects and reporting capabilities.
"""

from .analyser import Analyser
from .reporter import Reporter, SubReport, ReportFormatter
from .fasta_reader import FastaReader

__all__ = [
    'Analyser',
    'Reporter', 
    'SubReport',
    'ReportFormatter',
    "FastaReader",
]

__version__ = "0.1.0"
