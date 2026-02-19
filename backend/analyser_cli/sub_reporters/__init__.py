"""
Sub-reporters module for SNP Analyser

This module contains specialized reporting classes that generate focused analysis reports
for different aspects of DNA sequence analysis and mutation effects.

Each sub-reporter inherits from the SubReport base class and provides specific
analysis capabilities for different types of genetic analysis.
"""

from .mutation_effects import MutationEffectsReport
from .physicochemical_properties import PhysicochemicalProperties

__all__ = [
    'MutationEffectsReport',
    'PhysicochemicalProperties',
]

__version__ = "0.1.0"
