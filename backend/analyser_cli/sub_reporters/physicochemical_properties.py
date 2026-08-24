from functools import partial
from typing import Dict, List

from rpy2.robjects import conversion, default_converter
from rpy2.robjects.packages import importr

from ..analyser import Analyser
from ..report_base import SubReport, ReportFormatter


class PhysicochemicalProperties(SubReport):
    """Sub-report for effects of mutation on physicochemical properties of peptide."""

    def __init__(self):
       self.peptides_r = importr("Peptides")
       self.properties = {
           "Isoelectric Point": self.peptides_r.pI,
           "Molecular Weight": self.peptides_r.mw,
           "Charge": self.peptides_r.charge,
           "Hydrophobicity": partial(self.peptides_r.hydrophobicity, scale="KyteDoolittle"),
           "Instability Index": self.peptides_r.instaIndex,
       }
    
    
    @property
    def name(self) -> str:
        return "Physicochemical Properties"
    
    def generate_content(self, analyser: Analyser, mutations: List[Dict]) -> str:
        """Generate content showing changes in PC properties."""
        if not mutations or analyser.mutated_sequence is None:
            return "No mutations applied."
        
        report_lines = []
        for key, value in self.properties.items():
            # FastAPI runs sync endpoints in worker threads; set converter context explicitly for rpy2.
            with conversion.localconverter(default_converter):
                original_sequence_value = self._round_property_value(list(value(str(analyser.input_sequence)))[0])
                mutated_sequence_value = self._round_property_value(list(value(str(analyser.mutated_sequence)))[0])
            if mutated_sequence_value > original_sequence_value:
                change_symbol = "↑"
            elif mutated_sequence_value < original_sequence_value:
                change_symbol = "↓"
            else:
                change_symbol = "-"
            report_lines.extend([
                ReportFormatter.format_key_value(
                    f"{key} Original Sequence", self._format_property_value(original_sequence_value)
                ),
                ReportFormatter.format_key_value(
                    f"{key} Mutated Sequence", self._format_property_value(mutated_sequence_value)
                ),
                ReportFormatter.format_key_value("Change", change_symbol),
                ReportFormatter.format_detail_separator(),
            ])
        return "\n".join(report_lines[:-1])  # Assure that the last separator is not included, for aestetic reasons

    def _round_property_value(self, value: float) -> float:
        return round(float(value), 3)

    def _format_property_value(self, value: float) -> str:
        return f"{value:.3f}"
