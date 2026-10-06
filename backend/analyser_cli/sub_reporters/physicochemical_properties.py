from functools import partial
from typing import Callable, Dict, List, Optional

from rpy2.rinterface_lib.embedded import RRuntimeError
from rpy2.robjects import conversion, default_converter
from rpy2.robjects.packages import importr

from ..analyser import Analyser
from ..report_base import SubReport, ReportFormatter


STOP_CODON_SYMBOL = "*"
UNAVAILABLE_VALUE = "n/a"


class PhysicochemicalProperties(SubReport):
    """Sub-report for effects of mutation on physicochemical properties of peptide."""

    def __init__(self):
       self.peptides_r = importr("Peptides")
       self.properties = {
           # Bjellqvist pK scale, the one ExPASy ProtParam uses, to keep both tools comparable.
           "Isoelectric Point": partial(self.peptides_r.pI, pKscale="Bjellqvist"),
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

        original_peptide = self._translated_peptide(analyser.amino_acid_translation)
        mutated_peptide = self._translated_peptide(analyser.mutated_amino_acid_translation)
        if not original_peptide:
            return "Original sequence does not encode any residue, no properties to compare."
        if not mutated_peptide:
            return "Mutated sequence stops at the first codon, no peptide is produced."

        report_lines = []
        for key, value in self.properties.items():
            original_sequence_value = self._property_value(value, original_peptide)
            mutated_sequence_value = self._property_value(value, mutated_peptide)
            if original_sequence_value is None or mutated_sequence_value is None:
                change_symbol = UNAVAILABLE_VALUE
            elif mutated_sequence_value > original_sequence_value:
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

    def _translated_peptide(self, translation) -> str:
        """Residues the ribosome actually incorporates, i.e. everything before the first stop codon.

        Properties must be computed on the amino acid translation, never on the DNA: the
        Peptides package would silently read A/T/G/C as Ala/Thr/Gly/Cys. The stop codon
        symbol is dropped because it is not a residue and Peptides cannot score it.
        """
        return str(translation or "").split(STOP_CODON_SYMBOL)[0]

    def _property_value(self, property_function: Callable, peptide: str) -> Optional[float]:
        """Score a peptide, or return None when the property is undefined for it.

        Some Peptides functions need a minimum length - instaIndex, for instance, works on
        dipeptides and fails below three residues - which a nonsense mutation can easily
        produce.
        """
        try:
            # FastAPI runs sync endpoints in worker threads; set converter context explicitly for rpy2.
            with conversion.localconverter(default_converter):
                return self._round_property_value(list(property_function(peptide))[0])
        except RRuntimeError:
            return None

    def _round_property_value(self, value: float) -> float:
        return round(float(value), 3)

    def _format_property_value(self, value: Optional[float]) -> str:
        if value is None:
            return UNAVAILABLE_VALUE
        return f"{value:.3f}"
