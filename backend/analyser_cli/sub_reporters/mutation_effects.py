from typing import Dict, List
from ..analyser import Analyser
from ..report_base import SubReport, ReportFormatter


class MutationEffectsReport(SubReport):
    """Sub-report for DNA mutations and their amino acid effects."""
    
    @property
    def name(self) -> str:
        return "Mutation Effects"
    
    def generate_content(self, analyser: Analyser, mutations: List[Dict]) -> str:
        """Generate content showing DNA mutations and amino acid changes."""
        if not mutations or analyser.mutated_sequence is None:
            return "No mutations applied."
        
        report_lines = [
            ReportFormatter.format_key_value("DNA Sequence Type", analyser.dna_type.title()),
            ReportFormatter.format_key_value("Original Sequence Length", f"{len(analyser.input_sequence)} bp"),
            ReportFormatter.format_key_value("Number of Mutations", str(len(mutations))),
            "",
            ReportFormatter.format_detail_separator(),
            "MUTATION DETAILS:",
            ReportFormatter.format_detail_separator()
        ]
        
        # Group mutations by codon (3-nucleotide groups)
        codon_effects = self._analyze_codon_effects(analyser, mutations)
        
        for i, mutation in enumerate(mutations, 1):
            pos = mutation["position"]
            original_base = str(analyser.input_sequence[pos - 1])
            new_base = mutation["base"].upper()
            
            report_lines.extend([
                "",
                f"Mutation #{i}:",
                ReportFormatter.format_key_value("Position", str(pos), indent=2),
                ReportFormatter.format_key_value("Change", f"{original_base} → {new_base}", indent=2),
                ReportFormatter.format_key_value("Codon Position", self._get_codon_position(pos), indent=2),
            ])
            
            # Add amino acid change information
            if pos in codon_effects:
                effect = codon_effects[pos]
                report_lines.extend([
                    ReportFormatter.format_key_value("Codon", f"{effect['original_codon']} → {effect['mutated_codon']}", indent=2),
                    ReportFormatter.format_key_value("Amino Acid", f"{effect['original_aa']} → {effect['mutated_aa']}", indent=2),
                    ReportFormatter.format_key_value("Effect Type", effect['effect_type'], indent=2)
                ])
        
        # Add summary
        mutated_length = len(analyser.mutated_amino_acid_translation) if analyser.mutated_amino_acid_translation else 0
        report_lines.extend([
            "",
            "",
            ReportFormatter.format_detail_separator(),
            "AMINO ACID TRANSLATION SUMMARY:",
            ReportFormatter.format_detail_separator(),
            ReportFormatter.format_key_value("Original", str(analyser.amino_acid_translation)),
            ReportFormatter.format_key_value("Mutated", str(analyser.mutated_amino_acid_translation)),
            "",
            ReportFormatter.format_key_value("Translation Length", f"{len(analyser.amino_acid_translation)} → {mutated_length} amino acids")
        ])
        
        return "\n".join(report_lines)
    
    def _get_codon_position(self, nucleotide_pos: int) -> str:
        """Determine codon number and position within codon."""
        codon_num = (nucleotide_pos - 1) // 3 + 1
        pos_in_codon = (nucleotide_pos - 1) % 3 + 1
        return f"Codon {codon_num}, Position {pos_in_codon}"
    
    def _analyze_codon_effects(self, analyser: Analyser, mutations: List[Dict]) -> Dict:
        """Analyze the effect of mutations on codons and amino acids."""
        effects = {}
        
        original_seq = str(analyser.input_sequence)
        mutated_seq = str(analyser.mutated_sequence)
        
        for mutation in mutations:
            pos = mutation["position"]
            codon_start = ((pos - 1) // 3) * 3
            codon_end = codon_start + 3
            
            if codon_end <= len(original_seq):
                original_codon = original_seq[codon_start:codon_end]
                mutated_codon = mutated_seq[codon_start:codon_end]
                
                # Get amino acid translations for these codons
                original_aa = self._translate_codon(original_codon, analyser.dna_type)
                mutated_aa = self._translate_codon(mutated_codon, analyser.dna_type)
                
                # Determine effect type
                if original_aa == mutated_aa:
                    effect_type = "Silent (synonymous)"
                elif mutated_aa == "*":
                    effect_type = "Nonsense (stop codon)"
                elif original_aa == "*":
                    effect_type = "Readthrough (stop codon lost)"
                else:
                    effect_type = "Missense (amino acid change)"
                
                effects[pos] = {
                    "original_codon": original_codon,
                    "mutated_codon": mutated_codon,
                    "original_aa": original_aa,
                    "mutated_aa": mutated_aa,
                    "effect_type": effect_type
                }
        
        return effects
    
    def _translate_codon(self, codon: str, dna_type: str) -> str:
        """Translate a single codon to amino acid."""
        from Bio.Seq import Seq
        from ..analyser import MITOCHONDRIAL_TABLE, STANDARD_TABLE
        
        if len(codon) != 3:
            return "?"
        
        table = MITOCHONDRIAL_TABLE if dna_type == "mitochondrial" else STANDARD_TABLE
        try:
            return str(Seq(codon).translate(table=table))
        except:
            return "?"
