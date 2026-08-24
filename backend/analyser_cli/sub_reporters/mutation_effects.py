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
            ReportFormatter.format_key_value("Genetic Code", analyser.genetic_code_description),
            ReportFormatter.format_key_value("Original Sequence Length", f"{len(analyser.input_sequence)} bp"),
            ReportFormatter.format_key_value("Number of Mutations", str(len(mutations))),
            "",
            ReportFormatter.format_detail_separator(),
            "MUTATION DETAILS:",
            ReportFormatter.format_detail_separator()
        ]
        
        codon_effects = self._analyze_codon_effects(analyser, mutations)
        
        for i, mutation in enumerate(mutations, 1):
            pos = mutation["position"]
            end_pos = mutation.get("end_position", pos)
            mut_type = mutation.get("type", "sub")
            position_display = str(pos) if pos == end_pos else f"{pos}_{end_pos}"

            if mut_type == "sub":
                original_base = mutation.get("ref") or str(analyser.input_sequence[pos - 1])
                new_base = mutation.get("alt", "?")
                change_str = f"{original_base} -> {new_base}"
            elif mut_type == "del":
                if pos == end_pos:
                    original_base = mutation.get("ref") or str(analyser.input_sequence[pos - 1])
                    change_str = f"Deletion of {original_base}"
                else:
                    deleted_count = end_pos - pos + 1
                    change_str = f"Deletion of {deleted_count} nucleotides"
            elif mut_type == "ins":
                inserted_seq = mutation.get("alt", "?")
                change_str = f"Insertion of {inserted_seq}"
            elif mut_type == "delins":
                inserted_seq = mutation.get("alt", "?")
                deleted_count = end_pos - pos + 1
                change_str = f"Replacement of {deleted_count} nucleotide(s) with {inserted_seq}"
            elif mut_type == "dup":
                duplicated_seq = mutation.get("alt", "?")
                change_str = f"Duplication of {duplicated_seq}"
            else:
                change_str = "Unknown mutation type"

            consequence = self._describe_indel_consequence(mutation)

            report_lines.extend([
                "",
                f"Mutation #{i} ({mut_type.upper()}):",
                ReportFormatter.format_key_value("Position", position_display, indent=2),
                ReportFormatter.format_key_value("Change", change_str, indent=2),
            ])

            if mut_type == "sub":
                report_lines.append(
                    ReportFormatter.format_key_value(
                        "Codon Position",
                        self._get_codon_position(pos, analyser.coding_offset),
                        indent=2,
                    )
                )
            elif consequence:
                report_lines.append(
                    ReportFormatter.format_key_value("Effect Type", consequence, indent=2)
                )
                report_lines.append(
                    ReportFormatter.format_key_value(
                        "Length Change", self._format_length_change(mutation), indent=2
                    )
                )

            if pos in codon_effects:
                effect = codon_effects[pos]
                report_lines.extend([
                    ReportFormatter.format_key_value(
                        "Codon",
                        f"{effect['original_codon']} -> {effect['mutated_codon']}",
                        indent=2,
                    ),
                    ReportFormatter.format_key_value(
                        "Amino Acid",
                        f"{effect['original_aa']} -> {effect['mutated_aa']}",
                        indent=2,
                    ),
                    ReportFormatter.format_key_value("Effect Type", effect['effect_type'], indent=2),
                ])

        mutated_length = len(analyser.mutated_amino_acid_translation) if analyser.mutated_amino_acid_translation else 0
        report_lines.extend([
            "",
            "",
            ReportFormatter.format_detail_separator(),
            "AMINO ACID TRANSLATION SUMMARY:",
            ReportFormatter.format_detail_separator(),
            ReportFormatter.format_key_value("Original", str(analyser.amino_acid_translation)),
            ReportFormatter.format_key_value("Mutated", str(analyser.mutated_amino_acid_translation), indent=1),
            "",
            ReportFormatter.format_key_value("Translation Length", f"{len(analyser.amino_acid_translation)} -> {mutated_length} amino acids")
        ])
        
        return "\n".join(report_lines)
    
    def _get_codon_position(self, nucleotide_pos: int, coding_offset: int = 0) -> str:
        """Determine codon number and position within codon."""
        offset_pos = nucleotide_pos - coding_offset
        if offset_pos < 1:
            return "5' UTR (before the first codon)"
        codon_num = (offset_pos - 1) // 3 + 1
        pos_in_codon = (offset_pos - 1) % 3 + 1
        return f"Codon {codon_num}, Position {pos_in_codon}"
    
    def _analyze_codon_effects(self, analyser: Analyser, mutations: List[Dict]) -> Dict:
        """Analyze the effect of mutations on codons and amino acids."""
        effects = {}
        
        original_seq = str(analyser.input_sequence)
        mutated_seq = str(analyser.mutated_sequence)
        coding_offset = analyser.coding_offset
        
        for mutation in mutations:
            if mutation.get("type", "sub") != "sub":
                continue

            pos = mutation["position"]
            if pos <= coding_offset:
                continue

            codon_start = ((pos - coding_offset - 1) // 3) * 3 + coding_offset
            codon_end = codon_start + 3
            
            if codon_end <= len(original_seq) and codon_end <= len(mutated_seq):
                original_codon = original_seq[codon_start:codon_end]
                mutated_codon = mutated_seq[codon_start:codon_end]
                
                original_aa = self._translate_codon(original_codon, analyser.codon_table)
                mutated_aa = self._translate_codon(mutated_codon, analyser.codon_table)
                
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
    
    def _translate_codon(self, codon: str, codon_table) -> str:
        """Translate a single codon to amino acid using the supplied codon table."""
        from Bio.Seq import Seq

        if len(codon) != 3:
            return "?"

        try:
            return str(Seq(codon).translate(table=codon_table))
        except:
            return "?"

    def _indel_length_change(self, mutation: Dict) -> int | None:
        mut_type = mutation.get("type", "sub")
        if mut_type == "sub":
            return None

        if mut_type == "ins":
            return len(mutation.get("alt", ""))

        if mut_type == "del":
            start = mutation["position"]
            end = mutation.get("end_position", start)
            return -(end - start + 1)

        if mut_type == "delins":
            start = mutation["position"]
            end = mutation.get("end_position", start)
            inserted = len(mutation.get("alt", ""))
            deleted = end - start + 1
            return inserted - deleted

        if mut_type == "dup":
            return len(mutation.get("alt", ""))

        return None

    def _describe_indel_consequence(self, mutation: Dict) -> str | None:
        length_change = self._indel_length_change(mutation)
        if length_change is None:
            return None
        if length_change % 3 == 0:
            return "In-frame indel"
        return "Frameshift mutation"

    def _format_length_change(self, mutation: Dict) -> str:
        length_change = self._indel_length_change(mutation)
        if length_change is None:
            return "0"
        sign = "+" if length_change > 0 else ""
        return f"{sign}{length_change} nt"
